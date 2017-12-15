#!/usr/bin/env python

"""
ChIP-seq pipeline
"""

from argparse import ArgumentParser
from functools import partial
import os
import re
import sys

import pandas as pd
import yaml

import pypiper
from pypiper import NGSTk, PipelineError, Stage
from pypiper.utils import \
	build_sample_paths, is_fastq, is_unzipped_fastq, is_gzipped_fastq, \
	parse_cores
from pep import AttributeDict, Sample
from const import CHIP_COMPARE_COLUMN, CHIP_MARK_COLUMN
from pipelines.exceptions import InvalidFiletypeException


__author__ = "Andre Rendeiro"
__copyright__ = "Copyright 2015, Andre Rendeiro"
__credits__ = []
__license__ = "GPL2"
__version__ = "0.3"
__maintainer__ = "Andre Rendeiro"
__email__ = "arendeiro@cemm.oeaw.ac.at"
__status__ = "Development"



# Default to 4 cores for each operation that supports core count specification.
parse_cores = partial(parse_cores, default=4)


# Allow the sample's 'ip' (or other attribute indicating the antiody target)
# to determine the peak calling mode, based on known characteristics of
# binding density patterns for different targets.
BROAD_MARKS = {
		"H3K9ME1", "H3K9ME2", "H3K9ME3",
		"H3K27ME1", "H3K27ME2", "H3K27ME3",
		"H3K36ME1", "H3K36ME2", "H3K36ME3",
		"H3K72ME1", "H3K72ME2", "H3K72ME3"}
HISTONE_CODES = ["H3", "H2A", "H2B", "H4"]



def add_cmdl_options(parser):
	"""
	Global options for pipeline.
	"""
	parser.add_argument(
		"--sample-yaml", dest="sample_config",
		help="Yaml config file with sample attributes; in addition to "
			"sample_name, this should define '{rt}', as 'single' or "
			"'paired'; 'ip', with the mark analyzed in a sample, and "
			"'{comparison}' with the name of a control sample (if the "
			"sample itself is not a control.)".format(
			rt="read_type", comparison=CHIP_COMPARE_COLUMN)
	)
	parser.add_argument(
		"--post-hoc-headcrop", type=int,
		help="Number of bases to trim from read starts after adapter "
			 "trimming itself")
	parser.add_argument(
		"--peak-caller", choices=["macs2", "spp"], default="macs2",
		help="Name of peak calling program.",
	)
	parser.add_argument(
		"--pvalue", type=float, default=0.00001, help="MACS2 p-value")
	parser.add_argument(
		"--qvalue", type=float, help="Q-value for peak calling")
	return parser



def merge_input(sample, pipeline_manager, ngstk):
	"""
	Handle multiple input files for a single sample by merging.

	Sometimes, a ChIP-seq sample will have a single input but other times we
	want to merge input files (e.g., replicates). A call to this function
	handles either case, acting only if the given sample has multiple input
	files defined. If so, this call designates the sample's 'unmapped'
	attribute (a filepath) as its 'data_source.'

	:param pep.Sample sample: the being processed, for which to
		merge input files if needed.
	:param pypiper.PipelineManager pipeline_manager: Handler of control flow
		and data management for a running pipeline.
	:param pypiper.NGSTk ngstk: collection of NGS utility functions
	"""

	pipeline_manager.timestamp("Merging input files")


	def ensure_merged():
		assert os.path.isfile(sample.merged_input_path), \
			"Upon completion of input merger step, sample's unmapped reads " \
			"file must exist."
		sample.data_source = sample.merged_input_path


	if len(sample.input_file_paths) < 2:
		ensure_merged()
		return

	# Validate filetype validity and homogeneity.
	first_input = sample.input_file_paths[0]  # Path to first input file.
	# Candidate functions for file type homogeneity check.
	file_type_funcs = [is_gzipped_fastq, is_unzipped_fastq,
					   lambda f: os.path.splitext(f)[1].lower() == ".bam",
					   lambda f: os.path.splitext(f)[1].lower() == ".sam"]
	# The boolean function with which to check all files is the first one
	# that evaluates to true when passed the first input file's path.
	for i, ft_func in enumerate(file_type_funcs):
		if ft_func(first_input):
			match = ft_func
			break
	else:
		raise InvalidFiletypeException(first_input)

	# See if we have any mismatches.
	mismatched = [f for f in sample.input_file_paths if not match(f)]
	if mismatched:
		raise ValueError("Not all input file types match that of '{}': {}".
						 format(first_input, mismatched))

	# Build merger command based on (homogeneous) type of input files.
	_, unmapped_ext = os.path.splitext(sample.unmapped)
	get_merge_cmd = ngstk.merge_fastq if is_fastq(
		first_input) else ngstk.merge_bams
	merge_target = sample.merged_input_path
	# TODO: note that this assumes R1 + R2 in same file for paired-end.
	cmd = get_merge_cmd(sample.input_file_paths, merge_target)
	pipeline_manager.run(cmd, target=merge_target, shell=True)
	ensure_merged()


def ensure_fastq(sample, pipeline_manager, ngstk):
	"""
	Convert a sequencing reads file to another format.

	:param pep.Sample sample: sample for which to convert reads file
	:param pypiper.PipelineManager pipeline_manager: execution manager
	:param pypiper.NGSTk ngstk: configured NGS processing framework;
		required if and only if conversion function is not provided.
	"""


	def link_fastq(real, link):
		if not os.path.isfile(sample.fastq):
			link = "ln -ns {} {}".format(real, link)
			print("Linking standard FASTQ path to the merged path.")
			pipeline_manager.run(link, target=sample.fastq, shell=True)
		else:
			print("Sample's single FASTQ exists.")


	# If we already have fastq, short-circuit and return early.
	if is_unzipped_fastq(sample.merged_input_path):
		link_fastq(real=sample.merged_input_path, link=sample.unmapped)
		link_fastq(real=sample.merged_input_path, link=sample.fastq)
		clean_files = []
	elif is_gzipped_fastq(sample.merged_input_path):
		assert not os.path.isfile(sample.unmapped), \
			"Decompressed output must not exist prior to decompression."
		unzip = "{} -d -c {} > {}".format(
			ngstk.ziptool, sample.merged_input_path, sample.unmapped)
		pipeline_manager.run(unzip, target=sample.unmapped, shell=True)
		link_fastq(real=sample.unmapped, link=sample.fastq)
		clean_files = []
	else:
		# Read type determines which paths to pass to the conversion call.
		if sample.paired:
			out_fq1 = sample.fastq1
			out_fq2 = sample.fastq2
			outfiles = [out_fq1, out_fq2]
			unpaired_fq = sample.fastq_unpaired
			clean_files = [sample.fastq1, sample.fastq2, sample.unpaired]
		else:
			out_fq1 = sample.fastq
			out_fq2 = None
			outfiles = [out_fq1]
			unpaired_fq = None
			clean_files = [sample.fastq]
		kwargs = {"input_bam": sample.data_source, "output_fastq": out_fq1,
				  "output_fastq2": out_fq2, "unpaired_fastq": unpaired_fq}

		# Convert BAM to FASTQ.
		pipeline_manager.timestamp("Converting from BAM to FASTQ")

		# Perform and validate the reads file format conversion.
		cmd = ngstk.bam2fastq(**kwargs)
		validate = ngstk.check_fastq(
			input_files=sample.input_file_paths,
			output_files=outfiles, paired_end=sample.paired)
		pipeline_manager.run(cmd, target=out_fq1, follow=validate,
							 shell=True)

	for f in clean_files:
		pipeline_manager.clean_add(f, conditional=True)
	sample.data_source = sample.unmapped


# TODO: why were we only reporting trimming stats for skewer, not trimmomatic?
def trim_reads(sample, pipeline_manager, ngstk,
			   cores=None, fastqc_folder="fastqc", post_hoc_headcrop=None):
	"""
	Perform read trimming.

	:param pep.Sample sample: sample for which to convert reads file
	:param pypiper.PipelineManager pipeline_manager: execution manager
	:param pypiper.NGSTk ngstk: configured NGS processing framework
	:param int | str cores: number of CPUs to allow for read trimming process
	:param str fastqc_folder: name of folder for fastqc output
	:param int | NoneType post_hoc_headcrop: number of bases to chop from
		read head after initial trimming, optional; if unspecified, no
		additional trimming is done; only compatible with skewer
	"""

	pipeline_manager.timestamp("Trimming adapters from sample")

	# Collect trimmer-agnostic keyword arguments.
	if sample.paired:
		in_fq1 = sample.fastq1
		in_fq2 = sample.fastq2
		out_fq1 = sample.trimmed1
		out_fq2 = sample.trimmed2
	else:
		in_fq1 = sample.fastq
		in_fq2 = None
		out_fq1 = sample.trimmed
		out_fq2 = None
	cores = parse_cores(cores, pipeline_manager)
	kwargs = {"input_fastq1": in_fq1, "input_fastq2": in_fq2,
			  "output_fastq1": out_fq1, "output_fastq2": out_fq2,
			  "cpus": cores,
			  "adapters": pipeline_manager.config.resources.adapters,
			  "log": sample.trimlog}

	# Collect trimmer-specific keyword arguments and files to clean.
	trimmer = pipeline_manager.config.parameters.trimmer
	if trimmer == "trimmomatic":
		build_trim_cmd = ngstk.trimmomatic
		out_fq1_unpaired = getattr(sample, "trimmed1_unpaired", None)
		out_fq2_unpaired = getattr(sample, "trimmed2_unpaired", None)
		trimmer_specific_kwargs = {
			"output_fastq1_unpaired": out_fq1_unpaired,
			"output_fastq2_unpaired": out_fq2_unpaired}
		clean_files = [sample.trimmed1, sample.trimmed1_unpaired,
					   sample.trimmed2, sample.trimmed2_unpaired] \
			if sample.paired else [sample.trimmed]
	elif trimmer == "skewer":
		build_trim_cmd = ngstk.skewer
		out_pre = os.path.join(sample.paths.unmapped, sample.sample_name)
		trimmer_specific_kwargs = {"output_prefix": out_pre}
		clean_files = [sample.trimmed1, sample.trimmed2] \
			if sample.paired else [sample.trimmed]
	else:
		raise ValueError("Unsupported read trimmer: {}".format(trimmer))

	# Create the command and the path to the output target.
	kwargs.update(trimmer_specific_kwargs)
	cmd = build_trim_cmd(**kwargs)
	if post_hoc_headcrop and trimmer == "trimmomatic":
		cmd += " HEADCROP:{}".format(post_hoc_headcrop)
	target = sample.trimmed1 if sample.paired else sample.trimmed

	# Run, clean, and report.
	run_fastqc = ngstk.check_trim(
		target, paired_end=sample.paired, trimmed_fastq_R2=out_fq2,
		fastqc_folder=fastqc_folder)
	pipeline_manager.run(cmd, target, follow=run_fastqc, shell=True)
	for f in clean_files:
		pipeline_manager.clean_add(f, conditional=True)
	trim_stat = _parse_qc_metrics(
		sample.trimlog, prefix="trim_", paired_end=sample.paired)
	_report_dict(pipeline_manager, trim_stat)


def align_reads(sample, pipeline_manager, ngstk, cores=None):
	"""
	Align sequencing reads.

	:param pep.Sample sample: sample for which to align reads
	:param pypiper.PipelineManager pipeline_manager: execution manager
	:param pypiper.NGSTk ngstk: configured NGS processing framework
	:param int | str cores: number of cores allowed to be used for alignment
	"""
	# Map
	pipeline_manager.timestamp("Mapping reads with Bowtie2")
	cmd = ngstk.bowtie2_map(
		input_fastq1=sample.trimmed1 if sample.paired else sample.trimmed,
		input_fastq2=sample.trimmed2 if sample.paired else None,
		output_bam=sample.mapped,
		log=sample.aln_rates,
		metrics=sample.aln_metrics,
		genome_index=getattr(pipeline_manager.config.resources.genomes,
							 sample.genome),
		max_insert=pipeline_manager.config.parameters.max_insert,
		cpus=parse_cores(cores, pipeline_manager)
	)
	pipeline_manager.run(cmd, sample.mapped, shell=True)
	mapstats = _parse_mapping_stats(sample.aln_rates,
								   paired_end=sample.paired)
	_report_dict(pipeline_manager, mapstats)


def filter_reads(sample, pipeline_manager, ngstk, cores=None):
	"""
	Filter sequencing reads.

	:param pep.Sample sample: sample for which to filter reads
	:param pypiper.PipelineManager pipeline_manager: execution manager
	:param pypiper.NGSTk ngstk: configured NGS processing framework
	:param int | str cores: number of cores allowed to be used for filtration
	"""
	# Filter reads
	pipeline_manager.timestamp("Filtering reads for quality")
	cmd = ngstk.filter_reads(
		input_bam=sample.mapped,
		output_bam=sample.filtered,
		metrics_file=sample.dups_metrics,
		paired=sample.paired,
		cpus=parse_cores(cores, pipeline_manager),
		Q=pipeline_manager.config.parameters.read_quality
	)
	pipeline_manager.run(cmd, sample.filtered, shell=True)
	filter_stats = _parse_sambamba_duplicate_file(sample.dups_metrics)
	_report_dict(pipeline_manager, filter_stats)


def post_align_fastqc(sample, pipeline_manager, ngstk):
	"""
	Post-alignment quality control.

	:param pep.Sample sample: sample for which to run QC
	:param pypiper.PipelineManager pipeline_manager: execution manager
	:param pypiper.NGSTk ngstk: configured NGS processing framework
	"""
	cmd = ngstk.fastqc(sample.filtered, sample.paths.fastqc_aligned)
	fastqc_name = "{}_fastqc-aligned.zip".format(sample.name)
	fastqc_target = os.path.join(sample.paths.fastqc_aligned, fastqc_name)
	pipeline_manager.run(cmd, target=fastqc_target, shell=True)


def index_bams(sample, pipeline_manager, ngstk):
	""" Index the alignment files. """
	pipeline_manager.timestamp("Indexing bamfiles with samtools")
	cmd = ngstk.index_bam(input_bam=sample.mapped)
	pipeline_manager.run(cmd, sample.mapped + ".bai", shell=True)
	cmd = ngstk.index_bam(input_bam=sample.filtered)
	pipeline_manager.run(cmd, sample.filtered + ".bai", shell=True)


def make_tracks(sample, pipeline_manager, ngstk):
	""" Create tracks for viewing in browser. """

	track_dir = os.path.dirname(sample.bigwig)
	if not os.path.exists(track_dir):
		os.makedirs(track_dir)

	# right now tracks are only made for bams without duplicates
	pipeline_manager.timestamp("Making bigWig tracks from bam file")
	cmd = _bam_to_bigwig(
		input_bam=sample.filtered,
		output_bigwig=sample.bigwig,
		genome_sizes=getattr(
			pipeline_manager.config.resources.chromosome_sizes,
			sample.genome),
		genome=sample.genome,
		tagmented=pipeline_manager.config.parameters.tagmented,
		# by default make extended tracks
		normalize=pipeline_manager.config.parameters.normalize_tracks,
		norm_factor=pipeline_manager.config.parameters.norm_factor
	)
	pipeline_manager.run(cmd, sample.bigwig, shell=True)


def compute_metrics(sample, pipeline_manager, ngstk, cores=None):
	"""
	Determine fragment size distribution, genome-wide coverage, and NSC/RSC.

	:param pep.Sample sample: the sample undergoing processing and
		being analyzed
	:param pypiper.PipelineManager pipeline_manager: overseer of resources
		and other settings for a pipeline
	:param pypiper.NGSTk ngstk: Suite of common NGS tools, configured with
		a pipeline manager
	:param int | str cores: number of cores available for use in this
		phase of the pipeline
	"""
	# Plot fragment distribution
	if sample.paired and not os.path.exists(sample.insertplot):
		pipeline_manager.timestamp("Plotting insert size distribution")
		ngstk.plot_atacseq_insert_sizes(
			bam=sample.filtered,
			plot=sample.insertplot,
			output_csv=sample.insertdata
		)

	# Count coverage genome-wide
	pipeline_manager.timestamp("Calculating genome-wide coverage")
	cmd = ngstk.genome_wide_coverage(
		input_bam=sample.filtered,
		genome_windows=getattr(
			pipeline_manager.config.resources.genome_windows,
			sample.genome),
		output=sample.coverage
	)
	pipeline_manager.run(cmd, sample.coverage, shell=True)

	# Calculate NSC, RSC
	pipeline_manager.timestamp("Assessing signal/noise in sample")
	cmd = ngstk.run_spp(
		input_bam=sample.filtered,
		output=sample.qc,
		plot=sample.qc_plot,
		cpus=parse_cores(cores, pipeline_manager)
	)
	pipeline_manager.run(cmd, sample.qc_plot, shell=True, nofail=True)
	nsc_rsc_stats = _parse_nsc_rsc(sample.qc)
	_report_dict(pipeline_manager, nsc_rsc_stats)


def wait_for_control(sample, pipeline_manager, ngstk):
	""" Waiting for a ChIP sample's corresponding input DNA to complete. """
	comparison = sample.background_sample_name
	# The pipeline will now wait for the comparison sample file to be completed
	path_block_file = sample.filtered.replace(sample.name, comparison)
	pipeline_manager._wait_for_file(path_block_file)


# TODO: need to receive comparison sample name and args from caller.
# TODO: spin off the plotting functionality into its own stage, or use a
# TODO (cont.): substage mechanism if implemented
def call_peaks(sample, pipeline_manager, ngstk, cores=None, caller=None,
			   **caller_kwargs):
	"""
	Call peaks.

	:param pep.Sample sample: the sample for which to call peaks
	:param pypiper.PipelineManager pipeline_manager: the manager for a
		pipeline instance
	:param pypiper.NGSTk ngstk: configured NGS functions framework
	:param str | int cores: number of processing cores available for use
	:param str caller: name of the peak caller to use
	"""

	comparison = sample.background_sample_name

	# Call peaks.
	broad_mode = sample.broad
	peaks_folder = sample.paths.peaks
	treatment_file = sample.filtered
	control_file = sample.filtered.replace(sample.name, comparison)

	if not os.path.exists(peaks_folder):
		os.makedirs(peaks_folder)
	# TODO: include the filepaths as caller-neutral positionals/keyword args
	# TODO (cont.) once NGSTK API is tweaked.

	peak_call_kwargs = {"output_dir": peaks_folder,
						"broad": broad_mode,
						"qvalue": caller_kwargs.get("qvalue")}

	# Allow SPP but lean heavily toward MACS2.
	caller = caller or getattr(pipeline_manager, "peak_caller", "macs2")
	if caller not in ["spp", "macs2"]:
		raise ValueError("Unsupported peak caller: {}".format(caller))
	if caller == "spp":
		cmd = ngstk.spp_call_peaks(
			treatment_bam=treatment_file, control_bam=control_file,
			treatment_name=sample.name, control_name=comparison,
			cpus=parse_cores(cores, pipeline_manager), **peak_call_kwargs)
	else:
		cmd = ngstk.macs2_call_peaks(
			treatment_bams=treatment_file, control_bams=control_file,
			sample_name=sample.name, pvalue=caller_kwargs.get("pvalue"),
			genome=sample.genome, paired=sample.paired, **peak_call_kwargs)

	pipeline_manager.run(cmd, target=sample.peaks, shell=True)
	num_peaks = _parse_peak_count(sample.peaks)
	pipeline_manager.report_result("peaks", num_peaks)

	if caller == "macs2" and not broad_mode:
		pipeline_manager.timestamp("Plotting MACS2 model")
		model_files_base = sample.name + "_model"

		# Create the command to run the model script.
		name_model_script = model_files_base + ".r"
		path_model_script = os.path.join(peaks_folder, name_model_script)
		exec_model_script = \
			"{} {}".format(pipeline_manager.config.tools.Rscript,
						   path_model_script)

		# Create the command to create and rename the model plot.
		plot_name = model_files_base + ".pdf"
		src_plot_path = os.path.join(os.getcwd(), plot_name)
		dst_plot_path = os.path.join(peaks_folder, plot_name)
		rename_model_plot = "mv {} {}".format(src_plot_path, dst_plot_path)

		# Run the model script and rename the model plot.
		pipeline_manager.run([exec_model_script, rename_model_plot],
							 target=dst_plot_path, shell=True, nofail=True)


def calc_frip(sample, pipeline_manager, ngstk, cores=None):
	"""
	Calculate the fraction of reads in called peaks (FRIP).

	:param pep.Sample sample: the sample for which to call peaks
	:param pypiper.PipelineManager pipeline_manager: the manager for a
		pipeline instance
	:param pypiper.NGSTk ngstk: configured NGS functions framework
	:param str | int cores: number of processing cores available for use, 
		optional; if unspecified, use a single core.
	"""
	pipeline_manager.timestamp(
		"Calculating fraction of reads in peaks (FRiP)")
	cmd = ngstk.calculate_frip(
		input_bam=sample.filtered,
		input_bed=sample.peaks,
		output=sample.frip,
		cpus=parse_cores(cores, pipeline_manager)
	)
	pipeline_manager.run(cmd, sample.frip, shell=True)
	reads_SE = float(pipeline_manager.get_stat("filtered_single_ends"))
	reads_PE = float(pipeline_manager.get_stat("filtered_paired_ends"))
	total = 0.5 * (reads_SE + reads_PE)
	if not total or pd.np.isnan(total):
		frip = pd.np.nan
	else:
		frip = _parse_frip(sample.frip, total)
	pipeline_manager.report_result("frip", frip)



class ChIPseqSample(Sample):
	"""
	Class to model ChIP-seq samples based on the generic Sample class.

	:param series: Collection of sample attributes.
	:type series: Mapping | pandas.Series

	:Example:

	# create Samples through a project object
	from pep import Project
	prj = Project("project_config.yaml")
	prj.add_sample_sheet()
	s0 = prj.samples[0]  # here's a Sample

	# create Samples through a SampleSheet object
	from pep import SampleSheet, Sample
	sheet = SampleSheet("project_sheet.csv")
	s1 = Sample(sheet.ix[0])  # here's a Sample too
	"""

	__library__ = "ChIP-seq"


	def __init__(self, series):
		super(ChIPseqSample, self).__init__(series)
		self.tagmented = False

		# Set broad/histone status that may later be modified given
		# context of a pipeline configuration file, handling null/missing mark.
		mark = getattr(self, CHIP_MARK_COLUMN, None)
		if mark is None:
			self.broad = False
			self.histone = False
		else:
			mark = mark.upper()
			self.broad = mark in BROAD_MARKS
			self.histone = any([mark.startswith(histone_code)
								for histone_code in HISTONE_CODES])


	@property
	def background_sample_name(self):
		"""
		Determine the name of the background sample for this sample.

		:return str | NoneType: name of sample used as background/control
			for this one, perhaps null if this sample itself represents a
			background/control.
		"""
		return getattr(self, CHIP_COMPARE_COLUMN)


	@property
	def is_control(self):
		"""
		Determine whether this sample is a ChIP background/control.

		:return bool: whether this sample is a ChIP background/control
		"""
		return self.background_sample_name in [None, "", "NA"]


	def set_file_paths(self, project=None):
		"""
		Sets the paths of all files for this sample.

		:param AttributeDict project: Project or its metadata, from the
			Project instance that generated this sample
		"""

		# Get paths container structure and any contents used by any Sample.
		super(ChIPseqSample, self).set_file_paths(project)

		# Files in the root of the sample dir
		self.trimlog = os.path.join(self.paths.sample_root, self.name + ".trimlog.txt")
		self.aln_rates = os.path.join(self.paths.sample_root, self.name + ".aln_rates.txt")
		self.aln_metrics = os.path.join(self.paths.sample_root, self.name + ".aln_metrics.txt")
		self.dups_metrics = os.path.join(self.paths.sample_root, self.name + ".dups_metrics.txt")

		# Unmapped: merged bam, fastq, trimmed fastq
		self.paths.unmapped = os.path.join(self.paths.sample_root, "unmapped")
		first_input = self.input_file_paths[0]
		if is_fastq(first_input):
			if self.paired:
				raise PipelineError(
					"For paired-end reads, only BAM input is supported.")
			unmapped_ext = ".fastq"
			merge_ext = ".fastq.gz" if is_gzipped_fastq(first_input) \
					else ".fastq"
		else:
			unmapped_ext = ".bam"
			merge_ext = ".bam"

		if len(self.input_file_paths) > 1:
			merged_name = self.name + merge_ext
			merged_path = os.path.join(self.paths.sample_root, merged_name)
		else:
			merged_path = self.input_file_paths[0]
		self.merged_input_path = merged_path
		unmapped_name = self.name + unmapped_ext
		self.unmapped = os.path.join(self.paths.unmapped, unmapped_name)

		self.paths.fastqc = os.path.join(self.paths.sample_root, "fastqc")
		self.paths.fastqc_aligned = os.path.join(
			self.paths.sample_root, "aligned_fastqc")

		self.fastq = os.path.join(self.paths.unmapped, self.name + ".fastq")
		self.fastq1 = os.path.join(self.paths.unmapped, self.name + ".1.fastq")
		self.fastq2 = os.path.join(self.paths.unmapped, self.name + ".2.fastq")
		self.fastq_unpaired = os.path.join(self.paths.unmapped, self.name + ".unpaired.fastq")
		self.trimmed = os.path.join(self.paths.unmapped, self.name + ".trimmed.fastq")
		self.trimmed1 = os.path.join(self.paths.unmapped, self.name + ".1.trimmed.fastq")
		self.trimmed2 = os.path.join(self.paths.unmapped, self.name + ".2.trimmed.fastq")
		self.trimmed1_unpaired = os.path.join(self.paths.unmapped, self.name + ".1_unpaired.trimmed.fastq")
		self.trimmed2_unpaired = os.path.join(self.paths.unmapped, self.name + ".2_unpaired.trimmed.fastq")

		# Mapped: mapped, duplicates marked, removed, reads shifted
		self.paths.mapped = os.path.join(self.paths.sample_root, "mapped")
		self.mapped = os.path.join(self.paths.mapped, self.name + ".trimmed.bowtie2.bam")
		self.filtered = os.path.join(self.paths.mapped, self.name + ".trimmed.bowtie2.filtered.bam")

		# Files in the root of the sample dir
		self.frip = os.path.join(self.paths.sample_root, self.name + "_FRiP.txt")

		# Coverage: read coverage in windows genome-wide
		self.paths.coverage = os.path.join(self.paths.sample_root, "coverage")
		self.coverage = os.path.join(self.paths.coverage, self.name + ".cov")

		self.insertplot = os.path.join(self.paths.sample_root, self.name + "_insertLengths.pdf")
		self.insertdata = os.path.join(self.paths.sample_root, self.name + "_insertLengths.csv")

		self.qc = os.path.join(self.paths.sample_root, self.name + "_qc.tsv")
		self.qc_plot = os.path.join(self.paths.sample_root, self.name + "_qc.pdf")

		bigwig_subfolder = "bigwig_{}".format(self.genome)
		bigwig_folder = os.path.join(
				self.prj.metadata.results_subdir, self.name, bigwig_subfolder)
		bigwig_file = "CHIP_{}.bw".format(self.name)
		self.bigwig = os.path.join(bigwig_folder, bigwig_file)

		# Peak-related paths
		self.paths.peaks = os.path.join(self.paths.sample_root, "peaks")
		type_peaks_file = "broadPeak" if self.broad else "narrowPeak"
		peaks_fname = "{}_peaks.{}".format(self.name, type_peaks_file)
		self.peaks = os.path.join(self.paths.peaks, peaks_fname)



class ChIPmentationSample(ChIPseqSample):
	"""
	Class to model ChIPmentation samples based on the ChIPseqSample class.

	:param series: Collection of sample attributes.
	:type series: Mapping | pandas.Series
	"""

	__library__ = "ChIPmentation"


	def __init__(self, series):
		super(ChIPmentationSample, self).__init__(series)
		self.tagmented = True



class ChipseqPipeline(pypiper.Pipeline):
	""" Definition of the ChIP-seq pipeline in terms of processing stages. """


	def __init__(self, sample, manager, peak_caller,
				 cores=None, post_hoc_headcrop=None, **caller_kwargs):
		"""
		Define the pipeline instance with a sample and manager.
		
		Optionally, a peak caller name, number of processing cores, and peak
		calling keyword arguments.

		:param pep.Sample sample: the sample being processed by
			the pipeline
		:param pypiper.manager.PipelineManager manager: the collection of
			monitoring and resource specifications for the pipeline
		:param str peak_caller: name of peak caller to use
		:param str | int cores: number of cores to use for each pipeline
			operation that supports core count specification, optional;
			if this isn't specified, the value associated with the manager
			or the default for this module will be used.
		:param int | NoneType post_hoc_headcrop: number of bases to trim from
			head of read starts after initial trimming; if unspecified, no
			additional read head trimming is done.
		:param dict caller_kwargs: extra keyword arguments for peak calling
		"""

		pipe_name, _ = os.path.splitext(os.path.split(__file__)[1])

		# Essential pipeline attributes
		self.sample = sample
		self.manager = manager

		# The pipeline manager contextualizes the NGSTk, providing the
		# configuration data as well as I/O and logging infrastructure,
		# and a framework within which to run relevant commands.
		self.ngstk = NGSTk(pm=manager)

		self.peak_caller = peak_caller
		self._kwargs_by_func = {
			"call_peaks": caller_kwargs,
			"trim_reads": {"post_hoc_headcrop": post_hoc_headcrop}}

		super(ChipseqPipeline, self).__init__(pipe_name, manager)

		# Pass cores to stages via manager.
		if cores is not None:
			self.manager.cores = cores

		print("Ensuring necessary paths for sample processing")
		build_sample_paths(sample)


	def stages(self):
		"""
		The processing stages/phases that define the ChIP-seq pipeline.
		Gq
		:return list[pypiper.Stage]: sequence of pipeline phases/stages
		"""
		always = [merge_input, ensure_fastq,
				  trim_reads, align_reads, filter_reads,
				  post_align_fastqc,
				  index_bams, make_tracks, compute_metrics]
		treatment_only = [wait_for_control, call_peaks, calc_frip]
		f_args = (self.sample, self.manager, self.ngstk)
		funcs = always if self.sample.is_control else always + treatment_only
		return [Stage(f, f_args, self._kwargs_by_func.get(f.__name__, {}))
				for f in funcs]


	def wrapup(self):
		""" Print completion message before stopping pipeline. """
		print("Finished processing sample {}.".format(self.sample.name))
		super(ChipseqPipeline, self).wrapup()



# TODO: remove and use the pypiper version once it supports normalization factor.
def _bam_to_bigwig(input_bam, output_bigwig, genome_sizes, genome,
				   tagmented=False, normalize=False, norm_factor=1000000):
	import os

	# TODO:
	# Adjust fragment length dependent on read size and real fragment size
	# (right now it assumes 50bp reads with 180bp fragments.)
	cmds = list()
	transient_file = os.path.abspath(re.sub("\.bigWig", "", output_bigwig))
	cmd1 = "bedtools bamtobed -i {0} |".format(input_bam)
	if not tagmented:
		cmd1 += " " + "bedtools slop -i stdin -g {0} -s -l 0 -r 130 |".format(
			genome_sizes)
		bedfile_bounds_script = os.path.join(
			os.path.dirname(__file__), "tools",
			"fix_bedfile_genome_boundaries.py")
		cmd1 += " {0} {1} |".format(bedfile_bounds_script, genome)
	cmd1 += " " + "genomeCoverageBed {0}-bg -g {1} -i stdin > {2}.cov".format(
		"-5 " if tagmented else "",
		genome_sizes,
		transient_file
	)
	cmds.append(cmd1)
	if normalize:
		cmds.append(
			"""awk 'NR==FNR{{sum+= $4; next}}{{ $4 = ($4 / sum) * {1}; print}}' {0}.cov {0}.cov | sort -k1,1 -k2,2n > {0}.normalized.cov""".format(
				transient_file, norm_factor))
	cmds.append("bedGraphToBigWig {0}{1}.cov {2} {3}".format(transient_file,
															 ".normalized" if normalize else "",
															 genome_sizes,
															 output_bigwig))
	# remove tmp files
	cmds.append(
		"if [[ -s {0}.cov ]]; then rm {0}.cov; fi".format(transient_file))
	if normalize:
		cmds.append(
			"if [[ -s {0}.normalized.cov ]]; then rm {0}.normalized.cov; fi".format(
				transient_file))
	cmds.append("chmod 755 {0}".format(output_bigwig))
	return cmds



def _report_dict(pipe, stats_dict):
	"""
	Convenience wrapper to report a collection of pipeline results.

	This writes a collection of key-value pairs to the central stats/results
	file associated with the pipeline manager provided.

	:param pipe: Pipeline manager with which to do the reporting
	:type pipe: pypiper.PipelineManager
	:param stats_dict: Collection of results, each mapped to a name
	:type stats_dict: Mapping[str, object]
	"""
	for key, value in stats_dict.items():
		pipe.report_result(key, value)



def _parse_qc_metrics(stats_file, prefix="", paired_end=True):
	"""
	Parse and write some QC metrics.

	:param stats_file: sambamba output file with duplicate statistics.
	:type stats_file: str
	:param prefix: A string to be used as prefix to the output dictionary keys.
	:type stats_file: str
	:param paired_end: whether paired-end sequencing was used to generate
		reads for a sample being processed
	:type paired_end: bool
	"""

	error_dict = {
		prefix + "surviving_perc": pd.np.nan,
		prefix + "short_perc": pd.np.nan,
		prefix + "empty_perc": pd.np.nan,
		prefix + "trimmed_perc": pd.np.nan,
		prefix + "untrimmed_perc": pd.np.nan,
		prefix + "trim_loss_perc": pd.np.nan}
	try:
		with open(stats_file) as handle:
			content = handle.readlines()  # list of strings per line
	except:
		return error_dict

	suf = "s" if not paired_end else " pairs"

	try:
		line = [i for i in range(len(content)) if
				"read{} processed; of these:".format(suf) in content[i]][0]
		total = int(re.sub("\D", "", re.sub("\(.*", "", content[line])))
		line = [i for i in range(len(content)) if
				"read{} available; of these:".format(suf) in content[i]][0]
		surviving = int(re.sub("\D", "", re.sub("\(.*", "", content[line])))
		line = [i for i in range(len(content)) if
				"short read{} filtered out after trimming by size control".format(
					suf) in content[i]][0]
		short = int(re.sub(" \(.*", "", content[line]).strip())
		line = [i for i in range(len(content)) if
				"empty read{} filtered out after trimming by size control".format(
					suf) in content[i]][0]
		empty = int(re.sub(" \(.*", "", content[line]).strip())
		line = [i for i in range(len(content)) if
				"trimmed read{} available after processing".format(suf) in
				content[i]][0]
		trimmed = int(re.sub(" \(.*", "", content[line]).strip())
		line = [i for i in range(len(content)) if
				"untrimmed read{} available after processing".format(suf) in
				content[i]][0]
		untrimmed = int(re.sub(" \(.*", "", content[line]).strip())
		return {
			prefix + "surviving_perc": (float(surviving) / total) * 100,
			prefix + "short_perc": (float(short) / total) * 100,
			prefix + "empty_perc": (float(empty) / total) * 100,
			prefix + "trimmed_perc": (float(trimmed) / total) * 100,
			prefix + "untrimmed_perc": (float(untrimmed) / total) * 100,
			prefix + "trim_loss_perc": ((total - float(
				surviving)) / total) * 100}
	except IndexError:
		return error_dict



def _parse_mapping_stats(stats_file, prefix="", paired_end=True):
	"""
	Determine and write to stats file various alignment metrics.

	:param stats_file: sambamba output file with duplicate statistics.
	:type stats_file: str
	:param prefix: A string to be used as prefix to the output dictionary keys.
	:type stats_file: str
	:param paired_end: whether to determine the alignment metrics according
		to the notion that the sequencing was done in paired-end fashion.
	:type paired_end: bool
	"""

	if not paired_end:
		error_dict = {
			prefix + "not_aligned_perc": pd.np.nan,
			prefix + "unique_aligned_perc": pd.np.nan,
			prefix + "multiple_aligned_perc": pd.np.nan,
			prefix + "perc_aligned": pd.np.nan}
	else:
		error_dict = {
			prefix + "paired_perc": pd.np.nan,
			prefix + "concordant_perc": pd.np.nan,
			prefix + "concordant_unique_perc": pd.np.nan,
			prefix + "concordant_multiple_perc": pd.np.nan,
			prefix + "not_aligned_or_discordant_perc": pd.np.nan,
			prefix + "not_aligned_perc": pd.np.nan,
			prefix + "unique_aligned_perc": pd.np.nan,
			prefix + "multiple_aligned_perc": pd.np.nan,
			prefix + "perc_aligned": pd.np.nan}

	try:
		with open(stats_file) as handle:
			content = handle.readlines()  # list of strings per line
	except:
		return error_dict

	if not paired_end:
		try:
			line = [i for i in range(len(content)) if
					"reads; of these:" in content[i]][0]
			total = int(re.sub("\D", "", re.sub("\(.*", "", content[line])))
			line = [i for i in range(len(content)) if
					"aligned 0 times" in content[i]][0]
			not_aligned_perc = float(
				re.search("\(.*%\)", content[line]).group()[1:-2])
			line = [i for i in range(len(content)) if
					" aligned exactly 1 time" in content[i]][0]
			unique_aligned_perc = float(
				re.search("\(.*%\)", content[line]).group()[1:-2])
			line = [i for i in range(len(content)) if
					" aligned >1 times" in content[i]][0]
			multiple_aligned_perc = float(
				re.search("\(.*%\)", content[line]).group()[1:-2])
			line = [i for i in range(len(content)) if
					"overall alignment rate" in content[i]][0]
			perc_aligned = float(re.sub("%.*", "", content[line]).strip())
			return {
				prefix + "not_aligned_perc": not_aligned_perc,
				prefix + "unique_aligned_perc": unique_aligned_perc,
				prefix + "multiple_aligned_perc": multiple_aligned_perc,
				prefix + "perc_aligned": perc_aligned}
		except IndexError:
			return error_dict

	if paired_end:
		try:
			line = [i for i in range(len(content)) if
					"reads; of these:" in content[i]][0]
			total = int(re.sub("\D", "", re.sub("\(.*", "", content[line])))
			line = [i for i in range(len(content)) if
					" were paired; of these:" in content[i]][0]
			paired_perc = float(
				re.search("\(.*%\)", content[line]).group()[1:-2])
			line = [i for i in range(len(content)) if
					"aligned concordantly 0 times" in content[i]][0]
			concordant_unaligned_perc = float(
				re.search("\(.*%\)", content[line]).group()[1:-2])
			line = [i for i in range(len(content)) if
					"aligned concordantly exactly 1 time" in content[i]][0]
			concordant_unique_perc = float(
				re.search("\(.*%\)", content[line]).group()[1:-2])
			line = [i for i in range(len(content)) if
					"aligned concordantly >1 times" in content[i]][0]
			concordant_multiple_perc = float(
				re.search("\(.*%\)", content[line]).group()[1:-2])
			line = [i for i in range(len(content)) if
					"mates make up the pairs; of these:" in content[i]][0]
			not_aligned_or_discordant = int(
				re.sub("\D", "", re.sub("\(.*", "", content[line])))
			d_fraction = (not_aligned_or_discordant / float(total))
			not_aligned_or_discordant_perc = d_fraction * 100
			line = [i for i in range(len(content)) if
					"aligned 0 times\n" in content[i]][0]
			not_aligned_perc = float(
				re.search("\(.*%\)", content[line]).group()[1:-2]) * d_fraction
			line = [i for i in range(len(content)) if
					" aligned exactly 1 time" in content[i]][0]
			unique_aligned_perc = float(
				re.search("\(.*%\)", content[line]).group()[1:-2]) * d_fraction
			line = [i for i in range(len(content)) if
					" aligned >1 times" in content[i]][0]
			multiple_aligned_perc = float(
				re.search("\(.*%\)", content[line]).group()[1:-2]) * d_fraction
			line = [i for i in range(len(content)) if
					"overall alignment rate" in content[i]][0]
			perc_aligned = float(re.sub("%.*", "", content[line]).strip())
			return {
				prefix + "paired_perc": paired_perc,
				prefix + "concordant_unaligned_perc": concordant_unaligned_perc,
				prefix + "concordant_unique_perc": concordant_unique_perc,
				prefix + "concordant_multiple_perc": concordant_multiple_perc,
				prefix + "not_aligned_or_discordant_perc": not_aligned_or_discordant_perc,
				prefix + "not_aligned_perc": not_aligned_perc,
				prefix + "unique_aligned_perc": unique_aligned_perc,
				prefix + "multiple_aligned_perc": multiple_aligned_perc,
				prefix + "perc_aligned": perc_aligned}
		except IndexError:
			return error_dict



def _parse_sambamba_duplicate_file(stats_file, prefix=""):
	"""
	Parses sambamba markdup output, returns series with values.

	:param stats_file: sambamba output file with duplicate statistics.
	:type stats_file: str
	:param prefix: A string to be used as prefix to the output dictionary keys.
	:type stats_file: str
	"""

	error_dict = {
		prefix + "filtered_single_ends": pd.np.nan,
		prefix + "filtered_paired_ends": pd.np.nan,
		prefix + "duplicate_percentage": pd.np.nan}
	try:
		with open(stats_file) as handle:
			content = handle.readlines()  # list of strings per line
	except:
		return error_dict

	try:

		# Note that each comprehension assumes exactly 1 match for the line to parse.

		# First, parse number of single-end reads.
		line = [i for i in range(len(content)) if
				"single ends (among them " in content[i]][0]
		single_ends = int(re.sub("\D", "", re.sub("\(.*", "", content[line])))

		# Then, parse number of paired-end reads.
		# Version of sambamba before 0.6.6, maybe 0.5?
		# line = [i for i in range(len(content)) if " end pairs...   done in " in content[i]][0]
		line = [i for i in range(len(content)) if
				"sorted" in content[i] and "end pairs" in content[i]][0]
		paired_ends = int(
			re.sub("\D", "", re.sub("\.\.\..*", "", content[line])))

		# Parse the duplicates.
		# Version of sambamba before 0.6.6, maybe 0.5?
		# line = [i for i in range(len(content)) if " duplicates, sorting the list...   done in " in content[i]][0]
		line = [i for i in range(len(content)) if
				"found" in content[i] and "duplicates" in content[i]][0]
		duplicates = int(
			re.sub("\D", "", re.sub("\.\.\..*", "", content[line])))

		return {
			prefix + "filtered_single_ends": single_ends,
			prefix + "filtered_paired_ends": paired_ends,
			prefix + "duplicate_percentage": (float(duplicates) / (
			single_ends + paired_ends * 2)) * 100}

	except IndexError:
		return error_dict



def _parse_peak_count(peak_file):
	"""
	Count the called peaks.

	:param str peak_file: Path to called peaks file.
	:return: int | np.nan: Called peak count; null (np.nan) if an exception
		occurs during peak counting.
	"""
	from subprocess import check_output

	try:
		return int(check_output(["wc", "-l", peak_file]).split(" ")[0])
	except:
		return pd.np.nan



def _parse_frip(frip_file, total_reads):
	"""
	Calculates the fraction of reads in peaks for a given sample.

	:param frip_file: A path to a file with the FRiP output.
	:type frip_file: str
	:param total_reads: Number of total reads (i.e., the denominator for the
		FRiP calculation)
	:type total_reads: int | float
	:rtype: float
	"""

	try:
		with open(frip_file, "r") as handle:
			content = handle.readlines()
	except:
		print("Failed to read FRIP file: '{}'".format(frip_file))
		return pd.np.nan

	if content[0].strip() == "":
		print("Empty first line of FRIP file: '{}'".format(frip_file))
		return pd.np.nan

	reads_in_peaks = int(re.sub("\D", "", content[0]))
	frip = reads_in_peaks / float(total_reads)

	return frip



def _parse_nsc_rsc(nsc_rsc_file):
	"""
	Parses the values of NSC and RSC from a stats file.

	:param nsc_rsc_file: A sting path to a file with the NSC and RSC output (generally a tsv file).
	:type nsc_rsc_file: str
	"""
	try:
		nsc_rsc = pd.read_csv(nsc_rsc_file, header=None, sep="\t")
		return {"NSC": nsc_rsc[8].squeeze(), "RSC": nsc_rsc[9].squeeze()}
	except:
		return {"NSC": pd.np.nan, "RSC": pd.np.nan}



def main():
	""" Run the pipeline. """

	# Parse command-line arguments
	parser = ArgumentParser(
		prog="chipseq-pipeline",
		description="ChIP-seq pipeline."
	)
	parser = add_cmdl_options(parser)

	print("Adding pypiper arguments")
	parser = pypiper.add_pypiper_args(parser, all_args=True)
	args = parser.parse_args()
	if args.sample_config is None:
		parser.print_help()
		return 1

	# Read in yaml configs
	sample = Sample(yaml.load(open(args.sample_config, "r")))
	# Create Sample object

	if sample.protocol == "ChIPmentation":
		sample = ChIPmentationSample(sample)
	else:
		sample = ChIPseqSample(sample)

	# Check if merged
	if len(sample.input_file_paths) > 1:
		sample.merged = True
	else:
		sample.merged = False
	sample.prj = AttributeDict(sample.prj)
	sample.paths = AttributeDict(sample.paths.__dict__)

	# Flag version of read type since it's a binary; handle case vagaries.
	try:
		sample.paired = (sample.read_type.lower() == "paired")
	except AttributeError:
		print("WARNING: non-string read_type: {} ({})".format(
			sample.read_type, type(sample.read_type)))
		sample.paired = False

	# Set file paths
	sample.set_file_paths(sample.prj)
	# sample.make_sample_dirs()  # should be fixed to check if values of paths are strings and paths indeed

	# Start Pypiper object
	# Best practice is to name the pipeline with the name of the script;
	# or put the name in the pipeline interface.
	pl_mgr = pypiper.PipelineManager(
		name="chipseq", outfolder=sample.paths.sample_root, args=args)

	# With the sample and the manager created, we're ready to run the pipeline.
	pipeline = ChipseqPipeline(
		sample, pl_mgr, peak_caller=args.peak_caller,
		pvalue=args.pvalue, qvalue=args.qvalue,
		post_hoc_headcrop=args.post_hoc_headcrop)
	pipeline.run(start=args.start,
				 stop_at=args.stop_at, stop_after=args.stop_after)



if __name__ == '__main__':
	try:
		sys.exit(main())
	except KeyboardInterrupt:
		print("Program canceled by user!")
		sys.exit(1)
