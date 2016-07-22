#!/usr/bin/env python

"""
ChIP-seq pipeline
"""

import os
import sys
from argparse import ArgumentParser
import yaml
import pypiper
from pypiper.ngstk import NGSTk
from looper.models import AttributeDict, Sample

import pandas as pd


__author__ = "Andre Rendeiro"
__copyright__ = "Copyright 2015, Andre Rendeiro"
__credits__ = []
__license__ = "GPL2"
__version__ = "0.3"
__maintainer__ = "Andre Rendeiro"
__email__ = "arendeiro@cemm.oeaw.ac.at"
__status__ = "Development"


class ChIPseqSample(Sample):
	"""
	Class to model ChIP-seq samples based on the generic Sample class (itself a pandas.Series).

	:param series: Pandas `Series` object.
	:type series: pandas.Series

	:Example:

	# create Samples through a project object
	from looper.models import Project
	prj = Project("project_config.yaml")
	prj.add_sample_sheet()
	s0 = prj.samples[0]  # here's a Sample

	# create Samples through a SampleSheet object
	from looper.models import SampleSheet, Sample
	sheet = SampleSheet("project_sheet.csv")
	s1 = Sample(sheet.ix[0])  # here's a Sample too
	"""
	__library__ = "ChIP-seq"

	def __init__(self, series):
		# Passed series must either be a pd.Series or a daughter class
		if not isinstance(series, pd.Series):
			raise TypeError("Provided object is not a pandas Series.")
		super(ChIPseqSample, self).__init__(series)

		self.tagmented = False

		# Get type of factor
		# TODO: get config file specifying broad/narrow factors
		# e.g. self.broad = True if self.ip in self.prj.config["broadfactors"] else False
		# NS: I wrapped this in a try block because this makes it require that
		# 'ip' be defined, which may not be the case (like for Input samples)

		try:
			self.broad = True if any([ip in self.ip.upper() for ip in ["H3K27me3", "H3K36me3"]]) else False
			self.histone = True if any([ip in self.ip.upper() for ip in ["H3", "H2A", "H2B", "H4"]]) else False
		except:
			pass

	def __repr__(self):
		return "ChIP-seq sample '%s'" % self.sample_name

	def set_file_paths(self):
		"""
		Sets the paths of all files for this sample.
		"""
		# Inherit paths from Sample by running Sample's set_file_paths()
		super(ChIPseqSample, self).set_file_paths()

		# Files in the root of the sample dir
		self.fastqc = os.path.join(self.paths.sample_root, self.sample_name + ".fastqc.zip")
		self.trimlog = os.path.join(self.paths.sample_root, self.sample_name + ".trimlog.txt")
		self.aln_rates = os.path.join(self.paths.sample_root, self.sample_name + ".aln_rates.txt")
		self.aln_metrics = os.path.join(self.paths.sample_root, self.sample_name + ".aln_metrics.txt")
		self.dups_metrics = os.path.join(self.paths.sample_root, self.sample_name + ".dups_metrics.txt")

		# Unmapped: merged bam, fastq, trimmed fastq
		self.paths.unmapped = os.path.join(self.paths.sample_root, "unmapped")
		self.unmapped = os.path.join(self.paths.unmapped, self.sample_name + ".bam")
		self.fastq = os.path.join(self.paths.unmapped, self.sample_name + ".fastq")
		self.fastq1 = os.path.join(self.paths.unmapped, self.sample_name + ".1.fastq")
		self.fastq2 = os.path.join(self.paths.unmapped, self.sample_name + ".2.fastq")
		self.fastq_unpaired = os.path.join(self.paths.unmapped, self.sample_name + ".unpaired.fastq")
		self.trimmed = os.path.join(self.paths.unmapped, self.sample_name + ".trimmed.fastq")
		self.trimmed1 = os.path.join(self.paths.unmapped, self.sample_name + ".1.trimmed.fastq")
		self.trimmed2 = os.path.join(self.paths.unmapped, self.sample_name + ".2.trimmed.fastq")
		self.trimmed1_unpaired = os.path.join(self.paths.unmapped, self.sample_name + ".1_unpaired.trimmed.fastq")
		self.trimmed2_unpaired = os.path.join(self.paths.unmapped, self.sample_name + ".2_unpaired.trimmed.fastq")

		# Mapped: mapped, duplicates marked, removed, reads shifted
		self.paths.mapped = os.path.join(self.paths.sample_root, "mapped")
		self.mapped = os.path.join(self.paths.mapped, self.sample_name + ".trimmed.bowtie2.bam")
		self.filtered = os.path.join(self.paths.mapped, self.sample_name + ".trimmed.bowtie2.filtered.bam")

		# Files in the root of the sample dir
		self.frip = os.path.join(self.paths.sample_root, self.sample_name + "_FRiP.txt")

		# Coverage: read coverage in windows genome-wide
		self.paths.coverage = os.path.join(self.paths.sample_root, "coverage")
		self.coverage = os.path.join(self.paths.coverage, self.sample_name + ".cov")

		self.insertplot = os.path.join(self.paths.sample_root, self.name + "_insertLengths.pdf")
		self.insertdata = os.path.join(self.paths.sample_root, self.name + "_insertLengths.csv")

		self.qc = os.path.join(self.paths.sample_root, self.sample_name + "_qc.tsv")
		self.qc_plot = os.path.join(self.paths.sample_root, self.sample_name + "_qc.pdf")

		# Peaks: peaks called and derivate files
		self.paths.peaks = os.path.join(self.paths.sample_root, "peaks")
		self.peaks = os.path.join(self.paths.peaks, self.sample_name + ("_peaks.narrowPeak" if not self.broad else "_peaks.broadPeak"))


class ChIPmentation(ChIPseqSample):
	"""
	Class to model ChIPmentation samples based on the ChIPseqSample class.

	:param series: Pandas `Series` object.
	:type series: pandas.Series
	"""
	__library__ = "ChIPmentation"

	def __init__(self, series):
		# Use _pd.Series object to have all sample attributes
		if not isinstance(series, pd.Series):
			raise TypeError("Provided object is not a pandas Series.")
		super(ChIPmentation, self).__init__(series)

		self.tagmented = True

	def __repr__(self):
		return "ChIPmentation sample '%s'" % self.sample_name

	def set_file_paths(self):
		super(ChIPmentation, self).set_file_paths()


def main():
	# Parse command-line arguments
	parser = ArgumentParser(
		prog="chipseq-pipeline",
		description="ChIP-seq pipeline."
	)
	parser = arg_parser(parser)
	parser = pypiper.add_pypiper_args(parser, all_args=True)
	args = parser.parse_args()

	# Read in yaml configs
	series = pd.Series(yaml.load(open(args.sample_config, "r")))
	# Create Sample object
	if series["library"] != "ChIPmentation":
		sample = ChIPseqSample(series)
	else:
		sample = ChIPmentation(series)

	# Check if merged
	if len(sample.data_path.split(" ")) > 1:
		sample.merged = True
	else:
		sample.merged = False
	sample.prj = AttributeDict(sample.prj)
	sample.paths = AttributeDict(sample.paths.__dict__)

	# Shorthand for read_type
	if sample.read_type == "paired":
		sample.paired = True
	else:
		sample.paired = False

	# Set file paths
	sample.set_file_paths()
	sample.make_sample_dirs()

	# Start Pypiper object
	# Best practice is to name the pipeline with the name of the script;
	# or put the name in the pipeline interface.
	pipe_manager = pypiper.PipelineManager(name="chipseq", outfolder=sample.paths.sample_root, args=args)

	# Start main function
	process(sample, pipe_manager, args)


def arg_parser(parser):
	"""
	Global options for pipeline.
	"""
	parser.add_argument(
		"-y", "--sample-yaml",
		dest="sample_config",
		help="Yaml config file with sample attributes.",
		type=str
	)
	parser.add_argument(
		"-p", "--peak-caller",
		dest="peak_caller",
		help="Peak caller algorithm.",
		default="macs2",
		type=str
	)
	return parser


def process(sample, pipe_manager, args):
	"""
	This takes unmapped Bam files and makes trimmed, aligned, duplicate marked
	and removed, indexed (and shifted if necessary) Bam files
	along with a UCSC browser track.
	"""
	print("Start processing ChIP-seq sample %s." % sample.name)

	for path in ["sample_root"] + sample.paths.__dict__.keys():
		if not os.path.exists(sample.paths[path]):
			try:
				os.mkdir(sample.paths[path])
			except OSError("Cannot create '%s' path: %s" % (path, sample.paths[path])):
				raise

	# Create NGSTk instance
	tk = NGSTk(pm=pipe_manager)

	# Merge Bam files if more than one technical replicate
	if len(sample.data_path.split(" ")) > 1:
		pipe_manager.timestamp("Merging bam files from replicates")
		cmd = tk.merge_bams(
			input_bams=sample.data_path.split(" "),  # this is a list of sample paths
			merged_bam=sample.unmapped
		)
		pipe_manager.run(cmd, sample.unmapped, shell=True)
		sample.data_path = sample.unmapped

	# Fastqc
	pipe_manager.timestamp("Measuring sample quality with Fastqc")
	cmd = tk.fastqc_rename(
		input_bam=sample.data_path,
		output_dir=sample.paths.sample_root,
		sample_name=sample.sample_name
	)
	pipe_manager.run(cmd, os.path.join(sample.paths.sample_root, sample.sample_name + "_fastqc.zip"), shell=True)

	# Convert bam to fastq
	pipe_manager.timestamp("Converting to Fastq format")
	cmd = tk.bam2fastq(
		inputBam=sample.data_path,
		outputFastq=sample.fastq1 if sample.paired else sample.fastq,
		outputFastq2=sample.fastq2 if sample.paired else None,
		unpairedFastq=sample.fastq_unpaired if sample.paired else None
	)
	pipe_manager.run(cmd, sample.fastq1 if sample.paired else sample.fastq, shell=True)
	if not sample.paired:
		pipe_manager.clean_add(sample.fastq, conditional=True)
	if sample.paired:
		pipe_manager.clean_add(sample.fastq1, conditional=True)
		pipe_manager.clean_add(sample.fastq2, conditional=True)
		pipe_manager.clean_add(sample.fastq_unpaired, conditional=True)

	# Trim reads
	pipe_manager.timestamp("Trimming adapters from sample")
	if pipe_manager.config.parameters.trimmer == "trimmomatic":
		cmd = tk.trimmomatic(
			inputFastq1=sample.fastq1 if sample.paired else sample.fastq,
			inputFastq2=sample.fastq2 if sample.paired else None,
			outputFastq1=sample.trimmed1 if sample.paired else sample.trimmed,
			outputFastq1unpaired=sample.trimmed1_unpaired if sample.paired else None,
			outputFastq2=sample.trimmed2 if sample.paired else None,
			outputFastq2unpaired=sample.trimmed2_unpaired if sample.paired else None,
			cpus=args.cores,
			adapters=pipe_manager.config.resources.adapters,
			log=sample.trimlog
		)
		pipe_manager.run(cmd, sample.trimmed1 if sample.paired else sample.trimmed, shell=True)
		if not sample.paired:
			pipe_manager.clean_add(sample.trimmed, conditional=True)
		else:
			pipe_manager.clean_add(sample.trimmed1, conditional=True)
			pipe_manager.clean_add(sample.trimmed1_unpaired, conditional=True)
			pipe_manager.clean_add(sample.trimmed2, conditional=True)
			pipe_manager.clean_add(sample.trimmed2_unpaired, conditional=True)

	elif pipe_manager.config.parameters.trimmer == "skewer":
		cmd = tk.skewer(
			inputFastq1=sample.fastq1 if sample.paired else sample.fastq,
			inputFastq2=sample.fastq2 if sample.paired else None,
			outputPrefix=os.path.join(sample.paths.unmapped, sample.sample_name),
			outputFastq1=sample.trimmed1 if sample.paired else sample.trimmed,
			outputFastq2=sample.trimmed2 if sample.paired else None,
			trimLog=sample.trimlog,
			cpus=args.cores,
			adapters=pipe_manager.config.resources.adapters
		)
		pipe_manager.run(cmd, sample.trimmed1 if sample.paired else sample.trimmed, shell=True)
		if not sample.paired:
			pipe_manager.clean_add(sample.trimmed, conditional=True)
		else:
			pipe_manager.clean_add(sample.trimmed1, conditional=True)
			pipe_manager.clean_add(sample.trimmed2, conditional=True)

	# Map
	pipe_manager.timestamp("Mapping reads with Bowtie2")
	cmd = tk.bowtie2Map(
		inputFastq1=sample.trimmed1 if sample.paired else sample.trimmed,
		inputFastq2=sample.trimmed2 if sample.paired else None,
		outputBam=sample.mapped,
		log=sample.aln_rates,
		metrics=sample.aln_metrics,
		genomeIndex=getattr(pipe_manager.config.resources.genomes, sample.genome),
		maxInsert=pipe_manager.config.parameters.max_insert,
		cpus=args.cores
	)
	pipe_manager.run(cmd, sample.mapped, shell=True)

	# Filter reads
	pipe_manager.timestamp("Filtering reads for quality")
	cmd = tk.filterReads(
		inputBam=sample.mapped,
		outputBam=sample.filtered,
		metricsFile=sample.dups_metrics,
		paired=sample.paired,
		cpus=args.cores,
		Q=pipe_manager.config.parameters.read_quality
	)
	pipe_manager.run(cmd, sample.filtered, shell=True)

	# Index bams
	pipe_manager.timestamp("Indexing bamfiles with samtools")
	cmd = tk.indexBam(inputBam=sample.mapped)
	pipe_manager.run(cmd, sample.mapped + ".bai", shell=True)
	cmd = tk.indexBam(inputBam=sample.filtered)
	pipe_manager.run(cmd, sample.filtered + ".bai", shell=True)

	track_dir = os.path.dirname(sample.bigwig)
	if not os.path.exists(track_dir):
		os.makedirs(track_dir)

	# Make tracks
	# right now tracks are only made for bams without duplicates
	pipe_manager.timestamp("Making bigWig tracks from bam file")
	cmd = tk.bamToBigWig(
		inputBam=sample.filtered,
		outputBigWig=sample.bigwig,
		genomeSizes=getattr(pipe_manager.config.resources.chromosome_sizes, sample.genome),
		genome=sample.genome,
		tagmented=False,  # by default make extended tracks
		normalize=True
	)
	pipe_manager.run(cmd, sample.bigwig, shell=True)

	# Plot fragment distribution
	if sample.paired and not os.path.exists(sample.insertplot):
		pipe_manager.timestamp("Plotting insert size distribution")
		tk.plot_atacseq_insert_sizes(
			bam=sample.filtered,
			plot=sample.insertplot,
			output_csv=sample.insertdata
		)

	# Count coverage genome-wide
	pipe_manager.timestamp("Calculating genome-wide coverage")
	cmd = tk.genomeWideCoverage(
		inputBam=sample.filtered,
		genomeWindows=getattr(pipe_manager.config.resources.genome_windows, sample.genome),
		output=sample.coverage
	)
	pipe_manager.run(cmd, sample.coverage, shell=True)

	# Calculate NSC, RSC
	pipe_manager.timestamp("Assessing signal/noise in sample")
	cmd = tk.peakTools(
		inputBam=sample.filtered,
		output=sample.qc,
		plot=sample.qc_plot,
		cpus=args.cores
	)
	pipe_manager.run(cmd, sample.qc_plot, shell=True, nofail=True)

	# If sample does not have "ctrl" attribute, finish processing it.
	if not hasattr(sample, "compare_sample"):
		pipe_manager.stop_pipeline()
		print("Finished processing sample %s." % sample.name)
		return

	pipe_manager.wait_for_file(sample.filtered.replace(sample.name, sample.compare_sample))

	if args.peak_caller == "macs2":
		pipe_manager.timestamp("Calling peaks with MACS2")
		# make dir for output (macs fails if it does not exist)
		if not os.path.exists(sample.paths.peaks):
			os.makedirs(sample.paths.peaks)

		# For point-source factors use default settings
		# For broad factors use broad settings
		cmd = tk.macs2CallPeaks(
			treatmentBam=sample.filtered,
			controlBam=sample.filtered.replace(sample.name, sample.compare_sample),
			outputDir=sample.paths.peaks,
			sampleName=sample.name,
			genome=sample.genome,
			broad=True if sample.broad else False
		)
		pipe_manager.run(cmd, sample.peaks, shell=True)

		pipe_manager.timestamp("Ploting MACS2 model")
		cmd = tk.macs2PlotModel(
			sampleName=sample.name,
			outputDir=os.path.join(sample.paths.peaks, sample.name)
		)
		pipe_manager.run(cmd, os.path.join(sample.paths.peaks, sample.name, sample.name + "_model.pdf"), shell=True, nofail=True)
	elif args.peak_caller == "spp":
		pipe_manager.timestamp("Calling peaks with spp")
		# For point-source factors use default settings
		# For broad factors use broad settings
		cmd = tk.sppCallPeaks(
			treatmentBam=sample.filtered,
			controlBam=sample.filtered.replace(sample.name, sample.compare_sample),
			treatmentName=sample.name,
			controlName=sample.compare_sample,
			outputDir=os.path.join(sample.paths.peaks, sample.name),
			broad=True if sample.broad else False,
			cpus=args.cpus
		)
		pipe_manager.run(cmd, sample.peaks, shell=True)

	# Calculate fraction of reads in peaks (FRiP)
	pipe_manager.timestamp("Calculating fraction of reads in peaks (FRiP)")
	cmd = tk.calculateFRiP(
		inputBam=sample.filtered,
		inputBed=sample.peaks,
		output=sample.frip
	)
	pipe_manager.run(cmd, sample.frip, shell=True)

	print("Finished processing sample %s." % sample.name)
	pipe_manager.stop_pipeline()


def get_track_colour(sample, config):
	"""
	Get a colour for a genome browser track based on the IP.
	"""
	import random

	if not hasattr(config, "track_colours"):
		return "0,0,0"
	else:
		if hasattr(sample, "ip"):
			if sample.ip in config["track_colours"].__dict__.keys():
				sample.track_colour = config["track_colours"][sample.ip.upper()]
			else:
				if sample.library in ["ATAC", "ATACSEQ", "ATAC-SEQ"]:
					sample.track_colour = config["track_colours"]["ATAC"]
				elif sample.library in ["DNASE", "DNASESEQ", "DNASE-SEQ"]:
					sample.track_colour = config["track_colours"]["DNASE"]
				else:
					sample.track_colour = random.sample(config["track_colours"].__dict__.values(), 1)[0]  # pick one randomly
		else:
			sample.track_colour = random.sample(config["track_colours"].__dict__.values(), 1)[0]  # pick one randomly


if __name__ == '__main__':
	try:
		sys.exit(main())
	except KeyboardInterrupt:
		print("Program canceled by user!")
		sys.exit(1)
