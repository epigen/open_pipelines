#!/usr/bin/env python

"""
STARR-seq pipeline
"""

import os
import sys
from argparse import ArgumentParser
import yaml
import pypiper
from pypiper.ngstk import NGSTk
from pep import AttributeDict, Sample


__author__ = "Andre Rendeiro"
__copyright__ = "Copyright 2015, Andre Rendeiro"
__credits__ = []
__license__ = "GPL2"
__version__ = "0.1"
__maintainer__ = "Andre Rendeiro"
__email__ = "arendeiro@cemm.oeaw.ac.at"
__status__ = "Development"



class STARRseqSample(Sample):
	"""
	Class to model STARR-seq samples based on the ChIPseqSample class.

	:param series: Collection of sample attributes.
	:type series: Mapping | pandas.core.series.Series
	"""
	__library__ = "STARR-seq"

	def __init__(self, series):
		super(STARRseqSample, self).__init__(series)
		self.tagmented = False
		self.make_sample_dirs()

	def __repr__(self):
		return "STARR-seq sample '%s'" % self.sample_name

	def set_file_paths(self, project=None):
		"""
		Sets the paths of all files for this sample.
		"""
		# Inherit paths from Sample by running Sample's set_file_paths()
		super(STARRseqSample, self).set_file_paths(project)

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
		# this will create additional bam files with reads shifted
		self.filteredshifted = os.path.join(self.paths.mapped, self.name + ".trimmed.bowtie2.filtered.shifted.bam")

		# Files in the root of the sample dir
		self.frip = os.path.join(self.paths.sample_root, self.name + "_FRiP.txt")

		# Mapped: mapped, duplicates marked, removed, reads shifted
		# this will create additional bam files with reads shifted
		self.filteredshifted = os.path.join(self.paths.mapped, self.name + ".trimmed.bowtie2.filtered.shifted.bam")

		# Coverage: read coverage in windows genome-wide
		self.paths.coverage = os.path.join(self.paths.sample_root, "coverage")
		self.coverage = os.path.join(self.paths.coverage, self.name + ".cov")

		self.insertplot = os.path.join(self.paths.sample_root, self.name + "_insertLengths.pdf")
		self.insertdata = os.path.join(self.paths.sample_root, self.name + "_insertLengths.csv")
		self.qc = os.path.join(self.paths.sample_root, self.name + "_qc.tsv")
		self.qc_plot = os.path.join(self.paths.sample_root, self.name + "_qc.pdf")

		# Peaks: peaks called and derivate files
		self.paths.peaks = os.path.join(self.paths.sample_root, "peaks")
		self.peaks = os.path.join(self.paths.peaks, self.name + "_peaks.narrowPeak")
		self.filteredPeaks = os.path.join(self.paths.peaks, self.name + "_peaks.filtered.bed")


def main():
	# Parse command-line arguments
	parser = ArgumentParser(
		prog="starrseq-pipeline",
		description="STARR-seq pipeline."
	)
	parser = arg_parser(parser)
	parser = pypiper.add_pypiper_args(parser, all_args=True)
	args = parser.parse_args()
	if args.sample_config is None:
		parser.print_help()
		return 1

	# Read in yaml config and create Sample object
	sample = STARRseqSample(yaml.load(open(args.sample_config, "r")))

	# Check if merged
	if len(sample.data_source.split(" ")) > 1:
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
	pipe_manager = pypiper.PipelineManager(name="starrseq", outfolder=sample.paths.sample_root, args=args)

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
	return parser


def process(sample, pipe_manager, args):
	"""
	This takes unmapped Bam files and makes trimmed, aligned, duplicate marked
	and removed, indexed Bam files along with a UCSC browser track.
	Peaks are called and filtered.
	"""

	print("Start processing STARR-seq sample %s." % sample.sample_name)

	for path in ["sample_root"] + sample.paths.__dict__.keys():
		if not os.path.exists(sample.paths[path]):
			try:
				os.mkdir(sample.paths[path])
			except OSError("Cannot create '%s' path: %s" % (path, sample.paths[path])):
				raise

	# Create NGSTk instance
	tk = NGSTk(pm=pipe_manager)

	# Merge Bam files if more than one technical replicate
	if len(sample.data_source.split(" ")) > 1:
		pipe_manager.timestamp("Merging bam files from replicates")
		cmd = tk.mergeBams(
			input_bams=sample.data_source.split(" "),  # this is a list of sample paths
			output_bam=sample.unmapped
		)
		pipe_manager.run(cmd, sample.unmapped, shell=True)
		sample.data_source = sample.unmapped

	# Fastqc
	pipe_manager.timestamp("Measuring sample quality with Fastqc")
	cmd = tk.fastqc(sample.data_source, output_dir=sample.paths.sample_root)
	pipe_manager.run(cmd, os.path.join(sample.paths.sample_root, sample.sample_name + "_fastqc.zip"), shell=True)

	# Convert bam to fastq
	pipe_manager.timestamp("Converting to Fastq format")
	cmd = tk.bam2fastq(
		input_bam=sample.data_source,
		output_fastq=sample.fastq1 if sample.paired else sample.fastq,
		output_fastq2=sample.fastq2 if sample.paired else None,
		unpaired_fastq=sample.fastq_unpaired if sample.paired else None
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
	if pipe_manager.parameters.trimmer == "trimmomatic":
		cmd = tk.trimmomatic(
			input_fastq1=sample.fastq1 if sample.paired else sample.fastq,
			input_fastq2=sample.fastq2 if sample.paired else None,
			output_fastq1=sample.trimmed1 if sample.paired else sample.trimmed,
			output_fastq1_unpaired=sample.trimmed1_unpaired if sample.paired else None,
			output_fastq2=sample.trimmed2 if sample.paired else None,
			output_fastq2_unpaired=sample.trimmed2_unpaired if sample.paired else None,
			cpus=args.cores,
			adapters=pipe_manager.resources.adapters,
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

	elif pipe_manager.parameters.trimmer == "skewer":
		cmd = tk.skewer(
			input_fastq1=sample.fastq1 if sample.paired else sample.fastq,
			input_fastq2=sample.fastq2 if sample.paired else None,
			output_prefix=os.path.join(sample.paths.unmapped, sample.sample_name),
			output_fastq1=sample.trimmed1 if sample.paired else sample.trimmed,
			output_fastq2=sample.trimmed2 if sample.paired else None,
			log=sample.trimlog,
			cpus=args.cores,
			adapters=pipe_manager.resources.adapters
		)
		pipe_manager.run(cmd, sample.trimmed1 if sample.paired else sample.trimmed, shell=True)
		if not sample.paired:
			pipe_manager.clean_add(sample.trimmed, conditional=True)
		else:
			pipe_manager.clean_add(sample.trimmed1, conditional=True)
			pipe_manager.clean_add(sample.trimmed2, conditional=True)

	# Map
	pipe_manager.timestamp("Mapping reads with Bowtie2")
	cmd = tk.bowtie2_map(
		input_fastq1=sample.trimmed1 if sample.paired else sample.trimmed,
		input_fastq2=sample.trimmed2 if sample.paired else None,
		output_bam=sample.mapped,
		log=sample.aln_rates,
		metrics=sample.aln_metrics,
		genome_index=getattr(pipe_manager.resources.genomes, sample.genome),
		max_insert=pipe_manager.parameters.max_insert,
		cpus=args.cores
	)
	pipe_manager.run(cmd, sample.mapped, shell=True)

	# Filter reads
	pipe_manager.timestamp("Filtering reads for quality")
	cmd = tk.filter_reads(
		input_bam=sample.mapped,
		output_bam=sample.filtered,
		metrics_file=sample.dups_metrics,
		paired=sample.paired,
		cpus=args.cores,
		Q=pipe_manager.parameters.read_quality
	)
	pipe_manager.run(cmd, sample.filtered, shell=True)

	# Index bams
	pipe_manager.timestamp("Indexing bamfiles with samtools")
	cmd = tk.index_bam(input_bam=sample.mapped)
	pipe_manager.run(cmd, sample.mapped + ".bai", shell=True)
	cmd = tk.index_bam(input_bam=sample.filtered)
	pipe_manager.run(cmd, sample.filtered + ".bai", shell=True)

	# Make tracks
	# right now tracks are only made for bams without duplicates
	pipe_manager.timestamp("Making bigWig tracks from bam file")
	cmd = tk.bam_to_bigwig(
		input_bam=sample.filtered,
		output_bigwig=sample.bigwig,
		genome_sizes=getattr(pipe_manager.resources.chromosome_sizes, sample.genome),
		genome=sample.genome,
		tagmented=False,  # by default make extended tracks
		normalize=True
	)
	pipe_manager.run(cmd, sample.bigwig, shell=True)

	# Plot fragment distribution
	if sample.paired and not os.path.exists(sample.insertplot):
		pipe_manager.timestamp("Plotting insert size distribution")
		tk.plotInsertSizesFit(
			bam=sample.filtered,
			plot=sample.insertplot,
			outputCSV=sample.insertdata
		)

	# Count coverage genome-wide
	pipe_manager.timestamp("Calculating genome-wide coverage")
	cmd = tk.genome_wide_coverage(
		input_bam=sample.filtered,
		genome_windows=getattr(pipe_manager.resources.genome_windows, sample.genome),
		output=sample.coverage
	)
	pipe_manager.run(cmd, sample.coverage, shell=True)

	# Calculate NSC, RSC
	pipe_manager.timestamp("Assessing signal/noise in sample")
	cmd = tk.run_spp(
		input_bam=sample.filtered,
		output=sample.qc,
		plot=sample.qc_plot,
		cpus=args.cores
	)
	pipe_manager.run(cmd, sample.qc_plot, shell=True, nofail=True)

	# Call peaks
	pipe_manager.timestamp("Calling peaks with MACS2")
	# make dir for output (macs fails if it does not exist)
	if not os.path.exists(sample.paths.peaks):
		os.makedirs(sample.paths.peaks)

	cmd = tk.macs2_call_peaks_atacseq(
		treatment_bam=sample.filtered,
		output_dir=sample.paths.peaks,
		sample_name=sample.sample_name,
		genome=sample.genome
	)
	pipe_manager.run(cmd, sample.peaks, shell=True)

	# Calculate fraction of reads in peaks (FRiP)
	pipe_manager.timestamp("Calculating fraction of reads in peaks (FRiP)")
	cmd = tk.calculate_frip(
		input_bam=sample.filtered,
		input_bed=sample.peaks,
		output=sample.frip
	)
	pipe_manager.run(cmd, sample.frip, shell=True)

	print("Finished processing sample %s." % sample.sample_name)
	pipe_manager.stop_pipeline()



if __name__ == '__main__':
	try:
		sys.exit(main())
	except KeyboardInterrupt:
		print("Program canceled by user!")
		sys.exit(1)
