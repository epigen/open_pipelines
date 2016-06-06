#!/usr/bin/env python

"""
ATAC-seq pipeline
"""

import sys
from argparse import ArgumentParser
import yaml
import pypiper
import os
import pandas as pd

from looper.models import AttributeDict, Sample

from pipelines import toolkit as tk

__author__ = "Andre Rendeiro"
__copyright__ = "Copyright 2015, Andre Rendeiro"
__credits__ = []
__license__ = "GPL2"
__version__ = "0.2"
__maintainer__ = "Andre Rendeiro"
__email__ = "arendeiro@cemm.oeaw.ac.at"
__status__ = "Development"


class ATACseqSample(Sample):
	"""
	Class to model ATAC-seq samples based on the ChIPseqSample class.

	:param series: Pandas `Series` object.
	:type series: pandas.Series
	"""
	__library__ = "ATAC-seq"

	def __init__(self, series):

		# Use pd.Series object to have all sample attributes
		if not isinstance(series, pd.Series):
			raise TypeError("Provided object is not a pandas Series.")
		super(ATACseqSample, self).__init__(series)

		self.tagmented = True

		self.make_sample_dirs()

	def __repr__(self):
		return "ATAC-seq sample '%s'" % self.sample_name

	def set_file_paths(self):
		"""
		Sets the paths of all files for this sample.
		"""
		# Inherit paths from Sample by running Sample's set_file_paths()
		super(ATACseqSample, self).set_file_paths()

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
		self.fastqUnpaired = os.path.join(self.paths.unmapped, self.sample_name + ".unpaired.fastq")
		self.trimmed = os.path.join(self.paths.unmapped, self.sample_name + ".trimmed.fastq")
		self.trimmed1 = os.path.join(self.paths.unmapped, self.sample_name + ".1.trimmed.fastq")
		self.trimmed2 = os.path.join(self.paths.unmapped, self.sample_name + ".2.trimmed.fastq")
		self.trimmed1Unpaired = os.path.join(self.paths.unmapped, self.sample_name + ".1_unpaired.trimmed.fastq")
		self.trimmed2Unpaired = os.path.join(self.paths.unmapped, self.sample_name + ".2_unpaired.trimmed.fastq")

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


class DNaseSample(ATACseqSample):
	"""
	Class to model DNase-seq samples based on the ChIPseqSample class.

	:param series: Pandas `Series` object.
	:type series: pandas.Series
	"""
	__library__ = "DNase-seq"

	def __init__(self, series):

		# Use pd.Series object to have all sample attributes
		if not isinstance(series, pd.Series):
			raise TypeError("Provided object is not a pandas Series.")
		super(DNaseSample, self).__init__(series)

	def __repr__(self):
		return "DNase-seq sample '%s'" % self.sample_name

	def set_file_paths(self):
		super(DNaseSample, self).set_file_paths()


def main():
	# Parse command-line arguments
	parser = ArgumentParser(
		prog="atacseq-pipeline",
		description="ATAC-seq pipeline."
	)
	parser = arg_parser(parser)
	parser = pypiper.add_pypiper_args(parser, all_args=True)
	args = parser.parse_args()

	# Read in yaml configs
	series = pd.Series(yaml.load(open(args.sample_config, "r")))
	# Create Sample object
	if series["library"] != "DNase-seq":
		sample = ATACseqSample(series)
	else:
		sample = DNaseSample(series)
	# Set file paths
	sample.set_file_paths()
	sample.make_sample_dirs()

	pipeline_config = AttributeDict(yaml.load(open(os.path.join(os.path.dirname(__file__), args.config_file), "r")))

	# Start main function
	process(sample, pipeline_config, args)

	# # Remove sample config
	# if not args.dry_run:
	# 	os.system("rm %s" % args.sample_config)


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


def process(sample, pipeline_config, args):
	"""
	This takes unmapped Bam files and makes trimmed, aligned, duplicate marked
	and removed, indexed, shifted Bam files along with a UCSC browser track.
	Peaks are called and filtered.
	"""

	print("Start processing ATAC-seq sample %s." % sample.sample_name)

	for path in ["sample_root"] + sample.paths.__dict__.keys():
		if not os.path.exists(sample.paths[path]):
			try:
				os.mkdir(sample.paths[path])
			except OSError("Cannot create '%s' path: %s" % (path, sample.paths[path])):
				raise

	# Start Pypiper object
	pipe = pypiper.PipelineManager("pipe", sample.paths.sample_root, args=args)

	# Merge Bam files if more than one technical replicate
	if len(sample.data_path.split(" ")) > 1:
		pipe.timestamp("Merging bam files from replicates")
		cmd = tk.mergeBams(
			inputBams=sample.data_path.split(" "),  # this is a list of sample paths
			outputBam=sample.unmapped
		)
		pipe.run(cmd, sample.unmapped, shell=True)
		sample.data_path = sample.unmapped

	# Fastqc
	pipe.timestamp("Measuring sample quality with Fastqc")
	cmd = tk.fastqc(
		inputBam=sample.data_path,
		outputDir=sample.paths.sample_root,
		sampleName=sample.sample_name
	)
	pipe.run(cmd, os.path.join(sample.paths.sample_root, sample.sample_name + "_fastqc.zip"), shell=True)

	# Convert bam to fastq
	pipe.timestamp("Converting to Fastq format")
	cmd = tk.bam2fastq(
		inputBam=sample.data_path,
		outputFastq=sample.fastq1 if sample.paired else sample.fastq,
		outputFastq2=sample.fastq2 if sample.paired else None,
		unpairedFastq=sample.fastqUnpaired if sample.paired else None
	)
	pipe.run(cmd, sample.fastq1 if sample.paired else sample.fastq, shell=True)
	if not sample.paired:
		pipe.clean_add(sample.fastq, conditional=True)
	if sample.paired:
		pipe.clean_add(sample.fastq1, conditional=True)
		pipe.clean_add(sample.fastq2, conditional=True)
		pipe.clean_add(sample.fastqUnpaired, conditional=True)

	# Trim reads
	pipe.timestamp("Trimming adapters from sample")
	if pipeline_config.parameters.trimmer == "trimmomatic":
		cmd = tk.trimmomatic(
			inputFastq1=sample.fastq1 if sample.paired else sample.fastq,
			inputFastq2=sample.fastq2 if sample.paired else None,
			outputFastq1=sample.trimmed1 if sample.paired else sample.trimmed,
			outputFastq1unpaired=sample.trimmed1Unpaired if sample.paired else None,
			outputFastq2=sample.trimmed2 if sample.paired else None,
			outputFastq2unpaired=sample.trimmed2Unpaired if sample.paired else None,
			cpus=args.cores,
			adapters=pipeline_config.resources.adapters,
			log=sample.trimlog
		)
		pipe.run(cmd, sample.trimmed1 if sample.paired else sample.trimmed, shell=True)
		if not sample.paired:
			pipe.clean_add(sample.trimmed, conditional=True)
		else:
			pipe.clean_add(sample.trimmed1, conditional=True)
			pipe.clean_add(sample.trimmed1Unpaired, conditional=True)
			pipe.clean_add(sample.trimmed2, conditional=True)
			pipe.clean_add(sample.trimmed2Unpaired, conditional=True)

	elif pipeline_config.parameters.trimmer == "skewer":
		cmd = tk.skewer(
			inputFastq1=sample.fastq1 if sample.paired else sample.fastq,
			inputFastq2=sample.fastq2 if sample.paired else None,
			outputPrefix=os.path.join(sample.paths.unmapped, sample.sample_name),
			outputFastq1=sample.trimmed1 if sample.paired else sample.trimmed,
			outputFastq2=sample.trimmed2 if sample.paired else None,
			trimLog=sample.trimlog,
			cpus=args.cores,
			adapters=pipeline_config.resources.adapters
		)
		pipe.run(cmd, sample.trimmed1 if sample.paired else sample.trimmed, shell=True)
		if not sample.paired:
			pipe.clean_add(sample.trimmed, conditional=True)
		else:
			pipe.clean_add(sample.trimmed1, conditional=True)
			pipe.clean_add(sample.trimmed2, conditional=True)

	# Map
	pipe.timestamp("Mapping reads with Bowtie2")
	cmd = tk.bowtie2Map(
		inputFastq1=sample.trimmed1 if sample.paired else sample.trimmed,
		inputFastq2=sample.trimmed2 if sample.paired else None,
		outputBam=sample.mapped,
		log=sample.aln_rates,
		metrics=sample.aln_metrics,
		genomeIndex=getattr(pipeline_config.resources.genomes, sample.genome),
		maxInsert=pipeline_config.parameters.max_insert,
		cpus=args.cores
	)
	pipe.run(cmd, sample.mapped, shell=True)

	# Filter reads
	pipe.timestamp("Filtering reads for quality")
	cmd = tk.filterReads(
		inputBam=sample.mapped,
		outputBam=sample.filtered,
		metricsFile=sample.dups_metrics,
		paired=sample.paired,
		cpus=args.cores,
		Q=pipeline_config.parameters.read_quality
	)
	pipe.run(cmd, sample.filtered, shell=True)

	# Shift reads
	if sample.tagmented:
		pipe.timestamp("Shifting reads of tagmented sample")
		cmd = tk.shiftReads(
			inputBam=sample.filtered,
			genome=sample.genome,
			outputBam=sample.filteredshifted
		)
		pipe.run(cmd, sample.filteredshifted, shell=True)

	# Index bams
	pipe.timestamp("Indexing bamfiles with samtools")
	cmd = tk.indexBam(inputBam=sample.mapped)
	pipe.run(cmd, sample.mapped + ".bai", shell=True)
	cmd = tk.indexBam(inputBam=sample.filtered)
	pipe.run(cmd, sample.filtered + ".bai", shell=True)
	if sample.tagmented:
		cmd = tk.indexBam(inputBam=sample.filteredshifted)
		pipe.run(cmd, sample.filteredshifted + ".bai", shell=True)

	# Make tracks
	# right now tracks are only made for bams without duplicates
	pipe.timestamp("Making bigWig tracks from bam file")
	cmd = tk.bamToBigWig(
		inputBam=sample.filtered,
		outputBigWig=sample.bigwig,
		genomeSizes=getattr(pipeline_config.resources.chromosome_sizes, sample.genome),
		genome=sample.genome,
		tagmented=False,  # by default make extended tracks
		normalize=True
	)
	pipe.run(cmd, sample.bigwig, shell=True)
	cmd = tk.addTrackToHub(
		sampleName=sample.sample_name,
		trackURL=sample.track_url,
		trackHub=os.path.join(os.path.dirname(sample.bigwig), "trackHub_{0}.txt".format(sample.genome)),
		colour=get_track_colour(sample, pipeline_config)
	)
	pipe.run(cmd, lock_name=sample.sample_name + "addToTrackHub", shell=True)
	# tk.linkToTrackHub(
	# 	trackHubURL="/".join([prj.config["url"], prj.name, "trackHub_{0}.txt".format(sample.genome)]),
	# 	fileName=os.path.join(prj.dirs.root, "ucsc_tracks_{0}.html".format(sample.genome)),
	# 	genome=sample.genome
	# )

	# Plot fragment distribution
	if sample.paired and not os.path.exists(sample.insertplot):
		pipe.timestamp("Plotting insert size distribution")
		tk.plotInsertSizesFit(
			bam=sample.filtered,
			plot=sample.insertplot,
			outputCSV=sample.insertdata
		)

	# Count coverage genome-wide
	pipe.timestamp("Calculating genome-wide coverage")
	cmd = tk.genomeWideCoverage(
		inputBam=sample.filtered,
		genomeWindows=getattr(pipeline_config.resources.genome_windows, sample.genome),
		output=sample.coverage
	)
	pipe.run(cmd, sample.coverage, shell=True)

	# Calculate NSC, RSC
	pipe.timestamp("Assessing signal/noise in sample")
	cmd = tk.peakTools(
		inputBam=sample.filtered,
		output=sample.qc,
		plot=sample.qc_plot,
		cpus=args.cores
	)
	pipe.run(cmd, sample.qc_plot, shell=True, nofail=True)

	# Call peaks
	pipe.timestamp("Calling peaks with MACS2")
	# make dir for output (macs fails if it does not exist)
	if not os.path.exists(sample.paths.peaks):
		os.makedirs(sample.paths.peaks)

	cmd = tk.macs2CallPeaksATACSeq(
		treatmentBam=sample.filtered,
		outputDir=sample.paths.peaks,
		sampleName=sample.sample_name,
		genome=sample.genome
	)
	pipe.run(cmd, sample.peaks, shell=True)

	# # Filter peaks based on mappability regions
	# pipe.timestamp("Filtering peaks in low mappability regions")

	# # get closest read length of sample to available mappability read lengths
	# closestLength = min(pipeline_config.resources["alignability"][sample.genome].keys(), key=lambda x:abs(x - sample.readLength))

	# cmd = tk.filterPeaksMappability(
	#     peaks=sample.peaks,
	#     alignability=pipeline_config.resources["alignability"][sample.genome][closestLength],
	#     filteredPeaks=sample.filteredPeaks
	# )
	# pipe.run(cmd, sample.filteredPeaks, shell=True)

	# Calculate fraction of reads in peaks (FRiP)
	pipe.timestamp("Calculating fraction of reads in peaks (FRiP)")
	cmd = tk.calculateFRiP(
		inputBam=sample.filtered,
		inputBed=sample.peaks,
		output=sample.frip
	)
	pipe.run(cmd, sample.frip, shell=True)

	pipe.stop_pipeline()
	print("Finished processing sample %s." % sample.sample_name)


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
				sample.track_colour = config["track_colours"][sample.ip]
			else:
				if sample.library in ["ATAC", "ATACSEQ", "ATAC-SEQ"]:
					sample.track_colour = config["track_colours"]["ATAC"]
				elif sample.library in ["DNASE", "DNASESEQ", "DNASE-SEQ"]:
					sample.track_colour = config["track_colours"]["DNASE"]
				else:
					sample.track_colour = random.sample(config["track_colours"], 1)[0]  # pick one randomly
		else:
			sample.track_colour = random.sample(config["track_colours"], 1)[0]  # pick one randomly


if __name__ == '__main__':
	try:
		sys.exit(main())
	except KeyboardInterrupt:
		print("Program canceled by user!")
		sys.exit(1)
