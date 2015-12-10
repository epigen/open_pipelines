#!/usr/bin/env python

"""
STARR-seq pipeline
"""

import sys
from argparse import ArgumentParser
import yaml
import pypiper
import os

try:
	from pipelines.models import AttributeDict
	from pipelines import toolkit as tk
except:
	sys.path.append(os.path.join(os.path.dirname(__file__), "pipelines"))
	from models import AttributeDict
	import toolkit as tk


__author__ = "Andre Rendeiro"
__copyright__ = "Copyright 2015, Andre Rendeiro"
__credits__ = []
__license__ = "GPL2"
__version__ = "0.1"
__maintainer__ = "Andre Rendeiro"
__email__ = "arendeiro@cemm.oeaw.ac.at"
__status__ = "Development"


def main():
	# Parse command-line arguments
	parser = ArgumentParser(
		prog="starrseq-pipeline",
		description="STARR-seq pipeline."
	)
	parser = pypiper.add_pypiper_args(parser, all_args=True)
	parser = arg_parser(parser)
	args = parser.parse_args()

	# Read in yaml configs
	sample = AttributeDict(**yaml.load(open(args.sample_config, "r")))
	pipeline_config = AttributeDict(**yaml.load(open(os.path.join(os.path.dirname(__file__), "starrseq.yaml"), "r")))

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
	and removed, indexed Bam files along with a UCSC browser track.
	Peaks are called and filtered.
	"""

	print("Start processing STARR-seq sample %s." % sample.sample_name)

	for path in sample.paths.__dict__.keys():
		if not os.path.exists(path):
			try:
				os.mkdir(path)
			except OSError("Cannot create path: %s" % path):
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

	# Index bams
	pipe.timestamp("Indexing bamfiles with samtools")
	cmd = tk.indexBam(inputBam=sample.mapped)
	pipe.run(cmd, sample.mapped + ".bai", shell=True)
	cmd = tk.indexBam(inputBam=sample.filtered)
	pipe.run(cmd, sample.filtered + ".bai", shell=True)

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
		inputBam=sample.filteredshifted,
		genomeWindows=getattr(pipeline_config.resources.genome_windows, sample.genome),
		output=sample.coverage
	)
	pipe.run(cmd, sample.coverage, shell=True)

	# Calculate NSC, RSC
	pipe.timestamp("Assessing signal/noise in sample")
	cmd = tk.peakTools(
		inputBam=sample.filteredshifted,
		output=sample.qc,
		plot=sample.qcPlot,
		cpus=args.cores
	)
	pipe.run(cmd, sample.qcPlot, shell=True, nofail=True)

	# Call peaks
	pipe.timestamp("Calling peaks with MACS2")
	# make dir for output (macs fails if it does not exist)
	if not os.path.exists(sample.paths.peaks):
		os.makedirs(sample.paths.peaks)

	cmd = tk.macs2CallPeaksSTARRSeq(
		treatmentBam=sample.filteredshifted,
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
		inputBam=sample.filteredshifted,
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
		sample.track_colour = random.sample(config["colour_gradient"], 1)[0]  # pick one randomly


if __name__ == '__main__':
	try:
		sys.exit(main())
	except KeyboardInterrupt:
		print("Program canceled by user!")
		sys.exit(1)
