#!/usr/bin/env python

"""
QUANT-seq pipeline
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
__version__ = "0.2"
__maintainer__ = "Andre Rendeiro"
__email__ = "arendeiro@cemm.oeaw.ac.at"
__status__ = "Development"


def main():
	# Parse command-line arguments
	parser = ArgumentParser(
		prog="quantseq-pipeline",
		description="QUANT-seq pipeline."
	)
	parser = arg_parser(parser)
	parser = pypiper.add_pypiper_args(parser, all_args=True)
	args = parser.parse_args()

	# Read in yaml configs
	sample = AttributeDict(yaml.load(open(args.sample_config, "r")))
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

	print("Start processing QUANT-seq sample %s." % sample.sample_name)

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

	# With kallisto from unmapped reads
	pipe.timestamp("Quantifying read counts with kallisto")
	cmd = tk.kallisto(
		inputFastq=sample.trimmed1 if sample.paired else sample.trimmed,
		inputFastq2=sample.trimmed1 if sample.paired else None,
		outputDir=sample.paths.quant,
		outputBam=sample.pseudomapped,
		transcriptomeIndex=pipeline_config["resources"]["genome_index"][sample.genome],
		cpus=args.cpus
	)
	pipe.call_lock(cmd, sample.kallistoQuant, shell=True, nofail=True)

	pipe.stop_pipeline()
	print("Finished processing sample %s." % sample.sample_name)


if __name__ == '__main__':
	try:
		sys.exit(main())
	except KeyboardInterrupt:
		print("Program canceled by user!")
		sys.exit(1)
