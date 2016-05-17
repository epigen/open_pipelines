#!/usr/bin/env python

"""
QUANT-seq pipeline
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


class QuantseqSample(Sample):
	"""
	Class to model Quant-seq samples based on the generic Sample class (itself a pandas.Series).

	:param series: Pandas `Series` object.
	:type series: pandas.Series

	:Example:

	from pipelines import Project, SampleSheet, QuantseqSample
	prj = Project("ngs")
	sheet = SampleSheet("/projects/example/sheet.csv", prj)
	s1 = QuantseqSample(sheet.ix[0])
	"""
	def __init__(self, series):

		# Passed series must either be a pd.Series or a daugther class
		if not isinstance(series, pd.Series):
			raise TypeError("Provided object is not a pandas Series.")
		super(QuantseqSample, self).__init__(series)

	def __repr__(self):
		return "Quant-seq sample '%s'" % self.sample_name

	def set_file_paths(self):
		"""
		Sets the paths of all files for this sample.
		"""
		# Inherit paths from Sample by running Sample's set_file_paths()
		super(QuantseqSample, self).set_file_paths()

		# Files in the root of the sample dir
		self.fastqc = os.path.join(self.paths.sample_root, self.sample_name + ".fastqc.zip")
		self.trimlog = os.path.join(self.paths.sample_root, self.sample_name + ".trimlog.txt")

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

		# kallisto pseudoalignments
		self.pseudomapped = os.path.join(self.paths.mapped, self.name + ".pseudoalignment.bam")

		# RNA quantification
		self.paths.quant = os.path.join(self.paths.sample_root, "quantification")
		self.kallistoQuant = os.path.join(self.paths.quant, "abundance.tsv")


def main():
	# Parse command-line arguments
	parser = ArgumentParser(
		prog="quantseq-pipeline",
		description="QUANT-seq pipeline."
	)
	parser = arg_parser(parser)
	parser = pypiper.add_pypiper_args(parser, all_args=True)
	args = parser.parse_args()

	# Read in yaml config and create Sample object
	sample = QuantseqSample(pd.Series(yaml.load(open(args.sample_config, "r"))))
	# Set file paths
	sample.set_file_paths()

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
		transcriptomeIndex=pipeline_config["resources"]["genome_index"][sample.transcriptome],
		cpus=args.cores
	)
	pipe.run(cmd, sample.kallistoQuant, shell=True, nofail=True)

	pipe.stop_pipeline()
	print("Finished processing sample %s." % sample.sample_name)


if __name__ == '__main__':
	try:
		sys.exit(main())
	except KeyboardInterrupt:
		print("Program canceled by user!")
		sys.exit(1)
