#!/usr/bin/env python

"""
QUANT-seq pipeline
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
__version__ = "0.2"
__maintainer__ = "Andre Rendeiro"
__email__ = "arendeiro@cemm.oeaw.ac.at"
__status__ = "Development"



class QuantseqSample(Sample):
	"""
	Class to model Quant-seq samples based on the generic Sample class.

	:Example:

	from pipelines import Project, SampleSheet, QuantseqSample
	prj = Project("ngs")
	sheet = SampleSheet("/projects/example/sheet.csv", prj)
	s1 = QuantseqSample(sheet.ix[0])
	"""
	__library__ = "Quant-seq"

	def __repr__(self):
		return "Quant-seq sample '%s'" % self.sample_name

	def set_file_paths(self, project=None):
		"""
		Sets the paths of all files for this sample.
		"""
		# Inherit paths from Sample by running Sample's set_file_paths()
		super(QuantseqSample, self).set_file_paths(project)

		# Files in the root of the sample dir
		self.fastqc = os.path.join(self.paths.sample_root, self.sample_name + ".fastqc.zip")
		self.trimlog = os.path.join(self.paths.sample_root, self.sample_name + ".trimlog.txt")

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

		# kallisto pseudoalignments
		self.paths.mapped = os.path.join(self.paths.sample_root, "mapped")
		self.pseudomapped = os.path.join(self.paths.mapped, self.name + ".pseudoalignment.bam")

		# RNA quantification
		self.paths.quant = os.path.join(self.paths.sample_root, "quantification")
		self.kallisto_quant = os.path.join(self.paths.quant, "abundance.tsv")


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
	sample = QuantseqSample(yaml.load(open(args.sample_config, "r")))

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
	return parser


def process(sample, pipe_manager, args):
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
	cmd = tk.fastqc(sample.data_source, sample.paths.sample_root)
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

	# With kallisto from unmapped reads
	pipe_manager.timestamp("Quantifying read counts with kallisto")
	cmd = tk.kallisto(
		input_fastq=sample.trimmed1 if sample.paired else sample.trimmed,
		input_fastq2=sample.trimmed1 if sample.paired else None,
		output_dir=sample.paths.quant,
		output_bam=sample.pseudomapped,
		transcriptome_index=pipe_manager.resources.genome_index[sample.transcriptome],
		cpus=args.cores
	)
	pipe_manager.run(cmd, sample.kallisto_quant, shell=True, nofail=True)

	print("Finished processing sample %s." % sample.sample_name)
	pipe_manager.stop_pipeline()


if __name__ == '__main__':
	try:
		sys.exit(main())
	except KeyboardInterrupt:
		print("Program canceled by user!")
		sys.exit(1)
