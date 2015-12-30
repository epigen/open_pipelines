#!/usr/bin/env python

"""
Example pipeline
"""

import sys
import pypiper
from pypiper.ngstk import NGSTk
from argparse import ArgumentParser
import os
import yaml


__author__ = "Nathan Sheffield"
__credits__ = []
__license__ = "GPL3"
__version__ = "0.1"
__maintainer__ = "Nathan Sheffield"
__email__ = "nathan@code.databio.org"
__status__ = "Development"


def main():
	# Parse command-line arguments
	parser = ArgumentParser(
		prog="example-pipeline",
		description="Example pipeline."
	)

	# First, add arguments from Pypiper
	# this adds options including:
	# -R: Recover mode to overwrite locks
	# -D: Dirty mode to make suppress cleaning intermediate files
	parser = pypiper.add_pypiper_args(parser, all_args=True)
	# Add remaining pipeline-specific arguments
	parser = add_parser_arguments(parser)

	# Parse provided command-line arguments
	args = parser.parse_args()

	# Start sample processing function
	process(args)


def add_parser_arguments(parser):
	"""
	Add Pipeline options to argument parser.
	"""
	# Default config
	default_config = os.path.splitext(os.path.basename(__file__))[0] + ".yaml"

	parser.add_argument(
		"-c", "--config", dest="config_file", type=str,
		help="pipeline config file in YAML format; relative paths are considered \
		relative to the pipeline script. defaults to " + default_config,
		required=False, default=default_config, metavar="CONFIG_FILE")

	parser.add_argument(
		"-i", "--input", dest="input", type=str, nargs="+",
		help="one or more input files (required)",
		required=True, metavar="INPUT_FILES")
	# input was previously called unmapped_bam

	parser.add_argument(
		"-o", "--output_parent", dest="output_parent", type=str,
		help="parent output directory of the project (required). The sample_name \
		argument will be appended to this folder for output",
		required=True, metavar="PARENT_OUTPUT_FOLDER")
	# output_parent was previously called project_root

	parser.add_argument(
		"-s", "--sample_name", dest="sample_name", type=str,
		help="unique name for output subfolder and files (required)",
		required=True, metavar="SAMPLE_NAME")

	parser.add_argument(
		"-p", "--cores", dest="cores", type=str,
		help="number of cores to use for parallel processes",
		required=False, default=1, metavar="NUMBER_OF_CORES")

	return parser


def process(args):
	"""
	This pipeline does this and that... (describe)
	"""

	# Create PipelineManager object
	pipeline_outfolder = os.path.abspath(os.path.join(args.output_parent, args.sample_name))
	pipe = pypiper.PipelineManager(name="ExamplePipeline", outfolder=pipeline_outfolder, args=args)

	# Read YAML config file
	if not os.path.isabs(args.config_file):
		# Set the path to an absolute path, relative to pipeline script
		default_config_abs = os.path.join(os.path.dirname(__file__), args.config_file)
		if os.path.isfile(default_config_abs):
			args.config_file = default_config_abs
		else:
			args.config_file = None

		with open(args.config_file, 'r') as config_file:
			config = yaml.load(config_file)

	# Start pipeline
	pipe.timestamp("### Starting to process sample %s." % args.sample_name)

	# Add pipeline commands here

	# Example command
	# "fastqc --noextract -o <output_dir> <input>.bam" which produces "<input>.zip"
	# Build command - 2 options:
	# 1) build a command string yourself:
	cmd = " ".join([config["tools"]["fastqc"], "--noextract", "-o", args.output_parent, args.input])

	# 2) use pre-built commands from pypiper.ngstk:
	tk = NGSTk(config["tools"])
	cmd = tk.fastqc(args.input, args.output_parent)

	# build path to the file produced by the command
	output = os.path.join(args.output_parent, args.sample_name) + "_fastqc.zip"

	# run that command with pypiper, specifying the expected output file
	pipe.run(cmd, output)

	# Terminate
	pipe.timestamp("### Finished processing sample %s." % args.sample_name)
	pipe.stop_pipeline()

if __name__ == '__main__':
	try:
		sys.exit(main())
	except KeyboardInterrupt:
		print("Program canceled by user!")
		sys.exit(1)
