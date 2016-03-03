#!/usr/bin/env python

"""
Hi-C pipeline
"""

__author__ = "Charles Dietz"
__email__ = "cdietz@cemm.oeaw.ac.at"
__credits__ = [""]
__license__ = "GPL3"
__version__ = "0.1"
__status__ = "Development"

from argparse import ArgumentParser
import os, re
import sys
import subprocess
import yaml
import pypiper

parser = ArgumentParser(description='Pipeline')

# First, add arguments from Pypiper, including
# 1. pypiper options, 2. looper connections, 3. common options,
# using the all_args=True flag (you can customize this).
# Adds options including; for details, see:
# http://github.com/epigen/pypiper/command_line_args.md
parser = pypiper.add_pypiper_args(parser, all_args=True)

# Add any pipeline-specific arguments
#parser.add_argument('-t', '--trimgalore', dest='trimmomatic', action="store_false", default=True,
#	help='Use trimgalore instead of trimmomatic?')

args = parser.parse_args()

if args.single_or_paired == "paired":
	args.paired_end = True
else:
	args.paired_end = False

# Merging
################################################################################
# If 2 input files are given, then these are to be merged.
# Must be done here to initialize the sample name correctly

merge = False
if len(args.input) > 1:
	merge = True
	if args.sample_name == "default":
		args.sample_name = "merged"
else:
	if args.sample_name == "default":
		# Default sample name is derived from the input file
		args.sample_name = os.path.splitext(os.path.basename(args.input[0]))[0]

# Create a PipelineManager object and start the pipeline
pm = pypiper.PipelineManager(name = "HiC", outfolder = os.path.abspath(os.path.join(args.output_parent, args.sample_name)), args = args)

# Set up a few additional paths not in the config file
pm.config.tools.scripts_dir = os.path.join(os.path.dirname(os.path.realpath(__file__)), "tools")
pm.config.resources.ref_genome_fasta = os.path.join(pm.config.resources.genomes, args.genome_assembly, args.genome_assembly + ".fa")
pm.config.resources.chrom_sizes = os.path.join(pm.config.resources.genomes, args.genome_assembly, args.genome_assembly + ".chromSizes")
pm.config.resources.genomes_split = os.path.join(pm.config.resources.resources, "genomes_split")

pm.config.parameters.pipeline_outfolder = os.path.abspath(os.path.join(args.output_parent, args.sample_name))

print(pm.config)
tools = pm.config.tools  # Convenience alias
param = pm.config.parameters
resources = pm.config.resources

# Create a ngstk object
myngstk = pypiper.NGSTk(args.config_file)

myngstk.make_sure_path_exists(os.path.join(param.pipeline_outfolder, "unmapped_bam"))

if merge:
	# args.input is a list if merge is true;
	if not all([x.endswith(".bam") for x in args.input]):
		raise NotImplementedError("Currently we can only merge bam inputs")
	merge_folder = os.path.join(param.pipeline_outfolder, "unmapped_bam")
	sample_merged_bam = args.sample_name + ".merged.bam"
	output_merge = os.path.join(merge_folder, sample_merged_bam)
	cmd = myngstk.merge_bams(args.input, output_merge)

	pm.run(cmd, output_merge)
	args.input = os.path.join(param.pipeline_outfolder, "unmapped_bam", sample_merged_bam)  #update unmapped bam reference
	local_unmapped_bam_abs = os.path.join(param.pipeline_outfolder, "unmapped_bam", sample_merged_bam)
	input_ext = ".bam"
else:
	# Link the file into the unmapped_bam directory
	args.input = args.input[0]

	if not os.path.isabs(args.input):
		args.input = os.path.abspath(args.input)
		print args.input

	if args.input.endswith(".bam"):
		input_ext = ".bam"
	elif args.input.endswith(".fastq.gz"):
		input_ext = ".fastq.gz"
	else:
		raise NotImplementedError("This pipeline can only deal with .bam or .fastq.gz files")

	local_unmapped_bam_abs = os.path.join(param.pipeline_outfolder, "unmapped_bam", args.sample_name + input_ext)
	pm.callprint("ln -sf " + args.input + " " + local_unmapped_bam_abs, shell=True)

#check for file exists:
if not os.path.exists(local_unmapped_bam_abs):
	raise Exception(local_unmapped_bam_abs + " is not a file")

# Record file size of input file

cmd = "stat -Lc '%s' " + local_unmapped_bam_abs
input_size = pm.checkprint(cmd)
input_size = float(input_size.replace("'",""))

pm.report_result("File_mb", round((input_size/1024)/1024,2))
pm.report_result("Read_type",args.single_or_paired)
pm.report_result("Genome",args.genome_assembly)

# Fastq conversion can run out of heap space with the default java memory
# parameter for large input files.
#if input_size > 10000:
#	myngstk.set_java_mem("16g")

# Fastq conversion
################################################################################
pm.timestamp("### Fastq conversion: ")
fastq_folder = os.path.join(param.pipeline_outfolder, "00_fastq")
out_fastq_pre = os.path.join(fastq_folder, args.sample_name)
unaligned_fastq_R1 = out_fastq_pre + "_R1.fastq"
unaligned_fastq_R2 = out_fastq_pre + "_R2.fastq"

def check_fastq():
	raw_reads = myngstk.count_reads(local_unmapped_bam_abs, args.paired_end)
	pm.report_result("Raw_reads", str(raw_reads))
	fastq_reads = myngstk.count_reads(unaligned_fastq_R1, args.paired_end)
	pm.report_result("Fastq_reads", fastq_reads)
	fail_filter_reads = myngstk.count_fail_reads(local_unmapped_bam_abs, args.paired_end)
	pf_reads = int(raw_reads) - int(fail_filter_reads)
	pm.report_result("PF_reads", str(pf_reads))
	if fastq_reads != int(raw_reads):
		raise Exception("Fastq conversion error? Number of reads doesn't match unaligned bam")

if input_ext ==".bam":
	print("Found bam file")
	cmd = myngstk.bam_to_fastq(local_unmapped_bam_abs, out_fastq_pre, args.paired_end)
	pm.run(cmd, unaligned_fastq_R1, follow=check_fastq)
elif input_ext == ".fastq.gz":
	print("Found gz fastq file")
	cmd = "gunzip -c " + local_unmapped_bam_abs + " > " + unaligned_fastq_R1
	myngstk.make_sure_path_exists(fastq_folder)
	pm.run(cmd, unaligned_fastq_R1, shell=True, follow=lambda:
		pm.report_result("Fastq_reads",  myngstk.count_reads(unaligned_fastq_R1, args.paired_end)))


# HiCup Truncater
################################################################################
pm.timestamp("### HiCup Truncater: ")
trunc_folder = os.path.join(param.pipeline_outfolder, "01_trunc_" + args.genome_assembly)
myngstk.make_sure_path_exists(trunc_folder)

trunc_fastq_R1 = os.path.join(trunc_folder, args.sample_name) + "_R1.trunc.fastq"
trunc_fastq_R2 = os.path.join(trunc_folder, args.sample_name) + "_R2.trunc.fastq"

cmd = tools.perl + " " + tools.hicup_truncater
cmd += " --outdir " + trunc_folder
#cmd += " --sequence " + str(param.hicup_truncater.sequences)
cmd += " --re1 ATATCGCGG^CCGCGATAT "
cmd += " --threads " + str(param.hicup_truncater.threads)
#cmd += " --zip "
cmd += " " + unaligned_fastq_R1 + " " + unaligned_fastq_R2

pm.run(cmd, trunc_fastq_R1, shell=True)

# Clean up big intermediate files:
#pm.clean_add(os.path.join(trunc_folder, "*.fastq"))


# HiCup Mapper
################################################################################
pm.timestamp("### HiCup Mapper: ")
mapper_folder = os.path.join(param.pipeline_outfolder, "02_mapper_" + args.genome_assembly)
myngstk.make_sure_path_exists(mapper_folder)

paired_sam = os.path.join(mapper_folder, args.sample_name) + "_R1_2.pair.sam"

cmd = tools.perl + " " + tools.hicup_mapper
cmd += " --outdir " + mapper_folder
cmd += " --format Sanger "
#cmd += " --bowtie " + tools.bowtie
#cmd += " --index " + os.path.join(pm.config.resources.genomes, args.genome_assembly, indexed_bowtie, args.genome_assembly)
cmd += " --bowtie2 " + tools.bowtie2
cmd += " --index " + os.path.join(pm.config.resources.genomes, args.genome_assembly, "indexed_bowtie2", args.genome_assembly)
cmd += " --threads " + str(param.hicup_mapper.threads)
#cmd += " --zip "
cmd += " " + trunc_fastq_R1 + " " + trunc_fastq_R2

pm.run(cmd, paired_sam, shell=True)

# Clean up big intermediate files:
#pm.clean_add(os.path.join(mapper_folder, "*.fastq"))


# HiCup Deduplicator
################################################################################
pm.timestamp("### HiCup Deduplicator: ")
dedup_folder = os.path.join(param.pipeline_outfolder, "03_dedup_" + args.genome_assembly)
myngstk.make_sure_path_exists(dedup_folder)

cmd = tools.perl + " " + tools.hicup_deduplicator
cmd += " --outdir " + dedup_folder
cmd += " --threads " + str(param.hicup_deduplicator.threads)
#cmd += " --zip "
cmd += " " + paired_sam

pm.run(cmd, os.path.join(dedup_folder, args.sample_name) + "_R1_2.pair.dedup.sam", shell=True)

# Clean up big intermediate files:
#pm.clean_add(os.path.join(dedup_folder, "*.fastq"))


# Cleanup
################################################################################
pm.stop_pipeline()
