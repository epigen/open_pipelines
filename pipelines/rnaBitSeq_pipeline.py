#!/usr/bin/python2.7
"""
RNA BitSeq pipeline
documentation.
"""

from argparse import ArgumentParser
import os
import os.path
import sys
from subprocess import call
import subprocess
import re

from datetime import datetime




# Argument Parsing
# #######################################################################################
parser = ArgumentParser(description='Pypiper arguments.')

parser.add_argument('-P', '--pypiper', dest='pypiper_dir',
					default=os.path.dirname(os.path.dirname(os.path.dirname(os.path.realpath(__file__)))) + "/pypiper",
					type=str)
# Set up a pointer to the pypiper code you want to use:
# Just grab the single pypiper arg, and add it to the path here; leaving all other args for later parsing.
args = parser.parse_known_args()[0]
os.sys.path.insert(0, args.pypiper_dir)

parser.add_argument('-i', '--unmapped-bam',
					default="/fhgfs/groups/lab_bock/jklughammer/projects/otherProjects/CORE-seq/titration/CORE/unmapped_bam/BSF_0131_C5FD6ACXX_8__CORE_K562_500_1_sub.bam",
					nargs="+", dest='unmapped_bam', help="Input unmapped bam file(s))")

parser.add_argument('-s', '--sample-name', default="default",
					dest='sample_name', type=str, help='Sample Name')

parser.add_argument('-r', '--project-root', default="/fhgfs/groups/lab_bock/shared/COREseq/results_pipeline/",
					dest='project_root', type=str,
					help='Directory in which the project will reside. Default=/fhgfs/groups/lab_bock/shared/COREseq/results_pipeline/')

parser.add_argument('-g', '--genome', default="hg19_cdna",
					dest='genome_assembly', type=str, help='Genome Assembly')
parser.add_argument('-e', '--ercc', default="ERCC92",
					dest='ERCC_assembly', type=str, help='ERCC Assembly')
parser.add_argument('-f', dest='filter', action='store_false', default=True)

parser.add_argument('-p', '--paired-end', dest='paired_end', action='store_true', default=True, help='Paired End Mode')
parser.add_argument('-q', '--single-end', dest='paired_end', action='store_false', help='Single End Mode')
parser.add_argument('-C', '--no-checks', dest='no_check', action='store_true',
						default=False, help='Skip sanity checks')

# Core-seq as optional parameter
parser.add_argument('-cs', '--core-seq', default=False, dest='coreseq', action='store_true', help='CORE-seq Mode')

# Pipeline recovery interface/options
from pypiper import Pypiper
from pypiper import ngstk

parser = Pypiper.add_pypiper_args(parser)
args = parser.parse_args()



# Merging
########################################################################################
# If 2 unmapped bam files are given, then these are to be merged.
# Must be done here to initialize the sample name correctly
merge = False
if (len(args.unmapped_bam) > 1):
	merge = True
	if (args.sample_name == "default"):
		args.sample_name = "merged";
else:
	if (args.sample_name == "default"):
		args.sample_name = os.path.splitext(os.path.basename(args.unmapped_bam[0]))[0]

# Set up environment path variables
########################################################################################
# Set up an container class to hold paths
class Container:
	pass

paths = Container()
paths.scripts_dir = os.path.dirname(os.path.realpath(__file__))

# import the yaml config (work in progress)...
import yaml
pipelines_config_file = os.path.join(os.path.dirname(paths.scripts_dir), ".pipelines_config.yaml")
config = yaml.load(open(pipelines_config_file, 'r'))

# Resources
paths.resources_dir = config["paths"]["resources"]
paths.adapter_file = os.path.join(paths.resources_dir, "adapters", "epignome_adapters_2_add.fa")
paths.ref_genome = os.path.join(paths.resources_dir, "genomes")
paths.ref_genome_fasta = os.path.join(paths.resources_dir, "genomes", args.genome_assembly, args.genome_assembly + ".fa")
paths.ref_ERCC_fasta = os.path.join(paths.resources_dir, "genomes", args.ERCC_assembly, args.ERCC_assembly + ".fa")
paths.chrom_sizes = os.path.join(paths.resources_dir, "genomes", args.genome_assembly, args.genome_assembly + ".chromSizes")
paths.bowtie_indexed_genome = os.path.join(paths.resources_dir, "genomes", args.genome_assembly, "indexed_bowtie1", args.genome_assembly)
paths.bowtie_indexed_ERCC = os.path.join(paths.resources_dir, "genomes", args.ERCC_assembly, "indexed_bowtie1", args.ERCC_assembly)


# Tools
paths.trimmomatic_jar = "/cm/shared/apps/trimmomatic/0.32/trimmomatic-0.32-epignome.jar"
paths.bowtie1 = "/cm/shared/apps/bowtie/1.1.1/bin/bowtie"
paths.bowtie2 = "/cm/shared/apps/bowtie/2.2.3/bin"
paths.picard_dir = os.path.join(paths.resources_dir, "tools/picard-tools-1.100")
paths.bed2bigBed = os.path.join(paths.resources_dir, "tools", "bedToBigBed")
paths.bed2bigWig = os.path.join(paths.resources_dir, "tools", "bedGraphToBigWig")

# Output
paths.pipeline_outfolder = os.path.join(args.project_root + args.sample_name + "/")


# Run some initial setting and logging code to start the pipeline.
# Create a Pypiper object, and start the pipeline (runs initial setting and logging code to begin)
mypiper = Pypiper(name="rnaBitSeq", outfolder=paths.pipeline_outfolder, args=args)


print("N input bams:\t\t" + str(len(args.unmapped_bam)))
print("Sample name:\t\t" + args.sample_name)

sample_merged_bam = args.sample_name + ".merged.bam"
mypiper.make_sure_path_exists(paths.pipeline_outfolder + "unmapped_bam/")

if merge and not os.path.isfile(sample_merged_bam):
	print("Multiple unmapped bams found; merge requested")
	input_bams = args.unmapped_bam;
	print("input bams: " + str(input_bams))
	merge_folder = os.path.join(paths.pipeline_outfolder, "unmapped_bam/")
	input_string = " INPUT=" + " INPUT=".join(input_bams)
	output_merge = os.path.join(merge_folder, sample_merged_bam)
	cmd = "java -jar " + os.path.join(paths.picard_dir, "MergeSamFiles.jar")
	cmd += input_string
	cmd += " OUTPUT=" + output_merge
	cmd += " ASSUME_SORTED=TRUE"
	cmd += " CREATE_INDEX=TRUE"

	mypiper.call_lock(cmd, output_merge)
	args.unmapped_bam = paths.pipeline_outfolder+"unmapped_bam/" + sample_merged_bam  #update unmapped bam reference
	local_unmapped_bam = paths.pipeline_outfolder+"unmapped_bam/" + sample_merged_bam
else:
	# Link the file into the unmapped_bam directory
	print("Single unmapped bam found; no merge required")
	print("Unmapped bam:\t\t" + str(args.unmapped_bam[0]))
	args.unmapped_bam = args.unmapped_bam[0]
	local_unmapped_bam = paths.pipeline_outfolder+"unmapped_bam/"+args.sample_name+".bam"
	mypiper.callprint("ln -sf " + args.unmapped_bam + " " + local_unmapped_bam, shell=True)


print("Input Unmapped Bam: " + args.unmapped_bam)
#check for file exists:
if not os.path.isfile(local_unmapped_bam):
	print local_unmapped_bam + "is not a file"

# Fastq conversion
########################################################################################
# Uses ngstk module.

mypiper.timestamp("### Fastq conversion: ")
out_fastq_pre = os.path.join(paths.pipeline_outfolder, "fastq/", args.sample_name)
cmd = ngstk.bam_to_fastq(local_unmapped_bam, out_fastq_pre, args.paired_end, paths)
mypiper.call_lock(cmd, out_fastq_pre + "_R1.fastq")

mypiper.clean_add(out_fastq_pre + "*.fastq", conditional=True)

# Sanity checks:
if not args.no_check:
	raw_reads = ngstk.count_reads(local_unmapped_bam)
	mypiper.report_result("Raw_reads", str(raw_reads))
	fastq_reads = ngstk.count_reads(out_fastq_pre + "_R1.fastq", paired_end=args.paired_end)
	mypiper.report_result("Fastq_reads", fastq_reads)
	if (fastq_reads != int(raw_reads)):
		raise Exception("Fastq conversion error? Size doesn't match unaligned bam")


# Adapter trimming
########################################################################################
mypiper.timestamp("### Adapter trimming: ")

cmd = "java -jar "+ paths.trimmomatic_jar

if not args.paired_end:
	cmd += " SE -phred33 -threads 30"
	cmd += " -trimlog " + os.path.join(paths.pipeline_outfolder, "fastq/") + "trimlog.log "
	cmd += out_fastq_pre + "_R1.fastq "
	cmd += out_fastq_pre + "_R1_trimmed.fastq "

else:
	cmd += " PE -phred33 -threads 30"
	cmd += " -trimlog " + os.path.join(paths.pipeline_outfolder, "fastq/") + "trimlog.log "
	cmd += out_fastq_pre + "_R1.fastq "
	cmd += out_fastq_pre + "_R2.fastq "
	cmd += out_fastq_pre + "_R1_trimmed.fastq "
	cmd += out_fastq_pre + "_R1_unpaired.fastq "
	cmd += out_fastq_pre + "_R2_trimmed.fastq "
	cmd += out_fastq_pre + "_R2_unpaired.fastq "

# for Core-seq, trim off the first 6bp and the bit adjacent to identified adapter sequences:
if args.coreseq:
	cmd += " HEADCROP:6 ILLUMINACLIP:" + paths.adapter_file + ":2:10:4:1:true:epignome:5 SLIDINGWINDOW:4:1 MAXINFO:16:0.40 MINLEN:25"
# otherwise just look for normal adapters:
else:
	cmd += " ILLUMINACLIP:" + paths.adapter_file + ":2:10:4:1:true SLIDINGWINDOW:4:1 MAXINFO:16:0.40 MINLEN:21"

mypiper.call_lock(cmd, out_fastq_pre + "_R1_trimmed.fastq")
trimmed_fastq = out_fastq_pre + "_R1_trimmed.fastq"

if not args.no_check:
	x = ngstk.count_reads(trimmed_fastq,args.paired_end)
	mypiper.report_result("Trimmed_reads", x)

# RNA BitSeq pipeline.
########################################################################################
mypiper.timestamp("### Bowtie1 alignment: ")
bowtie1_folder = paths.pipeline_outfolder + "/bowtie1_" + args.genome_assembly + "/"
mypiper.make_sure_path_exists(bowtie1_folder)
out_bowtie1 = bowtie1_folder + args.sample_name + ".aln.sam"

if not args.paired_end:
	cmd = paths.bowtie1
	cmd += " -q -p 6 -a -m 100 --sam "
	cmd += paths.bowtie_indexed_genome + " "
	cmd += out_fastq_pre + "_R1_trimmed.fastq"
	cmd += " " + out_bowtie1
else:
	cmd = paths.bowtie1
	cmd += " -q -p 6 -a -m 100 --minins 0 --maxins 5000 --fr --sam --chunkmbs 200 "    # also checked --rf (1% aln) and --ff (0% aln) --fr(8% aln)
	cmd += paths.bowtie_indexed_genome
	cmd += " -1 " + out_fastq_pre + "_R1_trimmed.fastq"
	cmd += " -2 " + out_fastq_pre + "_R2_trimmed.fastq"
	cmd += " " + out_bowtie1

mypiper.call_lock(cmd, out_bowtie1)

if not args.no_check:
	x = ngstk.count_unique_mapped_reads(out_bowtie1, args.paired_end)
	mypiper.report_result("Aligned_reads", x)

mypiper.timestamp("### Raw: SAM to BAM conversion and sorting: ")

if args.filter:
	cmd = ngstk.sam_conversions(out_bowtie1,False)
	mypiper.call_lock(cmd,  re.sub(".sam$" , "_sorted.bam",out_bowtie1),shell=True)
else:
	cmd = ngstk.sam_conversions(out_bowtie1,True)
	mypiper.call_lock(cmd,  re.sub(".sam$" , "_sorted.depth",out_bowtie1),shell=True)

mypiper.clean_add(out_bowtie1, conditional=False)
mypiper.clean_add(re.sub(".sam$" , ".bam", out_bowtie1), conditional=False)


if not args.filter:
	mypiper.timestamp("### MarkDuplicates: ")

	aligned_file = re.sub(".sam$" , "_sorted.bam",out_bowtie1)
	out_file = re.sub(".sam$" , "_dedup.bam",out_bowtie1)
	metrics_file = re.sub(".sam$" , "_dedup.metrics",out_bowtie1)
	cmd = ngstk.markDuplicates(paths, aligned_file, out_file, metrics_file)
	mypiper.call_lock(cmd, out_file)

	if not args.no_check:
		x = ngstk.count_unique_mapped_reads(out_file, args.paired_end)
		mypiper.report_result("Deduplicated_reads", x)

if args.filter:
	mypiper.timestamp("### Aligned read filtering: ")
	out_sam_filter = bowtie1_folder + args.sample_name + ".aln.filt.sam"
	skipped_sam = out_sam_filter.replace(".filt." , ".skipped.")
	headerLines = subprocess.check_output("samtools view -SH " + out_bowtie1 + "|wc -l", shell=True).strip()
	cmd = "python " + paths.scripts_dir + "/bisulfiteReadFiltering_forRNA.py"
	cmd += " --infile=" + out_bowtie1
	cmd += " --outfile=" + out_sam_filter
	cmd += " --skipped=" + skipped_sam
	cmd += " --skipHeaderLines=" + headerLines
	cmd += " --genome=" + args.genome_assembly
	cmd += " --genomeDir=" + paths.ref_genome
	cmd += " --minNonCpgSites=3"
	cmd += " --minConversionRate=0.9"
	cmd += " --maxConversionRate=0.1"
	cmd += " -r"


	if args.paired_end:
		cmd = cmd + " --pairedEnd"

	mypiper.call_lock(cmd, out_sam_filter)

	if not args.no_check:
		x = ngstk.count_unique_mapped_reads(out_sam_filter, args.paired_end)
		mypiper.report_result("Filtered_reads", x)


#	Join skipped reads back in for bitseq.
#	mypiper.timestamp("### Set flag in skipped reads to 4 (unmapped): ")
#	joined_sam = out_sam_filter.replace(".filt." , ".filtFlag.")
#	skipped_sam = out_sam_filter.replace(".filt." , ".skipped.")
#	cmd = "samtools view -hS " + out_sam_filter + " > " + joined_sam + "\n"
#	cmd += "awk -v skip=" + headerLines + " -v OFS=\"\\t\" '{if (NR>skip){$2=4;$3=\"*\";$4=0;$5=0;$6=\"*\";$7=\"*\";$8=0;$9=0; print}}' " + skipped_sam
#	cmd += " >>" + joined_sam
#
#	mypiper.call_lock(cmd, joined_sam , shell=True)
#
#	if not args.no_check:
#		x = ngstk.count_unique_reads(joined_sam, args.paired_end)
#		mypiper.report_result("Joined_reads", x)
#
#	mypiper.timestamp("### Filtered: SAM to BAM conversion, sorting and depth calculation: ")
#	cmd = ngstk.sam_conversions(joined_sam)
#	mypiper.call_lock(cmd, joined_sam.replace(".sam" , "_sorted.depth"),shell=True)
#
#	mypiper.timestamp("### Skipped: SAM to BAM conversion and sorting: ")
#	cmd = ngstk.sam_conversions(skipped_sam,False)
#	mypiper.call_lock(cmd, skipped_sam.replace(".sam" , "_sorted.bam"),shell=True)
#
#	mypiper.timestamp("### MarkDuplicates: ")
#
#	aligned_file = joined_sam.replace(".sam" , "_sorted.bam")
#	out_file = joined_sam.replace(".sam" , "_dedup.bam")
#	metrics_file = joined_sam.replace(".sam" , "_dedup.metrics")
#	cmd = ngstk.markDuplicates(paths, aligned_file, out_file, metrics_file)
#	mypiper.call_lock(cmd, out_file)

	
	mypiper.timestamp("### Filtered: SAM to BAM conversion, sorting and depth calculation: ")
	cmd = ngstk.sam_conversions(out_sam_filter)
	mypiper.call_lock(cmd, re.sub(".sam$" , "_sorted.depth",out_sam_filter),shell=True)


	mypiper.timestamp("### Skipped: SAM to BAM conversion and sorting: ")
	cmd = ngstk.sam_conversions(skipped_sam,False)
	mypiper.call_lock(cmd, re.sub(".sam$", "_sorted.bam", skipped_sam),shell=True)
	
	mypiper.clean_add(skipped_sam, conditional=False)
	mypiper.clean_add(re.sub(".sam$" , ".bam", skipped_sam), conditional=False)
	mypiper.clean_add(out_sam_filter, conditional=False)
	mypiper.clean_add(re.sub(".sam$" , ".bam", out_sam_filter), conditional=False)



	mypiper.timestamp("### MarkDuplicates: ")
	
	aligned_file = re.sub(".sam$" , "_sorted.bam",out_sam_filter)
	out_file = re.sub(".sam$" , "_dedup.bam",out_sam_filter)
	metrics_file = re.sub(".sam$" , "_dedup.metrics",out_sam_filter)
	cmd = ngstk.markDuplicates(paths, aligned_file, out_file, metrics_file)
	mypiper.call_lock(cmd, out_file)
	
	if not args.no_check:
		x = ngstk.count_unique_mapped_reads(out_file, args.paired_end)
		mypiper.report_result("Deduplicated_reads", x)



# BitSeq
########################################################################################
mypiper.timestamp("### Expression analysis (BitSeq): ")

bitSeq_dir = bowtie1_folder + "/bitSeq"
mypiper.make_sure_path_exists(bitSeq_dir)
out_bitSeq = bitSeq_dir + "/" + args.sample_name + ".counts"

if args.filter:
	cmd = "Rscript " + paths.scripts_dir + "/bitSeq_parallel.R " + " " + out_sam_filter + " " + bitSeq_dir + " " + paths.ref_genome_fasta
else:
	cmd = "Rscript " + paths.scripts_dir + "/bitSeq_parallel.R " + " " + out_bowtie1 + " " + bitSeq_dir + " " + paths.ref_genome_fasta

mypiper.call_lock(cmd, out_bitSeq)


# ERCC Spike-in alignment
########################################################################################

mypiper.timestamp("### ERCC: Convert unmapped reads into fastq files: ")

unmappable_bam = re.sub(".sam$","_unmappable",out_bowtie1)
cmd = "samtools view -hbS -f4 " + out_bowtie1 + " > " + unmappable_bam + ".bam"
mypiper.call_lock(cmd, unmappable_bam, shell=True)


cmd = ngstk.bam_to_fastq(unmappable_bam + ".bam", unmappable_bam, args.paired_end, paths)
mypiper.call_lock(cmd, unmappable_bam + "_R1.fastq")

# Sanity checks:
if not args.no_check:
	raw_reads = ngstk.count_reads(unmappable_bam + ".bam",args.paired_end)
	mypiper.report_result("ERCC_raw_reads", str(raw_reads))
	fastq_reads = ngstk.count_reads(unmappable_bam + "_R1.fastq", paired_end=args.paired_end)
	mypiper.report_result("ERCC_fastq_reads", fastq_reads)
	if (fastq_reads != int(raw_reads)):
		raise Exception("Fastq conversion error? Size doesn't match unaligned bam")



mypiper.timestamp("### ERCC: Bowtie1 alignment: ")
bowtie1_folder = paths.pipeline_outfolder + "/bowtie1_" + args.ERCC_assembly + "/"
mypiper.make_sure_path_exists(bowtie1_folder)
out_bowtie1 = bowtie1_folder + args.sample_name + "_ERCC.aln.sam"

if not args.paired_end:
	cmd = paths.bowtie1
	cmd += " -q -p 6 -a -m 100 --sam "
	cmd += paths.bowtie_indexed_ERCC + " "
	cmd += unmappable_bam + "_R1.fastq"
	cmd += " " + out_bowtie1
else:
	cmd = paths.bowtie1
	cmd += " -q -p 6 -a -m 100 --minins 0 --maxins 5000 --fr --sam --chunkmbs 200 "
	cmd += paths.bowtie_indexed_ERCC
	cmd += " -1 " + unmappable_bam + "_R1.fastq"
	cmd += " -2 " + unmappable_bam + "_R2.fastq"
	cmd += " " + out_bowtie1

mypiper.call_lock(cmd, out_bowtie1)

if not args.no_check:
	x = ngstk.count_unique_mapped_reads(out_bowtie1, args.paired_end)
	mypiper.report_result("ERCC_aligned_reads", x)


mypiper.timestamp("### ERCC: SAM to BAM conversion, sorting and depth calculation: ")
cmd = ngstk.sam_conversions(out_bowtie1)
mypiper.call_lock(cmd, re.sub(".sam$" , "_sorted.depth", out_bowtie1), shell=True)

mypiper.clean_add(out_bowtie1, conditional=False)
mypiper.clean_add(re.sub(".sam$" , ".bam", out_bowtie1), conditional=False)
mypiper.clean_add(unmappable_bam + "*.fastq", conditional=False)

# BitSeq
########################################################################################
mypiper.timestamp("### ERCC: Expression analysis (BitSeq): ")

bitSeq_dir = bowtie1_folder + "/bitSeq"
mypiper.make_sure_path_exists(bitSeq_dir)
out_bitSeq = bitSeq_dir + "/" + re.sub(".aln.sam$" , ".counts",out_bowtie1)

cmd = "Rscript " + paths.scripts_dir + "/bitSeq_parallel.R " + " " + out_bowtie1 + " " + bitSeq_dir + " " + paths.ref_ERCC_fasta

mypiper.call_lock(cmd, out_bitSeq)



# Cleanup
########################################################################################


# remove temporary marker file:
mypiper.stop_pipeline()





