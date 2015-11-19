#!/usr/bin/env python
"""
RRBS pipeline

"""

from argparse import ArgumentParser
import os, re
import os.path
import sys
import subprocess
import yaml

# REMARK: The following slurm modules are needed:
# module load bsmap (even though we use path.bsmap as well, some BSMAP scripts need to be in path, potential source of bugs TODO)

# Argument Parsing
################################################################################
parser = ArgumentParser(description='Pypiper arguments.')

parser.add_argument('-P', '--pypiper', dest='pypiper_dir',
					default=os.path.dirname(os.path.dirname(os.path.dirname(os.path.realpath(__file__)))) + "/pypiper",
					type=str)
# Set up a pointer to the pypiper code you want to use:
# Just grab the single pypiper arg, and add it to the path here; leaving all other args for later parsing.
args = parser.parse_known_args()[0]
os.sys.path.insert(0, args.pypiper_dir)
from pypiper import Pypiper
from pypiper import ngstk
# Add Pypiper arguments to the arguments list and then reparse.
parser = Pypiper.add_pypiper_args(parser)

# Pipeline options
## Generic options for all pipelines
parser.add_argument('-i', '--unmapped-bam', default="default.bam", nargs="+", dest='unmapped_bam',
					help="Input unmapped bam file(s))") # here a default value doesn't make sense
parser.add_argument('-s', '--sample-name', default="default", dest='sample_name', type=str,
					help='Sample Name') # default means deduction from filename, except .bam extension
parser.add_argument('-r', '--project-root', default="", dest='project_root', type=str,
					help='Directory in which the project will reside.')
parser.add_argument('-g', '--genome', default="hg19", dest='genome_assembly', type=str, help='Genome Assembly')
parser.add_argument('-C', '--no-checks', dest='no_check', action='store_true',
						default=False, help='Skip sanity checks')
parser.add_argument('-p', '--paired-end', dest='paired_end', action='store_true', default=False, help='Paired End Mode')
parser.add_argument('-q', '--single-end', dest='paired_end', action='store_false', help='Single End Mode')
# AS: RRBS as used by us is mainly single-end... maybe it should be the default? Will break with other pipeline's defaults?
# AS: I got the code from Angelo to run it in paired end mode

## Pipeline-specific options
parser.add_argument('-t', '--trimgalore', dest='trimmomatic', action="store_false", default=True, help='Use trimgalore instead of trimmomatic?')
args = parser.parse_args()

# Merging
################################################################################
# If 2 unmapped bam files are given, then these are to be merged.
# Must be done here to initialize the sample name correctly

merge = False
if len(args.unmapped_bam) > 1:
	merge = True
	if args.sample_name == "default":
		args.sample_name = "merged"
else:
	if args.sample_name == "default":
		args.sample_name = os.path.splitext(os.path.basename(args.unmapped_bam[0]))[0]

# Set up environment path variables
################################################################################
# Set up an container class to hold paths
class Container:
	pass

paths = Container()
paths.scripts_dir = os.path.dirname(os.path.realpath(__file__))

# import the yaml config (work in progress)...
pipelines_config_file = os.path.join(os.path.dirname(paths.scripts_dir), ".pipelines_config.yaml")
config = yaml.load(open(pipelines_config_file, 'r'))

# Resources
paths.resources_dir = config["paths"]["resources"]
paths.adapter_file = os.path.join(paths.resources_dir, "adapters", "epignome_adapters_2_add.fa")
paths.rrbs_adapter_file = os.path.join(paths.resources_dir, "adapters", "RRBS_adapters.fa")
paths.ref_genome = os.path.join(paths.resources_dir, "genomes")
paths.ref_genome_fasta = os.path.join(paths.resources_dir, "genomes", args.genome_assembly, args.genome_assembly + ".fa")
paths.chrom_sizes = os.path.join(paths.resources_dir, "genomes", args.genome_assembly, args.genome_assembly + ".chromSizes")

paths.bismark_indexed_genome = os.path.join(paths.resources_dir, "genomes", args.genome_assembly, "indexed_bismark_bt2")
paths.bismark_spikein_genome = os.path.join(paths.resources_dir, "genomes", "meth_spikein_k1_k3", "indexed_bismark_bt1")

# Tools
paths.picard_dir = os.path.join(paths.resources_dir, "tools/picard-tools-1.100")
paths.trimmomatic_jar = "/cm/shared/apps/trimmomatic/0.32/trimmomatic-0.32-epignome.jar"
paths.bowtie1 = "/cm/shared/apps/bowtie/1.1.1/bin"
paths.bowtie2 = "/cm/shared/apps/bowtie/2.2.3/bin"
paths.bed2bigBed = os.path.join(paths.resources_dir, "tools", "bedToBigBed")
paths.bed2bigWig = os.path.join(paths.resources_dir, "tools", "bedGraphToBigWig")
paths.bismark = "/cm/shared/apps/bismark/0.12.2/bismark"
paths.deduplicate_bismark = "/cm/shared/apps/bismark/0.12.2/deduplicate_bismark"
#paths.bsmap = "/cm/shared/apps/bsmap/2.74/bsmap"   # module load bsmap/2.74

#Output
paths.pipeline_outfolder_abs = os.path.abspath(os.path.join(args.project_root, args.sample_name))

#Biseq paths:
paths.biseq_tools_dir = os.path.join(paths.resources_dir, "tools")
paths.genomes_split = os.path.join(paths.resources_dir, "genomes_split")

# Create a Pypiper object, and start the pipeline (runs initial setting and logging code to begin)

mypiper = Pypiper(name="RRBS", outfolder=paths.pipeline_outfolder_abs, args=args)

ngstk.make_sure_path_exists(os.path.join(paths.pipeline_outfolder_abs, "unmapped_bam"))

if merge:
	raise NotImplementedError("Sample merging currently not implemented for RRBS")  # TODO AS: merge currently deactivated for RRBS

	# inactive code
	input_bams = args.unmapped_bam
	merge_folder = os.path.join(paths.pipeline_outfolder_abs, "unmapped_bam")
	sample_merged_bam = args.sample_name + ".merged.bam"
	output_merge = os.path.join(merge_folder, sample_merged_bam)
	cmd = ngstk.merge_bams(input_bams, output_merge, paths)

	mypiper.call_lock(cmd, output_merge)
	args.unmapped_bam = os.path.join(paths.pipeline_outfolder_abs, "unmapped_bam", sample_merged_bam)  #update unmapped bam reference
	local_unmapped_bam_abs = os.path.join(paths.pipeline_outfolder_abs, "unmapped_bam", sample_merged_bam)

else:
	# Link the file into the unmapped_bam directory
#	print("Single unmapped bam found; no merge required")
#	print("Unmapped bam:\t\t" + str(args.unmapped_bam[0]))
	args.unmapped_bam = args.unmapped_bam[0]

	if not os.path.isabs(args.unmapped_bam):
		args.unmapped_bam = os.path.abspath(args.unmapped_bam)
		print args.unmapped_bam

	local_unmapped_bam_abs = os.path.join(paths.pipeline_outfolder_abs, "unmapped_bam", args.sample_name + ".bam")
	mypiper.callprint("ln -sf " + args.unmapped_bam + " " + local_unmapped_bam_abs, shell=True)

#check for file exists:
if not os.path.exists(local_unmapped_bam_abs):
	raise Exception(local_unmapped_bam_abs + " is not a file")


# Fastq conversion
################################################################################
mypiper.timestamp("### Fastq conversion: ")
fastq_folder = os.path.join(paths.pipeline_outfolder_abs, "fastq")
out_fastq_pre = os.path.join(fastq_folder, args.sample_name)
cmd = ngstk.bam_to_fastq(local_unmapped_bam_abs, out_fastq_pre, args.paired_end, paths)

unaligned_fastq = out_fastq_pre + "_R1.fastq"
mypiper.call_lock(cmd, unaligned_fastq)
mypiper.clean_add(out_fastq_pre + "*.fastq", conditional=True)

# Sanity checks:
if not args.no_check:
	raw_reads = ngstk.count_reads(local_unmapped_bam_abs, args.paired_end)
	mypiper.report_result("Raw_reads", str(raw_reads))
	fastq_reads = ngstk.count_reads(unaligned_fastq, args.paired_end)
	mypiper.report_result("Fastq_reads", fastq_reads)
	if fastq_reads != int(raw_reads):
		raise Exception("Fastq conversion error? Number of reads doesn't match unaligned bam")


# Adapter trimming (Trimmomatic)
################################################################################
mypiper.timestamp("### Adapter trimming: ")

if args.trimmomatic:
	trimmed_fastq = out_fastq_pre + "_R1_trimmed.fq"

	# REMARK AS: instead of trim_galore we try to use Trimmomatic for now
	# - we are more compatible with the other pipelines
	# - better code base, not a python wrapper of a perl script (as trim_galore)
	# - rrbs-mode not needed because biseq has the same functionality

	# REMARK NS:
	# The -Xmx4000m restricts heap memory allowed to java, and is necessary 
	#  to prevent java from allocating lots of memory willy-nilly 
	# if it's on a machine with lots of memory, which can lead
	# to jobs getting killed by a resource manager. By default, java will
	# use more memory on systems that have more memory, leading to node-dependent
	# killing effects that are hard to trace.
	if not args.paired_end:
		cmd = "java -Xmx4000m -jar  " + paths.trimmomatic_jar + " SE -phred33 -threads 30"
		cmd += " -trimlog " + os.path.join(fastq_folder, "trimlog.log") + " "
		cmd += out_fastq_pre + "_R1.fastq "
		cmd += out_fastq_pre + "_R1_trimmed.fq "
		cmd += "ILLUMINACLIP:" + paths.rrbs_adapter_file + ":2:40:7 SLIDINGWINDOW:4:20 MAXINFO:20:0.60 MINLEN:20"
		#updated by AS TO RRBS mode on 2015-05-06 according to 2014 optimizations

	else:
		raise NotImplementedError("Trimmomatic for PE RRBS not implemented")

		# inactive code from here for now
		cmd = "java -Xmx4000m -jar " + paths.trimmomatic_jar + " PE"
		cmd += " -phred33"
		cmd += " -threads 30"
		cmd += " -trimlog " + os.path.join(fastq_folder, "trimlog.log") + " "
		cmd += out_fastq_pre + "_R1.fastq "
		cmd += out_fastq_pre + "_R2.fastq "
		cmd += out_fastq_pre + "_R1_trimmed.fq "
		cmd += out_fastq_pre + "_R1_unpaired.fq "
		cmd += out_fastq_pre + "_R2_trimmed.fq "
		cmd += out_fastq_pre + "_R2_unpaired.fq "
		cmd += "ILLUMINACLIP:" + paths.adapter_file + ":2:10:4:1:true:epignome:5 SLIDINGWINDOW:4:1 MAXINFO:16:0.40 MINLEN:25"

else: # use trim_galore
	# Trim galore requires biopython, cutadapt modules. RSeQC as well (maybe?)
	#   --- $trim_galore -q $q --phred33 -a $a --stringency $s -e $e --length $l --output_dir $output_dir $input_fastq
	if args.paired_end:
		raise NotImplementedError("TrimGalore for PE RRBS not implemented")
	input_fastq = out_fastq_pre + "_R1.fastq "

	# With trimgalore, the output file is predetermined.
	trimmed_fastq = out_fastq_pre + "_R1_trimmed.fq"

	output_dir=fastq_folder

	#Adapter
	a="AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"

	cmd = "trim_galore"
	cmd += " -q 20" 	#quality trimming
	cmd += " --phred33 -a " + a
	cmd += " --stringency 1"		#stringency: Overlap with adapter sequence required to trim a sequence
	cmd += " -e 0.1" 	#Maximum allowed error rate
	cmd += " --length 16"	#Minimum Read length
	# by unchangeable default Trimmomatic discards reads of lenth 0 (produced by ILLUMINACLIP):
	cmd += " --output_dir " + output_dir + " " + input_fastq

# Trimming command has been constructed, using either trimming options.
# The code to run it is the same either way:

mypiper.call_lock(cmd, trimmed_fastq)   # TODO AS: maybe another lock_name?
if not args.no_check:
	x = ngstk.count_reads(trimmed_fastq, args.paired_end)
	mypiper.report_result("Trimmed_reads", x)

# RRBS alignment with BSMAP.
################################################################################
mypiper.timestamp("### BSMAP alignment: ")
bsmap_folder = os.path.join(paths.pipeline_outfolder_abs, "bsmap_" + args.genome_assembly)  # e.g. bsmap_hg19
ngstk.make_sure_path_exists(bsmap_folder)
# no tmp folder needed for BSMAP alignment

out_bsmap = os.path.join(bsmap_folder, args.sample_name + ".bam")
#original filename in bash version of pipeline:
#bsmap_aligned_bam=$align_dir/$sample.all.aligned.bsmap.$mismatches.mism.$allow_multimapper.r.bam


# REMARK NS: In previous versions of the pipeline, TRIMGALORE used a .fq
# extension, while trimmomatic used .fastq.
# I updated it so that the trimmmomatic path also outputs a .fq file, so this
# doesn't have to vary based on trimmer.
if not args.paired_end:
	cmd = "bsmap -a " + out_fastq_pre + "_R1_trimmed.fq"    # REMARK: bsmap needs to be in PATH (or load module)
	cmd += " -d " + paths.ref_genome_fasta
	cmd += " -o " + out_bsmap
	cmd += " -D C-CGG"  # " -D C-CGG "  # set to " " to disable!
	cmd += " -w 100"    # -w 100 is default (according to log output, help says up to 1000 locations are searched)
	cmd += " -v 0.08"
	cmd += " -r 1"      # -r <int>  = multimapper enabled (1) or disabled (0)
	cmd += " -p 4"      # -p <int> = num processors
	cmd += " -n 0"      # -n <int> = mapping strand information (2 strands or 4 strands)
	cmd += " -S 1"      # -S <int> = set a fixed seed for random number generator for reproducibility
	cmd += " -f 5"      # DEFAULT -f <int> = filter low-quality reads containing >n Ns, default=5
	#cmd += " -q 0"
	#cmd += " -A ''"    # -A <str>: Adapter trimming (is done in the step therefore disabled)
	cmd += " -u "       # report unmapped reads (into same bam file)

	# TODO AS: More configurability needed
	# possibility to enable / disable multimappers
	# possibility to enable / disable RRBS mode (only align reads starting with a certain restriction site)

else:
	raise NotImplementedError("Alignment for paired-end RRBS data currently not implemented!")

mypiper.call_lock(cmd, out_bsmap)

if not args.no_check:
	# BSMap apparently stores all the reads (mapped and unmapped) in
	# its output bam; to count aligned reads, then, we have to use
	# a -F4 flag (with count_mapped_reads instead of count_reads).
	x = ngstk.count_mapped_reads(out_bsmap, args.paired_end)
	mypiper.report_result("Aligned_reads", x)
	# In addition, BSMap can (if instructed by parameters) randomly assign
	# multimapping reads. It's useful to know how many in the final bam were such.
	x = ngstk.count_multimapping_reads(out_bsmap, args.paired_end)
	mypiper.report_result("Multimap_reads", x)


# The WGBS pipeline does some steps here which are not relevant for RRBS:
# mypiper.timestamp("### PCR duplicate removal: ")  # because the RRBS stacks of reads would be negatively affected
# mypiper.timestamp("### Aligned read filtering: ") #
# mypiper.timestamp("### Methylation calling (bismark extractor): ")    # we use biseq-methcalling instead

# Clean up big intermediate files:
mypiper.clean_add(bsmap_folder + "/*.fastq") # initial unmapped fastqs

# Run biseq-methcalling:
################################################################################
mypiper.timestamp("### biseqMethCalling: ")

# Python Software Requirements for biseqMethCalling
# REMARK AS: all packages are available via "easy_install --user <lib>"
# pip is also a possibility if available (currently not on CeMM infrastructure)
#
# Direct links just in case:
# - biopython: wget https://pypi.python.org/pypi/biopython or wget http://biopython.org/DIST/biopython-1.63.zip
# - bitarray: wget https://pypi.python.org/packages/source/b/bitarray/bitarray-0.8.1.tar.gz
# - guppy: wget https://pypi.python.org/packages/source/g/guppy/guppy-0.1.10.tar.gz
# - pysam: wget https://code.google.com/p/pysam/downloads/detail?name=pysam-0.7.5.tar.gz

biseq_output_path = os.path.join(paths.pipeline_outfolder_abs, "biseqMethcalling_" + args.genome_assembly)
biseq_output_path_web = os.path.join(biseq_output_path, "web")
biseq_output_path_temp = os.path.join(biseq_output_path, "temp")

ngstk.make_sure_path_exists (biseq_output_path)

cmd = "python -u " + os.path.join(paths.scripts_dir, "biseqMethCalling.py")
cmd += " --sampleName=" + args.sample_name
cmd += " --alignmentFile=" + out_bsmap      # this is the absolute path to the bsmap aligned bam file
cmd += " --methodPrefix=RRBS"
cmd += " --rrbsMode"
cmd += " --checkRestriction"
cmd += " --minFragmentLength=20"
cmd += " --maxFragmentLength=1000"
cmd += " --pfStatus=All"
cmd += " --maxMismatches=0.1"
cmd += " --baseQualityScoreC=20"
cmd += " --baseQualityScoreNextToC=10"
cmd += " --laneSpecificStatistics"
cmd += " --bigBedFormat"
cmd += " --deleteTemp"
cmd += " --toolsDir=" + paths.biseq_tools_dir
cmd += " --outputDir=" + biseq_output_path
cmd += " --webOutputDir=" + biseq_output_path_web
cmd += " --tempDir=" + biseq_output_path_temp
cmd += " --timeDelay=0"
cmd += " --genomeFraction=50"
cmd += " --smartWindows=250000"
cmd += " --maxProcesses=4"
cmd += " --genomeDir=" + paths.genomes_split
cmd += " --inGenome=" + args.genome_assembly
cmd += " --outGenome=" + args.genome_assembly
# TODO AS: Investigate what happens with biseq in the case of paired-end data

# The dog genome has 38 chromosomes (plus one X chromosome). It's probably best to check here for these rarely used
# reference genomes:
# The default value for includedChromosomes is chr1-30, X, Y, Z (sufficient for human and mouse genomes)
# REMARK NS: This is a hack to account for the way biseqMethCalling restricts to
# default chroms. THis should be fixed in biseq in the future, but for now, this
# lets us run dog samples using the default pipeline. hack!
if args.genome_assembly == "canFam3":
	cmd += ' --includedChromosomes="chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,' \
	       'chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chr23,chr24,chr25,chr26,chr27,chr28,chr29,chr30,chrX,' \
	       'chrY,chrZ,chr31,chr32,chr33,chr34,chr35,chr36,chr37,chr38"'

# Deactivated options:
#cmd += " --appendStatisticsOutput=" + stat_output  # TODO AS: I disable this option for now. This is an analysis-global file where every biseq run writes to
#stat_output = os.path.join(biseq_output_path, "RRBS_biseq_statistics.txt")  # general stats file independent of sample

biseq_finished_helper = os.path.join(biseq_output_path, "biseq.completed")
cmd2 = "touch " + biseq_finished_helper

mypiper.call_lock([cmd, cmd2], target=biseq_finished_helper)

# Now parse some results for pypiper result reporting.
cmd = "python " + os.path.join(paths.scripts_dir, "tsv_parser.py")
cmd += " -i " + biseq_output_path + "/RRBS_statistics_" + args.sample_name + ".txt"
cmd += " -c uniqueSeqMotifCount"

x = mypiper.checkprint(cmd, shell=True)
mypiper.report_result("Unique_CpGs", x)

################################################################################
mypiper.timestamp("### Make bigbed: ")
# REMARK AS: Make bigwig uses a bismark output file. For RRBS we don't have the bismark cov file
# (essentially a bedgraph file) which the tool bed2bigWig would need
# REMARK AS: UCSC tracks are generated by biseq-methcalling


# First, convert the bed format into the bigBed input style.
# This is how biseq did it, but it's actually unnecessary; instead we can just go straight off the output file.
# Left command here for posterity.
# awk '{ printf "%s\t%s\t%s\t\047%s%[\04720\047]\047\t%s\t%s\n", $1, $2, $3, $5/10, $5, $6 }' RRBS_cpgMethylation_01_2276TU.bed > f


# bigbed conversion input file is the biseq methylation calls output file
biseq_methcall_file = biseq_output_path + "/RRBS_cpgMethylation_" + args.sample_name + ".bed"

bigbed_output_path = os.path.join(paths.pipeline_outfolder_abs, "bigbed_" + args.genome_assembly)
bigwig_output_path = os.path.join(paths.pipeline_outfolder_abs, "bigwig_" + args.genome_assembly)


# bedToBigBed RRBS_cpgMethylation_01_2276TU.bed ~/linkto/resources/genomes/hg19/hg19.chromSizes RRBS_cpgMethylation_test2.bb

ngstk.make_sure_path_exists (bigbed_output_path)
ngstk.make_sure_path_exists (bigwig_output_path)
bigbed_output_file = bigbed_output_path + "/RRBS_" + args.sample_name + ".bb"
out_bedGraph = bigwig_output_path + "/RRBS_" + args.sample_name + ".bedGraph"
out_bigwig = bigwig_output_path + "/RRBS_" + args.sample_name + ".bw"


cmd = paths.bed2bigBed
cmd += " " + biseq_methcall_file
cmd += " " + paths.chrom_sizes
cmd += " " + bigbed_output_file

# REMARK NS: As of June 2015, IGV will load bigBed files for methylation
# in a unique format if the *filename contains  "RRBS_cpgMethylation" -- see
# https://github.com/igvteam/igv/blob/master/src/org/broad/igv/methyl/MethylTrack.java
# This is obviously not ideal, but I will create a link with this filename
# to the original file (even for WGBS tracks) so that you could load these into 
# IGV if you want:

filename_hack_link_file = bigbed_output_path + "/RRBS_cpgMethylation_" + args.sample_name + ".bb"
cmd2 = "ln -sf " + bigbed_output_file + " " + filename_hack_link_file

mypiper.call_lock([cmd, cmd2], bigbed_output_file)

# Let's also make bigwigs:

# First convert to bedGraph
# hard coding tabs doesn't seem to work:
#cmd = "awk '{ printf \\\"%s\t%s\t%s\t%s\n\\\", $1, $2, $3, $5/10 }'"
cmd = "awk -v OFS='\t' '{ print $1, $2, $3, $5/10 }'"
cmd += " " + biseq_methcall_file
cmd += " > " + out_bedGraph

cmd2 = paths.bed2bigWig
cmd2 += " " + out_bedGraph
cmd2 += " " + paths.chrom_sizes
cmd2 += " " + out_bigwig

mypiper.call_lock([cmd, cmd2], out_bigwig, shell=True)

################################################################################
# Calculate neighbor methylation matching
mypiper.timestamp("### Neighbor Methylation Matching: ")
nmm_output_dir = os.path.join(paths.pipeline_outfolder_abs, "nmm_" + args.genome_assembly)
ngstk.make_sure_path_exists (nmm_output_dir)
nmm_outfile=os.path.join(nmm_output_dir, args.sample_name + ".nmm.bed")

cmd = "python -u " + os.path.join(paths.scripts_dir, "methylMatch.py")
cmd += " --inFile=" + out_bsmap      # this is the absolute path to the bsmap aligned bam file
cmd += " --methFile=" + biseq_methcall_file
cmd += " --outFile=" + nmm_outfile
cmd += " --cores=4"
cmd += " -q"

mypiper.call_lock(cmd, nmm_outfile)

################################################################################
mypiper.timestamp("### Bismark spike-in alignment: ")
# currently using bowtie1 instead of bowtie2

# get unaligned reads out of BSMAP bam
bsmap_unalignable_bam = os.path.join(bsmap_folder, args.sample_name + "_unalignable.bam")
mypiper.call_lock("samtools view -bh -f 0x4 "+out_bsmap+" > " + bsmap_unalignable_bam, bsmap_unalignable_bam, shell=True)

# convert BAM to fastq
bsmap_fastq_unalignable_pre = os.path.join(bsmap_folder, args.sample_name + "_unalignable")
bsmap_fastq_unalignable = bsmap_fastq_unalignable_pre  + "_R1.fastq"
cmd = ngstk.bam_to_fastq(bsmap_unalignable_bam, bsmap_fastq_unalignable_pre, args.paired_end, paths)
mypiper.call_lock(cmd, bsmap_fastq_unalignable)

# actual spike-in analysis
spikein_folder = os.path.join(paths.pipeline_outfolder_abs, "bismark_spikein")
ngstk.make_sure_path_exists(spikein_folder)
spikein_temp = os.path.join(spikein_folder, "bismark_temp")
ngstk.make_sure_path_exists(spikein_temp)
out_spikein_base = args.sample_name + ".spikein.aln"

if not args.paired_end:
	out_spikein = os.path.join(spikein_folder, out_spikein_base + ".bam")
	cmd = paths.bismark + " " + paths.bismark_spikein_genome + " "
	cmd += bsmap_fastq_unalignable_pre + "_R1.fastq"
	cmd += " --bam --unmapped"
	cmd += " --path_to_bowtie " + paths.bowtie1
	#	cmd += " --bowtie2"
	cmd += " --temp_dir " + spikein_temp
	cmd += " --output_dir " + spikein_folder
	cmd += " --basename=" + out_spikein_base
#	cmd += " -p 4"

else:
	raise NotImplementedError("Spike-in not implemented for paired-end RRBS")

	# inaccessible code from here
	out_spikein = os.path.join(spikein_folder, out_spikein_base + "_pe.bam")
	cmd = paths.bismark + " " + paths.bismark_spikein_genome
	cmd += " --1 " + unmapped_reads_pre + "_unmapped_reads_1.fq"    # filenames not updated as in SE case
	cmd += " --2 " + unmapped_reads_pre + "_unmapped_reads_2.fq"
	#cmd += " --1 " + unmapped_reads_pre + "_R1_trimmed.fastq_unmapped_reads_1.fq"
	#cmd += " --2 " + unmapped_reads_pre + "_R2_trimmed.fastq_unmapped_reads_2.fq"
	cmd += " --bam --unmapped"
	cmd += " --path_to_bowtie " + paths.bowtie1
	#	cmd += " --bowtie2"
	cmd += " --temp_dir " + spikein_temp
	cmd += " --output_dir " + spikein_folder
	cmd += " --minins 0 --maxins 5000"
	cmd += " --basename=" + out_spikein_base
#	cmd += " -p 4"

mypiper.call_lock(cmd, out_spikein, nofail=True)

# Clean up the unmapped file which is copied from the parent
# bismark folder to here:
mypiper.clean_add(spikein_folder + "/*.fq", conditional=False)
mypiper.clean_add(spikein_temp) # For some reason, the temp folder is not deleted.



################################################################################
mypiper.timestamp("### PCR duplicate removal (Spike-in): ")
# Bismark's deduplication forces output naming, how annoying.
#out_spikein_dedup = spikein_folder + args.sample_name + ".spikein.aln.deduplicated.bam"
out_spikein_dedup = re.sub(r'.bam$', '.deduplicated.bam', out_spikein)
if not args.paired_end:
	cmd = paths.deduplicate_bismark + " --single "    # TODO: needs module load bismark or absolute path to this tool
	cmd += out_spikein
	cmd += " --bam"
else:
	cmd = paths.deduplicate_bismark + " --paired "
	cmd += out_spikein
	cmd += " --bam"

out_spikein_sorted = out_spikein_dedup.replace('.deduplicated.bam', '.deduplicated.sorted')
cmd2 = "samtools sort " + out_spikein_dedup + " " + out_spikein_sorted
cmd3 = "samtools index " + out_spikein_sorted + ".bam"
cmd4 = "rm " + out_spikein_dedup
mypiper.call_lock([cmd, cmd2, cmd3, cmd4], out_spikein_sorted + ".bam.bai", nofail=True)

# Spike-in methylation calling
################################################################################
mypiper.timestamp("### Methylation calling (testxmz) Spike-in: ")

cmd1 = "python -u " + os.path.join(paths.scripts_dir, "testxmz.py")
cmd1 += " " + out_spikein_sorted + ".bam" + " " + "K1_unmethylated"
cmd1 += " >> " + mypiper.pipeline_stats_file
cmd2 = cmd1.replace("K1_unmethylated", "K3_methylated")
mypiper.callprint(cmd1, shell=True, nofail=True)
mypiper.callprint(cmd2, shell=True, nofail=True)

# TODO: check if results file exists already...


# PDR calculation:
################################################################################

mypiper.timestamp("### PDR (Partial Disordered Methylation) analysis")

pdr_output_dir = os.path.join(paths.pipeline_outfolder_abs, "pdr_" + args.genome_assembly)
ngstk.make_sure_path_exists (pdr_output_dir)

# convert aligned bam to sam

pdr_in_samfile = os.path.join(pdr_output_dir, args.sample_name + ".aligned.sam") # gets deleted after, see some lines below
mypiper.call_lock("samtools view " + out_bsmap + " > " + pdr_in_samfile, pdr_in_samfile, shell=True)

# PDR calculation:
#
# output files:
pdr_bedfile=os.path.join(pdr_output_dir, args.sample_name + ".pdr.bed")

produce_sam = True  # TODO AS: make this an option somewhere
concordsam=os.path.join(pdr_output_dir, args.sample_name + ".concordant.sam")
discordsam=os.path.join(pdr_output_dir, args.sample_name + ".discordant.sam")

# command::
cmd1 = "python -u " + os.path.join(paths.scripts_dir, "bisulfiteReadConcordanceAnalysis.py")
cmd1 += " --infile=" + pdr_in_samfile
cmd1 += " --outfile=" + pdr_bedfile
cmd1 += " --skipHeaderLines=0"
cmd1 += " --genome=" + args.genome_assembly
cmd1 += " --genomeDir=" + paths.ref_genome
cmd1 += " --minNonCpgSites=3"   # These two parameters are not relevant for PDR analysis
cmd1 += " --minConversionRate=0.9"

if produce_sam == True:
	cmd1 += " --concordantOutfile=" + concordsam
	cmd1 += " --discordantOutfile=" + discordsam
	#TODO: perhaps convert them to bam *cough*

#call:
mypiper.call_lock(cmd1, pdr_bedfile)

# delete huge input SAM file
mypiper.clean_add(pdr_in_samfile)

# Final sorting and indexing
################################################################################
# create sorted and indexed BAM files for visualization and analysis
# bsmap already outputs a sorted and indexed bam file

# Cleanup
################################################################################
mypiper.stop_pipeline()
