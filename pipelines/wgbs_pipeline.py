#!/usr/bin/python2.7
"""
WGBS pipeline
documentation.
"""
from argparse import ArgumentParser
import os, re
import os.path
import sys
from subprocess import call
import subprocess
from datetime import datetime
import yaml # In process of moving config to yaml
# Argument Parsing
################################################################################
parser = ArgumentParser(description='Pypiper arguments.')
# Set up a pointer to the pypiper code you want to use:
# Just grab the single pypiper arg, and add it to the path here; leaving all other args for later parsing.
parser.add_argument('-P', '--pypiper', dest='pypiper_dir',
					default=os.path.dirname(os.path.dirname(os.path.dirname(os.path.realpath(__file__)))) + "/pypiper",
					type=str)

args = parser.parse_known_args()[0]
os.sys.path.insert(0, args.pypiper_dir)
from pypiper import Pypiper
from pypiper import ngstk
# Add Pypiper arguments to the arguments list
parser = Pypiper.add_pypiper_args(parser)

parser.add_argument('-i', '--unmapped-bam',
					default="",
					nargs="+", dest='unmapped_bam', help="Input unmapped bam file(s))")

parser.add_argument('-s', '--sample-name', default="default",
					dest='sample_name', type=str, help='Sample Name')

parser.add_argument('-r', '--project-root', default="/data/groups/lab_bock/shared/",
					dest='project_root', type=str,
					help='Directory in which the project will reside. Default=/data/groups/lab_bock/shared/')

parser.add_argument('-g', '--genome', default="hg19",
					dest='genome_assembly', type=str, help='Genome Assembly')
parser.add_argument('-C', '--no-checks', dest='no_check', action='store_true', default=False, help='Skip sanity checks')
parser.add_argument('-p', '--paired-end', dest='paired_end', action='store_true', default=True, help='Paired End Mode')
parser.add_argument('-q', '--single-end', dest='paired_end', action='store_false', help='Single End Mode')
args = parser.parse_args()

# Merging
################################################################################
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
paths.ref_genome = os.path.join(paths.resources_dir, "genomes")
paths.bismark_indexed_genome = os.path.join(paths.resources_dir, "genomes", args.genome_assembly, "indexed_bismark_bt2")
paths.bismark_spikein_genome = os.path.join(paths.resources_dir, "genomes", "meth_spikein_k1_k3", "indexed_bismark_bt1")
paths.chrom_sizes = os.path.join(paths.resources_dir, "genomes", args.genome_assembly, args.genome_assembly + ".chromSizes")

# Tools
paths.picard_dir = os.path.join(paths.resources_dir, "tools/picard-tools-1.100")
paths.trimmomatic_jar = "/cm/shared/apps/trimmomatic/0.32/trimmomatic-0.32-epignome.jar"
paths.bowtie1 = "/cm/shared/apps/bowtie/1.1.1/bin"
paths.bowtie2 = "/cm/shared/apps/bowtie/2.2.3/bin"
paths.bed2bigBed = os.path.join(paths.resources_dir, "tools", "bedToBigBed")
paths.bed2bigWig = os.path.join(paths.resources_dir, "tools", "bedGraphToBigWig")

# Output
paths.pipeline_outfolder = os.path.join(args.project_root + args.sample_name + "/")


# Create a Pypiper object, and start the pipeline
# (this runs initial setting and logging code to begin)

mypiper = Pypiper(name="WGBS", outfolder=paths.pipeline_outfolder, args=args)

print("N input bams:\t\t" + str(len(args.unmapped_bam)))
print("Sample name:\t\t" + args.sample_name)
ngstk.make_sure_path_exists(paths.pipeline_outfolder + "unmapped_bam/")

if merge:
	input_bams = args.unmapped_bam
	merge_folder = os.path.join(paths.pipeline_outfolder, "unmapped_bam/")
	sample_merged_bam = args.sample_name + ".merged.bam"
	output_merge = os.path.join(merge_folder, sample_merged_bam)
	cmd =  ngstk.merge_bams(input_bams, output_merge, paths)

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
	raise Exception(local_unmapped_bam + " is not a file")


# Fastq conversion
################################################################################
mypiper.timestamp("### Fastq conversion: ")
out_fastq_pre = os.path.join(paths.pipeline_outfolder, "fastq/", args.sample_name)
cmd = ngstk.bam_to_fastq(local_unmapped_bam, out_fastq_pre, args.paired_end, paths)
mypiper.call_lock(cmd, out_fastq_pre + "_R1.fastq")
mypiper.clean_add(out_fastq_pre + "*.fastq", conditional=True)

# Sanity checks:
if not args.no_check:
	raw_reads = ngstk.count_reads(local_unmapped_bam, args.paired_end)
	mypiper.report_result("Raw_reads", str(raw_reads))
	fastq_reads = ngstk.count_reads(out_fastq_pre + "_R1.fastq", paired_end=args.paired_end)
	mypiper.report_result("Fastq_reads", fastq_reads)
	if (fastq_reads != int(raw_reads)):
		raise Exception("Fastq conversion error? Size doesn't match unaligned bam")

# Adapter trimming
################################################################################
mypiper.timestamp("### Adapter trimming: ")

if not args.paired_end:
	cmd = "java -Xmx4000m -jar  " + paths.trimmomatic_jar + " SE -phred33 -threads 30"
	cmd += " -trimlog " + os.path.join(paths.pipeline_outfolder, "fastq/") + "trimlog.log "
	cmd += out_fastq_pre + "_R1.fastq "
	cmd += out_fastq_pre + "_R1_trimmed.fastq "
	cmd += " HEADCROP:6 ILLUMINACLIP:" + paths.adapter_file + ":2:10:4:1:true:epignome:5 SLIDINGWINDOW:4:1 MAXINFO:16:0.40 MINLEN:25"

else:
	cmd = "java -Xmx4000m -jar " + paths.trimmomatic_jar + " PE -phred33 -threads 30"
	cmd += " -trimlog " + os.path.join(paths.pipeline_outfolder, "fastq/") + "trimlog.log "
	cmd += out_fastq_pre + "_R1.fastq "
	cmd += out_fastq_pre + "_R2.fastq "
	cmd += out_fastq_pre + "_R1_trimmed.fastq "
	cmd += out_fastq_pre + "_R1_unpaired.fastq "
	cmd += out_fastq_pre + "_R2_trimmed.fastq "
	cmd += out_fastq_pre + "_R2_unpaired.fastq "
	cmd += " HEADCROP:6 ILLUMINACLIP:" + paths.adapter_file + ":2:10:4:1:true:epignome:5 SLIDINGWINDOW:4:1 MAXINFO:16:0.40 MINLEN:25"

mypiper.call_lock(cmd, out_fastq_pre + "_R1_trimmed.fastq")
trimmed_fastq = out_fastq_pre + "_R1_trimmed.fastq"

if not args.no_check:
	x = ngstk.count_reads(trimmed_fastq,args.paired_end)
	mypiper.report_result("Trimmed_reads", x)

# WGBS alignment with bismark.
################################################################################
mypiper.timestamp("### Bismark alignment: ")
bismark_folder = paths.pipeline_outfolder + "/bismark_" + args.genome_assembly + "/"
ngstk.make_sure_path_exists(bismark_folder)
bismark_temp = bismark_folder + "/" + "bismark_temp"
ngstk.make_sure_path_exists(bismark_temp)

if not args.paired_end:
#	out_bismark_temp = bismark_folder + args.sample_name + "_R1_trimmed.fastq_bismark_bt2.bam"
	out_bismark = bismark_folder + args.sample_name + ".bam"
	cmd = "bismark " + paths.bismark_indexed_genome + " "
	cmd += out_fastq_pre + "_R1_trimmed.fastq"
	cmd += " --bam --unmapped"
	cmd += " --path_to_bowtie " + paths.bowtie2
	cmd += " --bowtie2"
	cmd += " --temp_dir " + bismark_temp
	cmd += " --output_dir " + bismark_folder
	cmd += " -p 6 " # Number of processors
	cmd += " --basename=" +args.sample_name
else:
	#out_bismark = bismark_folder + args.sample_name + "_R1_trimmed.fastq_bismark_bt2_pe.bam"
	out_bismark = bismark_folder + args.sample_name + "_pe.bam"
	cmd = "bismark " + paths.bismark_indexed_genome
	cmd += " --1 " + out_fastq_pre + "_R1_trimmed.fastq"
	cmd += " --2 " + out_fastq_pre + "_R2_trimmed.fastq"
	cmd += " --bam --unmapped"
	cmd += " --path_to_bowtie " + paths.bowtie2
	cmd += " --bowtie2"
	cmd += " --temp_dir " + bismark_temp
	cmd += " --output_dir " + bismark_folder
	cmd += " --minins 0 --maxins 5000"
	cmd += " -p 6 " # Number of processors
	cmd += " --basename="  + args.sample_name
mypiper.call_lock(cmd, out_bismark)

# Rename annoying bismark outfile.
#cmd += "; mv " + out_bismark_temp + " " + out_bismark  # mv annoying bismark output

if not args.no_check:
	x = ngstk.count_reads(out_bismark, args.paired_end)
	mypiper.report_result("Aligned_reads", x)


mypiper.timestamp("### PCR duplicate removal: ")
# Bismark's deduplication forces output naming, how annoying.
out_dedup = bismark_folder + args.sample_name + "_pe.deduplicated.bam"
out_dedup = re.sub(r'.bam$', '.deduplicated.bam', out_bismark)
#out_dedup = bismark_folder + args.sample_name + "_R1_trimmed.fastq_bismark_bt2_pe.deduplicated.bam"

if not args.paired_end:
	cmd = "deduplicate_bismark --single "
	cmd += out_bismark
	cmd += " --bam"
else:
	cmd = "deduplicate_bismark --paired "
	cmd += out_bismark
	cmd += " --bam"

mypiper.call_lock(cmd, out_dedup)

if not args.no_check:
	deduplicated_reads = ngstk.count_reads(out_dedup, args.paired_end)
	mypiper.report_result("Deduplicated_reads", deduplicated_reads)


mypiper.timestamp("### Aligned read filtering: ")

# convert bam file into sam file and sort again to
# compensate for a sorting issue of "deduplicate_bismark"
sam_temp = bismark_folder + "/" + "sam_temp"
ngstk.make_sure_path_exists(sam_temp)
out_sam = bismark_folder + args.sample_name + ".aln.deduplicated.sam"
cmd = "samtools sort -n -o " + out_dedup + " " + out_dedup.replace(".bam", "_sorted") + " | samtools view -h - >" + out_sam

mypiper.call_lock(cmd, out_sam, shell=True)

if not args.no_check:
	#sorted file same size as presorted?
	sorted_reads = ngstk.count_reads(out_sam, paired_end=args.paired_end)
	if sorted_reads != deduplicated_reads:
		raise Exception("Sorted size doesn't match deduplicated size.")

out_sam_filter = bismark_folder + args.sample_name + ".aln.dedup.filt.sam"

headerLines = subprocess.check_output("samtools view -SH " + out_sam + "|wc -l", shell=True).strip()
cmd = "python " + paths.scripts_dir + "/bisulfiteReadFiltering_forRNA.py"
cmd += " --infile=" + out_sam
cmd += " --outfile=" + out_sam_filter
cmd += " --skipHeaderLines=" + headerLines
cmd += " --genome=" + args.genome_assembly
cmd += " --genomeDir=" + paths.ref_genome
cmd += " --minNonCpgSites=3"
cmd += " --minConversionRate=0.9"


if args.paired_end:
	cmd = cmd + " --pairedEnd"

mypiper.call_lock(cmd, out_sam_filter)

if not args.no_check:
	x = ngstk.count_reads(out_sam_filter, args.paired_end)
	mypiper.report_result("Filtered_reads", x)


# Clean up all intermediates
mypiper.clean_add(out_bismark) # initial mapped bam file
mypiper.clean_add(bismark_folder + "/*.fq") # initial unmapped fastqs
mypiper.clean_add(out_dedup) # deduplicated bam file
mypiper.clean_add(out_sam) # dedup conversion to sam
mypiper.clean_add(out_sam_filter) # after filtering


# Methylation extractor
################################################################################
# REMARK NS:
# Bismark methylation extractor produces various outpus, but unfortunately none
# are great. The default "coverage" (.bismark.cov) file is thus:
# chr	start	stop	meth	methylated	unmethylated
# chr17	4890653	4890653	100	1	0
# chr17	5334751	5334751	100	1	0
# This output lacks strand information, so you don't know if the coordinate is
# pointing to a C or G on the + strand unles you look it up in the reference genome.
# The "cytosine_report" file has all the info, but includes an entry for every
# CpG, covered or not:
# chr17	3000204	+	0	0	CG	CGT
# chr17	3000205	-	0	0	CG	CGA
# chr17	4890653	-	1	0	CG	CGA
# Solution: Use the cytosine_report file, and filter out any uncovered reads.

mypiper.timestamp("### Methylation calling (bismark extractor): ")

extract_dir = bismark_folder + "extractor"
ngstk.make_sure_path_exists(extract_dir)
#out_extractor = extract_dir + "/meth-output.txt"
out_extractor = extract_dir + "/" + re.sub(r'.sam$', '.bismark.cov', os.path.basename(out_sam_filter))
out_cpg_report = re.sub(r'.bismark.cov$', '.CpG_report.txt', out_extractor)


if not args.paired_end:
	cmd = "bismark_methylation_extractor --single-end"
	cmd += " --report"
	cmd += " --bedGraph"
	cmd += " --merge_non_CpG"
	cmd += " --cytosine_report"
	cmd += " --genome_folder " + paths.bismark_indexed_genome
	cmd += " --gzip"
	cmd += " --output " + extract_dir
	cmd += " " + out_sam_filter #input file
	#cmd += " > " + extract_dir + "/meth-output.txt" #not necessary since it's in the report.txt
else:
	cmd = "bismark_methylation_extractor --paired-end  --no_overlap"
	cmd += " --report"
	cmd += " --bedGraph"
	cmd += " --merge_non_CpG"
	cmd += " --cytosine_report"
	cmd += " --genome_folder " + paths.bismark_indexed_genome
	cmd += " --gzip"
	cmd += " --output " + extract_dir
	cmd += " " + out_sam_filter #input file
	#cmd += " > " + extract_dir + "/meth-output.txt" #not necessary since it's in the report.txt


mypiper.call_lock(cmd,  out_cpg_report)

# TODO: make these boolean flags options to the pipeline
keep_bismark_report = True
keep_non_standard_chromosomes = False
adjust_minus_strand = True

# prepare outputs:
out_cpg_report_filt = re.sub(r'.CpG_report.txt$', '.CpG_report_filt.txt', out_cpg_report)
out_cpg_report_filt_cov = re.sub(r'.CpG_report.txt$', '.CpG_report_filt.cov', out_cpg_report)

# remove uncovered regions:
cmd = "awk '{ if ($4+$5 > 0) print; }'"
cmd += " " + out_cpg_report
cmd += " > " + out_cpg_report_filt
mypiper.call_lock(cmd,  out_cpg_report_filt, shell=True)

# convert the bismark report to the simpler coverage format and adjust the coordinates
# of CpG's on the reverse strand while doing so (by substracting 1 from the start):
cmd = str(os.path.dirname(os.path.dirname(os.path.realpath(__file__)))) + "/bin/convertBismarkReport.R --formats=cov,min,gibberish --noCovFilter" # disable coverage filter, because we have already used `awk` to achieve this result
if keep_non_standard_chromosomes:
	cmd += " --noChromFilter" 
if not adjust_minus_strand:
	cmd += " --noAdjustMinusStrand" 
cmd += " -i " + out_cpg_report_filt
mypiper.call_lock(cmd,  out_cpg_report_filt_cov)


# tidy up:
if not keep_bismark_report:
	mypiper.clean_add(out_cpg_report_filt) 
	





# Make bigwig
################################################################################
mypiper.timestamp("### Make bigwig: ")

bedGraph = out_extractor.replace(".bismark.cov",".bedGraph")
out_bigwig = bedGraph.replace(".bedGraph", ".bw")
cmd = paths.bed2bigWig + " " + bedGraph + " " + paths.chrom_sizes
cmd += " " + out_bigwig

mypiper.call_lock(cmd, out_bigwig, shell=False)

#Make tracks to view in UCSC genome browser
################################################################################
#Makes bigBed (not necessary for this pipeline)
'''
mypiper.timestamp("### remove type line from bedGrap ")
bedGraph = out_extractor.replace(".bismark.cov",".bedGraph")
cmd = "sed '1d' " + bedGraph + " > " + bedGraph.replace(".bedGraph", ".bed")
mypiper.call_lock(cmd,  bedGraph.replace(".bedGraph", ".bed"),shell=True)
mypiper.timestamp("### make bigBed ")
cmd = paths.bed2bigBed + " " + bedGraph.replace(".bedGraph", ".bed") + " " + paths.chrom_sizes
cmd	+= " " + bedGraph.replace(".bedGraph", ".bb")
mypiper.call_lock(cmd,  bedGraph.replace(".bedGraph", ".bb"),shell=False)
'''

# Spike-in alignment
################################################################################
# currently using bowtie1 instead of bowtie2
mypiper.timestamp("### Bismark spike-in alignment: ")
spikein_folder = paths.pipeline_outfolder + "/bismark_spikein/"
ngstk.make_sure_path_exists(spikein_folder)
spikein_temp = spikein_folder + "/" + "bismark_temp"
ngstk.make_sure_path_exists(spikein_temp)
out_spikein_base = args.sample_name + ".spikein.aln"

#out_spikein = spikein_folder + args.sample_name + "_R1_trimmed.fastq_unmapped_reads_1.fq_bismark_pe.bam"

unmapped_reads_pre = bismark_folder + args.sample_name

if not args.paired_end:
	out_spikein = spikein_folder + out_spikein_base + ".bam"
	cmd = "bismark " + paths.bismark_spikein_genome + " "
	cmd += unmapped_reads_pre + "_unmapped_reads.fq"
	cmd += " --bam --unmapped"
	cmd += " --path_to_bowtie " + paths.bowtie1
#	cmd += " --bowtie2"
	cmd += " --temp_dir " + spikein_temp
	cmd += " --output_dir " + spikein_folder
	cmd += " --basename="  + out_spikein_base
#	cmd += " -p 4"
else:
	out_spikein = spikein_folder + out_spikein_base + "_pe.bam"
	cmd = "bismark " + paths.bismark_spikein_genome
	cmd += " --1 " + unmapped_reads_pre + "_unmapped_reads_1.fq"
	cmd += " --2 " + unmapped_reads_pre + "_unmapped_reads_2.fq"
	#cmd += " --1 " + unmapped_reads_pre + "_R1_trimmed.fastq_unmapped_reads_1.fq"
	#cmd += " --2 " + unmapped_reads_pre + "_R2_trimmed.fastq_unmapped_reads_2.fq"
	cmd += " --bam --unmapped"
	cmd += " --path_to_bowtie " + paths.bowtie1
#	cmd += " --bowtie2"
	cmd += " --temp_dir " + spikein_temp
	cmd += " --output_dir " + spikein_folder
	cmd += " --minins 0 --maxins 5000"
	cmd += " --basename="  + out_spikein_base
#	cmd += " -p 4"

mypiper.call_lock(cmd, out_spikein, nofail=True)
# Clean up the unmapped file which is copied from the parent
# bismark folder to here:
mypiper.clean_add(spikein_folder + "/*.fq", conditional=False)
mypiper.clean_add(spikein_temp) # For some reason, the temp folder is not deleted.


mypiper.timestamp("### PCR duplicate removal (Spike-in): ")
# Bismark's deduplication forces output naming, how annoying.
#out_spikein_dedup = spikein_folder + args.sample_name + ".spikein.aln.deduplicated.bam"
out_spikein_dedup = re.sub(r'.bam$', '.deduplicated.bam', out_spikein)
if not args.paired_end:
	cmd = "deduplicate_bismark --single "
	cmd += out_spikein
	cmd += " --bam"
else:
	cmd = "deduplicate_bismark --paired "
	cmd += out_spikein
	cmd += " --bam"

out_spikein_sorted = out_spikein_dedup.replace('.deduplicated.bam', '.deduplicated.sorted')
cmd2 = "samtools sort " + out_spikein_dedup + " " + out_spikein_sorted
cmd3 = "samtools index " + out_spikein_sorted + ".bam"
cmd4 = "rm " + out_spikein_dedup
mypiper.call_lock([cmd, cmd2, cmd3, cmd4], out_spikein_sorted +".bam.bai", nofail=True)

# Spike-in methylation calling
################################################################################
mypiper.timestamp("### Methylation calling (testxmz) Spike-in: ")

cmd1 = "python " + paths.scripts_dir + "/testxmz.py"
cmd1 += " " + out_spikein_sorted + ".bam" + " " + "K1_unmethylated"
cmd1 += " >> " + mypiper.pipeline_stats_file
cmd2 = cmd1.replace("K1_unmethylated", "K3_methylated")
mypiper.callprint(cmd1, shell=True)
mypiper.callprint(cmd2, shell=True)

# Methylation extractor
################################################################################
# Not necessary for spikein.... because the above does it.
'''
# choose your poison:
extract_from = out_spikein_dedup
#extract_from = out_spikein

extract_dir = spikein_folder + "extractor"
ngstk.make_sure_path_exists(extract_dir)
#out_extractor = extract_dir + args.sample_name + "/meth-output.txt"
out_extractor = extract_dir + "/" + re.sub(r'.bam$', '.bismark.cov', os.path.basename(extract_from))

if not args.paired_end:
	cmd = "bismark_methylation_extractor --single-end"
	cmd += " --report"
	cmd += " --bedGraph"
	cmd += " --merge_non_CpG"
	cmd += " --cytosine_report"
	cmd += " --genome_folder " + paths.bismark_spikein_genome
	cmd += " --gzip"
	cmd += " --output " + extract_dir
	cmd += " " + extract_from #input file
	#cmd += " > " + extract_dir + "/meth-output.txt" #not necessary since it's in the report.txt
else:
	cmd = "bismark_methylation_extractor --paired-end  --no_overlap"
	cmd += " --report"
	cmd += " --bedGraph"
	cmd += " --merge_non_CpG"
	cmd += " --cytosine_report"
	cmd += " --genome_folder " + paths.bismark_spikein_genome
	cmd += " --gzip"
	cmd += " --output " + extract_dir
	cmd += " " + extract_from #input file
	#cmd += " > " + extract_dir + "/meth-output.txt" #not necessary since it's in the report.txt

mypiper.call_lock(cmd, out_extractor)
'''

# Final sorting and indexing
################################################################################
# create sorted and indexed BAM files for visualization and analysis
mypiper.timestamp("### Final sorting and indexing: ")


#out_header = bismark_folder + args.sample_name + ".reheader.bam"
out_final = bismark_folder + args.sample_name + ".final.bam"

#cmd = "java -jar " + paths.picard_dir + "/ReplaceSamHeader.jar"
#cmd += " I=" + out_sam_filter
#cmd += " HEADER=" + out_dedup
#cmd += " O=" + out_header

#mypiper.call_lock(cmd, =None)


# Sort
cmd = "java -Xmx4000m -jar" 
# This sort can run out of temp space on big jobs; this puts the temp to a 
# local spot.
cmd += " -Djava.io.tmpdir=`pwd`/tmp"
cmd += " " + paths.picard_dir + "/SortSam.jar"
cmd +=" I=" + out_sam_filter
cmd +=" O=" + out_final
cmd +=" SORT_ORDER=coordinate VALIDATION_STRINGENCY=SILENT"
cmd +=" CREATE_INDEX=true"
mypiper.call_lock(['mkdir -p tmp', cmd], out_final, lock_name="final_sorting")

# Cleanup
################################################################################
# remove temporary folders
mypiper.clean_add(bismark_temp)
mypiper.clean_add(sam_temp)
mypiper.stop_pipeline()
