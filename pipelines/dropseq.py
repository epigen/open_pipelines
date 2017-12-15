#!/usr/bin/env python

"""
Drop-seq pipeline
"""

import os
import sys
from argparse import ArgumentParser

import pypiper
import yaml
from pep import AttributeDict


__author__ = "Andre Rendeiro"
__email__ = "arendeiro@cemm.oeaw.ac.at"
__credits__ = []
__license__ = "GPL2"
__version__ = "0.2"
__status__ = "Development"


def main():
	# Parse command-line arguments
	parser = ArgumentParser(
		prog="dropseq-pipeline",
		description="Drop-seq pipeline."
	)
	parser = arg_parser(parser)
	parser = pypiper.add_pypiper_args(parser, all_args=True)
	args = parser.parse_args()
	if args.sample_config is None:
		parser.print_help()
		return 1

	# Read in yaml configs
	sample = AttributeDict(yaml.load(open(args.sample_config, "r")))
	pipeline_config = AttributeDict(yaml.load(open(os.path.join(os.path.dirname(__file__), args.config_file), "r")))

	# Start main function
	process(sample, pipeline_config, args)


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
	parser.add_argument(
		"--no-stats",
		dest="debug",
		help="Do not output statistics. If not selected, statistics using pypiper report_result() function"
		" will be output in 'stats.tsv'",
		action="store_false"
	)
	return parser


def merge_bam_files(input_bams, output_bam, args, pipe, tmpdir):
	cmd = 'java -Djava.io.tmpdir={}'.format(tmpdir)
	cmd += " -Xmx{}g".format(int(args.mem) / 1000)
	cmd += " -jar " + pipe.config.tools.piccard_jar + " MergeSamFiles"
	cmd += " USE_THREADING=TRUE"
	cmd += " SORT_ORDER=queryname"
	cmd += " " + (" ".join(["INPUT=%s"] * len(input_bams))) % tuple(input_bams)
	cmd += " OUTPUT={0}".format(output_bam)
	return cmd


def report_flagstat(pipe, bam_file, prefix=""):
	"""
	Report flagstat statistics on a bam file.
	"""
	from collections import OrderedDict
	import re

	def flagstat(bam_file):
		pipe.run("sambamba flagstat {} > {}".format(bam_file, bam_file + ".flagstat"), bam_file + ".flagstat", shell=True)
		return bam_file + ".flagstat"

	def parse_flagstat(flagstat_file):
		stats = OrderedDict()
		with open(flagstat_file, "r") as handle:
			lines = handle.readlines()
		for line in lines:
			key = re.sub(r" \(.*\)", "", " ".join(line.strip().split(" ")[3:]))
			value = line.strip().split(" ")[0]
			stats[prefix + " " + key] = value
		return stats

	for key, value in parse_flagstat(flagstat(bam_file)).items():
		pipe.report_result(key, value)


def report_star_log(pipe, star_log, prefix=""):
	"""
	Report statistics on the STAR alignment step of Drop-seq.
	"""
	from collections import OrderedDict

	def parse_star_log(star_log):
		stats = OrderedDict()
		with open(star_log, "r") as handle:
			lines = handle.readlines()
		for line in lines:
			if len(line.split("\t")) == 2:
				key = line.strip().split("\t")[0].replace(" |", "")
				value = line.strip().split("\t")[1].replace("%", "")
				stats[prefix + " " + key] = value
		return stats

	for key, value in parse_star_log(star_log).items():
		pipe.report_result(key, value)


def report_bead_synthesis(pipe, synthesis_log, prefix=""):
	"""
	Report statistics on potential errors in bead synthesis.
	"""
	def parse_bead_synthesis_log(synthesis_log):
		import pandas as pd
		log = pd.read_csv(synthesis_log, sep="\t", comment="#").ix[0].astype(float)
		log["synthesis_error %"] = (1 - (log['NO_ERROR'] / log['NUM_BEADS'])) * 100
		log["fixed_first_base %"] = (log["FIXED_FIRST_BASE"] / log["NUM_BEADS"]) * 100.
		log["other_error_count %"] = (log["OTHER_ERROR_COUNT"] / log["NUM_BEADS"]) * 100.
		log["single_umi_error %"] = (log["SINGLE_UMI_ERROR"] / log["NUM_BEADS"]) * 100.
		log["primer_match %"] = (log["PRIMER_MATCH"] / log["NUM_BEADS"]) * 100.
		return log.to_dict()

	for key, value in parse_bead_synthesis_log(synthesis_log).items():
		pipe.report_result(key, value)


def report_digital_expression(pipe, expression_file, prefix=""):
	"""
	Report statistics on the digital expression output from Drop-seq across several thresholds
	(minimal expression value to call a gene "expressed", number of genes per cell and number of cells per gene).
	"""
	import pandas as pd
	from collections import OrderedDict

	def parse_digital_expression(expression_file):

		exp = pd.read_csv(expression_file, sep="\t").set_index("GENE")
		stats = OrderedDict()
		stats[prefix + " number_genes"] = exp.shape[0]
		stats[prefix + " number_cells"] = exp.shape[1]

		# reads per cell
		reads_per_cell = exp.sum(axis=0)
		stats[prefix + " total_used_reads"] = reads_per_cell.sum()
		# % mitochondrial
		mito_per_cell = (exp.ix[exp.index[exp.index.str.contains("^MT-")]].sum(axis=0) / reads_per_cell) * 100
		# % ribosomal proteins
		ribo_per_cell = (exp.ix[exp.index[exp.index.str.contains("^RP")]].sum(axis=0) / reads_per_cell) * 100

		# for each, add to stats table distribution values (mean, 5th, 25th, 50th, 75th, 95th percentiles)
		for metric in ["reads_per_cell", "mito_per_cell", "ribo_per_cell"]:
			stats[prefix + " {}:mean".format(metric)] = eval(metric).mean()
			stats[prefix + " {}:median".format(metric)] = eval(metric).median()
			stats[prefix + " {}:std".format(metric)] = eval(metric).std()

		for reads_to_coverage in [1, 2, 3]:  # , 4, 5
			cells_per_gene = exp.apply(lambda x: (x >= reads_to_coverage).sum(), axis=1)
			genes_per_cell = exp.apply(lambda x: (x >= reads_to_coverage).sum(), axis=0)

			stats[prefix + " {}reads_to_coverage_:genes_per_cell:mean".format(reads_to_coverage)] = genes_per_cell.mean()
			stats[prefix + " {}reads_to_coverage_:genes_per_cell:median".format(reads_to_coverage)] = genes_per_cell.median()
			stats[prefix + " {}reads_to_coverage_:genes_per_cell:std".format(reads_to_coverage)] = genes_per_cell.std()

			for c in [0.1, 0.5, 1, 2, 5, 10]:  # percentage of cells in which gene is covered
				for g in [100, 500]:  # how many genes does a cell has (dependent on the reads_to_coverage parameter)
					matrix = exp[cells_per_gene > (exp.shape[1] / (100 / float(c)))][exp.columns[genes_per_cell >= g]]
					stats[prefix + " {}reads_to_coverage_{}%cellswithgene_{}genespercell_:genes".format(reads_to_coverage, c, g)] = matrix.shape[0]
					stats[prefix + " {}reads_to_coverage_{}%cellswithgene_{}genespercell_:cells".format(reads_to_coverage, c, g)] = matrix.shape[1]
		return stats

	for key, value in parse_digital_expression(expression_file).items():
		pipe.report_result(key, value)


def report_umi_stats(pipe, umi_file, prefix=""):
	"""
	Report statistics on the digital expression output from Drop-seq across several thresholds
	(minimal expression value to call a gene "expressed", number of genes per cell and number of cells per gene).
	"""
	import pandas as pd
	from collections import OrderedDict

	def parse_digital_expression(umi_file):
		# UMI duplication stats
		umi = pd.read_csv(umi_file, sep="\t")
		stats = OrderedDict()
		# percent unique UMIs and UMI duplication
		perc_uni_umi = ((umi["Num_Obs"] == 1).sum() / float(umi.shape[0])) * 100
		stats[" percent_unique_umis"] = perc_uni_umi

		return stats

	for key, value in parse_digital_expression(umi_file).items():
		pipe.report_result(key, value)


def process(sample, pipeline_config, args):
	"""
	"""
	print("Start processing Drop-seq sample %s." % sample.sample_name)

	for path in ["sample_root"] + sample.paths.__dict__.keys():
		try:
			exists = os.path.exists(sample.paths[path])
		except TypeError:
			continue
		if not exists:
			try:
				os.mkdir(sample.paths[path])
			except OSError("Cannot create '%s' path: %s" % (path, sample.paths[path])):
				raise

	# Start Pypiper object
	pipe = pypiper.PipelineManager("dropseq", sample.paths.sample_root, args=args)

	# Set up a few handy shorthand variables
	dropseq_root = pipe.config.tools.dropseq_tools_root
	output_dir = sample.paths.sample_root

	# Merge Bam files if more than one technical replicate
	if len(sample.data_source.split(" ")) > 1:
		pipe.timestamp("## Merging bam files from replicates")
		cmd = merge_bam_files(
			input_bams=sample.data_source.split(" "),  # this is a list of sample paths
			output_bam=os.path.join(output_dir, "unaligned_merged.bam"),
			args=args, pipe=pipe,
			tmpdir=output_dir
		)
		pipe.run(cmd, os.path.join(output_dir, "unaligned_merged.bam"))
		pipe.clean_add(os.path.join(output_dir, "unaligned_merged.bam"), manual=True)

		input_file = os.path.join(output_dir, "unaligned_merged.bam")
	else:
		input_file = sample.data_source

	# Copy the input file if it is not writable
	# (the first step requires the file to be writable which is silly)
	if not os.access(input_file, os.W_OK):
		pipe.timestamp("## Copying input file to output directory")
		cmd = "cp {} {}".format(input_file, os.path.join(output_dir, "input_file.bam"))
		pipe.run(cmd, os.path.join(output_dir, "input_file.bam"))
		cmd = "chmod 664 {}".format(os.path.join(output_dir, "input_file.bam"))
		pipe.run(cmd, os.path.join(output_dir, "input_file.bam_chmod"))
		pipe.clean_add(os.path.join(output_dir, "input_file.bam"), manual=False)
		input_file = os.path.join(output_dir, "input_file.bam")

	os.environ['TMP_DIR'] = output_dir

	if args.debug:
		report_flagstat(pipe, os.path.join(output_dir, input_file), prefix="input_file")

	# Stage 1: pre-alignment tag and trim
	# Tag with cell barcode
	pipe.timestamp("## Tagging BAM file with cell barcode")
	cmd = os.path.join(dropseq_root, "TagBamWithReadSequenceExtended")
	cmd += " TMP_DIR=" + output_dir
	cmd += " SUMMARY=" + os.path.join(output_dir, "unaligned_tagged_Cellular.bam_summary.txt")
	cmd += " BASE_RANGE={}".format(pipe.config.parameters.cell_barcode_bases)
	cmd += " BASE_QUALITY={}".format(pipe.config.parameters.min_base_quality)
	cmd += " BARCODED_READ=1 DISCARD_READ=false TAG_NAME=XC NUM_BASES_BELOW_QUALITY={}".format(pipe.config.parameters.min_bases_below_quality)
	cmd += " INPUT=" + input_file
	cmd += " OUTPUT=" + os.path.join(output_dir, "unaligned_tagged_Cell.bam")
	pipe.run(cmd, os.path.join(output_dir, "unaligned_tagged_Cell.bam"))
	pipe.clean_add(os.path.join(output_dir, "unaligned_tagged_Cell.bam"), manual=True)

	# Tag with molecule barcode
	pipe.timestamp("## Tagging BAM file with molecule barcode (UMI)")
	cmd = os.path.join(dropseq_root, "TagBamWithReadSequenceExtended")
	cmd += " TMP_DIR=" + output_dir
	cmd += " SUMMARY=" + os.path.join(output_dir, "unaligned_tagged_Molecular.bam_summary.txt")
	cmd += " BASE_RANGE={}".format(pipe.config.parameters.umi_barcode_bases)
	cmd += " BASE_QUALITY={}".format(pipe.config.parameters.min_base_quality)
	cmd += " BARCODED_READ=1 DISCARD_READ=true TAG_NAME=XM NUM_BASES_BELOW_QUALITY={}".format(pipe.config.parameters.min_bases_below_quality)
	cmd += " INPUT=" + os.path.join(output_dir, "unaligned_tagged_Cell.bam")
	cmd += " OUTPUT=" + os.path.join(output_dir, "unaligned_tagged_CellMolecular.bam")
	pipe.run(cmd, os.path.join(output_dir, "unaligned_tagged_CellMolecular.bam"))
	pipe.clean_add(os.path.join(output_dir, "unaligned_tagged_CellMolecular.bam"), manual=True)

	# Filter bam
	pipe.timestamp("## Filtering BAM file")
	cmd = os.path.join(dropseq_root, "FilterBAM")
	cmd += " TAG_REJECT=XQ"
	cmd += " INPUT=" + os.path.join(output_dir, "unaligned_tagged_CellMolecular.bam")
	cmd += " OUTPUT=" + os.path.join(output_dir, "unaligned_tagged_filtered.bam")
	pipe.run(cmd, os.path.join(output_dir, "unaligned_tagged_filtered.bam"))
	pipe.clean_add(os.path.join(output_dir, "unaligned_tagged_filtered.bam"), manual=True)

	if args.debug:
		report_flagstat(pipe, os.path.join(output_dir, "unaligned_tagged_filtered.bam"), prefix="FilterBAM")

	# Trim starting sequence
	pipe.timestamp("## Triming starting sequence")
	cmd = os.path.join(dropseq_root, "TrimStartingSequence")
	cmd += " SEQUENCE={}".format(pipe.config.parameters.trim_sequence)
	cmd += " MISMATCHES=0 NUM_BASES={}".format(pipe.config.parameters.trim_sequence_length)
	cmd += " OUTPUT_SUMMARY=" + os.path.join(output_dir, "adapter_trimming_report.txt")
	cmd += " INPUT=" + os.path.join(output_dir, "unaligned_tagged_filtered.bam")
	cmd += " OUTPUT=" + os.path.join(output_dir, "unaligned_tagged_trimmed_smart.bam")
	pipe.run(cmd, os.path.join(output_dir, "unaligned_tagged_trimmed_smart.bam"))
	pipe.clean_add(os.path.join(output_dir, "unaligned_tagged_trimmed_smart.bam"), manual=True)

	if args.debug:
		report_flagstat(pipe, os.path.join(output_dir, "unaligned_tagged_trimmed_smart.bam"), prefix="TrimStartingSequence")

	# Trim polyA tail
	pipe.timestamp("## Trimming polyA tail")
	cmd = os.path.join(dropseq_root, "PolyATrimmer")
	cmd += " MISMATCHES=0 NUM_BASES={}".format(pipe.config.parameters.polya_size)
	cmd += " OUTPUT_SUMMARY=" + os.path.join(output_dir, "polyA_trimming_report.txt")
	cmd += " INPUT=" + os.path.join(output_dir, "unaligned_tagged_trimmed_smart.bam")
	cmd += " OUTPUT=" + os.path.join(output_dir, "unaligned_mc_tagged_polyA_filtered.bam")
	pipe.run(cmd, os.path.join(output_dir, "unaligned_mc_tagged_polyA_filtered.bam"))
	pipe.clean_add(os.path.join(output_dir, "unaligned_mc_tagged_polyA_filtered.bam"), manual=True)

	if args.debug:
		report_flagstat(pipe, os.path.join(output_dir, "unaligned_mc_tagged_polyA_filtered.bam"), prefix="PolyATrimmer")

	# Stage 2: alignment
	# Convert to fastq
	pipe.timestamp("## Converting to Fastq")
	cmd = "java -Xmx{}g -jar {} SamToFastq".format(int(args.mem) / 1000, pipe.config.tools.piccard_jar)
	cmd += " INPUT=" + os.path.join(output_dir, "unaligned_mc_tagged_polyA_filtered.bam")
	cmd += " FASTQ=" + os.path.join(output_dir, "unaligned_mc_tagged_polyA_filtered.fastq")
	pipe.run(cmd, os.path.join(output_dir, "unaligned_mc_tagged_polyA_filtered.fastq"))
	pipe.clean_add(os.path.join(output_dir, "unaligned_mc_tagged_polyA_filtered.fastq"), manual=True)

	# Align reads
	pipe.timestamp("## Aligning reads with STAR")
	cmd = pipe.config.tools.star
	cmd += " --genomeDir {}".format(getattr(pipe.config.resources.star_index, sample.genome))
	cmd += " --runThreadN {}".format(args.cores)
	cmd += " --outFileNamePrefix " + os.path.join(output_dir, "star.")
	cmd += " --readFilesIn " + os.path.join(output_dir, "unaligned_mc_tagged_polyA_filtered.fastq")
	pipe.run(cmd, os.path.join(output_dir, "star.Aligned.out.sam"))
	pipe.clean_add(os.path.join(output_dir, "star.Aligned.out.sam"), manual=True)

	if args.debug:
		report_star_log(pipe, os.path.join(output_dir, "star.Log.final.out"), prefix="STAR")

	# Stage 3: sort aligned reads (STAR does not necessarily emit reads in the same order as the input)
	pipe.timestamp("## Sorting aligned BAM file")
	cmd = "java -Dsamjdk.buffer_size=131072 -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Xmx{}g".format(int(args.mem) / 1000)
	cmd += " -jar {} SortSam".format(pipe.config.tools.piccard_jar)
	cmd += " INPUT=" + os.path.join(output_dir, "star.Aligned.out.sam")
	cmd += " OUTPUT=" + os.path.join(output_dir, "aligned.sorted.bam")
	cmd += " SORT_ORDER=queryname"
	cmd += " TMP_DIR=" + output_dir
	pipe.run(cmd, os.path.join(output_dir, "aligned.sorted.bam"))
	pipe.clean_add(os.path.join(output_dir, "aligned.sorted.bam"), manual=True)

	# Stage 4: merge and tag aligned reads
	# Merge
	pipe.timestamp("## Merging aligned with unaligned reads")
	cmd = "java -Djava.io.tmpdir={} -Xmx{}g -jar {} MergeBamAlignment".format(output_dir, int(args.mem) / 1000, pipe.config.tools.piccard_jar)
	cmd += " REFERENCE_SEQUENCE={}".format(getattr(pipe.config.resources.genome, sample.genome))
	cmd += " UNMAPPED_BAM=" + os.path.join(output_dir, "unaligned_mc_tagged_polyA_filtered.bam")
	cmd += " ALIGNED_BAM=" + os.path.join(output_dir, "aligned.sorted.bam")
	cmd += " INCLUDE_SECONDARY_ALIGNMENTS=false"
	cmd += " ALIGNED_READS_ONLY=false"
	cmd += " PAIRED_RUN=false"
	cmd += " OUTPUT=" + os.path.join(output_dir, "merged.bam")
	pipe.run(cmd, os.path.join(output_dir, "merged.bam"))
	pipe.clean_add(os.path.join(output_dir, "merged.bam"), manual=True)

	if args.debug:
		report_flagstat(pipe, os.path.join(output_dir, "merged.bam"), prefix="MergeBamAlignment")

	# Tag reads with exon
	pipe.timestamp("## Tagging reads with exon")
	cmd = os.path.join(dropseq_root, "TagReadWithGeneExon")
	cmd += " OUTPUT=" + os.path.join(output_dir, "star_gene_exon_tagged.bam")
	cmd += " ANNOTATIONS_FILE={}".format(getattr(pipe.config.resources.refflat, sample.genome))
	cmd += " TAG=GE CREATE_INDEX=true"
	cmd += " INPUT=" + os.path.join(output_dir, "merged.bam")
	pipe.run(cmd, os.path.join(output_dir, "star_gene_exon_tagged.bam"))

	if args.debug:
		report_flagstat(pipe, os.path.join(output_dir, "star_gene_exon_tagged.bam"), prefix="TagReadWithGeneExon")

	# QC time!

	if pipe.config.parameters.repair_barcodes:
		# Detect and fix bead synthesis errors
		pipe.timestamp("## Reporting and fixing bead synthesis errors")
		cmd = os.path.join(dropseq_root, "DetectBeadSynthesisErrors")
		cmd += " INPUT=" + os.path.join(output_dir, "star_gene_exon_tagged.bam")
		cmd += " OUTPUT=" + os.path.join(output_dir, "star_gene_exon_tagged.clean.bam")
		cmd += " OUTPUT_STATS=" + os.path.join(output_dir, "synthesis_statistics.txt")
		cmd += " SUMMARY=" + os.path.join(output_dir, "synthesis_statistics.summary.txt")
		cmd += " NUM_BARCODES={}".format(pipe.config.parameters.number_seq_error_barcodes_check)
		cmd += " PRIMER_SEQUENCE={}".format(pipe.config.parameters.bead_primer_sequence)
		cmd += " EDIT_DISTANCE={}".format(pipe.config.parameters.distance_to_bead_primer_seq)
		cmd += " MAX_NUM_ERRORS={}".format(pipe.config.parameters.max_number_barcode_bases_to_repair)
		cmd += " TMP_DIR=" + output_dir
		pipe.run(cmd, os.path.join(output_dir, "star_gene_exon_tagged.clean.bam"))

		if args.debug:
			report_bead_synthesis(pipe, os.path.join(output_dir, "synthesis_statistics.summary.txt"), prefix="DetectBeadSynthesisErrors")

		bam_file = os.path.join(output_dir, "star_gene_exon_tagged.clean.bam")
	else:
		bam_file = os.path.join(output_dir, "star_gene_exon_tagged.bam")

	# Distribution of read quality
	# cell barcode
	pipe.timestamp("## Read quality in cell barcodes")
	cmd = os.path.join(dropseq_root, "GatherReadQualityMetrics")
	cmd += " INPUT=" + bam_file
	cmd += " OUTPUT=" + os.path.join(output_dir, "quality_distribution.cell_barcode.txt")
	cmd += " TAG=XC"
	pipe.run(cmd, os.path.join(output_dir, "quality_distribution.cell_barcode.txt"))
	# UMI
	pipe.timestamp("## Read quality in molecule barcodes")
	cmd = os.path.join(dropseq_root, "GatherReadQualityMetrics")
	cmd += " INPUT=" + bam_file
	cmd += " OUTPUT=" + os.path.join(output_dir, "quality_distribution.mol_barcode.txt")
	cmd += " TAG=XM"
	pipe.run(cmd, os.path.join(output_dir, "quality_distribution.mol_barcode.txt"))

	# Distribution of bases in reads
	# cell barcode
	pipe.timestamp("## Distribution of bases in cell barcodes")
	cmd = os.path.join(dropseq_root, "BaseDistributionAtReadPosition")
	cmd += " INPUT=" + bam_file
	cmd += " OUTPUT=" + os.path.join(output_dir, "base_distribution.cell_barcode.txt")
	cmd += " TAG=XC"
	pipe.run(cmd, os.path.join(output_dir, "base_distribution.cell_barcode.txt"))
	# UMI
	pipe.timestamp("## Distribution of bases in molecule barcodes (UMI)")
	cmd = os.path.join(dropseq_root, "BaseDistributionAtReadPosition")
	cmd += " INPUT=" + bam_file
	cmd += " OUTPUT=" + os.path.join(output_dir, "base_distribution.mol_barcode.txt")
	cmd += " TAG=XM"
	pipe.run(cmd, os.path.join(output_dir, "base_distribution.mol_barcode.txt"))

	# Expression time!

	# Reads per cell summary
	pipe.timestamp("## Reporting summary of reads per cell")
	cmd = os.path.join(dropseq_root, "BAMTagHistogram")
	cmd += " INPUT=" + bam_file
	cmd += " OUTPUT=" + os.path.join(output_dir, "cell_readcounts.txt")
	cmd += " FILTER_PCR_DUPLICATES=true"
	cmd += " TAG=XC"
	pipe.run(cmd, os.path.join(output_dir, "cell_readcounts.txt"))

	if args.debug:
		report_flagstat(pipe, bam_file, prefix="BAMTagHistogram")

	# Perform digital gene expression analysis selecting all cells that have at least minGenes genes covered
	for n_genes in pipe.config.parameters.min_genes_per_cell:

		pipe.timestamp("## Perform digital gene expression analysis for cells with at least {} genes covered".format(n_genes))
		cmd = os.path.join(dropseq_root, "DigitalExpression")
		cmd += " -m {}g".format(int(args.mem) / 1000)
		cmd += " TMP_DIR=" + output_dir
		cmd += " INPUT=" + bam_file
		cmd += " OUTPUT=" + os.path.join(output_dir, "digital_expression.{}genes.tsv".format(n_genes))
		cmd += " SUMMARY=" + os.path.join(output_dir, "digital_expression.summary.{}genes.tsv".format(n_genes))
		cmd += " MIN_NUM_GENES_PER_CELL={}".format(n_genes)
		pipe.run(cmd, os.path.join(output_dir, "digital_expression.{}genes.tsv".format(n_genes)), nofail=True)

		if args.debug:
			if os.path.exists(os.path.join(output_dir, "digital_expression.{}genes.tsv".format(n_genes))):
				try:
					print("Reporting digital expression for cells with at least {} genes covered".format(n_genes))
					report_digital_expression(pipe, os.path.join(output_dir, "digital_expression.{}genes.tsv".format(n_genes)), prefix="DigitalExpression_{}genes".format(n_genes))
				except IOError:
					print("Digital expression for cells with at least {} genes covered could not be open.".format(n_genes))

	# Report how often the same UMI is found per cell per gene --> estimate of PCR duplicates
	for n_genes in pipe.config.parameters.min_genes_per_cell:
		pipe.timestamp("## Report UMI count per cell per gene for cells with at least {} genes covered".format(n_genes))
		cmd = os.path.join(dropseq_root, "GatherMolecularBarcodeDistributionByGene")
		cmd += " -m {}g".format(int(args.mem) / 1000)
		cmd += " TMP_DIR=" + output_dir
		cmd += " INPUT=" + bam_file
		cmd += " OUTPUT=" + os.path.join(output_dir, "cell_umi_barcodes.{}genes.tsv".format(n_genes))
		cmd += " MIN_NUM_GENES_PER_CELL={}".format(n_genes)
		pipe.run(cmd, os.path.join(output_dir, "cell_umi_barcodes.{}genes.tsv".format(n_genes)))

	print("Finished processing sample %s." % sample.sample_name)
	pipe.stop_pipeline()


if __name__ == '__main__':
	try:
		sys.exit(main())
	except KeyboardInterrupt:
		print("Program canceled by user!")
		sys.exit(1)
