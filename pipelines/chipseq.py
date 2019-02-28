#!/usr/bin/env python

"""
ChIP-seq pipeline
"""

from argparse import ArgumentParser
import os
import re
import sys

import pandas as pd
import yaml

import pypiper
from pypiper.ngstk import NGSTk
from looper.models import AttributeDict, Sample
from const import CHIP_COMPARE_COLUMN, CHIP_MARK_COLUMN


__author__ = "Andre Rendeiro"
__copyright__ = "Copyright 2015, Andre Rendeiro"
__credits__ = []
__license__ = "GPL2"
__version__ = "0.3"
__maintainer__ = "Andre Rendeiro"
__email__ = "arendeiro@cemm.oeaw.ac.at"
__status__ = "Development"



BROAD_MARKS = {
		"H3K9ME1", "H3K9ME2", "H3K9ME3",
		"H3K27ME1", "H3K27ME2", "H3K27ME3",
		"H3K36ME1", "H3K36ME2", "H3K36ME3",
		"H3K72ME1", "H3K72ME2", "H3K72ME3"}
HISTONE_CODES = ["H3", "H2A", "H2B", "H4"]



class ChIPseqSample(Sample):
	"""
	Class to model ChIP-seq samples based on the generic Sample class (itself a pandas.Series).

	:param series: Pandas `Series` object.
	:type series: pandas.Series

	:Example:

	# create Samples through a project object
	from looper.models import Project
	prj = Project("project_config.yaml")
	prj.add_sample_sheet()
	s0 = prj.samples[0]  # here's a Sample

	# create Samples through a SampleSheet object
	from looper.models import SampleSheet, Sample
	sheet = SampleSheet("project_sheet.csv")
	s1 = Sample(sheet.ix[0])  # here's a Sample too
	"""
	__library__ = "ChIP-seq"

	def __init__(self, series):
		super(ChIPseqSample, self).__init__(series)
		self.tagmented = False

		# Set broad/histone status that may later be modified given
		# context of a pipeline configuration file, handling null/missing mark.
		mark = getattr(self, CHIP_MARK_COLUMN, None)
		if mark is None:
			self.broad = False
			self.histone = False
		else:
			mark = mark.upper()
			self.broad = mark in BROAD_MARKS
			self.histone = any([mark.startswith(histone_code)
								for histone_code in HISTONE_CODES])

	def __repr__(self):
		return "ChIP-seq sample '%s'" % self.sample_name

	def set_file_paths(self, project):
		"""
		Sets the paths of all files for this sample.
		"""

		# Get paths container structure and any contents used by any Sample.
		super(ChIPseqSample, self).set_file_paths(project)

		# Files in the root of the sample dir
		self.fastqc = os.path.join(self.paths.sample_root, self.sample_name + ".fastqc.zip")
		self.trimlog = os.path.join(self.paths.sample_root, self.sample_name + ".trimlog.txt")
		self.aln_rates = os.path.join(self.paths.sample_root, self.sample_name + ".aln_rates.txt")
		self.aln_metrics = os.path.join(self.paths.sample_root, self.sample_name + ".aln_metrics.txt")
		self.dups_metrics = os.path.join(self.paths.sample_root, self.sample_name + ".dups_metrics.txt")

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

		# Mapped: mapped, duplicates marked, removed, reads shifted
		self.paths.mapped = os.path.join(self.paths.sample_root, "mapped")
		self.mapped = os.path.join(self.paths.mapped, self.sample_name + ".trimmed.bowtie2.bam")
		self.filtered = os.path.join(self.paths.mapped, self.sample_name + ".trimmed.bowtie2.filtered.bam")

		# Files in the root of the sample dir
		self.frip = os.path.join(self.paths.sample_root, self.sample_name + "_FRiP.txt")

		# Coverage: read coverage in windows genome-wide
		self.paths.coverage = os.path.join(self.paths.sample_root, "coverage")
		self.coverage = os.path.join(self.paths.coverage, self.sample_name + ".cov")

		self.insertplot = os.path.join(self.paths.sample_root, self.name + "_insertLengths.pdf")
		self.insertdata = os.path.join(self.paths.sample_root, self.name + "_insertLengths.csv")

		self.qc = os.path.join(self.paths.sample_root, self.sample_name + "_qc.tsv")
		self.qc_plot = os.path.join(self.paths.sample_root, self.sample_name + "_qc.pdf")

		bigwig_subfolder = "bigwig_{}".format(self.genome)
		bigwig_folder = os.path.join(
				self.prj.metadata.results_subdir, self.name, bigwig_subfolder)
		bigwig_file = "CHIP_{}.bw".format(self.sample_name)
		self.bigwig = os.path.join(bigwig_folder, bigwig_file)

		# Peaks: peaks called and derivate files
		self.paths.peaks = os.path.join(self.paths.sample_root, "peaks")
		self.peaks = os.path.join(self.paths.peaks, self.sample_name + ("_peaks.narrowPeak" if not self.broad else "_peaks.broadPeak"))


class ChIPmentation(ChIPseqSample):
	"""
	Class to model ChIPmentation samples based on the ChIPseqSample class.

	:param series: Pandas `Series` object.
	:type series: pandas.Series
	"""
	__library__ = "ChIPmentation"

	def __init__(self, series):
		super(ChIPmentation, self).__init__(series)
		self.tagmented = True

	def __repr__(self):
		return "ChIPmentation sample '%s'" % self.sample_name

	def set_file_paths(self, project):
		super(ChIPmentation, self).set_file_paths(project)


# TODO: remove and use the pypiper version once it supports normalization factor.
def bam_to_bigwig(input_bam, output_bigwig, genome_sizes, genome, tagmented=False, normalize=False, norm_factor=1000000):
	import os
	# TODO:
	# Adjust fragment length dependent on read size and real fragment size
	# (right now it assumes 50bp reads with 180bp fragments.)
	cmds = list()
	transient_file = os.path.abspath(re.sub("\.bigWig", "", output_bigwig))
	cmd1 = "bedtools bamtobed -i {0} |".format(input_bam)
	if not tagmented:
		cmd1 += " " + "bedtools slop -i stdin -g {0} -s -l 0 -r 130 |".format(genome_sizes)
		bedfile_bounds_script = os.path.join(
			os.path.dirname(__file__), "tools",
			"fix_bedfile_genome_boundaries.py")
		cmd1 += " {0} {1} |".format(bedfile_bounds_script, genome)
	cmd1 += " " + "genomeCoverageBed {0}-bg -g {1} -i stdin > {2}.cov".format(
		"-5 " if tagmented else "",
		genome_sizes,
		transient_file
	)
	cmds.append(cmd1)
	if normalize:
		cmds.append("""awk 'NR==FNR{{sum+= $4; next}}{{ $4 = ($4 / sum) * {1}; print}}' {0}.cov {0}.cov | sort -k1,1 -k2,2n > {0}.normalized.cov""".format(transient_file, norm_factor))
	cmds.append("bedGraphToBigWig {0}{1}.cov {2} {3}".format(transient_file, ".normalized" if normalize else "", genome_sizes, output_bigwig))
	# remove tmp files
	cmds.append("if [[ -s {0}.cov ]]; then rm {0}.cov; fi".format(transient_file))
	if normalize:
		cmds.append("if [[ -s {0}.normalized.cov ]]; then rm {0}.normalized.cov; fi".format(transient_file))
	cmds.append("chmod 755 {0}".format(output_bigwig))
	return cmds


def report_dict(pipe, stats_dict):
	"""
	Convenience wrapper to report a collection of pipeline results.

	This writes a collection of key-value pairs to the central stats/results
	file associated with the pipeline manager provided.

	:param pipe: Pipeline manager with which to do the reporting
	:type pipe: pypiper.PipelineManager
	:param stats_dict: Collection of results, each mapped to a name
	:type stats_dict: Mapping[str, object]
	"""
	for key, value in stats_dict.items():
		pipe.report_result(key, value)


def parse_fastqc(fastqc_zip, prefix=""):
	"""
	Fetch a handful of fastqc metrics from its output file.

	Count of total reads passing filtration, poor quality reads, read length,
	and GC content are the metrics parsed and returned.

	:param fastqc_zip: Path to the zipfile created by fastqc
	:type fastqc_zip: str
	:param prefix: Prefix for name of each metric
	:type prefix: str
	:rtype: Mapping[str, object]
	"""
	import StringIO
	import zipfile

	stat_names = [prefix + stat for stat in
				  ["total_pass_filter_reads", "poor_quality",
				   "read_length", "GC_perc"]]

	def error_dict():
		return dict((sn, pd.np.nan) for sn in stat_names)

	# Read the zipfile from fastqc.
	try:
		zfile = zipfile.ZipFile(fastqc_zip)
		content = StringIO.StringIO(zfile.read(os.path.join(zfile.filelist[0].filename, "fastqc_data.txt"))).readlines()
	except:
		return error_dict()

	# Parse the unzipped fastqc data, fetching the desired metrics.
	try:
		line = [i for i in range(len(content)) if "Total Sequences" in content[i]][0]
		total = int(re.sub("\D", "", re.sub("\(.*", "", content[line])))
		line = [i for i in range(len(content)) if "Sequences flagged as poor quality" in content[i]][0]
		poor_quality = int(re.sub("\D", "", re.sub("\(.*", "", content[line])))
		line = [i for i in range(len(content)) if "Sequence length	" in content[i]][0]
		seq_len = int(re.sub("\D", "", re.sub(" \(.*", "", content[line]).strip()))
		line = [i for i in range(len(content)) if "%GC" in content[i]][0]
		gc_perc = int(re.sub("\D", "", re.sub(" \(.*", "", content[line]).strip()))
		values = [total, 100 * float(poor_quality) / total, seq_len, gc_perc]
		return dict(zip(stat_names, values))

	except IndexError:
		return error_dict()


def parse_trim_stats(stats_file, prefix="", paired_end=True):
	"""
	:param stats_file: sambamba output file with duplicate statistics.
	:type stats_file: str
	:param prefix: A string to be used as prefix to the output dictionary keys.
	:type stats_file: str
	"""

	error_dict = {
		prefix + "surviving_perc": pd.np.nan,
		prefix + "short_perc": pd.np.nan,
		prefix + "empty_perc": pd.np.nan,
		prefix + "trimmed_perc": pd.np.nan,
		prefix + "untrimmed_perc": pd.np.nan,
		prefix + "trim_loss_perc": pd.np.nan}
	try:
		with open(stats_file) as handle:
			content = handle.readlines()  # list of strings per line
	except:
		return error_dict

	suf = "s" if not paired_end else " pairs"

	try:
		line = [i for i in range(len(content)) if "read{} processed; of these:".format(suf) in content[i]][0]
		total = int(re.sub("\D", "", re.sub("\(.*", "", content[line])))
		line = [i for i in range(len(content)) if "read{} available; of these:".format(suf) in content[i]][0]
		surviving = int(re.sub("\D", "", re.sub("\(.*", "", content[line])))
		line = [i for i in range(len(content)) if "short read{} filtered out after trimming by size control".format(suf) in content[i]][0]
		short = int(re.sub(" \(.*", "", content[line]).strip())
		line = [i for i in range(len(content)) if "empty read{} filtered out after trimming by size control".format(suf) in content[i]][0]
		empty = int(re.sub(" \(.*", "", content[line]).strip())
		line = [i for i in range(len(content)) if "trimmed read{} available after processing".format(suf) in content[i]][0]
		trimmed = int(re.sub(" \(.*", "", content[line]).strip())
		line = [i for i in range(len(content)) if "untrimmed read{} available after processing".format(suf) in content[i]][0]
		untrimmed = int(re.sub(" \(.*", "", content[line]).strip())
		return {
			prefix + "surviving_perc": (float(surviving) / total) * 100,
			prefix + "short_perc": (float(short) / total) * 100,
			prefix + "empty_perc": (float(empty) / total) * 100,
			prefix + "trimmed_perc": (float(trimmed) / total) * 100,
			prefix + "untrimmed_perc": (float(untrimmed) / total) * 100,
			prefix + "trim_loss_perc": ((total - float(surviving)) / total) * 100}
	except IndexError:
		return error_dict


def parse_mapping_stats(stats_file, prefix="", paired_end=True):
	"""
	:param stats_file: sambamba output file with duplicate statistics.
	:type stats_file: str
	:param prefix: A string to be used as prefix to the output dictionary keys.
	:type stats_file: str
	"""

	if not paired_end:
		error_dict = {
			prefix + "not_aligned_perc": pd.np.nan,
			prefix + "unique_aligned_perc": pd.np.nan,
			prefix + "multiple_aligned_perc": pd.np.nan,
			prefix + "perc_aligned": pd.np.nan}
	else:
		error_dict = {
			prefix + "paired_perc": pd.np.nan,
			prefix + "concordant_perc": pd.np.nan,
			prefix + "concordant_unique_perc": pd.np.nan,
			prefix + "concordant_multiple_perc": pd.np.nan,
			prefix + "not_aligned_or_discordant_perc": pd.np.nan,
			prefix + "not_aligned_perc": pd.np.nan,
			prefix + "unique_aligned_perc": pd.np.nan,
			prefix + "multiple_aligned_perc": pd.np.nan,
			prefix + "perc_aligned": pd.np.nan}

	try:
		with open(stats_file) as handle:
			content = handle.readlines()  # list of strings per line
	except:
		return error_dict

	if not paired_end:
		try:
			line = [i for i in range(len(content)) if "reads; of these:" in content[i]][0]
			total = int(re.sub("\D", "", re.sub("\(.*", "", content[line])))
			line = [i for i in range(len(content)) if "aligned 0 times" in content[i]][0]
			not_aligned_perc = float(re.search("\(.*%\)", content[line]).group()[1:-2])
			line = [i for i in range(len(content)) if " aligned exactly 1 time" in content[i]][0]
			unique_aligned_perc = float(re.search("\(.*%\)", content[line]).group()[1:-2])
			line = [i for i in range(len(content)) if " aligned >1 times" in content[i]][0]
			multiple_aligned_perc = float(re.search("\(.*%\)", content[line]).group()[1:-2])
			line = [i for i in range(len(content)) if "overall alignment rate" in content[i]][0]
			perc_aligned = float(re.sub("%.*", "", content[line]).strip())
			return {
				prefix + "not_aligned_perc": not_aligned_perc,
				prefix + "unique_aligned_perc": unique_aligned_perc,
				prefix + "multiple_aligned_perc": multiple_aligned_perc,
				prefix + "perc_aligned": perc_aligned}
		except IndexError:
			return error_dict

	if paired_end:
		try:
			line = [i for i in range(len(content)) if "reads; of these:" in content[i]][0]
			total = int(re.sub("\D", "", re.sub("\(.*", "", content[line])))
			line = [i for i in range(len(content)) if " were paired; of these:" in content[i]][0]
			paired_perc = float(re.search("\(.*%\)", content[line]).group()[1:-2])
			line = [i for i in range(len(content)) if "aligned concordantly 0 times" in content[i]][0]
			concordant_unaligned_perc = float(re.search("\(.*%\)", content[line]).group()[1:-2])
			line = [i for i in range(len(content)) if "aligned concordantly exactly 1 time" in content[i]][0]
			concordant_unique_perc = float(re.search("\(.*%\)", content[line]).group()[1:-2])
			line = [i for i in range(len(content)) if "aligned concordantly >1 times" in content[i]][0]
			concordant_multiple_perc = float(re.search("\(.*%\)", content[line]).group()[1:-2])
			line = [i for i in range(len(content)) if "mates make up the pairs; of these:" in content[i]][0]
			not_aligned_or_discordant = int(re.sub("\D", "", re.sub("\(.*", "", content[line])))
			d_fraction = (not_aligned_or_discordant / float(total))
			not_aligned_or_discordant_perc = d_fraction * 100
			line = [i for i in range(len(content)) if "aligned 0 times\n" in content[i]][0]
			not_aligned_perc = float(re.search("\(.*%\)", content[line]).group()[1:-2]) * d_fraction
			line = [i for i in range(len(content)) if " aligned exactly 1 time" in content[i]][0]
			unique_aligned_perc = float(re.search("\(.*%\)", content[line]).group()[1:-2]) * d_fraction
			line = [i for i in range(len(content)) if " aligned >1 times" in content[i]][0]
			multiple_aligned_perc = float(re.search("\(.*%\)", content[line]).group()[1:-2]) * d_fraction
			line = [i for i in range(len(content)) if "overall alignment rate" in content[i]][0]
			perc_aligned = float(re.sub("%.*", "", content[line]).strip())
			return {
				prefix + "paired_perc": paired_perc,
				prefix + "concordant_unaligned_perc": concordant_unaligned_perc,
				prefix + "concordant_unique_perc": concordant_unique_perc,
				prefix + "concordant_multiple_perc": concordant_multiple_perc,
				prefix + "not_aligned_or_discordant_perc": not_aligned_or_discordant_perc,
				prefix + "not_aligned_perc": not_aligned_perc,
				prefix + "unique_aligned_perc": unique_aligned_perc,
				prefix + "multiple_aligned_perc": multiple_aligned_perc,
				prefix + "perc_aligned": perc_aligned}
		except IndexError:
			return error_dict


def parse_duplicate_stats(stats_file, prefix=""):
	"""
	Parses sambamba markdup output, returns series with values.

	:param stats_file: sambamba output file with duplicate statistics.
	:type stats_file: str
	:param prefix: A string to be used as prefix to the output dictionary keys.
	:type stats_file: str
	"""

	error_dict = {
		prefix + "filtered_single_ends": pd.np.nan,
		prefix + "filtered_paired_ends": pd.np.nan,
		prefix + "duplicate_percentage": pd.np.nan}
	try:
		with open(stats_file) as handle:
			content = handle.readlines()  # list of strings per line
	except:
		return error_dict

	try:
		line = [i for i in range(len(content)) if "single ends (among them " in content[i]][0]
		single_ends = int(re.sub("\D", "", re.sub("\(.*", "", content[line])))
		line = [i for i in range(len(content)) if " end pairs...   done in " in content[i]][0]
		paired_ends = int(re.sub("\D", "", re.sub("\.\.\..*", "", content[line])))
		line = [i for i in range(len(content)) if " duplicates, sorting the list...   done in " in content[i]][0]
		duplicates = int(re.sub("\D", "", re.sub("\.\.\..*", "", content[line])))
		return {
			prefix + "filtered_single_ends": single_ends,
			prefix + "filtered_paired_ends": paired_ends,
			prefix + "duplicate_percentage": (float(duplicates) / (single_ends + paired_ends * 2)) * 100}
	except IndexError:
		return error_dict


def parse_peak_number(peak_file):
    from subprocess import check_output
    try:
        return {"peaks": int(check_output(["wc", "-l", peak_file]).split(" ")[0])}
    except:
        return {"peaks": pd.np.nan}


def calculate_frip(input_bam, input_bed, output, cpus=4):
    return ("sambamba view -t {0} -c  -L {1}  {2} > {3}"
            .format(cpus, input_bed, input_bam, output))


def parse_frip(frip_file, total_reads):
	"""
	Calculates the fraction of reads in peaks for a given sample.

	:param frip_file: A path to a file with the FRiP output.
	:type frip_file: str
	:param total_reads: Number of total reads (i.e., the denominator for the
		FRiP calculation)
	:type total_reads: int | float
	:rtype: float
	"""

	try:
		with open(frip_file, "r") as handle:
			content = handle.readlines()
	except:
		return pd.np.nan
	if content[0].strip() == "":
		return pd.np.nan

	reads_in_peaks = int(re.sub("\D", "", content[0]))
	frip = reads_in_peaks / float(total_reads)
	return frip


def parse_nsc_rsc(nsc_rsc_file):
    """
    Parses the values of NSC and RSC from a stats file.

    :param nsc_rsc_file: A sting path to a file with the NSC and RSC output (generally a tsv file).
    :type nsc_rsc_file: str
    """
    try:
        nsc_rsc = pd.read_csv(nsc_rsc_file, header=None, sep="\t")
        return {"NSC": nsc_rsc[8].squeeze(), "RSC": nsc_rsc[9].squeeze()}
    except:
        return {"NSC": pd.np.nan, "RSC": pd.np.nan}


def main():
	# Parse command-line arguments
	parser = ArgumentParser(
		prog="chipseq-pipeline",
		description="ChIP-seq pipeline."
	)
	parser = arg_parser(parser)
	parser = pypiper.add_pypiper_args(parser, groups=["all"])
	args = parser.parse_args()
	if args.sample_config is None:
		parser.print_help()
		return 1

	# Read in yaml configs
	series = pd.Series(yaml.load(open(args.sample_config, "r")))
	sample = Sample(series)
	# Create Sample object
	if sample.protocol == "ChIPmentation":
		sample = ChIPmentation(sample)
	else:
		sample = ChIPseqSample(sample)

	# Check if merged
	if len(sample.input_file_paths) > 1:
		sample.merged = True
	else:
		sample.merged = False
	sample.prj = AttributeDict(sample.prj)
	sample.paths = AttributeDict(sample.paths.__dict__)

	# Flag version of read type since it's a binary; handle case vagaries.
	read_type = sample.read_type
	try:
		sample.paired = (read_type.lower() == "paired")
	except AttributeError:
		print("WARNING: non-string read_type: {} ({})".format(
			read_type, type(read_type)))
		sample.paired = False

	# Set file paths
	sample.set_file_paths(sample.prj)
	# sample.make_sample_dirs()  # should be fixed to check if values of paths are strings and paths indeed

	# Start Pypiper object
	# Best practice is to name the pipeline with the name of the script;
	# or put the name in the pipeline interface.
	pipe_manager = pypiper.PipelineManager(name="chipseq", outfolder=sample.paths.sample_root, args=args)

	# Start main function
	process(sample, pipe_manager, args)

    if genome not in ['hg19', 'hg38', 'mm10', 'mm9']:
        print("Genome assembly is not known. Using size of human genome. Beware.")

def arg_parser(parser):
	"""
	Global options for pipeline.
	"""
	parser.add_argument(
		"-y", "--sample-yaml",
		dest="sample_config",
		help="Yaml config file with sample attributes; in addition to "
			"sample_name, this should define '{rt}', as 'single' or "
			"'paired'; 'ip', with the mark analyzed in a sample, and "
			"'{comparison}' with the name of a control sample (if the "
			"sample itself is not a control.)".format(rt="read_type", comparison=CHIP_COMPARE_COLUMN)
	)
	parser.add_argument(
		"-p", "--peak-caller",
		dest="peak_caller",
		choices=["macs2", "spp"],
		help="Peak caller algorithm.",
		default="macs2"
	)
	parser.add_argument("--pvalue", type=float, default=0.001, help="MACS2 p-value")
	parser.add_argument("--qvalue", type=float, help="Q-value for peak calling")
	return parser


def process(sample, pipe_manager, args):
	"""
	This takes unmapped Bam files and makes trimmed, aligned, duplicate marked
	and removed, indexed (and shifted if necessary) Bam files
	along with a UCSC browser track.

	:param Sample sample: individual Sample object to process
	:param pypiper.PipelineManager pipe_manager: PipelineManager to use during
		Sample processing
	:param argparse.Namespace args: binding between command-line option and
		argument, for specifying values various pipeline parameters
	"""
	print("Start processing ChIP-seq sample %s." % sample.name)

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

	# Create NGSTk instance
	tk = NGSTk(pm=pipe_manager)

	# Merge Bam files if more than one technical replicate
	if len(sample.input_file_paths) > 1:
		pipe_manager.timestamp("Merging bam files from replicates")
		cmd = tk.merge_bams(
			input_bams=sample.input_file_paths, merged_bam=sample.unmapped)
		pipe_manager.run(cmd, sample.unmapped, shell=True)
		sample.data_source = sample.unmapped

	# Fastqc
	pipe_manager.timestamp("Measuring sample quality with Fastqc")
	cmd = tk.fastqc_rename(
		input_bam=sample.data_source,
		output_dir=sample.paths.sample_root,
		sample_name=sample.sample_name
	)
	fastqc_target = os.path.join(sample.paths.sample_root, sample.sample_name + "_fastqc.zip")
	pipe_manager.run(cmd, fastqc_target, shell=True)
	report_dict(pipe_manager, parse_fastqc(fastqc_target, prefix="fastqc_"))

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
	if pipe_manager.config.parameters.trimmer == "trimmomatic":
		cmd = tk.trimmomatic(
			input_fastq1=sample.fastq1 if sample.paired else sample.fastq,
			input_fastq2=sample.fastq2 if sample.paired else None,
			output_fastq1=sample.trimmed1 if sample.paired else sample.trimmed,
			output_fastq1_unpaired=sample.trimmed1_unpaired if sample.paired else None,
			output_fastq2=sample.trimmed2 if sample.paired else None,
			output_fastq2_unpaired=sample.trimmed2_unpaired if sample.paired else None,
			cpus=args.cores,
			adapters=pipe_manager.config.resources.adapters,
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

	elif pipe_manager.config.parameters.trimmer == "skewer":
		cmd = tk.skewer(
			input_fastq1=sample.fastq1 if sample.paired else sample.fastq,
			input_fastq2=sample.fastq2 if sample.paired else None,
			output_prefix=os.path.join(sample.paths.unmapped, sample.sample_name),
			output_fastq1=sample.trimmed1 if sample.paired else sample.trimmed,
			output_fastq2=sample.trimmed2 if sample.paired else None,
			trim_log=sample.trimlog,
			cpus=args.cores,
			adapters=pipe_manager.config.resources.adapters
		)
		pipe_manager.run(cmd, sample.trimmed1 if sample.paired else sample.trimmed, shell=True)
		if not sample.paired:
			pipe_manager.clean_add(sample.trimmed, conditional=True)
		else:
			pipe_manager.clean_add(sample.trimmed1, conditional=True)
			pipe_manager.clean_add(sample.trimmed2, conditional=True)
		report_dict(pipe_manager, parse_trim_stats(sample.trimlog, prefix="trim_", paired_end=sample.paired))

	# Map
	pipe_manager.timestamp("Mapping reads with Bowtie2")
	cmd = tk.bowtie2_map(
		input_fastq1=sample.trimmed1 if sample.paired else sample.trimmed,
		input_fastq2=sample.trimmed2 if sample.paired else None,
		output_bam=sample.mapped,
		log=sample.aln_rates,
		metrics=sample.aln_metrics,
		genome_index=getattr(pipe_manager.config.resources.genomes, sample.genome),
		max_insert=pipe_manager.config.parameters.max_insert,
		cpus=args.cores
	)
	pipe_manager.run(cmd, sample.mapped, shell=True)
	report_dict(pipe_manager, parse_mapping_stats(sample.aln_rates, paired_end=sample.paired))

	# Filter reads
	pipe_manager.timestamp("Filtering reads for quality")
	cmd = tk.filter_reads(
		input_bam=sample.mapped,
		output_bam=sample.filtered,
		metrics_file=sample.dups_metrics,
		paired=sample.paired,
		cpus=args.cores,
		Q=pipe_manager.config.parameters.read_quality
	)
	pipe_manager.run(cmd, sample.filtered, shell=True)
	report_dict(pipe_manager, parse_duplicate_stats(sample.dups_metrics))

	# Index bams
	pipe_manager.timestamp("Indexing bamfiles with samtools")
	cmd = tk.index_bam(input_bam=sample.mapped)
	pipe_manager.run(cmd, sample.mapped + ".bai", shell=True)
	cmd = tk.index_bam(input_bam=sample.filtered)
	pipe_manager.run(cmd, sample.filtered + ".bai", shell=True)

	track_dir = os.path.dirname(sample.bigwig)
	if not os.path.exists(track_dir):
		os.makedirs(track_dir)

	# Make tracks
	# right now tracks are only made for bams without duplicates
	pipe_manager.timestamp("Making bigWig tracks from bam file")
	cmd = bam_to_bigwig(
		input_bam=sample.filtered,
		output_bigwig=sample.bigwig,
		genome_sizes=getattr(pipe_manager.config.resources.chromosome_sizes, sample.genome),
		genome=sample.genome,
		tagmented=pipe_manager.config.parameters.tagmented,  # by default make extended tracks
		normalize=pipe_manager.config.parameters.normalize_tracks,
		norm_factor=pipe_manager.config.parameters.norm_factor
	)
	pipe_manager.run(cmd, sample.bigwig, shell=True)

	# Plot fragment distribution
	if sample.paired and not os.path.exists(sample.insertplot):
		pipe_manager.timestamp("Plotting insert size distribution")
		tk.plot_atacseq_insert_sizes(
			bam=sample.filtered,
			plot=sample.insertplot,
			output_csv=sample.insertdata
		)

	# Count coverage genome-wide
	pipe_manager.timestamp("Calculating genome-wide coverage")
	cmd = tk.genome_wide_coverage(
		input_bam=sample.filtered,
		genome_windows=getattr(pipe_manager.config.resources.genome_windows, sample.genome),
		output=sample.coverage
	)
	pipe_manager.run(cmd, sample.coverage, shell=True)

	# Calculate NSC, RSC
	pipe_manager.timestamp("Assessing signal/noise in sample")
	cmd = tk.run_spp(
		input_bam=sample.filtered,
		output=sample.qc,
		plot=sample.qc_plot,
		cpus=args.cores
	)
	pipe_manager.run(cmd, sample.qc_plot, shell=True, nofail=True)
	report_dict(pipe_manager, parse_nsc_rsc(sample.qc))

	# If the sample is a control, we're finished.
	# The type/value for the comparison Sample in this case should be either
	# absent or a null-indicative/-suggestive value.
	comparison = getattr(sample, CHIP_COMPARE_COLUMN, None)
	if comparison in [None, "", "NA"]:
		pipe_manager.stop_pipeline()
		print("Finished processing sample {}".format(sample.name))
		return

	# The pipeline will now wait for the comparison sample file to be completed
	pipe_manager._wait_for_file(sample.filtered.replace(sample.name, comparison))

	# Call peaks.
	broad_mode = sample.broad
	peaks_folder = sample.paths.peaks
	treatment_file = sample.filtered
	control_file = sample.filtered.replace(sample.name, comparison)
	if not os.path.exists(peaks_folder):
		os.makedirs(peaks_folder)
	# TODO: include the filepaths as caller-neutral positionals/keyword args
	# TODO (cont.) once NGSTK API is tweaked.
	peak_call_kwargs = {
		"output_dir": peaks_folder,
		"broad": broad_mode,
		"qvalue": args.qvalue
	}
	if args.peak_caller == "macs2":
		cmd = tk.macs2_call_peaks(
				treatment_bams=treatment_file, control_bams=control_file,
                sample_name=sample.name, pvalue=args.pvalue,
                genome=sample.genome, paired=sample.paired, **peak_call_kwargs)
	else:
		cmd = tk.spp_call_peaks(
				treatment_bam=treatment_file, control_bam=control_file,
				treatment_name=sample.name, control_name=comparison,
				cpus=args.cpus, **peak_call_kwargs)
	pipe_manager.run(cmd, target=sample.peaks, shell=True)
	report_dict(pipe_manager, parse_peak_number(sample.peaks))

	# Do plotting as desired.
	if args.peak_caller == "macs2" and not broad_mode:
		pipe_manager.timestamp("Plotting MACS2 model")
		model_files_base = sample.name + "_model"

		# Create the command to run the model script.
		name_model_script = model_files_base + ".r"
		path_model_script = os.path.join(peaks_folder, name_model_script)
		exec_model_script = \
				"{} {}".format(pipe_manager.config.tools.Rscript, path_model_script)

		# Create the command to create and rename the model plot.
		plot_name = model_files_base + ".pdf"
		src_plot_path = os.path.join(os.getcwd(), plot_name)
		dst_plot_path = os.path.join(peaks_folder, plot_name)
		rename_model_plot = "mv {} {}".format(src_plot_path, dst_plot_path)

		# Run the model script and rename the model plot.
		pipe_manager.run([exec_model_script, rename_model_plot],
						 target=dst_plot_path, shell=True, nofail=True)

	# Calculate fraction of reads in peaks (FRiP)
	pipe_manager.timestamp("Calculating fraction of reads in peaks (FRiP)")
	cmd = tk.calculate_frip(
		input_bam=sample.filtered,
		input_bed=sample.peaks,
		output=sample.frip,
		cpus=args.cores
	)
	pipe_manager.run(cmd, sample.frip, shell=True)
	reads_SE = float(pipe_manager.get_stat("filtered_single_ends"))
	reads_PE = float(pipe_manager.get_stat("filtered_paired_ends"))
	total = 0.5 * (reads_SE + reads_PE)
	frip = parse_frip(sample.frip, total)
	pipe_manager.report_result("frip", frip)

	print("Finished processing sample %s." % sample.name)
	pipe_manager.stop_pipeline()



if __name__ == '__main__':
    try:
        sys.exit(main())
    except KeyboardInterrupt:
        print("Program canceled by user!")
        sys.exit(1)
