#!/usr/bin/env python

"""
RNA-seq pipeline
"""

import os
import sys
from argparse import ArgumentParser
import yaml
import pypiper
from pypiper.ngstk import NGSTk
from looper.models import AttributeDict, Sample

import pandas as pd


__author__ = "Andre Rendeiro"
__copyright__ = "Copyright 2018, Andre Rendeiro"
__credits__ = []
__license__ = "GPL2"
__version__ = "0.1"
__maintainer__ = "Andre Rendeiro"
__email__ = "arendeiro@cemm.oeaw.ac.at"
__status__ = "Development"


class RNASeqSample(Sample):
    """
    Class to model RNA-seq samples based on the ChIPseqSample class.

    :param pandas.Series series: Pandas `Series` object.
    """
    __library__ = "RNA-seq"

    def __init__(self, series):

        # Use pd.Series object to have all sample attributes
        if not isinstance(series, pd.Series):
            raise TypeError("Provided object is not a pandas Series.")
        super(RNASeqSample, self).__init__(series)

    def __repr__(self):
        return "RNA-seq sample '%s'" % self.sample_name

    def set_file_paths(self):
        """
        Sets the paths of all files for this sample.
        """
        # Inherit paths from Sample by running Sample's set_file_paths()
        super(RNASeqSample, self).set_file_paths()

        # Files in the root of the sample dir
        self.fastqc = os.path.join(
            self.paths.sample_root, self.sample_name + ".fastqc.zip")
        self.trimlog = os.path.join(
            self.paths.sample_root, self.sample_name + ".trimlog.txt")
        self.dups_metrics = os.path.join(
            self.paths.sample_root, self.sample_name + ".dups_metrics.txt")

        # Unmapped: merged bam, fastq, trimmed fastq
        self.paths.unmapped = os.path.join(self.paths.sample_root, "unmapped")
        self.unmapped = os.path.join(
            self.paths.unmapped, self.sample_name + ".bam")
        self.fastq = os.path.join(
            self.paths.unmapped, self.sample_name + ".fastq")
        self.fastq1 = os.path.join(
            self.paths.unmapped, self.sample_name + ".1.fastq")
        self.fastq2 = os.path.join(
            self.paths.unmapped, self.sample_name + ".2.fastq")
        self.fastq_unpaired = os.path.join(
            self.paths.unmapped, self.sample_name + ".unpaired.fastq")
        self.trimmed = os.path.join(
            self.paths.unmapped, self.sample_name + ".trimmed.fastq")
        self.trimmed1 = os.path.join(
            self.paths.unmapped, self.sample_name + ".1.trimmed.fastq")
        self.trimmed2 = os.path.join(
            self.paths.unmapped, self.sample_name + ".2.trimmed.fastq")
        self.trimmed1_unpaired = os.path.join(
            self.paths.unmapped, self.sample_name + ".1_unpaired.trimmed.fastq")
        self.trimmed2_unpaired = os.path.join(
            self.paths.unmapped, self.sample_name + ".2_unpaired.trimmed.fastq")

        # Expression quantification
        self.paths.kallisto = os.path.join(self.paths.sample_root, "kallisto")
        self.kallisto_output_dir = self.paths.kallisto
        self.kallisto_quantification = os.path.join(self.kallisto_output_dir, "abundance.tsv")


def kallisto(
        fastq_files, kallisto_index, read_type, output_dir,
        threads=4, bootstrap_number=100, fragment_size=300, fragment_std=200):
    """
    Quantify gene expression using Kallisto given reads in FASTQ format and an indexed transcriptome assembly.

    :param list fastq_files: List of paths to FASTQ files to be quantified.
    :param str kallisto_index: Path to Kallisto index.
    :param str read_type: Either 'single' or 'paired' for single-end or paired-end mode, respectively.
    :param str output_dir: Path to output_directory.
    :param int threads: Number of CPU threads to use.
    :param int bootstrap_number: Number of sampling bootstraps to perform.
    :param int fragment_size: Mean fragment size of library. Only required if `read_type` is 'single'.
    :param int fragment_std: Standard deviation of libray fragment size. Only required if `read_type` is 'single'.
    """
    cmd = """kallisto quant -t {} -b {} --index {}""".format(threads, bootstrap_number, kallisto_index)
    if read_type == "single":
        cmd += """ --single -l {} -s {}""".format(fragment_size, fragment_std)
    cmd += """ --output-dir {} {}""".format(output_dir, " ".join(fastq_files))

    return cmd


def report_dict(pipe, stats_dict):
    for key, value in stats_dict.items():
        pipe.report_result(key, value)


def parse_fastqc(fastqc_zip, prefix=""):
    """
    """
    import StringIO
    import zipfile
    import re

    error_dict = {
        prefix + "total_pass_filter_reads": pd.np.nan,
        prefix + "poor_quality": pd.np.nan,
        prefix + "read_length": pd.np.nan,
        prefix + "GC_perc": pd.np.nan}

    try:
        zfile = zipfile.ZipFile(fastqc_zip)
        content = StringIO.StringIO(zfile.read(os.path.join(
            zfile.filelist[0].filename, "fastqc_data.txt"))).readlines()
    except:
        return error_dict
    try:
        line = [i for i in range(len(content))
                if "Total Sequences" in content[i]][0]
        total = int(re.sub("\D", "", re.sub("\(.*", "", content[line])))
        line = [i for i in range(
            len(content)) if "Sequences flagged as poor quality" in content[i]][0]
        poor_quality = int(re.sub("\D", "", re.sub("\(.*", "", content[line])))
        line = [i for i in range(len(content))
                if "Sequence length	" in content[i]][0]
        seq_len = int(re.sub("\D", "", re.sub(
            " \(.*", "", content[line]).strip()))
        line = [i for i in range(len(content)) if "%GC" in content[i]][0]
        gc_perc = int(re.sub("\D", "", re.sub(
            " \(.*", "", content[line]).strip()))
        return {
            prefix + "total_pass_filter_reads": total,
            prefix + "poor_quality_perc": (float(poor_quality) / total) * 100,
            prefix + "read_length": seq_len,
            prefix + "GC_perc": gc_perc}
    except IndexError:
        return error_dict


def parse_trim_stats(stats_file, prefix="", paired_end=True):
    """
    :param stats_file: sambamba output file with duplicate statistics.
    :type stats_file: str
    :param prefix: A string to be used as prefix to the output dictionary keys.
    :type stats_file: str
    """
    import re

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
        line = [i for i in range(len(content)) if "read{} processed; of these:".format(
            suf) in content[i]][0]
        total = int(re.sub("\D", "", re.sub("\(.*", "", content[line])))
        line = [i for i in range(len(content)) if "read{} available; of these:".format(
            suf) in content[i]][0]
        surviving = int(re.sub("\D", "", re.sub("\(.*", "", content[line])))
        line = [i for i in range(len(
            content)) if "short read{} filtered out after trimming by size control".format(suf) in content[i]][0]
        short = int(re.sub(" \(.*", "", content[line]).strip())
        line = [i for i in range(len(
            content)) if "empty read{} filtered out after trimming by size control".format(suf) in content[i]][0]
        empty = int(re.sub(" \(.*", "", content[line]).strip())
        line = [i for i in range(len(
            content)) if "trimmed read{} available after processing".format(suf) in content[i]][0]
        trimmed = int(re.sub(" \(.*", "", content[line]).strip())
        line = [i for i in range(len(
            content)) if "untrimmed read{} available after processing".format(suf) in content[i]][0]
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


def parse_kallisto_stats(abundance):
    import numpy as np

    stats = dict()
    df = pd.read_table(abundance, sep="\t")
    stats['transcripts'] = df.shape[0]
    stats['zero-count_transcripts'] = (df['est_counts'] == 0).sum()
    stats['non-zero-count_transcripts'] = (df['est_counts'] > 0).sum()
    log_tpm = np.log2(1 + df['tpm'])
    stats['log2tpm_mean'] = log_tpm.mean()
    stats['log2tpm_median'] = log_tpm.median()
    p_log_tpm = np.log2(1 + df['tpm'].where(lambda x: x > 0)).dropna()
    stats['non-zero_log2tpm_mean'] = p_log_tpm.mean()
    stats['non-zero_log2tpm_median'] = p_log_tpm.median()
    try:
        import scipy
        stats['log2tpm_iqr'] = scipy.stats.iqr(log_tpm)
        stats['non-zero_log2tpm_iqr'] = scipy.stats.iqr(p_log_tpm)
    except ImportError:
        stats['log2tpm_iqr'] = np.nan
        stats['non-zero_log2tpm_iqr'] = np.nan

    return stats


def main():
    # Parse command-line arguments
    parser = ArgumentParser(
        prog="rnaseq-pipeline",
        description="RNA-seq pipeline.")
    parser = arg_parser(parser)
    parser = pypiper.add_pypiper_args(parser, groups=["all"])
    args = parser.parse_args()

    # Read in yaml configs
    sample = RNASeqSample(pd.Series(yaml.load(open(args.sample_config, "r"))))

    # Check if merged
    if len(sample.data_path.split(" ")) > 1:
        sample.merged = True
    else:
        sample.merged = False
    sample.prj = AttributeDict(sample.prj)
    sample.paths = AttributeDict(sample.paths.__dict__)

    # Check read type if not provided
    if not hasattr(sample, "ngs_inputs"):
        sample.ngs_inputs = [sample.data_source]
    if not hasattr(sample, "read_type"):
        sample.set_read_type()

    # Shorthand for read_type
    if sample.read_type == "paired":
        sample.paired = True
    else:
        sample.paired = False

    # Set file paths
    sample.set_file_paths()
    # sample.make_sample_dirs()  # should be fixed to check if values of paths are strings and paths indeed

    # Start Pypiper object
    # Best practice is to name the pipeline with the name of the script;
    # or put the name in the pipeline interface.
    pipe_manager = pypiper.PipelineManager(
        name="rnaseq", outfolder=sample.paths.sample_root, args=args)
    pipe_manager.config.tools.scripts_dir = os.path.join(
        os.path.dirname(os.path.realpath(__file__)), "tools")

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
        type=str)

    return parser


def process(sample, pipe_manager, args):
    """
    This takes unmapped Bam files and makes trimmed, aligned, duplicate marked
    and removed, indexed, shifted Bam files along with a UCSC browser track.
    Peaks are called and filtered.
    """
    print("Start processing RNA-seq sample %s." % sample.sample_name)

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
    if len(sample.data_path.split(" ")) > 1:
        pipe_manager.timestamp("Merging bam files from replicates")
        cmd = tk.merge_bams(
            # this is a list of sample paths
            input_bams=sample.data_path.split(" "),
            merged_bam=sample.unmapped
        )
        pipe_manager.run(cmd, sample.unmapped, shell=True)
        sample.data_path = sample.unmapped

    # Fastqc
    pipe_manager.timestamp("Measuring sample quality with Fastqc")
    cmd = tk.fastqc_rename(
        input_bam=sample.data_path,
        output_dir=sample.paths.sample_root,
        sample_name=sample.sample_name
    )
    pipe_manager.run(cmd, os.path.join(sample.paths.sample_root,
                                       sample.sample_name + "_fastqc.zip"), shell=True)
    report_dict(pipe_manager, parse_fastqc(os.path.join(
        sample.paths.sample_root, sample.sample_name + "_fastqc.zip"), prefix="fastqc_"))

    # Convert bam to fastq
    pipe_manager.timestamp("Converting to Fastq format")
    cmd = tk.bam2fastq(
        inputBam=sample.data_path,
        outputFastq=sample.fastq1 if sample.paired else sample.fastq,
        outputFastq2=sample.fastq2 if sample.paired else None,
        unpairedFastq=sample.fastq_unpaired if sample.paired else None
    )
    pipe_manager.run(
        cmd, sample.fastq1 if sample.paired else sample.fastq, shell=True)
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
            inputFastq1=sample.fastq1 if sample.paired else sample.fastq,
            inputFastq2=sample.fastq2 if sample.paired else None,
            outputFastq1=sample.trimmed1 if sample.paired else sample.trimmed,
            outputFastq1unpaired=sample.trimmed1_unpaired if sample.paired else None,
            outputFastq2=sample.trimmed2 if sample.paired else None,
            outputFastq2unpaired=sample.trimmed2_unpaired if sample.paired else None,
            cpus=args.cores,
            adapters=pipe_manager.config.resources.adapters,
            log=sample.trimlog
        )
        pipe_manager.run(
            cmd, sample.trimmed1 if sample.paired else sample.trimmed, shell=True)
        if not sample.paired:
            pipe_manager.clean_add(sample.trimmed, conditional=True)
        else:
            pipe_manager.clean_add(sample.trimmed1, conditional=True)
            pipe_manager.clean_add(sample.trimmed1_unpaired, conditional=True)
            pipe_manager.clean_add(sample.trimmed2, conditional=True)
            pipe_manager.clean_add(sample.trimmed2_unpaired, conditional=True)

    elif pipe_manager.config.parameters.trimmer == "skewer":
        cmd = tk.skewer(
            inputFastq1=sample.fastq1 if sample.paired else sample.fastq,
            inputFastq2=sample.fastq2 if sample.paired else None,
            outputPrefix=os.path.join(
                sample.paths.unmapped, sample.sample_name),
            outputFastq1=sample.trimmed1 if sample.paired else sample.trimmed,
            outputFastq2=sample.trimmed2 if sample.paired else None,
            trimLog=sample.trimlog,
            cpus=args.cores,
            adapters=pipe_manager.config.resources.adapters
        )
        pipe_manager.run(
            cmd, sample.trimmed1 if sample.paired else sample.trimmed, shell=True)
        if not sample.paired:
            pipe_manager.clean_add(sample.trimmed, conditional=True)
        else:
            pipe_manager.clean_add(sample.trimmed1, conditional=True)
            pipe_manager.clean_add(sample.trimmed2, conditional=True)

        report_dict(pipe_manager, parse_trim_stats(
            sample.trimlog, prefix="trim_", paired_end=sample.paired))

    # Quantify gene expression
    pipe_manager.timestamp("Quantifying expression with Kallisto")
    cmd = kallisto(
        fastq_files=[sample.trimmed1, sample.trimmed2 if sample.paired else sample.trimmed],
        kallisto_index=getattr(pipe_manager.config.resources.kallisto_index, sample.genome),
        read_type=sample.read_type,
        output_dir=sample.kallisto_output_dir,
        threads=args.cores,
        bootstrap_number=pipe_manager.config.parameters.bootstrap_number,
        fragment_size=pipe_manager.config.parameters.fragment_size,
        fragment_std=pipe_manager.config.parameters.fragment_std)
    pipe_manager.run(cmd, sample.kallisto_quantification, shell=True)
    report_dict(pipe_manager, parse_kallisto_stats(sample.kallisto_quantification))

    # Finish up
    print(pipe_manager.stats_dict)

    pipe_manager.stop_pipeline()
    print("Finished processing sample %s." % sample.sample_name)


if __name__ == '__main__':
    try:
        sys.exit(main())
    except KeyboardInterrupt:
        print("Program canceled by user!")
        sys.exit(1)
