#!/usr/bin/env python

"""
ATAC-seq pipeline
"""

import os
from os.path import join as pjoin
import sys
from argparse import ArgumentParser
import yaml
import pypiper
from pypiper.ngstk import NGSTk
from attmap import AttributeDict

import pandas as pd

# for TSS analysis
import pybedtools as bedtools
import matplotlib
import matplotlib.pyplot as plt

matplotlib.use("Agg")


__author__ = "Andre Rendeiro"
__copyright__ = "Copyright 2019, Andre Rendeiro"
__credits__ = []
__license__ = "GPL2"
__version__ = "0.4"
__maintainer__ = "Andre Rendeiro"
__email__ = "arendeiro@cemm.oeaw.ac.at"
__status__ = "Development"


class ATACseqSample:
    """
    Class to model ATAC-seq samples based on the ChIPseqSample class.

    :param series: Pandas `Series` object.
    :type series: pandas.Series
    """

    __library__ = "ATAC-seq"

    def __init__(self, series):

        # Use pd.Series object to have all sample attributes
        if not isinstance(series, pd.Series):
            raise TypeError("Provided object is not a pandas Series.")
        # super(ATACseqSample, self).__init__(series)

        self.tagmented = True
        for k, v in series.items():
            setattr(self, k, v)

    def __repr__(self):
        return "ATAC-seq sample '%s'" % self.sample_name

    def set_file_paths(self, project):
        """
        Sets the paths of all files for this sample.
        """
        # Inherit paths from Sample by running Sample's set_file_paths()
        super(ATACseqSample, self)  # .set_file_paths(project)

        # Files in the root of the sample dir
        prefix = pjoin(self.sample_root, self.sample_name)
        self.fastqc_initial_output = (
            pjoin(
                self.sample_root,
                os.path.splitext(os.path.basename(self.data_source))[0],
            )
            + "_fastqc.zip"
        )
        self.fastqc = prefix + ".fastqc.zip"
        self.trimlog = prefix + ".trimlog.txt"
        self.aln_rates = prefix + ".aln_rates.txt"
        self.aln_metrics = prefix + ".aln_metrics.txt"
        self.dups_metrics = prefix + ".dups_metrics.txt"

        # Unmapped: merged bam, fastq, trimmed fastq
        self.unmapped_dir = pjoin(self.sample_root, "unmapped")
        unmapped_p = pjoin(self.unmapped_dir, self.sample_name)
        self.unmapped = unmapped_p + ".bam"
        self.fastq = unmapped_p + ".fastq"
        self.fastq1 = unmapped_p + ".1.fastq"
        self.fastq2 = unmapped_p + ".2.fastq"
        self.fastq_unpaired = unmapped_p + ".unpaired.fastq"
        self.trimmed = unmapped_p + ".trimmed.fastq"
        self.trimmed1 = unmapped_p + ".1.trimmed.fastq"
        self.trimmed2 = unmapped_p + ".2.trimmed.fastq"
        self.trimmed1_unpaired = unmapped_p + ".1_unpaired.trimmed.fastq"
        self.trimmed2_unpaired = unmapped_p + ".2_unpaired.trimmed.fastq"

        # Mapped: mapped, duplicates marked, removed, reads shifted
        self.mapped_dir = pjoin(self.sample_root, "mapped")
        mapped_p = pjoin(self.mapped_dir, self.sample_name)
        self.mapped = mapped_p + ".trimmed.bowtie2.bam"
        self.filtered = mapped_p + ".trimmed.bowtie2.filtered.bam"
        # this will create additional bam files with reads shifted
        self.filteredshifted = (
            mapped_p + ".trimmed.bowtie2.filtered.shifted.bam"
        )

        # Files in the root of the sample dir
        self.frip = prefix + "_FRiP.txt"
        self.oracle_frip = prefix + "_oracle_FRiP.txt"

        # Coverage: read coverage in windows genome-wide
        self.coverage_dir = pjoin(self.sample_root, "coverage")
        self.coverage = pjoin(self.coverage_dir, self.sample_name + ".cov")

        # self.bigwig = pjoin(self.coverage, self.name + ".bigWig")
        self.insertplot = prefix + "_insertLengths.pdf"
        self.insertdata = prefix + "_insertLengths.csv"
        self.mitochondrial_stats = prefix + "_mitochondrial_stats.tsv"
        self.qc = prefix + "_qc.tsv"
        self.qc_plot = prefix + "_qc.pdf"
        self.tss_dir = pjoin(self.sample_root, "tss")
        tss_p = pjoin(self.tss_dir, self.sample_name)
        self.tss_plot = tss_p + "_TSS.svg"
        self.tss_hist = tss_p + "_TSS_histogram.csv"
        self.tss_lock = tss_p + "_TSS.complete"

        # Peaks: peaks called and derivate files
        self.peaks_dir = pjoin(self.sample_root, "peaks")
        peaks_p = pjoin(self.peaks_dir, self.sample_name)
        self.peaks = peaks_p + "_peaks.narrowPeak"
        self.summits = peaks_p + "_summits.bed"
        self.filteredPeaks = peaks_p + "_peaks.filtered.bed"


class DNaseSample(ATACseqSample):
    """
    Class to model DNase-seq samples based on the ChIPseqSample class.

    :param series: Pandas `Series` object.
    :type series: pandas.Series
    """

    __library__ = "DNase-seq"

    def __init__(self, series):

        # Use pd.Series object to have all sample attributes
        if not isinstance(series, pd.Series):
            raise TypeError("Provided object is not a pandas Series.")
        super(DNaseSample, self).__init__(series)

    def __repr__(self):
        return "DNase-seq sample '%s'" % self.sample_name

    def set_file_paths(self, project):
        super(DNaseSample, self).set_file_paths(project)


def main():
    # Parse command-line arguments
    parser = ArgumentParser(
        prog="atacseq-pipeline", description="ATAC-seq pipeline."
    )
    parser = arg_parser(parser)
    parser = pypiper.add_pypiper_args(
        parser, groups=["ngs", "looper", "resource", "pypiper"]
    )
    args = parser.parse_args()
    if args.sample_config is None or args.output_parent is None:
        parser.print_help()
        return 1

    # Read in yaml configs
    series = pd.Series(yaml.safe_load(open(args.sample_config, "r")))
    series["sample_root"] = args.output_parent
    print(series)
    # Create Sample object
    if series["protocol"] != "DNase-seq":
        sample = ATACseqSample(series)
    else:
        sample = DNaseSample(series)

    print(sample)
    # Check if merged
    if len(sample.data_source.split(" ")) > 1:
        sample.merged = True
    else:
        sample.merged = False
    sample.prj = AttributeDict(sample.prj)
    sample.paths = AttributeDict(sample.__dict__)

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
    sample.set_file_paths(sample.prj)

    # Start Pypiper object
    # Best practice is to name the pipeline with the name of the script;
    # or put the name in the pipeline interface.
    pipe_manager = pypiper.PipelineManager(
        name="atacseq", outfolder=sample.sample_root, args=args
    )
    pipe_manager.config.tools.scripts_dir = pjoin(
        os.path.dirname(os.path.realpath(__file__)), "tools"
    )

    # Start main function
    process(sample, pipe_manager, args)


def process(sample, pipe_manager, args):
    """
    This takes unmapped Bam files and makes trimmed, aligned, duplicate marked
    and removed, indexed, shifted Bam files along with a UCSC browser track.
    Peaks are called and filtered.
    """
    print("Start processing ATAC-seq sample %s." % sample.sample_name)

    # for path in ["sample_root"] + list(sample.__dict__.keys()):
    for path in [
        "sample_root",
        "unmapped_dir",
        "mapped_dir",
        "peaks_dir",
        "coverage_dir",
        "tss_dir",
    ]:
        p = getattr(sample, path)
        try:
            exists = os.path.exists(p)
        except TypeError:
            continue
        if not exists:
            msg = "Cannot create '{}' path: {}".format(path, p)
            try:
                os.mkdir(p)
            except OSError(msg):
                raise

    # Create NGSTk instance
    tk = NGSTk(pm=pipe_manager)

    # Merge Bam files if more than one technical replicate
    if len(sample.data_source.split(" ")) > 1:
        pipe_manager.timestamp("Merging bam files from replicates")
        cmd = tk.merge_bams(
            input_bams=sample.data_source.split(
                " "
            ),  # this is a list of sample paths
            merged_bam=sample.unmapped,
        )
        pipe_manager.run(cmd, sample.unmapped, shell=True)
        sample.data_source = sample.unmapped

    # Fastqc
    pipe_manager.timestamp("Measuring sample quality with Fastqc")
    if not os.path.exists(sample.fastqc):
        cmd = tk.fastqc(file=sample.data_source, output_dir=sample.sample_root)
        pipe_manager.run(cmd, sample.fastqc_initial_output, shell=False)
    # # rename output
    if os.path.exists(sample.fastqc_initial_output):
        os.rename(sample.fastqc_initial_output, sample.fastqc)
    report_dict(pipe_manager, parse_fastqc(sample.fastqc, prefix="fastqc_"))

    # Convert bam to fastq
    pipe_manager.timestamp("Converting to Fastq format")
    cmd = tk.bam2fastq(
        input_bam=sample.data_source,
        output_fastq=sample.fastq1 if sample.paired else sample.fastq,
        output_fastq2=sample.fastq2 if sample.paired else None,
        unpaired_fastq=sample.fastq_unpaired if sample.paired else None,
    )
    pipe_manager.run(
        cmd, sample.fastq1 if sample.paired else sample.fastq, shell=True
    )
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
            output_fastq1_unpaired=sample.trimmed1_unpaired
            if sample.paired
            else None,
            output_fastq2=sample.trimmed2 if sample.paired else None,
            output_fastq2_unpaired=sample.trimmed2_unpaired
            if sample.paired
            else None,
            cpus=args.cores,
            adapters=pipe_manager.config.resources.adapters,
            log=sample.trimlog,
        )
        pipe_manager.run(
            cmd,
            sample.trimmed1 if sample.paired else sample.trimmed,
            shell=True,
        )
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
            output_prefix=pjoin(sample.unmapped, sample.sample_name),
            output_fastq1=sample.trimmed1 if sample.paired else sample.trimmed,
            output_fastq2=sample.trimmed2 if sample.paired else None,
            log=sample.trimlog,
            cpus=args.cores,
            adapters=pipe_manager.config.resources.adapters,
        )
        pipe_manager.run(
            cmd,
            sample.trimmed1 if sample.paired else sample.trimmed,
            shell=True,
        )
        if not sample.paired:
            pipe_manager.clean_add(sample.trimmed, conditional=True)
        else:
            pipe_manager.clean_add(sample.trimmed1, conditional=True)
            pipe_manager.clean_add(sample.trimmed2, conditional=True)

        report_dict(
            pipe_manager,
            parse_trim_stats(
                sample.trimlog, prefix="trim_", paired_end=sample.paired
            ),
        )

    # Map
    pipe_manager.timestamp("Mapping reads with Bowtie2")
    cmd = tk.bowtie2_map(
        input_fastq1=sample.trimmed1 if sample.paired else sample.trimmed,
        input_fastq2=sample.trimmed2 if sample.paired else None,
        output_bam=sample.mapped,
        log=sample.aln_rates,
        metrics=sample.aln_metrics,
        genome_index=getattr(
            pipe_manager.config.resources.genome_index, sample.genome
        ),
        max_insert=pipe_manager.config.parameters.max_insert,
        cpus=args.cores,
    )
    pipe_manager.run(cmd, sample.mapped, shell=True)
    report_dict(
        pipe_manager,
        parse_mapping_stats(sample.aln_rates, paired_end=sample.paired),
    )

    # Get mitochondrial reads
    pipe_manager.timestamp("Getting mitochondrial stats")
    cmd = tk.get_mitochondrial_reads(
        bam_file=sample.mapped,
        output=sample.mitochondrial_stats,
        cpus=args.cores,
    )
    pipe_manager.run(cmd, sample.mitochondrial_stats, shell=True, nofail=True)
    report_dict(
        pipe_manager,
        parse_duplicate_stats(sample.mitochondrial_stats, prefix="MT_"),
    )

    # Filter reads
    pipe_manager.timestamp("Filtering reads for quality")
    cmd = tk.filter_reads(
        input_bam=sample.mapped,
        output_bam=sample.filtered,
        metrics_file=sample.dups_metrics,
        paired=sample.paired,
        cpus=args.cores,
        Q=pipe_manager.config.parameters.read_quality,
    )
    pipe_manager.run(cmd, sample.filtered, shell=True)
    report_dict(pipe_manager, parse_duplicate_stats(sample.dups_metrics))

    # Index bams
    pipe_manager.timestamp("Indexing bamfiles with samtools")
    cmd = tk.index_bam(input_bam=sample.mapped)
    pipe_manager.run(cmd, sample.mapped + ".bai", shell=True)
    cmd = tk.index_bam(input_bam=sample.filtered)
    pipe_manager.run(cmd, sample.filtered + ".bai", shell=True)

    # Shift reads
    if args.shift_reads:
        pipe_manager.timestamp("Shifting reads of tagmented sample")
        cmd = tk.shift_reads(
            input_bam=sample.filtered,
            genome=sample.genome,
            output_bam=sample.filteredshifted,
        )
        pipe_manager.run(cmd, sample.filteredshifted, shell=True)

        cmd = tk.index_bam(input_bam=sample.filteredshifted)
        pipe_manager.run(cmd, sample.filteredshifted + ".bai", shell=True)

    # Run TSS enrichment
    tss_enrichment = run_tss_analysis(
        sample=sample,
        bam_file=sample.filtered,
        chrom_file=getattr(
            pipe_manager.config.resources.chromosome_sizes, sample.genome
        ),
        tss_file=getattr(
            pipe_manager.config.resources.unique_tss, sample.genome
        ),
    )
    report_dict(pipe_manager, {"tss_enrichment": tss_enrichment})

    # Call peaks
    pipe_manager.timestamp("Calling peaks with MACS2")
    # make dir for output (macs fails if it does not exist)
    if not os.path.exists(sample.peaks):
        os.makedirs(sample.peaks)

    cmd = tk.macs2_call_peaks_atacseq(
        treatment_bam=sample.filtered,
        output_dir=sample.peaks,
        sample_name=sample.sample_name,
        genome=sample.genome,
    )
    pipe_manager.run(cmd, sample.peaks, shell=True)
    report_dict(pipe_manager, parse_peak_number(sample.peaks))

    # Calculate fraction of reads in peaks (FRiP)
    pipe_manager.timestamp("Calculating fraction of reads in peaks (FRiP)")
    cmd = tk.calculate_frip(
        input_bam=sample.filtered,
        input_bed=sample.peaks,
        output=sample.frip,
        cpus=args.cores,
    )
    pipe_manager.run(cmd, sample.frip, shell=True)
    total = float(pipe_manager.stats_dict["filtered_single_ends"]) + (
        float(pipe_manager.stats_dict["filtered_paired_ends"]) / 2.0
    )
    report_dict(pipe_manager, parse_frip(sample.frip, total))

    # on an oracle peak list
    if hasattr(
        pipe_manager.config.resources.oracle_peak_regions, sample.genome
    ):
        cmd = calculate_frip(
            input_bam=sample.filtered,
            input_bed=getattr(
                pipe_manager.config.resources.oracle_peak_regions, sample.genome
            ),
            output=sample.oracle_frip,
            cpus=args.cores,
        )
        pipe_manager.run(cmd, sample.oracle_frip, shell=True)
        report_dict(
            pipe_manager,
            parse_frip(sample.oracle_frip, total, prefix="oracle_"),
        )

    # Plot fragment distribution
    if sample.paired and not os.path.exists(sample.insertplot):
        pipe_manager.timestamp("Plotting insert size distribution")
        tk.plot_atacseq_insert_sizes(
            bam=sample.filtered,
            plot=sample.insertplot,
            output_csv=sample.insertdata,
        )

    # # Count coverage genome-wide
    # pipe_manager.timestamp("Calculating genome-wide coverage")
    # cmd = tk.genome_wide_coverage(
    #     input_bam=sample.filtered,
    #     genome_windows=getattr(pipe_manager.config.resources.genome_windows, sample.genome),
    #     output=sample.coverage)
    # pipe_manager.run(cmd, sample.coverage, shell=True)

    # Calculate NSC, RSC
    pipe_manager.timestamp("Assessing signal/noise in sample")
    cmd = tk.run_spp(
        input_bam=sample.filtered,
        output=sample.qc,
        plot=sample.qc_plot,
        cpus=args.cores,
    )
    pipe_manager.run(cmd, sample.qc_plot, shell=True, nofail=True)
    report_dict(pipe_manager, parse_nsc_rsc(sample.qc))

    # Make tracks
    track_dir = os.path.dirname(sample.bigwig)
    if not os.path.exists(track_dir):
        os.makedirs(track_dir)
    # right now tracks are only made for bams without duplicates
    pipe_manager.timestamp("Making bigWig tracks from BAM file")
    cmd = bam_to_bigwig(
        input_bam=sample.filtered,
        output_bigwig=sample.bigwig,
        genome=sample.genome,
        normalization_method="RPGC",
    )
    pipe_manager.run(cmd, sample.bigwig, shell=True)

    print(pipe_manager.stats_dict)

    pipe_manager.stop_pipeline()
    print("Finished processing sample %s." % sample.sample_name)


def arg_parser(parser):
    """
    Global options for pipeline.
    """
    parser.add_argument(
        "-y",
        "--sample-yaml",
        dest="sample_config",
        help="Yaml config file with sample attributes; in addition to "
        "sample_name, this should define '{rt}', as 'single' or "
        "'paired'".format(rt="read_type"),
    )
    parser.add_argument(
        "-p",
        "--peak-caller",
        dest="peak_caller",
        choices=["macs2", "spp"],
        help="Peak caller algorithm.",
        default="macs2",
    )
    parser.add_argument(
        "--shift-reads",
        dest="shift_reads",
        action="store_true",
        help="Whether to produce a file with 5' read positions shifted to their original position.",
    )
    parser.add_argument(
        "--pvalue", type=float, default=0.001, help="MACS2 p-value"
    )
    parser.add_argument("--qvalue", type=float, help="Q-value for peak calling")
    return parser


def report_dict(pipe, stats_dict):
    for key, value in stats_dict.items():
        pipe.report_result(key, value)


def parse_fastqc(fastqc_zip, prefix=""):
    """
    """
    # TODO: review the parsing of paired-end files
    import zipfile
    import re

    error_dict = {
        prefix + "total_pass_filter_reads": pd.np.nan,
        prefix + "poor_quality": pd.np.nan,
        prefix + "read_length": pd.np.nan,
        prefix + "GC_perc": pd.np.nan,
    }

    try:
        zfile = zipfile.ZipFile(fastqc_zip)
        content = (
            zfile.read(pjoin(zfile.filelist[0].filename, "fastqc_data.txt"))
            .decode()
            .split("\n")
        )
    except:
        return error_dict
    try:
        line = [
            i for i in range(len(content)) if "Total Sequences" in content[i]
        ][0]
        total = int(re.sub(r"\D", "", re.sub(r"\(.*", "", content[line])))
        line = [
            i
            for i in range(len(content))
            if "Sequences flagged as poor quality" in content[i]
        ][0]
        poor_quality = int(
            re.sub(r"\D", "", re.sub(r"\(.*", "", content[line]))
        )
        line = [
            i for i in range(len(content)) if "Sequence length" in content[i]
        ][0]
        seq_len = int(
            re.sub(r"\D", "", re.sub(r" \(.*", "", content[line]).strip())
        )
        line = [i for i in range(len(content)) if "%GC" in content[i]][0]
        gc_perc = int(
            re.sub(r"\D", "", re.sub(r" \(.*", "", content[line]).strip())
        )
        return {
            prefix + "total_pass_filter_reads": total,
            prefix + "poor_quality_perc": (float(poor_quality) / total) * 100,
            prefix + "read_length": seq_len,
            prefix + "GC_perc": gc_perc,
        }
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

    stats_dict = {
        prefix + "surviving_perc": pd.np.nan,
        prefix + "short_perc": pd.np.nan,
        prefix + "empty_perc": pd.np.nan,
        prefix + "trimmed_perc": pd.np.nan,
        prefix + "untrimmed_perc": pd.np.nan,
        prefix + "trim_loss_perc": pd.np.nan,
    }
    try:
        with open(stats_file) as handle:
            content = handle.readlines()  # list of strings per line
    except:
        return stats_dict

    suf = "s" if not paired_end else " pairs"

    try:
        line = [
            i
            for i in range(len(content))
            if "read{} processed; of these:".format(suf) in content[i]
        ][0]
        total = int(re.sub(r"\D", "", re.sub(r"\(.*", "", content[line])))
    except IndexError:
        return stats_dict
    try:
        line = [
            i
            for i in range(len(content))
            if "read{} available; of these:".format(suf) in content[i]
        ][0]
        surviving = int(re.sub(r"\D", "", re.sub(r"\(.*", "", content[line])))
        stats_dict[prefix + "surviving_perc"] = (float(surviving) / total) * 100
        stats_dict[prefix + "trim_loss_perc"] = (
            (total - float(surviving)) / total
        ) * 100
    except IndexError:
        pass
    try:
        line = [
            i
            for i in range(len(content))
            if "short read{} filtered out after trimming by size control".format(
                suf
            )
            in content[i]
        ][0]
        short = int(re.sub(r" \(.*", "", content[line]).strip())
        stats_dict[prefix + "short_perc"] = (float(short) / total) * 100
    except IndexError:
        pass
    try:
        line = [
            i
            for i in range(len(content))
            if "empty read{} filtered out after trimming by size control".format(
                suf
            )
            in content[i]
        ][0]
        empty = int(re.sub(r" \(.*", "", content[line]).strip())
        stats_dict[prefix + "empty_perc"] = (float(empty) / total) * 100
    except IndexError:
        pass
    try:
        line = [
            i
            for i in range(len(content))
            if "trimmed read{} available after processing".format(suf)
            in content[i]
        ][0]
        trimmed = int(re.sub(r" \(.*", "", content[line]).strip())
        stats_dict[prefix + "trimmed_perc"] = (float(trimmed) / total) * 100
    except IndexError:
        pass
    try:
        line = [
            i
            for i in range(len(content))
            if "untrimmed read{} available after processing".format(suf)
            in content[i]
        ][0]
        untrimmed = int(re.sub(r" \(.*", "", content[line]).strip())
        stats_dict[prefix + "untrimmed_perc"] = (float(untrimmed) / total) * 100
    except IndexError:
        pass
    return stats_dict


def parse_mapping_stats(stats_file, prefix="", paired_end=True):
    """
    :param stats_file: sambamba output file with duplicate statistics.
    :type stats_file: str
    :param prefix: A string to be used as prefix to the output dictionary keys.
    :type stats_file: str
    """
    import re

    if not paired_end:
        error_dict = {
            prefix + "not_aligned_perc": pd.np.nan,
            prefix + "unique_aligned_perc": pd.np.nan,
            prefix + "multiple_aligned_perc": pd.np.nan,
            prefix + "perc_aligned": pd.np.nan,
        }
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
            prefix + "perc_aligned": pd.np.nan,
        }

    try:
        with open(stats_file) as handle:
            content = handle.readlines()  # list of strings per line
    except:
        return error_dict

    if not paired_end:
        try:
            line = [
                i
                for i in range(len(content))
                if "reads; of these:" in content[i]
            ][0]
            total = int(re.sub(r"\D", "", re.sub(r"\(.*", "", content[line])))
            line = [
                i
                for i in range(len(content))
                if "aligned 0 times" in content[i]
            ][0]
            not_aligned_perc = float(
                re.search(r"\(.*%\)", content[line]).group()[1:-2]
            )
            line = [
                i
                for i in range(len(content))
                if " aligned exactly 1 time" in content[i]
            ][0]
            unique_aligned_perc = float(
                re.search(r"\(.*%\)", content[line]).group()[1:-2]
            )
            line = [
                i
                for i in range(len(content))
                if " aligned >1 times" in content[i]
            ][0]
            multiple_aligned_perc = float(
                re.search(r"\(.*%\)", content[line]).group()[1:-2]
            )
            line = [
                i
                for i in range(len(content))
                if "overall alignment rate" in content[i]
            ][0]
            perc_aligned = float(re.sub("%.*", "", content[line]).strip())
            return {
                prefix + "not_aligned_perc": not_aligned_perc,
                prefix + "unique_aligned_perc": unique_aligned_perc,
                prefix + "multiple_aligned_perc": multiple_aligned_perc,
                prefix + "perc_aligned": perc_aligned,
            }
        except IndexError:
            return error_dict

    if paired_end:
        try:
            line = [
                i
                for i in range(len(content))
                if "reads; of these:" in content[i]
            ][0]
            total = int(re.sub(r"\D", "", re.sub(r"\(.*", "", content[line])))
            line = [
                i
                for i in range(len(content))
                if " were paired; of these:" in content[i]
            ][0]
            paired_perc = float(
                re.search(r"\(.*%\)", content[line]).group()[1:-2]
            )
            line = [
                i
                for i in range(len(content))
                if "aligned concordantly 0 times" in content[i]
            ][0]
            concordant_unaligned_perc = float(
                re.search(r"\(.*%\)", content[line]).group()[1:-2]
            )
            line = [
                i
                for i in range(len(content))
                if "aligned concordantly exactly 1 time" in content[i]
            ][0]
            concordant_unique_perc = float(
                re.search(r"\(.*%\)", content[line]).group()[1:-2]
            )
            line = [
                i
                for i in range(len(content))
                if "aligned concordantly >1 times" in content[i]
            ][0]
            concordant_multiple_perc = float(
                re.search(r"\(.*%\)", content[line]).group()[1:-2]
            )
            line = [
                i
                for i in range(len(content))
                if "mates make up the pairs; of these:" in content[i]
            ][0]
            not_aligned_or_discordant = int(
                re.sub(r"\D", "", re.sub(r"\(.*", "", content[line]))
            )
            d_fraction = not_aligned_or_discordant / float(total)
            not_aligned_or_discordant_perc = d_fraction * 100
            line = [
                i
                for i in range(len(content))
                if "aligned 0 times\n" in content[i]
            ][0]
            not_aligned_perc = (
                float(re.search(r"\(.*%\)", content[line]).group()[1:-2])
                * d_fraction
            )
            line = [
                i
                for i in range(len(content))
                if " aligned exactly 1 time" in content[i]
            ][0]
            unique_aligned_perc = (
                float(re.search(r"\(.*%\)", content[line]).group()[1:-2])
                * d_fraction
            )
            line = [
                i
                for i in range(len(content))
                if " aligned >1 times" in content[i]
            ][0]
            multiple_aligned_perc = (
                float(re.search(r"\(.*%\)", content[line]).group()[1:-2])
                * d_fraction
            )
            line = [
                i
                for i in range(len(content))
                if "overall alignment rate" in content[i]
            ][0]
            perc_aligned = float(re.sub("%.*", "", content[line]).strip())
            return {
                prefix + "paired_perc": paired_perc,
                prefix + "concordant_unaligned_perc": concordant_unaligned_perc,
                prefix + "concordant_unique_perc": concordant_unique_perc,
                prefix + "concordant_multiple_perc": concordant_multiple_perc,
                prefix
                + "not_aligned_or_discordant_perc": not_aligned_or_discordant_perc,
                prefix + "not_aligned_perc": not_aligned_perc,
                prefix + "unique_aligned_perc": unique_aligned_perc,
                prefix + "multiple_aligned_perc": multiple_aligned_perc,
                prefix + "perc_aligned": perc_aligned,
            }
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
    import re

    error_dict = {
        prefix + "filtered_single_ends": pd.np.nan,
        prefix + "filtered_paired_ends": pd.np.nan,
        prefix + "duplicate_percentage": pd.np.nan,
    }
    try:
        with open(stats_file) as handle:
            content = handle.readlines()  # list of strings per line
    except:
        return error_dict

    try:
        line = [i for i in range(len(content)) if "single ends" in content[i]][
            0
        ]
        single_ends = int(re.sub(r"\D", "", re.sub(r"\(.*", "", content[line])))
        line = [i for i in range(len(content)) if " end pairs" in content[i]][0]
        paired_ends = int(
            re.sub(r"\D", "", re.sub(r"\.\.\..*", "", content[line]))
        )
        line = [i for i in range(len(content)) if " duplicates" in content[i]][
            0
        ]
        duplicates = int(
            re.sub(r"\D", "", re.sub(r"\.\.\..*", "", content[line]))
        )
        return {
            prefix + "filtered_single_ends": single_ends,
            prefix + "filtered_paired_ends": paired_ends,
            prefix
            + "duplicate_percentage": (
                float(duplicates) / (single_ends + paired_ends * 2)
            )
            * 100,
        }
    except IndexError:
        return error_dict


def parse_peak_number(peak_file, prefix=""):
    import subprocess

    try:
        return {
            prefix
            + "peaks": int(
                subprocess.check_output(["wc", "-l", peak_file])
                .decode()
                .strip()
                .split(" ")[0]
            )
        }
    except (TypeError, IndexError, subprocess.CalledProcessError):
        return {prefix + "peaks": pd.np.nan}


def calculate_frip(input_bam, input_bed, output, cpus=4):
    return "sambamba view -t {0} -c  -L {1}  {2} > {3}".format(
        cpus, input_bed, input_bam, output
    )


def parse_frip(frip_file, total_reads, prefix=""):
    """
    Calculates the fraction of reads in peaks for a given sample.

    :param frip_file: A sting path to a file with the FRiP output.
    :type frip_file: str
    :param total_reads: A Sample object with the "peaks" attribute.
    :type total_reads: int
    """
    import re

    error_dict = {prefix + "frip": pd.np.nan}
    try:
        with open(frip_file, "r") as handle:
            content = handle.readlines()
    except:
        return error_dict

    if content[0].strip() == "":
        return error_dict

    reads_in_peaks = int(re.sub(r"\D", "", content[0]))

    return {prefix + "frip": reads_in_peaks / float(total_reads)}


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


def run_tss_analysis(
    sample, bam_file, chrom_file, tss_file, read_length=50, slop_size=1000
):

    if not os.path.exists(bam_file):
        print("Bam File {} not found!".format(bam_file))
        return 1

    tss_bed = bedtools.BedTool(tss_file)
    tss_bed = tss_bed.slop(b=slop_size, g=chrom_file).sort(faidx=chrom_file)
    # sample = os.path.splitext(os.path.split(bam_file)[1])[0]

    alignments = (
        bedtools.BedTool(bam_file)
        .bam_to_bed()
        .shift(g=chrom_file, p=-read_length / 2, m=read_length / 2)
        .sort(faidx=chrom_file)
    )

    coverage = tss_bed.coverage(
        alignments, g=chrom_file, sorted=True, d=True, F=0.5
    )

    histogram = coverage.to_dataframe(
        names=["chrom", "start", "end", "gene", "X", "strand", "base", "count"],
        usecols=["base", "count", "strand"],
    )
    histogram.loc[histogram["strand"] == "+", "base"] = (
        histogram.loc[histogram["strand"] == "+", "base"] - slop_size - 1
    )
    histogram.loc[histogram["strand"] == "-", "base"] = -(
        histogram.loc[histogram["strand"] == "-", "base"] - slop_size - 1
    )
    histogram = (
        histogram[["base", "count"]]
        .sort_values(by=["base"])
        .groupby("base")
        .sum()
    )

    noise = (histogram[:100].sum() + histogram[-100:].sum()) / 200
    normalized_histogram = histogram / noise
    normalized_histogram.to_csv(sample.tss_hist)

    fig, ax = plt.subplots(1, 1)
    ax.plot(normalized_histogram.index, normalized_histogram, color="k")
    ax.axvline(0, linestyle=":", color="k")

    ax.set_xlabel("Distance from TSS (bp)")
    ax.set_ylabel("Normalized coverage")
    fig.savefig(sample.tss_plot)
    plt.close(fig)

    open(sample.tss_lock, "w")

    enr = normalized_histogram.max()["count"]
    print("TSS enrichment: {}".format(enr))
    return enr


def filter_peaks(peaks, exclude, filtered_peaks):
    return "bedtools intersect -v -wa -a {} -b {} > {}".format(
        peaks, exclude, filtered_peaks
    )


def bam_to_bigwig(
    input_bam, output_bigwig, genome, normalization_method="RPGC"
):
    from collections import defaultdict

    if genome not in ["hg19", "hg38", "mm10", "mm9"]:
        print(
            "Genome assembly is not known. Using size of human genome. Beware."
        )

    genome_size = defaultdict(lambda: 3300000000)
    for g in ["mm9", "mm10"]:
        genome_size[g] = 2800000000

    cmd = "bamCoverage --bam {bam_file} -o {bigwig}"
    cmd += (
        " -p max --binSize 10  --normalizeUsing {norm} --effectiveGenomeSize {genome_size} --extendReads 175"
        ""
    )
    cmd = cmd.format(
        bam_file=input_bam,
        bigwig=output_bigwig,
        norm=normalization_method,
        genome_size=genome_size[genome],
    )
    return cmd


if __name__ == "__main__":
    try:
        sys.exit(main())
    except KeyboardInterrupt:
        print("Program canceled by user!")
        sys.exit(1)
