#!/usr/bin/env python

"""
ATAC-seq pipeline
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
__copyright__ = "Copyright 2015, Andre Rendeiro"
__credits__ = []
__license__ = "GPL2"
__version__ = "0.2"
__maintainer__ = "Andre Rendeiro"
__email__ = "arendeiro@cemm.oeaw.ac.at"
__status__ = "Development"


class ATACseqSample(Sample):
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
        super(ATACseqSample, self).__init__(series)

        self.tagmented = True

    def __repr__(self):
        return "ATAC-seq sample '%s'" % self.sample_name

    def set_file_paths(self):
        """
        Sets the paths of all files for this sample.
        """
        # Inherit paths from Sample by running Sample's set_file_paths()
        super(ATACseqSample, self).set_file_paths()

        # Files in the root of the sample dir
        sample_root = self.paths.sample_root
        root = os.path.join(sample_root, self.sample_name)
        self.fastqc = root + ".fastqc.zip"
        self.trimlog = root + ".trimlog.txt"
        self.aln_rates = root + ".aln_rates.txt"
        self.aln_metrics = root + ".aln_metrics.txt"
        self.dups_metrics = root + ".dups_metrics.txt"

        # Unmapped: merged bam, fastq, trimmed fastq
        self.paths.unmapped = os.path.join(sample_root, "unmapped")
        unmapped = os.path.join(self.paths.unmapped, self.sample_name)
        self.unmapped = unmapped + ".bam"
        self.fastq = unmapped + ".fastq"
        self.fastqc = root + "_fastqc.zip"
        self.fastq1 = unmapped + "_R1.fastq"
        self.fastqc1 = root + "_R1_fastqc.zip"
        self.fastq2 = unmapped + "_R2.fastq"
        self.fastqc2 = root + "_R2_fastqc.zip"
        self.fastq_unpaired = unmapped + ".unpaired.fastq"
        self.trimmed = unmapped + ".trimmed.fastq"
        self.trimmed1 = unmapped + "_R1.trimmed.fastq"
        self.trimmed2 = unmapped + "_R2.trimmed.fastq"
        self.trimmed1_unpaired = unmapped + "_R1.trimmed.unpaired.fastq"
        self.trimmed2_unpaired = unmapped + "_R2.trimmed.unpaired.fastq"

        # Mapped: mapped, duplicates marked, removed, reads shifted
        self.paths.mapped = os.path.join(sample_root, "mapped")
        mapped = os.path.join(self.paths.mapped, self.sample_name)
        self.mapped = mapped + ".trimmed.bowtie2.bam"
        self.filtered = mapped + ".trimmed.bowtie2.filtered.bam"
        # this will create additional bam files with reads shifted
        self.filteredshifted = mapped + ".trimmed.bowtie2.filtered.shifted.bam"

        # Files in the root of the sample dir
        self.frip = root + "_FRiP.txt"
        self.oracle_frip = root + "_oracle_FRiP.txt"

        # Mapped: mapped, duplicates marked, removed, reads shifted
        # this will create additional bam files with reads shifted
        self.filteredshifted = mapped + ".trimmed.bowtie2.filtered.shifted.bam"

        # Coverage: read coverage in windows genome-wide
        self.paths.coverage = os.path.join(sample_root, "coverage")
        coverage = os.path.join(self.paths.coverage, self.sample_name)
        self.coverage = coverage + ".cov"
        self.bigwig = coverage + ".bigWig"

        self.insertplot = root + "_insertLengths.pdf"
        self.insertdata = root + "_insertLengths.csv"
        self.mitochondrial_stats = root + "_mitochondrial_stats.tsv"
        self.qc = root + "_qc.tsv"
        self.qc_plot = root + "_qc.pdf"

        # Peaks: peaks called and derivate files
        self.paths.peaks = os.path.join(sample_root, "peaks")
        peaks = os.path.join(self.paths.peaks, self.name)
        self.peaks = peaks + "_peaks.narrowPeak"
        self.summits = peaks + "_summits.bed"
        self.filtered_peaks = peaks + "_peaks.filtered.bed"


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

    def set_file_paths(self):
        super(DNaseSample, self).set_file_paths()


def main():
    # Parse command-line arguments
    parser = ArgumentParser(
        prog="python atacseq.py",
        description="\n".join([
            "ATAC-seq pipeline from Bock lab. See https://github.com/epigen/open_pipelines and " +
            "https://github.com/epigen/open_pipelines/blob/master/pipelines/atacseq.md for specific documentation. ",
            "Input options should be either a sample config file or 'input', 'sample_name' and 'output_parent'."])
    )
    parser = arg_parser(parser)
    parser = pypiper.add_pypiper_args(parser, groups=["all"])
    args = parser.parse_args()
    if (
        (args.sample_config is None) and
        (
            (args.input is None) or
            (args.sample_name is None) or
            (args.output_parent is None) or
            (args.genome_assembly is None) or
            (args.single_or_paired is None))):
        print("No sample config or (input files, sample name, genome assembly and output directory) provided!\n")
        parser.print_help()
        return 1

    # Read in yaml configs
    if (args.sample_config is not None):
        series = pd.Series(yaml.load(open(args.sample_config, "r")))
    else:
        if args.input2 is not None:
            if len(args.input) != len(args.input2):
                raise ValueError(
                    "Provided number of input files does not match!")
        else:
            args.input2 = []

        series = pd.Series({
            "sample_name": args.sample_name,
            "input": args.input,
            "input2": args.input2,
            "output_parent": args.output_parent,
            "genome": args.genome_assembly,
            "protocol": args.protocol})

    # looper 0.6/0.7 compatibility:
    if "protocol" in series.index:
        key = "protocol"
    elif "library" in series.index:
        key = "library"
    else:
        raise KeyError(
            "Sample does not contain either a 'protocol' or 'library' attribute!")

    # Create Sample object
    if series[key] != "DNase-seq":
        sample = ATACseqSample(series)
    else:
        sample = DNaseSample(series)

    # Standardize input files
    # overwrite config values with args if provided
    print(args)
    if args.sample_config is not None:
        if sample.read1 is not None:
            sample.read1_files = sample.read1.split(" ")
        if sample.read2 is not None:
            sample.read2_files = sample.read2.split(" ")
    if args.input is not None:
        sample.read1_files = args.input
    if args.input2 is not None:
        sample.read2_files = args.input2

    # Check if merged
    if not hasattr(sample, "merged"):
        sample.merged = len(sample.read1_files) > 1
    # Check read type if not provided
    if not hasattr(sample, "ngs_inputs"):
        sample.ngs_inputs = sample.read1_files + sample.read2_files
    if not hasattr(sample, "read_type"):
        try:
            sample.set_read_type(permissive=False)
        except NotImplementedError:
            if len(sample.read2_files) > 0:
                sample.read_type = "paired"
    print(getattr(sample, "read_type"))

    # sample.read_type = "single" if (sample.read2_files is None) else "paired"
    sample.paired = sample.read_type == "paired"

    # Set paths if no sample config
    if args.sample_config is None:
        sample.prj = AttributeDict({'metadata': AttributeDict()})
        sample.prj.metadata['results_subdir'] = os.path.abspath(args.output_parent)
        sample.paths = dict({"sample_root": os.path.abspath(os.path.join(args.output_parent, args.sample_name))})
    sample.paths = AttributeDict(sample.paths.__dict__)

    # Set file paths
    sample.set_file_paths()


    # Start Pypiper object
    # Best practice is to name the pipeline with the name of the script;
    # or put the name in the pipeline interface.
    pipe_manager = pypiper.PipelineManager(name="atacseq", outfolder=sample.paths.sample_root, args=args)
    pipe_manager.config.tools.scripts_dir = os.path.join(os.path.dirname(os.path.realpath(__file__)), "tools")


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
    choices = ["ATAC-seq", "DNase-seq"]
    parser.add_argument(
        "-p", "--protocol",
        dest="protocol",
        default=choices[0],
        choices=choices,
        help="Experimental protocol of the sample. " +
             "One of ['{}']. Default is '{}'.".format(
                "', '".join(choices), choices[0]),
        type=str
    )
    return parser


def process(sample, pipe_manager, args):
    """
    This takes unmapped Bam files and makes trimmed, aligned, duplicate marked
    and removed, indexed, shifted Bam files along with a UCSC browser track.
    Peaks are called and filtered.
    """
    import glob

    print("Start processing ATAC-seq sample %s." % sample.sample_name)

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

    # Merge input files if needed handling various input formats
    pipe_manager.timestamp("Merging or linking input files")
    local_input_files = tk.merge_or_link(
        input_args=[sample.read1_files, sample.read2_files],
        raw_folder=sample.paths.unmapped)

    # Convert to FASTQ
    pipe_manager.timestamp("Converting to FASTQ format")
    cmd, out_fastq_prefix, unaligned_fastq = tk.input_to_fastq(
        input_file=local_input_files,
        sample_name=sample.sample_name,
        paired_end=sample.paired,
        fastq_folder=sample.paths.unmapped)
    pipe_manager.run(cmd, unaligned_fastq)
    pipe_manager.clean_add(out_fastq_prefix + "*.fastq", conditional=True)

    # Run FASTQC
    for (prefix, fastq_file, fastqc_file) in zip(
        ["fastq_", "fastq_R1_", "fastq_R2_"],
        [out_fastq_prefix + ".fastq", out_fastq_prefix + "_R1.fastq", out_fastq_prefix + "_R2.fastq"],
        [sample.fastqc, sample.fastqc1, sample.fastqc2]):
        if not os.path.exists(fastq_file):
            continue

        pipe_manager.timestamp("Running FastQC in file: '{}'".format(fastq_file))
        cmd = tk.fastqc(fastq_file, sample.paths.sample_root)
        pipe_manager.run(cmd, fastqc_file)
        report_dict(pipe_manager, parse_fastqc(fastqc_file, prefix=prefix))

    # report joint metrics
    m = pd.Series(pipe_manager.stats_dict).astype(float)
    for metric in ["total_pass_filter_reads"]:
        pipe_manager.report_result("fastqc_total_" + metric, m[m.index.str.contains(metric)].sum())
    for metric in ["read_length", "GC_perc", "poor_quality_perc"]:
        pipe_manager.report_result("fastqc_total_" + metric, m[m.index.str.contains(metric)].mean())

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
            inputFastq1=sample.fastq1 if sample.paired else sample.fastq,
            inputFastq2=sample.fastq2 if sample.paired else None,
            outputPrefix=os.path.join(sample.paths.unmapped, sample.sample_name),
            outputFastq1=sample.trimmed1 if sample.paired else sample.trimmed,
            outputFastq2=sample.trimmed2 if sample.paired else None,
            trimLog=sample.trimlog,
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
    cmd = tk.bowtie2Map(
        inputFastq1=sample.trimmed1 if sample.paired else sample.trimmed,
        inputFastq2=sample.trimmed2 if sample.paired else None,
        outputBam=sample.mapped,
        log=sample.aln_rates,
        metrics=sample.aln_metrics,
        genomeIndex=getattr(pipe_manager.config.resources.genomes, sample.genome),
        maxInsert=pipe_manager.config.parameters.max_insert,
        cpus=args.cores
    )
    pipe_manager.run(cmd, sample.mapped, shell=True)
    report_dict(pipe_manager, parse_mapping_stats(sample.aln_rates, paired_end=sample.paired))

    # Get mitochondrial reads
    pipe_manager.timestamp("Getting mitochondrial stats")
    cmd = tk.get_mitochondrial_reads(
        bam_file=sample.mapped,
        output=sample.mitochondrial_stats,
        cpus=args.cores
    )
    pipe_manager.run(cmd, sample.mitochondrial_stats, shell=True, nofail=True)
    report_dict(pipe_manager, parse_duplicate_stats(sample.mitochondrial_stats, prefix="MT_"))

    # Filter reads
    pipe_manager.timestamp("Filtering reads for quality")
    cmd = tk.filterReads(
        inputBam=sample.mapped,
        outputBam=sample.filtered,
        metricsFile=sample.dups_metrics,
        paired=sample.paired,
        cpus=args.cores,
        Q=pipe_manager.config.parameters.read_quality
    )
    pipe_manager.run(cmd, sample.filtered, shell=True)
    report_dict(pipe_manager, parse_duplicate_stats(sample.dups_metrics))

    # Shift reads
    if sample.tagmented:
        pipe_manager.timestamp("Shifting reads of tagmented sample")
        cmd = tk.shiftReads(
            inputBam=sample.filtered,
            genome=sample.genome,
            outputBam=sample.filteredshifted
        )
        pipe_manager.run(cmd, sample.filteredshifted, shell=True)

    # Index bams
    pipe_manager.timestamp("Indexing bamfiles with samtools")
    cmd = tk.indexBam(inputBam=sample.mapped)
    pipe_manager.run(cmd, sample.mapped + ".bai", shell=True)
    cmd = tk.indexBam(inputBam=sample.filtered)
    pipe_manager.run(cmd, sample.filtered + ".bai", shell=True)
    if sample.tagmented:
        cmd = tk.indexBam(inputBam=sample.filteredshifted)
        pipe_manager.run(cmd, sample.filteredshifted + ".bai", shell=True)

    # Make tracks
    # right now tracks are only made for bams without duplicates
    pipe_manager.timestamp("Making bigWig tracks from bam file")
    cmd = bamToBigWig(
        inputBam=sample.filtered,
        outputBigWig=sample.bigwig,
        genomeSizes=getattr(pipe_manager.config.resources.chromosome_sizes, sample.genome),
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
        pipe_manager.report_figure("insert_sizes", sample.insertplot)

    # Count coverage genome-wide
    pipe_manager.timestamp("Calculating genome-wide coverage")
    cmd = tk.genomeWideCoverage(
        inputBam=sample.filtered,
        genomeWindows=getattr(pipe_manager.config.resources.genome_windows, sample.genome),
        output=sample.coverage
    )
    pipe_manager.run(cmd, sample.coverage, shell=True)

    # Calculate NSC, RSC
    pipe_manager.timestamp("Assessing signal/noise in sample")
    cmd = tk.peakTools(
        inputBam=sample.filtered,
        output=sample.qc,
        plot=sample.qc_plot,
        cpus=args.cores
    )
    pipe_manager.run(cmd, sample.qc_plot, shell=True, nofail=True)
    report_dict(pipe_manager, parse_nsc_rsc(sample.qc))
    pipe_manager.report_figure("cross_correlation", sample.qc_plot)

    # Call peaks
    pipe_manager.timestamp("Calling peaks with MACS2")
    # make dir for output (macs fails if it does not exist)
    if not os.path.exists(sample.paths.peaks):
        os.makedirs(sample.paths.peaks)

    cmd = tk.macs2CallPeaksATACSeq(
        treatmentBam=sample.filtered,
        outputDir=sample.paths.peaks,
        sampleName=sample.sample_name,
        genome=sample.genome
    )
    pipe_manager.run(cmd, sample.peaks, shell=True)
    report_dict(pipe_manager, parse_peak_number(sample.peaks))

    # Filter peaks
    if hasattr(pipe_manager.config.resources.blacklisted_regions, sample.genome):
        pipe_manager.timestamp("Filtering peaks from blacklisted regions")
        cmd = filter_peaks(
            peaks=sample.peaks,
            exclude=getattr(pipe_manager.config.resources.blacklisted_regions, sample.genome),
            filtered_peaks=sample.filtered_peaks
        )
        pipe_manager.run(cmd, sample.filtered_peaks, shell=True)
        report_dict(pipe_manager, parse_peak_number(sample.filtered_peaks, prefix="filtered_"))

    # Calculate fraction of reads in peaks (FRiP)
    pipe_manager.timestamp("Calculating fraction of reads in peaks (FRiP)")
    # on the sample's peaks
    cmd = tk.calculate_FRiP(
        inputBam=sample.filtered,
        inputBed=sample.peaks,
        output=sample.frip,
        cpus=args.cores
    )
    pipe_manager.run(cmd, sample.frip, shell=True)
    total = (float(pipe_manager.stats_dict["filtered_single_ends"]) + (float(pipe_manager.stats_dict["filtered_paired_ends"]) / 2.))
    report_dict(pipe_manager, parse_FRiP(sample.frip, total))

    # on an oracle peak list
    if hasattr(pipe_manager.config.resources.oracle_peak_regions, sample.genome):
        cmd = tk.calculate_FRiP(
            inputBam=sample.filtered,
            inputBed=getattr(pipe_manager.config.resources.oracle_peak_regions, sample.genome),
            output=sample.oracle_frip,
            cpus=args.cores
        )
        pipe_manager.run(cmd, sample.oracle_frip, shell=True)
        report_dict(pipe_manager, parse_FRiP(sample.oracle_frip, total, prefix="oracle_"))

    # Finish up
    print(pipe_manager.stats_dict)

    pipe_manager.stop_pipeline()
    print("Finished processing sample %s." % sample.sample_name)


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
        content = StringIO.StringIO(zfile.read(os.path.join(zfile.filelist[0].filename, "fastqc_data.txt"))).readlines()
    except:
        return error_dict
    try:
        line = [i for i in range(len(content)) if "Total Sequences\t" in content[i]][0]
        total = int(re.sub("\D", "", re.sub("\(.*", "", content[line])))
        line = [i for i in range(len(content)) if "Sequences flagged as poor quality\t" in content[i]][0]
        poor_quality = int(re.sub("\D", "", re.sub("\(.*", "", content[line])))
        line = [i for i in range(len(content)) if "Sequence length\t" in content[i]][0]
        seq_len = int(re.sub("\D", "", re.sub(" \(.*", "", content[line]).strip()))
        line = [i for i in range(len(content)) if "%GC\t" in content[i]][0]
        gc_perc = int(re.sub("\D", "", re.sub(" \(.*", "", content[line]).strip()))
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
    import re

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
    import re

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


def parse_peak_number(peak_file, prefix=""):
    from subprocess import check_output
    try:
        return {prefix + "peaks": int(check_output(["wc", "-l", peak_file]).split(" ")[0])}
    except:
        return {prefix + "peaks": pd.np.nan}


def parse_FRiP(frip_file, total_reads, prefix=""):
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

    reads_in_peaks = int(re.sub("\D", "", content[0]))

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


def bamToBigWig(inputBam, outputBigWig, genomeSizes, genome, tagmented=False, normalize=False, norm_factor=1000000):
    import os
    import re
    # TODO:
    # addjust fragment length dependent on read size and real fragment size
    # (right now it asssumes 50bp reads with 180bp fragments)
    cmds = list()
    transientFile = os.path.abspath(re.sub("\.bigWig", "", outputBigWig))
    cmd1 = "bedtools bamtobed -i {0} |".format(inputBam)
    if not tagmented:
        cmd1 += " " + "bedtools slop -i stdin -g {0} -s -l 0 -r 130 |".format(genomeSizes)
        cmd1 += " fix_bedfile_genome_boundaries.py {0} |".format(genome)
    cmd1 += " " + "genomeCoverageBed {0}-bg -g {1} -i stdin > {2}.cov".format(
        "-5 " if tagmented else "",
        genomeSizes,
        transientFile
    )
    cmds.append(cmd1)
    if normalize:
        cmds.append("""awk 'NR==FNR{{sum+= $4; next}}{{ $4 = ($4 / sum) * {1}; print}}' {0}.cov {0}.cov | sort -k1,1 -k2,2n > {0}.normalized.cov""".format(transientFile, norm_factor))
    cmds.append("bedGraphToBigWig {0}{1}.cov {2} {3}".format(transientFile, ".normalized" if normalize else "", genomeSizes, outputBigWig))
    # remove tmp files
    cmds.append("if [[ -s {0}.cov ]]; then rm {0}.cov; fi".format(transientFile))
    if normalize:
        cmds.append("if [[ -s {0}.normalized.cov ]]; then rm {0}.normalized.cov; fi".format(transientFile))
    cmds.append("chmod 755 {0}".format(outputBigWig))
    return cmds


def filter_peaks(peaks, exclude, filtered_peaks):
    return "bedtools intersect -v -wa -a {} -b {} > {}".format(
        peaks, exclude, filtered_peaks)


if __name__ == '__main__':
    try:
        sys.exit(main())
    except KeyboardInterrupt:
        print("Program canceled by user!")
        sys.exit(1)
