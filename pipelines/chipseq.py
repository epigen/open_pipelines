#!/usr/bin/env python

"""
ChIP-seq pipeline
"""

import os
import sys
from argparse import ArgumentParser
import yaml
import pypiper
from pypiper.ngstk import NGSTk
from peppy import AttributeDict, Sample

import pandas as pd


__author__ = "Andre Rendeiro"
__copyright__ = "Copyright 2015, Andre Rendeiro"
__credits__ = []
__license__ = "GPL2"
__version__ = "0.3"
__maintainer__ = "Andre Rendeiro"
__email__ = "arendeiro@cemm.oeaw.ac.at"
__status__ = "Development"


class ChIPseqSample(Sample):
    """
    Class to model ChIP-seq samples based on the generic Sample class.

    :param series: Pandas `Series` object.
    :type series: pandas.Series
    """
    __library__ = "ChIP-seq"

    def __init__(self, series):
        # Passed series must either be a pd.Series or a daughter class
        if not isinstance(series, pd.Series):
            raise TypeError("Provided object is not a pandas Series.")
        super(ChIPseqSample, self).__init__(series)

        self.tagmented = False

        # Get type of factor
        # TODO: get config file specifying broad/narrow factors
        # e.g. self.broad = True if self.ip in self.prj.config["broadfactors"] else False
        # NS: I wrapped this in a try block because this makes it require that
        # 'ip' be defined, which may not be the case (like for Input samples)

        try:
            self.broad = True if any([ip in self.ip.upper() for ip in ["H3K27me3", "H3K36me3"]]) else False
            self.histone = True if any([ip in self.ip.upper() for ip in ["H3", "H2A", "H2B", "H4"]]) else False
        except AttributeError:
            pass

    def __repr__(self):
        return "ChIP-seq sample '%s'" % self.sample_name

    def set_file_paths(self):
        """
        Sets the paths of all files for this sample.
        """
        # Inherit paths from Sample by running Sample's set_file_paths()
        super(ChIPseqSample, self).set_file_paths()

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

        # Peaks: peaks called and derivate files
        self.paths.peaks = os.path.join(self.paths.sample_root, "peaks")
        self.peaks = os.path.join(self.paths.peaks, self.sample_name + ("_peaks.narrowPeak" if not self.broad else "_peaks.broadPeak"))
        self.r_peak_model_file = os.path.join(self.paths.peaks, self.sample_name + "_model.r")


class ChIPmentation(ChIPseqSample):
    """
    Class to model ChIPmentation samples based on the ChIPseqSample class.

    :param series: Pandas `Series` object.
    :type series: pandas.Series
    """
    __library__ = "ChIPmentation"

    def __init__(self, series):
        # Use _pd.Series object to have all sample attributes
        if not isinstance(series, pd.Series):
            raise TypeError("Provided object is not a pandas Series.")
        super(ChIPmentation, self).__init__(series)

        self.tagmented = True

    def __repr__(self):
        return "ChIPmentation sample '%s'" % self.sample_name

    def set_file_paths(self):
        super(ChIPmentation, self).set_file_paths()


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
        line = [i for i in range(len(content)) if "Total Sequences" in content[i]][0]
        total = int(re.sub("\D", "", re.sub("\(.*", "", content[line])))
        line = [i for i in range(len(content)) if "Sequences flagged as poor quality" in content[i]][0]
        poor_quality = int(re.sub("\D", "", re.sub("\(.*", "", content[line])))
        line = [i for i in range(len(content)) if "Sequence length	" in content[i]][0]
        seq_len = int(re.sub("\D", "", re.sub(" \(.*", "", content[line]).strip()))
        line = [i for i in range(len(content)) if "%GC" in content[i]][0]
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

    stats_dict = {
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
        return stats_dict

    suf = "s" if not paired_end else " pairs"

    try:
        line = [i for i in range(len(content)) if "read{} processed; of these:".format(suf) in content[i]][0]
        total = int(re.sub("\D", "", re.sub("\(.*", "", content[line])))
    except IndexError:
        return stats_dict
    try:
        line = [i for i in range(len(content)) if "read{} available; of these:".format(suf) in content[i]][0]
        surviving = int(re.sub("\D", "", re.sub("\(.*", "", content[line])))
        stats_dict[prefix + "surviving_perc"] = (float(surviving) / total) * 100
        stats_dict[prefix + "trim_loss_perc"] = ((total - float(surviving)) / total) * 100
    except IndexError:
        pass
    try:
        line = [i for i in range(len(content)) if "short read{} filtered out after trimming by size control".format(suf) in content[i]][0]
        short = int(re.sub(" \(.*", "", content[line]).strip())
        stats_dict[prefix + "short_perc"] = (float(short) / total) * 100
    except IndexError:
        pass
    try:
        line = [i for i in range(len(content)) if "empty read{} filtered out after trimming by size control".format(suf) in content[i]][0]
        empty = int(re.sub(" \(.*", "", content[line]).strip())
        stats_dict[prefix + "empty_perc"] = (float(empty) / total) * 100
    except IndexError:
        pass
    try:
        line = [i for i in range(len(content)) if "trimmed read{} available after processing".format(suf) in content[i]][0]
        trimmed = int(re.sub(" \(.*", "", content[line]).strip())
        stats_dict[prefix + "trimmed_perc"] = (float(trimmed) / total) * 100
    except IndexError:
        pass
    try:
        line = [i for i in range(len(content)) if "untrimmed read{} available after processing".format(suf) in content[i]][0]
        untrimmed = int(re.sub(" \(.*", "", content[line]).strip())
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
        line = [i for i in range(len(content)) if "single ends" in content[i]][0]
        single_ends = int(re.sub("\D", "", re.sub("\(.*", "", content[line])))
        line = [i for i in range(len(content)) if " end pairs" in content[i]][0]
        paired_ends = int(re.sub("\D", "", re.sub("\.\.\..*", "", content[line])))
        line = [i for i in range(len(content)) if " duplicates" in content[i]][0]
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


def parse_FRiP(frip_file, total_reads):
    """
    Calculates the fraction of reads in peaks for a given sample.

    :param frip_file: A sting path to a file with the FRiP output.
    :type frip_file: str
    :param total_reads: A Sample object with the "peaks" attribute.
    :type total_reads: int
    """
    import re

    error_dict = {"frip": pd.np.nan}
    try:
        with open(frip_file, "r") as handle:
            content = handle.readlines()
    except:
        return error_dict

    if content[0].strip() == "":
        return error_dict

    reads_in_peaks = int(re.sub("\D", "", content[0]))

    return {"frip": reads_in_peaks / float(total_reads)}


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

    # looper 0.6/0.7 compatibility:
    if "protocol" in series.index:
        key = "protocol"
    elif "library" in series.index:
        key = "library"
    else:
        raise KeyError(
            "Sample does not contain either a 'protocol' or 'library' attribute!")

    # Create Sample object
    if series[key] != "ChIPmentation":
        sample = ChIPseqSample(series)
    else:
        sample = ChIPmentation(series)

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
    else:
        if sample.read_type not in ['single', 'paired']:
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
    pipe_manager = pypiper.PipelineManager(name="chipseq", outfolder=sample.paths.sample_root, args=args)

    # Start main function
    if not args.only_peaks:
        pipe_manager = process(sample, pipe_manager, args)
    else:
        print("Skipped processing sample '{}'.".format(sample.name))

    # If sample does not have "ctrl" attribute, finish processing it.
    if not hasattr(sample, "compare_sample"):
        pipe_manager.stop_pipeline()
        print("Finished processing sample '{}'.".format(sample.name))
        return

    # The pipeline will now wait for the comparison sample file to be completed
    pipe_manager._wait_for_file(sample.filtered.replace(sample.name, sample.compare_sample))

    # Start peak calling function
    call_peaks(sample, pipe_manager, args)


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
        "-p", "--only-call-peaks",
        dest="only_peaks",
        help="Skip preprocessing and only call peaks.",
        action="store_true",
    )
    return parser


def process(sample, pipe_manager, args):
    """
    This takes unmapped Bam files and makes trimmed, aligned, duplicate marked
    and removed, indexed (and shifted if necessary) Bam files
    along with a UCSC browser track.
    """
    print("Start processing ChIP-seq sample '{}'.".format(sample.name))

    for path in ["sample_root"] + list(sample.paths.__dict__.keys()):
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
            input_bams=sample.data_path.split(" "),  # this is a list of sample paths
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
    pipe_manager.run(cmd, os.path.join(sample.paths.sample_root, sample.sample_name + "_fastqc.zip"), shell=True)
    report_dict(pipe_manager, parse_fastqc(os.path.join(sample.paths.sample_root, sample.sample_name + "_fastqc.zip"), prefix="fastqc_"))

    # Convert bam to fastq
    pipe_manager.timestamp("Converting to Fastq format")
    cmd = tk.bam2fastq(
        inputBam=sample.data_path,
        outputFastq=sample.fastq1 if sample.paired else sample.fastq,
        outputFastq2=sample.fastq2 if sample.paired else None,
        unpairedFastq=sample.fastq_unpaired if sample.paired else None
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

    # Index bams
    pipe_manager.timestamp("Indexing bamfiles with samtools")
    cmd = tk.indexBam(inputBam=sample.mapped)
    pipe_manager.run(cmd, sample.mapped + ".bai", shell=True)
    cmd = tk.indexBam(inputBam=sample.filtered)
    pipe_manager.run(cmd, sample.filtered + ".bai", shell=True)

    track_dir = os.path.dirname(sample.bigwig)
    if not os.path.exists(track_dir):
        os.makedirs(track_dir)

    # Report total efficiency
    usable = (
        float(pipe_manager.stats_dict["filtered_single_ends"]) +
        (float(pipe_manager.stats_dict["filtered_paired_ends"]) / 2.))
    total = float(pipe_manager.stats_dict['fastqc_total_pass_filter_reads'])
    report_dict(
        pipe_manager,
        {"total_efficiency": (usable / total) * 100})

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

    print("Finished processing sample '{}'.".format(sample.name))
    return pipe_manager


def call_peaks(sample, pipe_manager, args):
    """
    Calls peaks for a sample give a control sample specified as a `compare_sample` attribute.
    """
    print("Start calling peaks for sample '{}'.".format(sample.name))

    for path in ["sample_root"] + list(sample.paths.__dict__.keys()):
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

    if pipe_manager.config.parameters.peak_caller == "macs2":
        pipe_manager.timestamp("Calling peaks with MACS2")
        # make dir for output (macs fails if it does not exist)
        if not os.path.exists(sample.paths.peaks):
            os.makedirs(sample.paths.peaks)

        # For point-source factors use default settings
        # For broad factors use broad settings
        cmd = tk.macs2CallPeaks(
            treatmentBams=sample.filtered,
            controlBams=sample.filtered.replace(sample.name, sample.compare_sample),
            outputDir=sample.paths.peaks,
            sampleName=sample.name,
            genome=sample.genome,
            broad=sample.broad,
            paired=sample.paired
        )
        pipe_manager.run(cmd, sample.peaks, shell=True)
        report_dict(pipe_manager, parse_peak_number(sample.peaks))

        if not sample.broad:
            pipe_manager.timestamp("Ploting MACS2 model")
            cmd = tk.macs2PlotModel(
                r_peak_model_file=sample.r_peak_model_file,
                sample_name=sample.name,
                output_dir=sample.paths.peaks
            )
            pipe_manager.run(cmd, os.path.join(sample.paths.peaks, sample.name + "_model.pdf"), shell=True, nofail=True)

    elif pipe_manager.config.parameters.peak_caller == "spp":
        pipe_manager.timestamp("Calling peaks with spp")
        # For point-source factors use default settings
        # For broad factors use broad settings
        cmd = tk.sppCallPeaks(
            treatmentBam=sample.filtered,
            controlBam=sample.filtered.replace(sample.name, sample.compare_sample),
            treatmentName=sample.name,
            controlName=sample.compare_sample,
            outputDir=os.path.join(sample.paths.peaks, sample.name),
            broad=sample.broad,
            cpus=args.cpus
        )
        pipe_manager.run(cmd, sample.peaks, shell=True)

    # Calculate fraction of reads in peaks (FRiP)
    pipe_manager.timestamp("Calculating fraction of reads in peaks (FRiP)")
    cmd = tk.calculate_FRiP(
        inputBam=sample.filtered,
        inputBed=sample.peaks,
        output=sample.frip,
        cpus=args.cores
    )
    pipe_manager.run(cmd, sample.frip, shell=True)
    total = (
        float(pipe_manager.stats_dict["filtered_single_ends"]) +
        (float(pipe_manager.stats_dict["filtered_paired_ends"]) / 2.))
    report_dict(pipe_manager, parse_FRiP(sample.frip, total))

    print("Finished calling peaks for sample '{}'.".format(sample.name))
    pipe_manager.stop_pipeline()


if __name__ == '__main__':
    try:
        sys.exit(main())
    except KeyboardInterrupt:
        print("Program canceled by user!")
        sys.exit(1)
