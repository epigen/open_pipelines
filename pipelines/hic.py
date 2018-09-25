#!/usr/bin/env python

"""
Hi-C pipeline
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
__version__ = "0.2"
__maintainer__ = "Andre Rendeiro"
__email__ = "arendeiro@cemm.oeaw.ac.at"
__status__ = "Development"


class HiCSample(Sample):
    """
    Class to model Hi-C samples based on the ChIPseqSample class.

    :param series: Pandas `Series` object.
    :type series: pandas.Series
    """
    __library__ = "Hi-C"

    def __init__(self, series):

        # Use pd.Series object to have all sample attributes
        if not isinstance(series, pd.Series):
            raise TypeError("Provided object is not a pandas Series.")
        super(HiCSample, self).__init__(series)

        self.tagmented = True

    def __repr__(self):
        return "Hi-C sample '%s'" % self.sample_name

    def set_file_paths(self):
        """
        Sets the paths of all files for this sample.
        """
        # Inherit paths from Sample by running Sample's set_file_paths()
        super(HiCSample, self).set_file_paths()

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

        # HiC-Pro
        self.paths.hicpro = os.path.join(self.paths.sample_root, "hic-pro")
        self.hicpro_config = os.path.join(self.paths.hicpro, "hic-pro.sample_config.txt")


class HiChIPSample(HiCSample):
    """
    Class to model HiChIP samples based on the ChIPseqSample class.

    :param series: Pandas `Series` object.
    :type series: pandas.Series
    """
    __library__ = "HiChIP"

    def __init__(self, series):

        # Use pd.Series object to have all sample attributes
        if not isinstance(series, pd.Series):
            raise TypeError("Provided object is not a pandas Series.")
        super(HiChIPSample, self).__init__(series)

    def __repr__(self):
        return "HiChIP sample '%s'" % self.sample_name

    def set_file_paths(self):
        super(HiChIPSample, self).set_file_paths()


def main():
    # Parse command-line arguments
    parser = ArgumentParser(
        prog="hic-pipeline",
        description="Hi-C pipeline."
    )
    parser = arg_parser(parser)
    parser = pypiper.add_pypiper_args(parser, groups=["all"])
    args = parser.parse_args()

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
    if series[key] != "HiChIP":
        sample = HiCSample(series)
    else:
        sample = HiChIPSample(series)

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
    pipe_manager = pypiper.PipelineManager(name="hic", outfolder=sample.paths.sample_root, args=args)
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
    parser.add_argument(
        "-s", "--serial",
        dest="serial",
        help="Whether HiCPro pipeline should run in a serial manner.",
        action="store_true"
    )
    return parser


def process(sample, pipe_manager, args):
    """
    This takes unmapped Bam files and makes trimmed, aligned, duplicate marked
    and removed, indexed, shifted Bam files along with a UCSC browser track.
    Peaks are called and filtered.
    """
    import textwrap

    print("Start processing Hi-C sample %s." % sample.sample_name)

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
            input_bams=sample.data_path.split(" "),  # this is a list of sample paths
            merged_bam=sample.unmapped
        )
        pipe_manager.run(cmd, sample.unmapped, shell=True)
        sample.data_path = sample.unmapped

    # Fastqc
    pipe_manager.timestamp("Measuring read quality with Fastqc")
    cmd = tk.fastqc_rename(
        input_bam=sample.data_path,
        output_dir=sample.paths.sample_root,
        sample_name=sample.sample_name
    )
    pipe_manager.run(
        cmd, os.path.join(sample.paths.sample_root, sample.sample_name + "_fastqc.zip"), shell=True)
    report_dict(pipe_manager, parse_fastqc(os.path.join(sample.paths.sample_root, sample.sample_name + "_fastqc.zip"), prefix="fastqc_"))

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

    # HiC-Pro pipeline
    # make dir with linked fastq files for HiC-Pro
    sample.paths.hicpro_input = os.path.join(sample.paths.unmapped, sample.name)
    if not os.path.exists(sample.paths.hicpro_input):
        os.makedirs(sample.paths.hicpro_input)

    fq1 = os.path.join(sample.paths.hicpro_input, sample.name + "_R1.fastq")
    if not os.path.exists(fq1):
        pipe_manager.run(
            "ln -s {} {}".format(sample.fastq1, fq1),
            target=os.path.join(sample.paths.hicpro_input, os.path.basename(sample.fastq1)))
    fq2 = os.path.join(sample.paths.hicpro_input, sample.name + "_R2.fastq")
    if not os.path.exists(fq2):
        pipe_manager.run(
            "ln -s {} {}".format(sample.fastq2, fq2),
            target=os.path.join(sample.paths.hicpro_input, os.path.basename(sample.fastq2)))

    # edit config
    hicpro_config = open(pipe_manager.config.parameters.hicpro_template_config, 'r').read()
    with open(sample.hicpro_config, 'w') as handle:
        handle.write(hicpro_config.replace("\nJOB_NAME = \n", "\nJOB_NAME = {}\n".format(sample.name)))

    # run
    sample.paths.hicpro_output = os.path.join(sample.paths.sample_root, "hic-pro_output")
    if args.serial:
        # run the whole HiC-Pro pipeline as once
        pipe_manager.run(
            """{} -i {} -o {} -c {}""".format(
            pipe_manager.config.tools.hicpro, sample.paths.hicpro_input,
            sample.paths.hicpro_output, sample.hicpro_config),
            target=os.path.join(
                    sample.paths.hicpro_output,
                    "hic_results", "data", sample.name,
                    sample.name + "_allValidPairs"))
    else:
        # run each step in sequence
        pipe_manager.run(
            "{} -s mapping -i {} -o {} -c {}".format(
                pipe_manager.config.tools.hicpro,
                sample.paths.unmapped,
                sample.paths.hicpro_output,
                sample.hicpro_config),
            target=os.path.join(
                sample.paths.hicpro_output,
                "bowtie_results", "bwt2_global", sample.name,
                sample.name + "_R2_{}.bwt2glob.bam".format(sample.genome)))

        pipe_manager.run(
            "{} -s proc_hic -i {} -o {} -c {}".format(
                pipe_manager.config.tools.hicpro,
                os.path.join(sample.paths.hicpro_output, "bowtie_results", "bwt2"),
                sample.paths.hicpro_output,
                sample.hicpro_config),
            target=os.path.join(
                sample.paths.hicpro_output,
                "bowtie_results", "bwt2", sample.name,
                sample.name + "_{}.bwt2pairs.bam".format(sample.genome)))

        pipe_manager.run(
            "{} -s quality_checks -i {} -o {} -c {}".format(
                pipe_manager.config.tools.hicpro,
                sample.paths.unmapped,
                sample.paths.hicpro_output,
                sample.hicpro_config),
            target=os.path.join(
                sample.paths.hicpro_output,
                "hic_results", "pic", sample.name,
                "plotMappingPairing_" + sample.name + ".pdf"), nofail=True)

        pipe_manager.run(
            "{} -s merge_persample -i {} -o {} -c {}".format(
                pipe_manager.config.tools.hicpro,
                os.path.join(sample.paths.hicpro_output, "hic_results", "data"),
                sample.paths.hicpro_output,
                sample.hicpro_config),
            target=os.path.join(
                sample.paths.hicpro_output,
                "hic_results", "data", sample.name,
                sample.name + "_allValidPairs.mergestat"))

        pipe_manager.run(
            "{} -s build_contact_maps -i {} -o {} -c {}".format(
                pipe_manager.config.tools.hicpro,
                os.path.join(sample.paths.hicpro_output, "hic_results", "data"),
                sample.paths.hicpro_output,
                sample.hicpro_config),
            target=os.path.join(
                sample.paths.hicpro_output,
                "hic_results", "matrix", sample.name,
                "raw", "1000", sample.name + "_1000.matrix"))

        pipe_manager.run(
            "{} -s ice_norm -i {} -o {} -c {}".format(
                pipe_manager.config.tools.hicpro,
                os.path.join(sample.paths.hicpro_output, "hic_results", "matrix", sample.name, "raw"),
                sample.paths.hicpro_output,
                sample.hicpro_config),
            target=os.path.join(
                sample.paths.hicpro_output,
                "hic_results", "matrix",
                "1000", "iced", "1000", "1000_1000_iced.matrix"))

    # Report stats
    stats = get_hicpro_stats(sample)
    report_dict(pipe_manager, stats.to_dict())

    ## Convertions

    ### HiC-Pro output to Juicebox ".hic"
    pipe_manager.run(
        "{} -i {} -g {} -j {} -r {} -o {}"
            .format(pipe_manager.config.tools.hicpro2juicebox,
                os.path.join(
                    sample.paths.hicpro_output,
                    "hic_results", "data", sample.name,
                    sample.name + "_allValidPairs"),
                pipe_manager.config.resources.chromosome_sizes[sample.genome],
                pipe_manager.config.tools.juicertools,
                pipe_manager.config.parameters.hicpro_restriction_fragments,
                sample.paths.hicpro_output),
        target=os.path.join(sample.paths.hicpro_output, sample.name + "_allValidPairs.hic"))

    ### make pairix indexed BEDPE
    pipe_manager.run(
        "awk -v OFS='\\t' '{{print $2,$3,$3+75,$5,$6,$6+75,\".\",\".\",$4,$7}}' {} | sort -k1,1V -k4,4V -k2,2n -k5,5n | bgzip -@ {} > {}".format(
            os.path.join(sample.paths.hicpro_output, "hic_results", "data", sample.name, sample.name + "_allValidPairs"),
            args.cores,
            os.path.join(sample.paths.hicpro_output, sample.name + "_allValidPairs.bed.gz")),
        target=os.path.join(sample.paths.hicpro_output, sample.name + "_allValidPairs.bed.gz"))
    pipe_manager.run(
        "pairix  -s 1 -d 4 -b 2 -e 3 -u 5 -v 6 {}".format(
            os.path.join(sample.paths.hicpro_output, sample.name + "_allValidPairs.bed.gz")),
        target=os.path.join(sample.paths.hicpro_output, sample.name + "_allValidPairs.bed.gz.px2"))

    ### make cool
    pipe_manager.run(
        "hic2cool {} {}".format(
            os.path.join(sample.paths.hicpro_output, sample.name + "_allValidPairs.hic"),
            os.path.join(sample.paths.hicpro_output, sample.name + "_allValidPairs.cool")),
        target=os.path.join(sample.paths.hicpro_output, sample.name + "_allValidPairs.multi.cool"))

    # add balanced normalizations to cooler file
    for resolution in [1, 5, 10, 25, 100, 250, 500, 1000]:
        pipe_manager.run(
            "cooler balance -p {} --blacklist {} {}::/resolutions/{}".format(
                args.cores,
                pipe_manager.config.resources.blacklisted_regions[sample.genome],
                os.path.join(sample.paths.hicpro_output, sample.name + "_allValidPairs.multi.cool"),
                resolution * 1000),
            lock_name="cooler.balance.{}kb".format(resolution), nofail=True)

    # Call peaks with MACS2
    ## TODO: optimize parameters further
    pipe_manager.run(
        "macs2 callpeak -t {} -f BEDPE --keep-dup auto --nomodel --extsize 147 -g hs -n {} --outdir {}".format(
            os.path.join(sample.paths.hicpro_output, sample.name + "_allValidPairs.bed.gz"),
            sample.name,
            os.path.join(sample.paths.hicpro_output, "hic_results", "peaks")),
        target=os.path.join(sample.paths.hicpro_output, "hic_results", "peaks", sample.name + "_peaks.narrowPeak"), nofail=True)

    # Call loops
    ### with cLoops
    if not os.path.exists(os.path.join(sample.paths.hicpro_output, "hic_results", "cLoops")):
        os.makedirs(os.path.join(sample.paths.hicpro_output, "hic_results", "cLoops"))
    pipe_manager.run(
        "cLoops -f {} -o {} ".format(
            os.path.join(sample.paths.hicpro_output, sample.name + "_allValidPairs.bed.gz"),
            os.path.join(sample.paths.hicpro_output, "hic_results", "cLoops", sample.name)
        ) +
        "-m 4 " +
        "-eps 5000,7500,10000 " +
        "-minPts 10,20,30,40,50 " +
        "-p {} ".format(args.cores) +
        "-w -j -s -hic",
        target=os.path.join(sample.paths.hicpro_output, "hic_results", "cLoops", sample.name + ".loop"), nofail=True)

    ### with hichipper
    #### make hichipper config file
    yaml = textwrap.dedent("""
    peaks:
     - {}
    resfrags:
     - {}
    hicpro_output:
     - {}""".format(
        os.path.join(sample.paths.hicpro_output, "hic_results", "peaks", sample.name + "_peaks.narrowPeak"),
        pipe_manager.config.resources.hicpro_restriction_fragments,
        os.path.join(sample.paths.hicpro_output)))

    if os.path.exists(os.path.join(sample.paths.sample_root, "hichipper")):
        import shutil
        shutil.rmtree(os.path.join(sample.paths.sample_root, "hichipper"))

    hichipper_config = os.path.join(sample.paths.sample_root, "hichipper_config.yaml")
    with open(hichipper_config, 'w') as handle:
        handle.write(yaml)
    #### run
    pipe_manager.run(  # TODO: I think this command has to be run from sample.paths.sample_root, needs testing
        "hichipper --out {} {}".format(
            os.path.join(sample.paths.sample_root, "hichipper"),
            hichipper_config),
        target=os.path.join(sample.paths.sample_root, "hichipper", sample.name + ".filt.intra.loop_counts.bedpe"), nofail=True)
    # or target to os.path.join(sample.paths.hicpro_output, "hic_results", "hichipper", "qcReport_make.html")

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
        line = [i for i in range(len(content)) if "Total Sequences" in content[i]][0]
        total = int(re.sub("\D", "", re.sub("\(.*", "", content[line])))
        line = [i for i in range(len(content)) if "Sequences flagged as poor quality" in content[i]][0]
        poor_quality = int(re.sub("\D", "", re.sub("\(.*", "", content[line])))
        line = [i for i in range(len(content)) if "Sequence length  " in content[i]][0]
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


def get_hicpro_stats(sample):
    stats = pd.Series()
    stats_files = [
        os.path.join(sample.paths.hicpro_output, 'bowtie_results', 'bwt2', sample.name, sample.name + '.mpairstat'),
        os.path.join(sample.paths.hicpro_output, 'hic_results', 'data', sample.name, sample.name + '.mRSstat'),
        os.path.join(sample.paths.hicpro_output, 'hic_results', 'data', sample.name, sample.name + '_allValidPairs.mergestat')]
    for stats_file in stats_files:
        try:
            s = pd.read_table(stats_file, header=None, comment="#", index_col=0)
            stats = stats.append(s[1])
        except IOError:
            continue

    trim_file = os.path.join(sample.paths.hicpro_output, 'logs', sample.name, 'readsTrimming.log')
    try:
        stats['trimmed_reads'] = int(os.popen('grep reads {}'.format(trim_file)).read().strip().split(": ")[1])
    except:
        pass
    all_valid_pairs = os.path.join(sample.paths.hicpro_output, "hic_results", "data", sample.name, sample.name + "_allValidPairs")
    try:
        stats['all_valid_pairs'] = int(os.popen('wc -l {}'.format(all_valid_pairs)).read().strip().split(" ")[0])
    except:
        pass

    return stats


if __name__ == '__main__':
    try:
        sys.exit(main())
    except KeyboardInterrupt:
        print("Program canceled by user!")
        sys.exit(1)
