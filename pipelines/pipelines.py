#!/usr/bin/env python

"""
pipelines
=========

Project management and Sample loop.
"""

from argparse import ArgumentParser
from .models import Project
from . import toolkit as tk
import cPickle as pickle
import os
import pandas as pd
import sys
import textwrap
import time

__author__ = "Andre Rendeiro"
__copyright__ = "Copyright 2015, Andre Rendeiro"
__credits__ = []
__license__ = "GPL2"
__version__ = "0.1"
__maintainer__ = "Andre Rendeiro"
__email__ = "arendeiro@cemm.oeaw.ac.at"
__status__ = "Development"


def main():
    # Parse command-line arguments
    parser = ArgumentParser(
        prog="pipelines",
        description="pipelines. Project management and sample loop."
    )
    parser = add_args(parser)

    # Parse
    args = parser.parse_args()

    # Start project
    prj = Project(args.project_name)
    prj.addSampleSheet(args.csv)

    # Start main function
    if args.stats:
        read_stats(prj)
    elif args.compare:
        compare()
    else:
        sample_loop(args, prj)

    # Exit
    print("Finished and exiting.")
    sys.exit(0)


def add_args(parser):
    """
    Options for project and pipelines.
    """
    # Project
    parser.add_argument(dest="project_name", help="Project name.", type=str)
    parser.add_argument(dest="csv", help="CSV file with sample annotation.", type=str)  # improvement: check project dirs for csv

    # Behaviour
    parser.add_argument("--stats", dest="stats", action="store_true",
                        help="Do not run pipelines, but gather stats on produced files.")
    parser.add_argument("--compare", dest="compare", action="store_true",
                        help="Do not loop through samples, but perform comparisons betweem them.")
    parser.add_argument("-r", "--rm-tmp", dest="rm_tmp", action="store_true",
                        help="Remove intermediary files. If not it will preserve all intermediary files. Default=False.")
    parser.add_argument("--dry-run", dest="dry_run", action="store_true",
                        help="Dry run. Assemble commands, but do not submit jobs to slurm. Default=False.")
    parser.add_argument("--no-checks", dest="checks", action="store_false",
                        help="Don't check file existence and integrity. Default=False.")

    # Pypiper
    parser.add_argument("--overwrite", dest="recover", action="store_true",
                        help="Overwrite existing files. Default=False.")
    parser.add_argument("--fresh-start", dest="fresh", action="store_true",
                        help="Start from beginning of pipeline. Default=False.")
    parser.add_argument("--manual-clean", dest="manual_clean", action="store_true",
                        help="Manually clean temporary files. Default=False.")

    # Slurm-related
    parser.add_argument("-c", "--cpus", default=4, dest="cpus",
                        help="Number of CPUs to use. Default is specified in the pipeline config file.", type=int)
    parser.add_argument("-m", "--mem-per-cpu", default=4000, dest="mem",
                        help="Memory per CPU to use. Default is specified in the pipeline config file.", type=int)
    parser.add_argument("-q", "--queue", default="shortq", dest="queue",
                        choices=["develop", "shortq", "mediumq", "longq"],
                        help="Queue to submit jobs to. Default is specified in the pipeline config file.", type=str)
    parser.add_argument("-t", "--time", default="10:00:00", dest="time",
                        help="Maximum time for jobs to run. Default is specified in the pipeline config file.", type=str)
    parser.add_argument("-u", "--user-mail", default="mail@example.com", dest="user_mail",
                        help="User email.", type=str)

    # Preprocessing: trimming, mapping, etc...
    parser.add_argument("--trimmer", default="skewer", choices=["trimmomatic", "skewer"],
                        dest="trimmer", help="Trimmer to use. Default=skewer.", type=str)
    parser.add_argument("-i", "--max-insert-size", default=2000,
                        dest="maxinsert",
                        help="Maximum allowed insert size allowed for paired end mates. Default=2000.",
                        type=int)
    parser.add_argument("-Q", "--quality", default=30,
                        dest="quality",
                        help="Minimum read quality to keep. Default=30.",
                        type=int)

    # Further downstream
    parser.add_argument("--window-size", default=1000, dest="windowsize",
                        help="Window size used for genome-wide correlations. Default=1000.",
                        type=int)
    parser.add_argument("--peak-caller", default="macs2", choices=["macs2", "spp"],
                        dest="peak_caller", help="Peak caller to use. Default=macs2.", type=str)
    parser.add_argument("--peak-window-width", default=2000,
                        dest="peak_window_width",
                        help="Width of window around peak motifs. Default=2000.",
                        type=int)

    return parser


def sample_loop(args, prj):
    """
    Loop through all samples and submit jobs to the pipeline under Slurm.

    :param args: Parsed ArgumentParser object.
    :type args: argparse.ArgumentParser
    :param prj: `Project` object.
    :type prj: pipelines.Project
    """

    print("Starting sample preprocessing into jobs.")

    # start pipeline
    run_name = "_".join([prj.name, time.strftime("%Y%m%d-%H%M%S")])

    # add track headers to track hubs
    for genome in pd.Series([s.genome for s in prj.samples]).unique():
        if not os.path.exists(os.path.join(prj.dirs.html, "trackHub_{0}.txt".format(genome))):
            with open(os.path.join(prj.dirs.html, "trackHub_{0}.txt".format(genome)), "w") as handle:
                handle.write("browser position {0}\n".format(prj.config["defaultposition"]))

    # Loop through samples, submit to corresponding job (preprocess, analyse)
    for sample in prj.samples:
        # get job_name
        job_name = "_".join([run_name, sample.name])

        # if unmappedBam is a list, add final "unmapped" attr to sample object
        if type(sample.unmappedBam) is list:
            sample.unmapped = os.path.join(sample.dirs.unmapped, sample.name + ".bam")

        # assemble command
        # slurm header
        job_code = tk.slurmHeader(
            jobName=job_name,
            output=os.path.join(prj.dirs.logs, job_name + ".slurm.log"),
            queue=args.queue,
            time=args.time,
            cpusPerTask=args.cpus,
            memPerCpu=args.mem,
            userMail=args.user_mail
        )

        sample_pickle = os.path.join(prj.dirs.pickles, job_name + ".pickle")
        # self reference the pickle file in its sample
        sample.pickle = sample_pickle

        # If sample has control attribute, get that sample and pair them
        if hasattr(sample, "controlname"):
            if type(sample.controlname) == str:
                # Assign the sample with that name to ctrl
                ctrl = [s for s in prj.samples if s.name == sample.controlname]
                # if there is only one record, use that as control
                if len(ctrl) == 1:
                    sample.ctrl = ctrl[0]
                else:
                    # if not, process sample anyway, but without a matched control
                    print("Provided control sample name does not exist or is ambiguous: %s" % sample.controlname)

        # save pickle with all objects (this time, 2nd element is a tuple!)
        pickle.dump((prj, sample, args), open(sample_pickle, "wb"))

        # Actual call to pipeline
        technique = sample.technique.upper()
        if technique in prj.config["techniques"]["chipseq"]:
            job_code += "chipseq-pipeline {0}\n".format(sample_pickle)
        elif technique in prj.config["techniques"]["cm"]:
            job_code += "chipseq-pipeline {0}\n".format(sample_pickle)
        elif technique in prj.config["techniques"]["atacseq"]:
            job_code += "atacseq-pipeline {0}\n".format(sample_pickle)
        elif technique in prj.config["techniques"]["dnase"]:
            job_code += "atacseq-pipeline {0}\n".format(sample_pickle)
        elif technique in prj.config["techniques"]["quantseq"]:
            job_code += "quantseq-pipeline {0}\n".format(sample_pickle)
        elif technique in prj.config["techniques"]["chemseq"]:
            job_code += "chipseq-pipeline {0}\n".format(sample_pickle)
        else:
            raise TypeError("Sample is not in known sample class.")

        # Slurm footer
        job_code += tk.slurmFooter()

        # Save code as executable
        job_file = os.path.join(prj.dirs.executables, job_name + ".sh")
        with open(job_file, 'w') as handle:
            handle.write(textwrap.dedent(job_code))

        # Submit to slurm
        if not args.dry_run:
            status = tk.slurmSubmitJob(job_file)

            if status != 0:
                print("Could not submit job '%s' to slurm." % job_file)
                sys.exit(1)
            print("Submitted job to slurm: '%s'" % job_name)

        # Create link to trackHub in project folder
        tk.linkToTrackHub(
            trackHubURL=os.path.join(prj.dirs.html, "trackHub_{0}.txt".format(sample.genome)),
            fileName=os.path.join(prj.dirs.root, "ucsc_tracks_{0}.html".format(sample.genome)),
            genome=sample.genome
        )

    # write original annotation sheet to project folder
    # add field for manual sample-control pairing
    prj.sheet.df.controlname = None
    prj.sheet.to_csv(os.path.join(prj.dirs.root, prj.name + ".annotation_sheet.csv"))

    print("Finished preprocessing")


def read_stats(prj):
    """
    Given an annotation sheet with replicates, gets number of reads mapped, duplicates, etc...

    :param prj: `Project` object.
    :type prj: pipelines.Project
    """

    print("Starting sample read stats.")

    bowtie_cols = ["readCount", "unpaired", "unaligned", "unique", "multiple", "alignmentRate"]

    samples = pd.DataFrame(index=["name"] + bowtie_cols + ["single-ends", "paired-ends", "duplicates", "NSC", "RSC", "qualityTag", "peakNumber", "FRiP"])

    for sample in prj.samples:
        sample = sample.asSeries()
        # Get alignment stats
        try:
            sample = sample.append(parse_bowtie_stats(sample.alnRates))
        except:
            print("Record with alignment rates is empty or not found for sample %s" % sample.name)

        # Get duplicate stats
        try:
            sample = sample.append(parse_duplicate_stats(sample.dupsMetrics))
        except:
            print("Record with duplicates is empty or not found for sample %s" % sample.name)

        # Get NSC and RSC
        try:
            sample = sample.append(parse_qc(sample.qc))
        except:
            print("Record with quality control is empty or not found for sample %s" % sample.name)

        # Count peak number (if peaks exist)
        if hasattr(sample, "peaks"):
            # and if sample has peaks
            if str(sample.peaks) != "nan":
                # if peak file exist
                if os.path.exists(sample.peaks):
                    sample = get_peak_number(sample)

        # Get FRiP from file (if exists) and add to sheet
        if hasattr(sample, "peaks"):
            # if sample has peaks
            if str(sample.peaks) != "nan":
                try:
                    sample = get_frip(sample)
                except:
                    print("Record with FRiP value is empty or not found for sample %s" % sample.name)
        samples[sample["name"]] = sample

    # write annotation sheet with statistics
    samples.T.to_csv(prj.sampleStats, index=False)

    print("Finished getting read statistics.")


def compare():
    raise NotImplementedError


def parse_bowtie_stats(stats_file):
    """
    Parses Bowtie2 stats file, returns series with values.

    :param stats_file: Bowtie2 output file with alignment statistics.
    :type stats_file: str
    """
    import re

    stats = pd.Series(index=["readCount", "unpaired", "unaligned", "unique", "multiple", "alignmentRate"])

    try:
        with open(stats_file) as handle:
            content = handle.readlines()  # list of strings per line
    except:
        return stats

    # total reads
    try:
        line = [i for i in range(len(content)) if " reads; of these:" in content[i]][0]
        stats["readCount"] = re.sub("\D.*", "", content[line])
        if 7 > len(content) > 2:
            line = [i for i in range(len(content)) if "were unpaired; of these:" in content[i]][0]
            stats["unpaired"] = re.sub("\D", "", re.sub("\(.*", "", content[line]))
        else:
            line = [i for i in range(len(content)) if "were paired; of these:" in content[i]][0]
            stats["unpaired"] = stats["readCount"] - int(re.sub("\D", "", re.sub("\(.*", "", content[line])))
        line = [i for i in range(len(content)) if "aligned 0 times" in content[i]][0]
        stats["unaligned"] = re.sub("\D", "", re.sub("\(.*", "", content[line]))
        line = [i for i in range(len(content)) if "aligned exactly 1 time" in content[i]][0]
        stats["unique"] = re.sub("\D", "", re.sub("\(.*", "", content[line]))
        line = [i for i in range(len(content)) if "aligned >1 times" in content[i]][0]
        stats["multiple"] = re.sub("\D", "", re.sub("\(.*", "", content[line]))
        line = [i for i in range(len(content)) if "overall alignment rate" in content[i]][0]
        stats["alignmentRate"] = re.sub("\%.*", "", content[line]).strip()
    except IndexError:
        pass
    return stats


def parse_duplicate_stats(stats_file):
    """
    Parses Bowtie2 stats file, returns series with values.

    :param stats_file: Bowtie2 output file with alignment statistics.
    :type stats_file: str
    """
    import re

    series = pd.Series()

    try:
        with open(stats_file) as handle:
            content = handle.readlines()  # list of strings per line
    except:
        return series

    try:
        line = [i for i in range(len(content)) if "single ends (among them " in content[i]][0]
        series["single-ends"] = re.sub("\D", "", re.sub("\(.*", "", content[line]))
        line = [i for i in range(len(content)) if " end pairs...   done in " in content[i]][0]
        series["paired-ends"] = re.sub("\D", "", re.sub("\.\.\..*", "", content[line]))
        line = [i for i in range(len(content)) if " duplicates, sorting the list...   done in " in content[i]][0]
        series["duplicates"] = re.sub("\D", "", re.sub("\.\.\..*", "", content[line]))
    except IndexError:
        pass
    return series


def parse_qc(qc_file):
    """
    Parses QC table produced by phantompeakqualtools (spp) and returns sample quality metrics.

    :param qc_file: phantompeakqualtools output file sample quality measurements.
    :type qc_file: str
    """
    series = pd.Series()
    try:
        with open(qc_file) as handle:
            line = handle.readlines()[0].strip().split("\t")  # list of strings per line
        series["NSC"] = line[-3]
        series["RSC"] = line[-2]
        series["qualityTag"] = line[-1]
    except:
        pass
    return series


def get_peak_number(sample):
    """
    Counts number of peaks from a sample's peak file.

    :param sample: A Sample object with the "peaks" attribute.
    :type name: pipelines.Sample
    """
    import subprocess
    import re

    proc = subprocess.Popen(["wc", "-l", sample.peaks], stdout=subprocess.PIPE)
    out, err = proc.communicate()
    sample["peakNumber"] = re.sub("\D.*", "", out)

    return sample


def get_frip(sample):
    """
    Calculates the fraction of reads in peaks for a given sample.

    :param sample: A Sample object with the "peaks" attribute.
    :type name: pipelines.Sample
    """
    import re

    with open(sample.frip, "r") as handle:
        content = handle.readlines()

    reads_in_peaks = int(re.sub("\D", "", content[0]))
    mapped_reads = sample["readCount"] - sample["unaligned"]

    return pd.Series(reads_in_peaks / mapped_reads, index="FRiP")


if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        print("Program canceled by user!")
        sys.exit(1)
