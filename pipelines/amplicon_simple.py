#!/usr/bin/env python

"""
Simple CRISPR amplicon pipeline
"""

import os
import sys
from argparse import ArgumentParser
import pandas as pd
import pypiper


__author__ = "Andre Rendeiro"
__copyright__ = "Copyright 2016, Andre Rendeiro"
__credits__ = []
__license__ = "GPL2"
__version__ = "0.4"
__maintainer__ = "Andre Rendeiro"
__email__ = "arendeiro@cemm.oeaw.ac.at"
__status__ = "Development"


def arg_parser(parser):
    """
    Global options for pipeline.
    """
    parser.add_argument(
        "-i",
        dest="input",
        help="Input Fastq file.",
        type=str
    )
    parser.add_argument(
        "-n",
        dest="sample_name",
        help="Sample name.",
        type=str
    )
    parser.add_argument(
        "-a", "--amplicon",
        dest="amplicon",
        help="Full amplicon sequence.",
        type=str
    )
    parser.add_argument(
        "-g", "--guide-rna",
        dest="guide_rna",
        help="Guide RNA sequence used to target the genome.",
        type=str
    )
    return parser


def count_sizes(fastq_file, amplicon, guide_rna, window=20, anchor_length=10):
    """
    Count indel sizes based on sequence identity (pattern-based method).
    """
    import itertools
    import re
    from collections import Counter

    editing_position = amplicon.index(guide_rna)

    guide_plus_window = (editing_position - window, editing_position + len(guide_rna) + window)

    left_anchor_start = guide_plus_window[0] - anchor_length
    right_anchor_end = guide_plus_window[1] + anchor_length

    left_anchor_end = left_anchor_start + anchor_length
    right_anchor_start = right_anchor_end - anchor_length

    a = amplicon[left_anchor_start: left_anchor_end]
    b = amplicon[right_anchor_start: right_anchor_end]

    pattern = a + "(.*)" + b

    window_size = (right_anchor_start - left_anchor_end)

    # Open file handle
    if fastq_file.endswith(".gz"):
        import gzip
        handle = gzip.open(fastq_file, 'r')
    else:
        handle = open(fastq_file, 'r')

    # Iterate over sequence lines, append matches
    seqs = itertools.islice(handle, 1, None, 4)
    lines = list()
    for seq in seqs:
        m = re.findall(pattern, seq.decode("ascii").strip())
        if len(m) != 0:
            lines += m
    handle.close()

    # return difference between window size and pattern match length
    return pd.Series(Counter([len(x) - window_size for x in lines]))


def main():
    # Parse command-line arguments
    parser = ArgumentParser(
        prog="amplicon-pipeline",
        description="Amplicon pipeline."
    )
    parser = arg_parser(parser)
    parser = pypiper.add_pypiper_args(parser, all_args=True)
    args = parser.parse_args()

    print("Processing sample {}.".format(args.sample_name))

    output_folder = os.path.abspath(os.path.join(args.output_parent, args.sample_name))

    # Create output directories if not existing
    for path in [args.output_parent, output_folder]:
        if not os.path.exists(path):
            try:
                os.mkdir(path)
            except OSError("Cannot create directory '{}'".format(path)):
                raise

    # Count length of pattern matches
    sizes = count_sizes(
        fastq_file=args.input,
        amplicon=args.amplicon,
        guide_rna=args.guide_rna)

    # Calculate efficiency
    efficiency = (sizes[sizes.index != 0].sum() / float(sizes.sum())) * 100

    # Save
    with open(os.path.join(output_folder, args.sample_name + ".efficiency.csv"), 'w') as handle:
        handle.write("{},{}\n".format(args.sample_name, efficiency))

    print("Sample {} has an editing efficiency of {}.".format(args.sample_name, efficiency))
    print("Finished processing sample {}.".format(args.sample_name))


if __name__ == '__main__':
    try:
        sys.exit(main())
    except KeyboardInterrupt:
        print("Program canceled by user!")
        sys.exit(1)
