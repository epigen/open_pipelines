#!/usr/bin/env python

import csv
import os
import sys

EXTENSION = "chromSizes"
DEFAULT_CHROMSIZES_LOCATION = "/data/groups/lab_bock/shared/resources/genomes"




def getChrSizes(chrmFile):
    """
    Reads tab-delimiter file with two rows describing the chromosomes and its lengths.
    Returns dictionary of chr:sizes.
    """
    with open(chrmFile, 'r') as f:
        chrmSizes = {}
        for line in enumerate(f):
            row = line[1].strip().split('\t')
            chrmSizes[str(row[0])] = int(row[1])
    return chrmSizes



def getChrFile(assembly):
    filename = "{}.{}".format(assembly, EXTENSION)
    default_filepath = os.path.join(
            DEFAULT_CHROMSIZES_LOCATION, assembly, filename)
    genomes_folder_path = os.getenv("GENOMES")
    try:
        env_based_filepath = os.path.join(
                genomes_folder_path, assembly, filename)
    except TypeError:
        return default_filepath
    return env_based_filepath if os.path.exists(env_based_filepath) \
            else default_filepath



genome = sys.argv[1]
chrSizesFilepath = getChrFile(genome)
chrms = getChrSizes(chrSizesFilepath)  # get size of chromosomes

wr = csv.writer(sys.stdout, delimiter='\t', lineterminator='\n')

for row in csv.reader(iter(sys.stdin.readline, ''), delimiter='\t'):
    chrm = row[0]
    start = int(row[1])
    end = int(row[2])

    if chrm in chrms.keys():  # skip weird chromosomes
        if start >= 1 and end <= chrms[chrm] and start < end:
            wr.writerow(row)
