#!/usr/bin/env python

"""
pipelines
=========

Pipelines for NGS

Workflow explained:
    - Project is created
    - Add Sample sheet to project (spawns next)
        - Samples are created and added to project (automatically)

In the process, stuff is checked:
    - project structure (created if not existing)
    - existance of csv sample sheet
    - existance of bam files from samples
    - read type of samples

:Example:

from pipelines import Project
prj = Project("ngs")
prj.addSampleSheet("sample_annotation.csv")
# that's it!

# explore!
prj.samples
prj.samples[0].mapped
prj.samples[3].nodupsshifted

prj.dirs.results
prj.sheet.to_csv(_os.path.join(prj.dirs.root, "sample_annotation.csv"))

# project options are read from the config file
# but can be changed on the fly:
prj = Project("ngs")
prj.config["mergetechnical"] = False
prj.addSampleSheet("sample_annotation.csv")

"""

import os as _os
import pandas as _pd
import yaml as _yaml


class Paths(object):
    """
    A class to hold paths as attributes.
    """
    pass


class Project(object):
    """
    A class to model a Project.

    :param name: Project name.
    :type name: str

    Kwargs (will overule specified in config):
    :param parent: Path to where the project structure will be created.
    :type parent: str
    :param parenthtml: Path to where the project structure will be created.
    :type parenthtml: str
    :param dry: If dry mode is activated, no directories will be created.
    :type dry: bool

    :Example:

    from pipelines import Project
    prj = Project("ngs")
    prj = Project("ngs2", parent="/projects", parenthtml="/public_html")
    """
    def __init__(self, name, dry=False, **kwargs):
        super(Project, self).__init__()
        # check it's a string
        self.name = name

        # Path structure
        self.dirs = Paths()

        # Read configuration file
        with open(_os.path.join(_os.path.expanduser("~"), ".pipelines_config.yaml"), 'r') as handle:
            self.config = _yaml.load(handle)

        # If kwargs were passed, overule paths specified in the config with new ones.
        # parse kwargs
        for key in ('parent', 'parenthtml'):
            if key in kwargs:
                # check they're strings
                setattr(self.dirs, key, kwargs[key])
            else:
                setattr(self.dirs, key, self.config["paths"][key])

        # flow
        self.setProjectDirs()
        if not dry:
            self.makeProjectDirs()
            self.setProjectPermissions()

        # samples
        self.samples = list()

    def __repr__(self):
        return "Project '%s'" % self.name

    def setProjectDirs(self):
        """
        Atributes directories for the project.
        """
        # check self.project_root exists and user has write access
        self.dirs.parent = _os.path.abspath(self.dirs.parent)
        if not _os.access(self.dirs.parent, _os.W_OK):
            raise IOError("%s does not exist, or user has no write access.\n\
            Use option '-r' '--project-root' to set a non-default project root path." % self.dirs.parent)

        # check self.html_root exists and user has write access
        self.dirs.parenthtml = _os.path.abspath(self.dirs.parenthtml)
        if not _os.access(self.dirs.parenthtml, _os.W_OK):
            raise IOError("%s does not exist, or user has no write access.\n\
            Use option '--html-root' to set a non-default html root path." % self.dirs.parenthtml)

        # Project directory structure
        self.dirs.root = _os.path.join(self.dirs.parent, self.name)
        # Runs
        self.dirs.runs = _os.path.join(self.dirs.root, "runs")
        self.dirs.pickles = _os.path.join(self.dirs.runs, "pickles")
        self.dirs.executables = _os.path.join(self.dirs.runs, "executables")
        self.dirs.logs = _os.path.join(self.dirs.runs, "logs")
        # Data
        self.dirs.data = _os.path.join(self.dirs.root, "data")
        # Results
        self.dirs.results = _os.path.join(self.dirs.root, "results")
        self.dirs.plots = _os.path.join(self.dirs.results, "plots")
        # Html structure
        self.dirs.html = _os.path.join(self.dirs.parenthtml, self.name)

        # Sample stats csv
        self.sampleStats = _os.path.join(self.dirs.root, self.name + ".sample_stats.csv")
        # Diffbind file
        self.diffBindCSV = _os.path.join(self.dirs.root, self.name + ".diffBind.csv")

    def makeProjectDirs(self):
        """
        Creates project directory structure if it doesn't exist.
        """
        for name, path in self.dirs.__dict__.items():
            if not _os.path.exists(path):
                _os.makedirs(path)

    def setProjectPermissions(self):
        """
        Makes the project's public_html folder executable.
        """
        for d in [self.dirs.parenthtml, self.dirs.html]:
            try:
                _os.chmod(d, 0755)
            except OSError:
                # logger.error("cannot change folder's mode: %s" % d)
                continue

    def addSampleSheet(self, csv):
        """
        Build a `SampleSheet` object from a csv file and
        add it and its samples to the project.

        :param csv: Path to csv file.
        :type csv: str
        """
        # Make SampleSheet object
        self.sheet = SampleSheet(csv)

        # pair project and sheet
        self.sheet.project = self

        # Generate sample objects from annotation sheet
        self.sheet.makeSamples()

        # Generate sample objects if merging options are on:
        if self.config["options"]["mergetechnical"] and hasattr(self.sheet.df, "technicalReplicate"):
            self.sheet.getBiologicalReplicates()
        if self.config["options"]["mergebiological"] and hasattr(self.sheet.df, "biologicalReplicate"):
            self.sheet.getMergedBiologicalReplicates()

        # Add samples to Project
        for sample in self.sheet.samples:
            # Check sample is from a supported genome
            if sample.genome not in self.config["genomes"]:
                raise TypeError("Sample's genome is not supported.")
            self.addSample(sample)
            sample.setFilePaths()
            sample.makeSampleDirs()

    def addSample(self, sample):
        """
        Adds a sample to the project's `samples`.
        """
        # Check sample is Sample object
        if not isinstance(sample, Sample):
            raise TypeError("Provided object is not a Sample object.")

        # Tie sample and project bilateraly
        sample.project = self
        # Append
        self.samples.append(sample)


class SampleSheet(object):
    """
    Class to model a sample annotation sheet.

    :param csv: Path to csv file.
    :type csv: str

    Kwargs (will overule specified in config):
    :param mergetechnical: Should technical replicates be merged to create biological replicate samples?
    :type mergetechnical: bool
    :param mergebiological: Should biological replicates be merged?
    :type mergebiological: bool

    :Example:

    from pipelines import Project, SampleSheet
    prj = Project("ngs")
    sheet = SampleSheet("/projects/example/sheet.csv")
    """
    def __init__(self, csv, **kwargs):

        super(SampleSheet, self).__init__()
        # TODO: checks on given args
        self.csv = csv

        # Read configuration file
        with open(_os.path.join(_os.path.expanduser("~"), ".pipelines_config.yaml"), 'r') as handle:
            self.config = _yaml.load(handle)

        # Sample merging options:
        # If kwargs were passed, overule options specified in the config with new ones.
        # parse kwargs
        self.opts = dict()
        for key in ('mergetechnical', 'mergebiological'):
            if key in kwargs:
                # check they're strings
                self.opts[key] = kwargs[key]
            else:
                self.opts[key] = self.config["options"][key]

        self.samples = list()
        self.checkSheet()

    def __repr__(self):
        if hasattr(self, "project"):
            return "SampleSheet for project '%s' with %i samples." % (self.project, len(self.df))
        else:
            return "SampleSheet with %i samples." % len(self.df)

    def checkSheet(self):
        """
        Check if csv file exists and has all required columns.
        """
        try:
            self.df = _pd.read_csv(self.csv)
        except IOError("Given csv file couldn't be read.") as e:
            raise e

        req = ["technique", "genome", "unmappedBam"]
        missing = [col for col in req if col not in self.df.columns]

        if len(missing) != 0:
            raise TypeError("Annotation sheet is missing columns: %s" % " ".join(missing))

    def makeSample(self, series):
        """
        Make a children of class Sample dependent on its "technique" attribute.

        :param series: Pandas `Series` object.
        :type series: pandas.Series
        :return: An object or class `Sample` or a child of that class.
        :rtype: pipelines.Sample
        """
        technique = series["technique"].upper()
        if technique in self.config["techniques"]["chipseq"]:
            return ChIPseqSample(series)
        elif technique in self.config["techniques"]["cm"]:
            return CMSample(series)
        elif technique in self.config["techniques"]["dnase"]:
            return DNaseSample(series)
        elif technique in self.config["techniques"]["atacseq"]:
            return ATACseqSample(series)
        elif technique in self.config["techniques"]["quantseq"]:
            return QuantseqSample(series)
        else:
            raise TypeError("Sample is not in known sample class.")
            # I might want to change this behaviour to return a generic sample
            # self.samples.append(Sample(self.df.ix[sample]))

    def makeSamples(self):
        """
        Creates samples from annotation sheet dependent on technique and adds them to the project.
        """
        for i in range(len(self.df)):
            self.samples.append(self.makeSample(self.df.ix[i]))

    def getBiologicalReplicates(self):
        """
        Produces biological replicate samples from merged technical replicates.
        """
        # copy columns list
        attributes = self.df.columns.tolist()[:]
        # ignore some fields in the annotation sheet
        attributes.pop(attributes.index("unmappedBam"))
        attributes.pop(attributes.index("technicalReplicate"))
        if "controlSampleName" in attributes:
            attributes.pop(attributes.index("controlSampleName"))

        # get technical replicates -> biological replicates
        for key, values in self.df.groupby(attributes).groups.items():
            rep = self.df.ix[values][attributes].reset_index(drop=True).ix[0]
            if len(values) > 1:
                rep["technicalReplicate"] = 0
                rep["unmappedBam"] = self.df.ix[values]["unmappedBam"].tolist()
                # append biological replicate to samples
                self.samples.append(self.makeSample(rep))

    def getMergedBiologicalReplicates(self):
        """
        Produces samples from merged biological replicates.
        """
        # copy columns list
        attributes = self.df.columns.tolist()[:]
        # ignore some fields in the annotation sheet
        attributes.pop(attributes.index("unmappedBam"))
        attributes.pop(attributes.index("technicalReplicate"))
        attributes.pop(attributes.index("biologicalReplicate"))
        if "controlSampleName" in attributes:
            attributes.pop(attributes.index("controlSampleName"))

        # get biological replicates -> merged biological replicates
        for key, values in self.df.groupby(attributes).groups.items():
            rep = self.df.ix[values][attributes].reset_index(drop=True).ix[0]
            if len(values) > 1:
                # check samples in group are from different biological replicates
                if len(self.df.ix[self.df.groupby(attributes).groups[key]]['biologicalReplicate'].unique()) > 1:
                    rep["biologicalReplicate"] = 0
                    rep["technicalReplicate"] = 0
                    rep["unmappedBam"] = self.df.ix[values]["unmappedBam"].tolist()
                    # append merged biological replicate to samples
                    self.samples.append(self.makeSample(rep))

    def asDataFrame(self, all=True):
        """
        Returns a `pandas.DataFrame` representation of self.
        """
        df = _pd.DataFrame([s.asSeries() for s in self.samples])

        # Filter some columns out
        if not all:
            columns = self.df.columns.tolist()
            if hasattr(df, "unmapped"):
                columns[columns.index("unmappedBam")] = "unmapped"
            df = df[["sampleName"] + columns + ["mapped"]]

        return df

    def to_csv(self, path, all=False):
        """
        Saves a csv annotation sheet from the samples.

        :param path: Path to csv file to be written.
        :type path: str
        :param all: If all sample attributes should be kept in the annotation sheet.
        :type all: bool

        :Example:

        from pipelines import SampleSheet
        sheet = SampleSheet("/projects/example/sheet.csv")
        sheet.to_csv("/projects/example/sheet2.csv")
        """
        df = self.asDataFrame(all=all)
        df.to_csv(path, index=False)


class Sample(object):
    """
    Class to model Samples basd on a pandas Series.

    :param series: Pandas `Series` object.
    :type series: pandas.Series

    :Example:

    from pipelines import Project, SampleSheet, Sample
    prj = Project("ngs")
    sheet = SampleSheet("/projects/example/sheet.csv", prj)
    s1 = Sample(sheet.ix[0])
    """
    # Originally, this object was inheriting from _pd.Series,
    # but complications with serializing and code maintenance
    # made me go back and implement it as a top-level object
    def __init__(self, series):
        from os.path import expanduser

        # Passed series must either be a pd.Series or a daugther class
        if not isinstance(series, _pd.Series):
            raise TypeError("Provided object is not a pandas Series.")
        super(Sample, self).__init__()

        # Set series attributes on self
        for key, value in series.to_dict().items():
            setattr(self, key, value)

        # Check if required attributes exist and are not empty
        self.checkValid()

        # Get name for sample:
        # this is a concatenation of all passed Series attributes except "unmappedBam"
        self.generateName()

        # Read configuration file
        with open(_os.path.join(expanduser("~"), ".pipelines_config.yaml"), 'r') as handle:
            self.config = _yaml.load(handle)

        # check if sample is to be analysed with cuts
        cuts = self.config["techniques"]["atacseq"] + self.config["techniques"]["cm"]
        self.tagmented = True if self.technique.upper() in cuts else False

        # Get track colour
        self.getTrackColour()

        # Get read type
        self.getReadType()

        # Sample dirs
        self.dirs = Paths()
        # Only when sample is added to project, can paths be added.
        # The SampleSheet object, after beeing assigned to a project, will
        # call Sample.setFilePaths()

    def __repr__(self):
        return "Sample '%s'" % self.name

    def checkValid(self):
        """
        Check if any of its important attributes is None.
        """
        req = ["technique", "genome", "unmappedBam"]

        if not all([hasattr(self, attr) for attr in req]):
            raise ValueError("Required columns for sample do not exist.")

        if any([attr == "nan" for attr in req]):
            raise ValueError("Required values for sample are empty.")

    def generateName(self):
        """
        Generates a name for the sample by joining some of its attribute strings.
        """
        self.name = "_".join(
            [str(self.__getattribute__(attr)) for attr in [
                "cellLine", "numberCells", "technique", "ip",
                "patient", "patientID", "sampleID", "treatment", "condition",
                "biologicalReplicate", "technicalReplicate",
                "experimentName", "genome"] if hasattr(self, attr) and str(self.__getattribute__(attr)) != "nan"]
        )
        self.sampleName = self.name

    def asSeries(self):
        """
        Returns a `pandas.Series` object with all the sample's attributes.
        """
        return _pd.Series(self.__dict__)

    def getReadType(self, n=10):
        """
        Gets the read type (single, paired) and length of bam file.
        Returns tuple of (readType=string, readLength=int).

        :param n: Number of reads to read to determine read type. Default=10.
        :type n: int
        """
        import subprocess as sp
        from collections import Counter

        # for samples with multiple original bams, get only first
        if type(self.unmappedBam) == list:
            bam = self.unmappedBam[0]
        else:
            bam = self.unmappedBam

        try:
            # view reads
            p = sp.Popen(['samtools', 'view', bam], stdout=sp.PIPE)

            # Count paired alignments
            paired = 0
            readLength = Counter()
            while n > 0:
                line = p.stdout.next().split("\t")
                flag = int(line[1])
                readLength[len(line[9])] += 1
                if 1 & flag:  # check decimal flag contains 1 (paired)
                    paired += 1
                n -= 1
            p.kill()
        except:
            raise IOError("Bam file does not exist or cannot be read: %s" % bam)

        # Get most abundant read length
        self.readLength = sorted(readLength)[-1]

        # If at least half is paired, consider paired end reads
        if paired > (n / 2):
            self.readType = "PE"
            self.paired = True
        else:
            self.readType = "SE"
            self.paired = False

    def setFilePaths(self):
        """
        Sets the paths of all files for this sample.
        """
        self.dirs.sampleRoot = _os.path.join(self.project.dirs.data, self.name)

        # Files in the root of the sample dir
        self.fastqc = self.dirs.sampleRoot
        self.trimlog = _os.path.join(self.dirs.sampleRoot, self.name + ".trimlog.txt")
        self.alnRates = _os.path.join(self.dirs.sampleRoot, self.name + ".alnRates.txt")
        self.alnMetrics = _os.path.join(self.dirs.sampleRoot, self.name + ".alnMetrics.txt")
        self.dupsMetrics = _os.path.join(self.dirs.sampleRoot, self.name + ".duplicates.txt")

        # Unmapped: merged bam, fastq, trimmed fastq
        self.dirs.unmapped = _os.path.join(self.dirs.sampleRoot, "unmapped")
        self.unmapped = _os.path.join(self.dirs.unmapped, self.name + ".bam")
        if self.readType == "SE":
            self.fastq = _os.path.join(self.dirs.unmapped, self.name + ".fastq")
        else:
            self.fastq1 = _os.path.join(self.dirs.unmapped, self.name + ".1.fastq")
            self.fastq2 = _os.path.join(self.dirs.unmapped, self.name + ".2.fastq")
            self.fastqUnpaired = _os.path.join(self.dirs.unmapped, self.name + ".unpaired.fastq")
        if self.readType == "SE":
            self.trimmed = _os.path.join(self.dirs.unmapped, self.name + ".trimmed.fastq")
        else:
            self.trimmed1 = _os.path.join(self.dirs.unmapped, self.name + ".1.trimmed.fastq")
            self.trimmed2 = _os.path.join(self.dirs.unmapped, self.name + ".2.trimmed.fastq")
            self.trimmed1Unpaired = _os.path.join(self.dirs.unmapped, self.name + ".1_unpaired.trimmed.fastq")
            self.trimmed2Unpaired = _os.path.join(self.dirs.unmapped, self.name + ".2_unpaired.trimmed.fastq")

        # Mapped: mapped, duplicates marked, removed, reads shifted
        self.dirs.mapped = _os.path.join(self.dirs.sampleRoot, "mapped")
        self.mapped = _os.path.join(self.dirs.mapped, self.name + ".trimmed.bowtie2.bam")
        self.filtered = _os.path.join(self.dirs.mapped, self.name + ".trimmed.bowtie2.filtered.bam")

        # Project's public_html folder
        self.bigwig = _os.path.join(self.project.dirs.html, self.name + ".bigWig")

        # Track url
        self.trackURL = "/".join([self.config["url"], self.project.name, self.name + ".bigWig"])

    def makeSampleDirs(self):
        """
        Creates sample directory structure if it doesn't exist.
        """
        for path in self.dirs.__dict__.values():
            if not _os.path.exists(path):
                _os.makedirs(path)

    def getTrackColour(self):
        """
        Get a colour for a genome browser track based on the IP.
        """
        # This is ChIP-centric, and therefore if no "ip" attrbute,
        # will just pick one color randomly from a gradient.
        import random

        if hasattr(self, "ip"):
            if self.ip in self.config["trackcolours"].keys():
                self.trackColour = self.config["trackcolours"][self.ip]
            else:
                if self.technique in ["ATAC", "ATACSEQ", "ATAC-SEQ"]:
                    self.trackColour = self.config["trackcolours"]["ATAC"]
                elif self.technique in ["DNASE", "DNASESEQ", "DNASE-SEQ"]:
                    self.trackColour = self.config["trackcolours"]["DNASE"]
                else:
                    self.trackColour = random.sample(self.config["colourgradient"], 1)[0]  # pick one randomly
        else:
            self.trackColour = random.sample(self.config["colourgradient"], 1)[0]  # pick one randomly


class ChIPseqSample(Sample):
    """
    Class to model ChIP-seq samples based on the generic Sample class (itself a pandas.Series).

    :param series: Pandas `Series` object.
    :type series: pandas.Series

    :Example:

    from pipelines import Project, SampleSheet, ChIPseqSample
    prj = Project("ngs")
    sheet = SampleSheet("/projects/example/sheet.csv", prj)
    s1 = ChIPseqSample(sheet.ix[0])
    """
    def __init__(self, series):

        # Passed series must either be a pd.Series or a daugther class
        if not isinstance(series, _pd.Series):
            raise TypeError("Provided object is not a pandas Series.")
        super(ChIPseqSample, self).__init__(series)

        # Get type of factor
        self.broad = True if self.technique in self.config["broadfactors"] else False
        self.histone = True if self.technique in self.config["histones"] else False

    def __repr__(self):
        return "ChIP-seq sample '%s'" % self.name

    def setFilePaths(self):
        """
        Sets the paths of all files for this sample.
        """
        # Inherit paths from Sample by running Sample's setFilePaths()
        super(ChIPseqSample, self).setFilePaths()

        # Files in the root of the sample dir
        self.frip = _os.path.join(self.dirs.sampleRoot, self.name + "_FRiP.txt")

        # Mapped: mapped, duplicates marked, removed, reads shifted
        # this will create additional bam files with reads shifted
        if self.tagmented:
            self.filteredshifted = _os.path.join(self.dirs.mapped, self.name + ".trimmed.bowtie2.filtered.shifted.bam")

        # Coverage: read coverage in windows genome-wide
        self.dirs.coverage = _os.path.join(self.dirs.sampleRoot, "coverage")
        self.coverage = _os.path.join(self.dirs.coverage, self.name + ".cov")

        self.qc = _os.path.join(self.dirs.sampleRoot, self.name + "_QC.tsv")
        self.qcPlot = _os.path.join(self.dirs.sampleRoot, self.name + "_QC.pdf")

        # Peaks: peaks called and derivate files
        self.dirs.peaks = _os.path.join(self.dirs.sampleRoot, "peaks")
        self.peaks = _os.path.join(self.dirs.peaks, self.name + ("_peaks.narrowPeak" if not self.broad else "_peaks.broadPeak"))
        self.peaksMotifCentered = _os.path.join(self.dirs.peaks, self.name + "_peaks.motifCentered.bed")
        self.peaksMotifAnnotated = _os.path.join(self.dirs.peaks, self.name + "_peaks.motifAnnotated.bed")

        # Motifs
        self.dirs.motifs = _os.path.join(self.dirs.sampleRoot, "motifs", self.name)


class CMSample(ChIPseqSample):
    """
    Class to model CM samples based on the ChIPseqSample class.

    :param series: Pandas `Series` object.
    :type series: pandas.Series
    """
    def __init__(self, series):

        # Use _pd.Series object to have all sample attributes
        if not isinstance(series, _pd.Series):
            raise TypeError("Provided object is not a pandas Series.")
        super(CMSample, self).__init__(series)

    def __repr__(self):
        return "CM sample '%s'" % self.name


class DNaseSample(ChIPseqSample):
    """
    Class to model DNase-seq samples based on the ChIPseqSample class.

    :param series: Pandas `Series` object.
    :type series: pandas.Series
    """
    def __init__(self, series):

        # Use _pd.Series object to have all sample attributes
        if not isinstance(series, _pd.Series):
            raise TypeError("Provided object is not a pandas Series.")
        super(DNaseSample, self).__init__(series)

    def __repr__(self):
        return "DNase-seq sample '%s'" % self.name


class ATACseqSample(ChIPseqSample):
    """
    Class to model ATAC-seq samples based on the ChIPseqSample class.

    :param series: Pandas `Series` object.
    :type series: pandas.Series
    """
    def __init__(self, series):

        # Use _pd.Series object to have all sample attributes
        if not isinstance(series, _pd.Series):
            raise TypeError("Provided object is not a pandas Series.")
        super(ATACseqSample, self).__init__(series)

    def __repr__(self):
        return "ATAC-seq sample '%s'" % self.name

    def setFilePaths(self):
        """
        Sets the paths of all files for this sample.
        """
        # Inherit paths from Sample by running Sample's setFilePaths()
        super(ChIPseqSample, self).setFilePaths()

        # Files in the root of the sample dir
        self.frip = _os.path.join(self.dirs.sampleRoot, self.name + "_FRiP.txt")

        # Mapped: mapped, duplicates marked, removed, reads shifted
        # this will create additional bam files with reads shifted
        self.filteredshifted = _os.path.join(self.dirs.mapped, self.name + ".trimmed.bowtie2.filtered.shifted.bam")

        # Coverage: read coverage in windows genome-wide
        self.dirs.coverage = _os.path.join(self.dirs.sampleRoot, "coverage")
        self.coverage = _os.path.join(self.dirs.coverage, self.name + ".cov")

        self.insertplot = _os.path.join(self.dirs.sampleRoot, self.name + "_insertLengths.pdf")
        self.insertdata = _os.path.join(self.dirs.sampleRoot, self.name + "_insertLengths.csv")
        self.qc = _os.path.join(self.dirs.sampleRoot, self.name + "_QC.tsv")
        self.qcPlot = _os.path.join(self.dirs.sampleRoot, self.name + "_QC.pdf")

        # Peaks: peaks called and derivate files
        self.dirs.peaks = _os.path.join(self.dirs.sampleRoot, "peaks")
        self.peaks = _os.path.join(self.dirs.peaks, self.name + "_peaks.broadPeak")
        self.filteredPeaks = _os.path.join(self.dirs.peaks, self.name + "_peaks.filtered.bed")


class QuantseqSample(Sample):
    """
    Class to model Quant-seq samples based on the generic Sample class (itself a pandas.Series).

    :param series: Pandas `Series` object.
    :type series: pandas.Series

    :Example:

    from pipelines import Project, SampleSheet, QuantseqSample
    prj = Project("ngs")
    sheet = SampleSheet("/projects/example/sheet.csv", prj)
    s1 = QuantseqSample(sheet.ix[0])
    """
    def __init__(self, series):

        # Passed series must either be a pd.Series or a daugther class
        if not isinstance(series, _pd.Series):
            raise TypeError("Provided object is not a pandas Series.")
        super(QuantseqSample, self).__init__(series)

    def __repr__(self):
        return "Quant-seq sample '%s'" % self.name

    def setFilePaths(self):
        """
        Sets the paths of all files for this sample.
        """
        # Inherit paths from Sample by running Sample's setFilePaths()
        super(QuantseqSample, self).setFilePaths()

        self.erccAlnRates = _os.path.join(self.dirs.sampleRoot, self.name + "_ercc.alnRates.txt")
        self.erccAlnMetrics = _os.path.join(self.dirs.sampleRoot, self.name + "_ercc.alnMetrics.txt")
        self.erccDupsMetrics = _os.path.join(self.dirs.sampleRoot, self.name + "_ercc.duplicates.txt")

        # Mapped: Tophat mapped, duplicates marked, removed
        self.dirs.mapped = _os.path.join(self.dirs.sampleRoot, "mapped")
        self.mapped = _os.path.join(self.dirs.mapped, "accepted_hits.bam")
        self.filtered = _os.path.join(self.dirs.mapped, self.name + ".trimmed.bowtie2.filtered.bam")
        # ercc alignments
        self.erccMapped = _os.path.join(self.dirs.mapped, self.name + "_ercc.bam")
        self.erccDups = _os.path.join(self.dirs.mapped, self.name + "_ercc.dups.bam")
        self.erccNodups = _os.path.join(self.dirs.mapped, self.name + "_ercc.nodups.bam")
        # kallisto pseudoalignments
        self.pseudomapped = _os.path.join(self.dirs.mapped, self.name + ".pseudoalignment.bam")

        # RNA quantification
        self.dirs.quant = _os.path.join(self.dirs.sampleRoot, "quantification")
        self.quant = _os.path.join(self.dirs.quant, "tophat-htseq_quantification.tsv")
        self.erccQuant = _os.path.join(self.dirs.quant, "tophat-htseq_quantification_ercc.tsv")
        self.kallistoQuant = _os.path.join(self.dirs.quant, "abundance.tsv")
