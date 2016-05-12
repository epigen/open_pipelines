from looper.models import Sample
import pandas as _pd
import os as _os
import yaml as _yaml
import warnings as _warnings


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

		# Passed series must either be a pd.Series or a daughter class
		if not isinstance(series, _pd.Series):
			raise TypeError("Provided object is not a pandas Series.")
		super(ChIPseqSample, self).__init__(series)

		# Get type of factor
		# TODO: get config file specifying broad/narrow factors
		# e.g. self.broad = True if self.ip in self.prj.config["broadfactors"] else False
		# NS: I wrapped this in a try block because this makes it require that
		# 'ip' be defined, which may not be the case (like for Input samples)

		try:
			self.broad = True if any([ip in self.ip.upper() for ip in ["H3K27me3", "H3K36me3"]]) else False
			self.histone = True if any([ip in self.ip.upper() for ip in ["H3", "H2A", "H2B", "H4"]]) else False
		except:
			pass
		self.make_sample_dirs()

	def __repr__(self):
		return "ChIP-seq sample '%s'" % self.sample_name

	def set_file_paths(self):
		"""
		Sets the paths of all files for this sample.
		"""
		# Inherit paths from Sample by running Sample's set_file_paths()
		super(ChIPseqSample, self).set_file_paths()

		# Files in the root of the sample dir
		self.fastqc = _os.path.join(self.paths.sample_root, self.sample_name + ".fastqc.zip")
		self.trimlog = _os.path.join(self.paths.sample_root, self.sample_name + ".trimlog.txt")
		self.aln_rates = _os.path.join(self.paths.sample_root, self.sample_name + ".aln_rates.txt")
		self.aln_metrics = _os.path.join(self.paths.sample_root, self.sample_name + ".aln_metrics.txt")
		self.dups_metrics = _os.path.join(self.paths.sample_root, self.sample_name + ".dups_metrics.txt")

		# Unmapped: merged bam, fastq, trimmed fastq
		self.paths.unmapped = _os.path.join(self.paths.sample_root, "unmapped")
		self.unmapped = _os.path.join(self.paths.unmapped, self.sample_name + ".bam")
		self.fastq = _os.path.join(self.paths.unmapped, self.sample_name + ".fastq")
		self.fastq1 = _os.path.join(self.paths.unmapped, self.sample_name + ".1.fastq")
		self.fastq2 = _os.path.join(self.paths.unmapped, self.sample_name + ".2.fastq")
		self.fastqUnpaired = _os.path.join(self.paths.unmapped, self.sample_name + ".unpaired.fastq")
		self.trimmed = _os.path.join(self.paths.unmapped, self.sample_name + ".trimmed.fastq")
		self.trimmed1 = _os.path.join(self.paths.unmapped, self.sample_name + ".1.trimmed.fastq")
		self.trimmed2 = _os.path.join(self.paths.unmapped, self.sample_name + ".2.trimmed.fastq")
		self.trimmed1Unpaired = _os.path.join(self.paths.unmapped, self.sample_name + ".1_unpaired.trimmed.fastq")
		self.trimmed2Unpaired = _os.path.join(self.paths.unmapped, self.sample_name + ".2_unpaired.trimmed.fastq")

		# Mapped: mapped, duplicates marked, removed, reads shifted
		self.paths.mapped = _os.path.join(self.paths.sample_root, "mapped_"+self.genome)
		self.mapped = _os.path.join(self.paths.mapped, self.sample_name + ".trimmed.bowtie2.bam")
		self.filtered = _os.path.join(self.paths.mapped, self.sample_name + ".trimmed.bowtie2.filtered.bam")

		# Files in the root of the sample dir
		self.frip = _os.path.join(self.paths.sample_root, self.sample_name + "_FRiP.txt")

		# Coverage: read coverage in windows genome-wide
		self.paths.coverage = _os.path.join(self.paths.sample_root, "coverage_" + self.genome)
		self.coverage = _os.path.join(self.paths.coverage, self.sample_name + ".cov")

		self.insertplot = _os.path.join(self.paths.sample_root, self.name + "_insertLengths.pdf")
		self.insertdata = _os.path.join(self.paths.sample_root, self.name + "_insertLengths.csv")

		self.qc = _os.path.join(self.paths.sample_root, self.sample_name + "_qc.tsv")
		self.qc_plot = _os.path.join(self.paths.sample_root, self.sample_name + "_qc.pdf")

		# Peaks: peaks called and derivate files
		self.paths.peaks = _os.path.join(self.paths.sample_root, "peaks_"+self.genome)
		self.peaks = _os.path.join(self.paths.peaks, self.sample_name + ("_peaks.narrowPeak" if not self.broad else "_peaks.broadPeak"))
		self.peaks_motif_centered = _os.path.join(self.paths.peaks, self.sample_name + "_peaks.motif_centered.bed")
		self.peaks_motif_annotated = _os.path.join(self.paths.peaks, self.sample_name + "_peaks._motif_annotated.bed")

		# Motifs
		self.paths.motifs = _os.path.join(self.paths.sample_root, "motifs")


