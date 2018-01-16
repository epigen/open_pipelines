# ATAC-seq pipeline

The [ATAC-seq pipeline](atacseq.py) processes ATAC-seq and DNAse-seq data.
It does adapter trimming, mapping, peak calling, and creates a bigwig file for visualization.
For each step statistics are collected and several library quality metrics are computed.


## Installation and configuration

### Prequisites

**Python packages**. This pipeline uses [pypiper](https://github.com/epigen/pypiper) and [looper](https://github.com/epigen/looper). You can do a user-specific install of these like this:

```
pip install --user https://github.com/epigen/pypiper/zipball/v0.6
pip install --user https://github.com/epigen/looper/zipball/v0.7.2
```

**Required executables**. You will need some common bioinformatics software installed. The list is specified in the pipeline configuration file ([atacseq.yaml](atacseq.yaml)) tools section.

**Static files**. This pipeline requires static files which are specified in the [pipeline configuration file](atacseq.yaml).


## Usage

 - Clone the pipeline repository: `git clone git@github.com:epigen/open_pipelines.git`;
 - Adapt the [pipeline configuration file](atacseq.yaml) to point to required software needed by the pipeline. All runnables are values under the "tools" section;
 - Create a sample annotation sheet containing the variables `sample_name`, `protocol`, and `organism`;
 - Create a project configuration file that points to the [pipeline interface file](../pipeline_interface.yaml) and the sample annotation sheet;
 - Run pipelines using looper `looper run project_config.yaml`.

More detailed instructions or creating a project configuration file and sample annotation sheet canbe found in the [Looper documentation](http://looper.readthedocs.io).


## Pipeline steps

### Merging input files

If given more than one BAM file as input, the pipeline will merge them before begining processing. The merged, unmapped inpu BAM file will be output in `$sample_name/unmapped`. This file is temporary and will be removed if the pipeline finishes successfully.

### Sequencing read quality control with FASTQC

[FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) is ran on the unaligned input BAM files for quality control.

An HTML report and accompaning zip file will be output in the root directory of each sample.

### Read trimming

Reads are trimmed for adapters prior to alignment. 

Adapter sequences to be trimmed can be specified in a FASTA file which is stated in the [pipeline configuration file](atacseq.yaml) under `resources: adapters`.

Two trimming programs can be selected: [trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic) and [skewer](https://github.com/relipmoc/skewer) in the [pipeline configuration file](atacseq.yaml) under `parameters: trimmer`. While rigorous benchmarking of both trimmers could be performed, the reason to use skewer is its considerable speed compared with trimmomatic and the fact that it is available as a binary executable rather than a Java jar.

These produce FASTQ files for each read pair and one file with unpaired reads, which are stored under `$sample_name/unmapped/`. These files are temporary and will be removed if the pipeline finishes sucessfully.

### Read alignment

Read alignment is done with [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml). This file is then sorted and indexed with [sambamba](http://lomereiter.github.io/sambamba/). Alignment is done with the `--very-sensitive` parameters and in paired-end mode, only fragments with at most 2000bp are allowed. An aligned BAM file is produced in sample output directory `$sample_name/mapped`.

### Read filtering

Reads are filtered quite stringently, specially in paired-end mode. This fits the many clinically-oriented projects in the Bock lab but may be adjusted. 

The minimum accepted mapping quality in 30 (phred scale) and in paired-end mode only concordant proper pairs are kept. In addition, optical duplicates are removed with [sambamba](http://lomereiter.github.io/sambamba/) and only nuclear genome reads are kept.

### Transposition event reconstruction

The location of Tn5 transposition events can be reconstituted by shifting the 5' position of the mapped reads by a known distance.
The pipeline will perform this step and create a BAM file with shifted read positions which can be used for several downstream applications (e.g. footprinting).
However, do not use this file for any other application.

### Peak calling

Peaks can be called with [MACS2](https://github.com/taoliu/MACS) or [spp](https://github.com/kundajelab/phantompeakqualtools) but since ATAC-seq peak calling with MACS2 works (surprisingly) very well and spp is not really a command-line program but a script, I strongly recommend using MACS2 only.

MACS2 is run with the `--no-model --extsize 147` parameters and outputs peak calls in the `$sample_name/peaks` directory in both narrowPeak and excel formats.
For genomes with annotation available, I strongly recommend filtering the peak calls using a blacklist file such as the ENCODE (https://sites.google.com/site/anshulkundaje/projects/blacklists) - this is performed if such a list is provided in the pipeline's configuration file.

### Quality control

#### bigWig file

A bigWig file with signal along the genome will be produced which can be used for visualizations in genome browsers. This requires some [UCSC tools](http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/) to be available, particularly `bedGraphToBigWig`.
This file can be normalized to the total library size if specified in the [pipeline configuration file](atacseq.yaml) under `parameters: normalize_tracks`. The normalization can take a normalization factor `parameters: norm_factor` under consideration which will make the scalling of the bigWig file comparable between samples. This can however create files with numerical scales which are impractical and therefore producing unnormalized bigWig files or making a posterior normalization with all samples jointly can sometimes be preferable.

#### Fragment distributions

Library fragment distribution can be an indicator of library quality and refects the abundance of nucleosome-free or nucleosome fragments. **For paired-end samples** a plot of fragment distributions will be produced and stored in the sample's root directory. A mixture of distributions will try to be fit to the observed distribution: an inverse exponential distribution modeling the nucleosome-free fragments and four gaussians modeling various nucleosomal distributions.

#### FRiP

One of the most used measures of signal-to-noise ratio in an ChIP-seq library is the fraction of reads in peaks (FRiP). This is simply the fraction of filtered reads that is overlapping peaks (from the own sample or from an oracle region list) over all filtered reads. For a more complete description and reasoning behind this metric see [Landt et al 2012](https://dx.doi.org/doi:10.1101/gr.136184.111). Files holding the number of reads overlapping peaks will be output in the sample's root directory and the actual FRiP value(s) will be reported in the statistics. This can also be calculated from an oracle region set (e.g. Ensembl regulatory build annotations) if provided in the pipeline's configuration file.

#### Cross-correlation enrichment

Two additional quality metric will be computed: the normalized strand cross-correlation (NSC) and relative strand cross-correlation (RSC). These metrics were defined as a measure of signal-to-noise ratio in ChIP-seq experiments and don't have the same meaning for ATAC-seq data (see [Landt et al 2012](https://dx.doi.org/doi:10.1101/gr.136184.111)). However, a passing quality ATAC-seq library should have a NSC higher than 1 and generally a higher RSC value is related with better quality samples.


## Collecting statistics from pipeline runs

You can easily collect statistics from all runs using looper: `lopper summarize project_config.yaml`


### Statistics

 - `fastqc_GC_perc`: GC percentage of all sequenced reads from FASTQC report.
 - `fastqc_read_length`: read length as determined from FASTQC report.
 - `fastqc_total_pass_filter_reads`: number of pass filter reads from FASTQC report.
 - `fastqc_poor_quality_perc`: number of poor quality reads from FASTQC report
 - `trim_short_perc`: percentage of reads dropped because of too short length after trimming
 - `trim_empty_perc`: percentage of reads dropped because empty after trimming
 - `trim_trim_loss_perc`: percentage of reads lost during trimming
 - `trim_surviving_perc`: percentage of reads surviving after trimming
 - `trim_trimmed_perc`: percentage of reads that were trimmed
 - `trim_untrimmed_perc`: percentage of reads that were untrimmed
 - `unique_aligned_perc`: percentage of reads that aligned uniquely
 - `perc_aligned`: percentage of reads that were aligned
 - `not_aligned_perc`: percentage of reads that were not aligned
 - `multiple_aligned_perc`: percentage of reads that aligned to multiple reference locations
 - `MT_filtered_paired_ends`: number of mitochondrial paired-end reads available after filtering
 - `MT_duplicate_percentage`: percentage of mitochondrial reads that are duplicate
 - `MT_filtered_single_ends`: number of mithochondrial single-end reads available after filtering
 - `filtered_paired_ends`: number of paired-end reads available after filtering
 - `duplicate_percentage`: percentage of duplicate reads
 - `filtered_single_ends`: number of single-end reads available after filtering
 - `NSC`: normalized strand cross-correlation
 - `RSC`: relative strand cross-correlation
 - `peaks`: number of called peaks
 - `filtered_peaks`: number of called peaks filtered from genome-specific blacklisted regions, if provided
 - `frip`: fraction of reads in peaks (FRiP) called from the sample
 - `oracle_frip`: fraction of reads in peaks (FRiP) from an external oracle region set
 - `Time`: pipeline run time 
 - `Success`: time pipeline run finished


## Sample quality, wet lab troubleshooting and advice


Here are a few comments on issues that can condition the quality of an ATAC-seq library by our experience (by rough order of importance):
 - A crucial factor for a successful library is the quality of the starting sample. ATAC-seq requires viable cells, so make sure to use only living cells by either sorting (good to increase purity too) or simply careful decantation - even if processing a fresh sample. Some cell types/lines are more sensitive to freezing/thawing and this should be done with great care.
 - One issue that can reduce the efficiency of the experiment considerably is the amount of mitochondrial DNA that is tagmented, specially in cell types with high mitochondrial content. This has been improved in the FastATAC protocol (Corces et al, 2016) by the use of digitonin as a permeabilizing agent as opposed to more harsher detergents which often lyse mitochondria too. High DNA mitochondrial content will not increase noise levels in the experiment but will reduce the amount of useful nuclear reads. Typical mitochondrial read fractions in the first protocol version were between 20-40% and with the FastATAC protocol have been dramatically reduced to 1-8%.
 - Another issue that sometimes occurs (specially for new users) is how to size select the libraries. This can often be the confluence point of several factors/problems: too small fragments can have origin in some adapter dimers or from genomic DNA that has been too tagmented. The first case is the most serious as these fragments are not alignable and easily outcompete larger fragments for sequencing - this can be easily spotted by assessing the quality of the raw reads with FASTQC for example. The later is not necessarily a problem since the smallest sequencable fragment from an ATAC-seq library is usually larger than 28-31bp which is mappable. It can however reduce the return of sequencing with large read lengths in paired-end mode. This can also be monitored by observing the abundance of Illumina Nextera sequences in the reads starting with basepair ~30 and progressively increasing. Another source of small fragments can be the tagmentation of naked DNA due to excessive cell lysis or death and is easy to spot on a bioanalyzer. A good library will contain an initial negative exponential distribution of nucleosome free fragments and a series of gaussian-like distributions of nucleossomal fragments (see Buenrostro et al Nat Methods 2013).
 - PCR overamplification of a library will result in a low complexity library (high redundancy) and should be avoided. Always perform an initial qPCR amplification of an aliquot of the library to estimate the number of cycles needed to amplify each individual library. While cumbersome (each library may need a specific number of cycles), this process is quite crucial.
 - Lastly, the type of Tn5 enzyme used can also be a source of biases in the data. The standard protocol uses the Tn5 enzyme from the Illumina Nextera kit, but one can also use Tn5 recombinant protein produced in-house. The later has to be loaded with adapters which adds a few extra steps to the protocol including an important step of removal of unloaded adapter sequences. If excess adapters are present this will greatly reduce the efficiency of the experiment. However, when comparing the quality (in quantitative but also qualitative terms) of the Illumina or in-house Tn5, we have detectable but only subtle differences, which indicated that both could be used for profiling accessible chromatin. However, avoid using the two types of enzyme in a given project to minimize bias just to be on the safe side. Monitor also the date of use of the transposase enzyme as its efficiency tends to decay with time.

## Contributing

Pull requests welcome. Active development should occur in the development branch.

## Contributors

* Andre Rendeiro, arendeiro@cemm.oeaw.ac.at
