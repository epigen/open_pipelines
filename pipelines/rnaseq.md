# RNA-seq pipeline

This simple [RNA-seq pipeline](rnaseq.py) processes most RNA-seq protocols.
It uses the ultra-fast pseudomapping-based kallisto to generate transcript abundance estimations from unmapped reads.


## Installation and configuration

### Prequisites

**Python packages**. This pipeline uses [pypiper](https://github.com/epigen/pypiper) and [looper](https://github.com/epigen/looper). You can do a user-specific install of these like this:

```
pip install --user https://github.com/epigen/pypiper/zipball/v0.6
pip install --user https://github.com/epigen/looper/zipball/v0.7.2
```

**Required executables**. You will need some common bioinformatics software installed. The list is specified in the pipeline configuration file ([rnaseq.yaml](rnaseq.yaml)) tools section.

**Static files**. This pipeline requires static files which are specified in the [pipeline configuration file](rnaseq.yaml).


## Usage

 - Clone the pipeline repository: `git clone git@github.com:epigen/open_pipelines.git`;
 - Adapt the [pipeline configuration file](rnaseq.yaml) to point to specific software if needed;
 - Create a sample annotation sheet containing the variables `sample_name`, `protocol`, and `organism`;
 - Create a project configuration file that points to the [pipeline interface file](../pipeline_interface.yaml) and the sample annotation sheet;
 - Run pipelines using looper `looper run project_config.yaml`.

More detailed instructions or creating a project configuration file and sample annotation sheet canbe found in the [Looper documentation](http://looper.readthedocs.io).

In the particular case of the RNA-seq pipeline, one special column in the annotation sheet can be used to pair samples for peak calling. Add a column named "compare_sample" containing the name ("sample_name" column) of the sample to use as background.


## Pipeline steps

### Merging input files

If given more than one BAM file as input, the pipeline will merge them before begining processing. The merged, unmapped inpu BAM file will be output in `$sample_name/unmapped`. This file is temporary and will be removed if the pipeline finishes successfully.

### Sequencing read quality control with FASTQC

[FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) is ran on the unaligned input BAM files for quality control.

An HTML report and accompaning zip file will be output in the root directory of each sample.

### Read trimming

Reads are trimmed for adapters prior to alignment.

Adapter sequences to be trimmed can be specified in a FASTA file which is stated in the [pipeline configuration file](rnaseq.yaml) under `resources: adapters`.

Two trimming programs can be selected: [trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic) and [skewer](https://github.com/relipmoc/skewer) in the [pipeline configuration file](rnaseq.yaml) under `parameters: trimmer`. While rigorous benchmarking of both trimmers could be performed, the reason to use skewer is its considerable speed compared with trimmomatic and the fact that it is available as a binary executable rather than a Java jar.

These produce FASTQ files for each read pair and one file with unpaired reads, which are stored under `$sample_name/unmapped/`. These files are temporary and will be removed if the pipeline finishes sucessfully.

### Expression quantification trimming

This pipeline uses [Kallisto](https://pachterlab.github.io/kallisto/) for transcript quantification without the need of alignment.

Kallisto needs a transcriptome index which should be specified in the [pipeline configuration file](rnaseq.yaml) under `resources: kallisto_index`. This can be easily created with the `kallisto index` command.

A TSV file containing estimation of transcript abundances is created under `$sample_name/kallisto/`.


## Collecting statistics from pipeline runs

You can easily collect statistics from all runs using looper: `lopper summarize project_config.yaml`


### Quality control and Statistics

Due to the minimal pipeline size/steps, statistics produced are limited.
Here are the reported statistics and their description:

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
 - `transcripts`: number of transcripts quantified
 - `zero-count_transcripts`: number of transcripts with estimated zero counts.
 - `non-zero-count_transcripts`: number of transcripts with more than zero estimated counts.
 - `log2tpm_mean`: Mean expression (log2(1 + TPM)) of all transcripts
 - `log2tpm_median`: Median expression (log2(1 + TPM)) of all transcripts
 - `log2tpm_iqr`: Interquantile range (IQR) of the expression (log2(1 + TPM)) of all transcripts
 - `non-zero_log2tpm_mean`: Mean expression (log2 TPM) of transcripts with >0 counts
 - `non-zero_log2tpm_median`: Mean expression (log2 TPM) of transcripts with >0 counts
 - `non-zero_log2tpm_iqr`: Interquantile range (IQR) of the expression (log2(1 + TPM)) of transcripts with >0 counts
 - `Time`: pipeline run time
 - `Success`: time pipeline run finished


## Contributing

Pull requests welcome. Active development should occur in the development branch.

## Contributors

* Andre Rendeiro, arendeiro@cemm.oeaw.ac.at
