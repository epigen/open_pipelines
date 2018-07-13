# Hi-C and Hi-ChIP pipeline

This [pipeline](hic.py) processes Hi-C and HiChIP data.
It takes raw unaligned reads and produces mapped valid interaction pairs, raw and normalized interaction matrices, and peak and loop calls.
Right now it won't report any performance statistics, but this will be implemented briefly.

> TODO: add statistic collection


## Installation and configuration

### Prequisites

**Python packages**. This pipeline uses [pypiper](https://github.com/epigen/pypiper) and [looper](https://github.com/epigen/looper). You can do a user-specific install of these like this:

```
pip install --user https://github.com/epigen/pypiper/zipball/v0.6
pip install --user https://github.com/epigen/looper/zipball/v0.7.2
```

currently these exact versions are needed.

**Required executables**. You will need some common bioinformatics software installed. The list is specified in the pipeline configuration file ([hic.yaml](hic.yaml)) tools section.


**Static files**. This pipeline requires static files which are specified in the [pipeline configuration file](hic.yaml).


## Usage


### Standalone usage

This pipeline can be used on its own. For that you need to pass several inputs to the command-line interface (CLI):

> TODO

### With Looper

 - Clone the pipeline repository: `git clone git@github.com:epigen/open_pipelines.git`;
 - Adapt the [pipeline configuration file](hic.yaml) to point to required software needed by the pipeline. All runnables are values under the "tools" section;
 - Create a sample annotation sheet containing the variables `sample_name`, `protocol`, `organism` and `read1` and `read2` (for paired FASTQ files if needed);
 - Create a project configuration file that points to the [pipeline interface file](../pipeline_interface.yaml) and the sample annotation sheet;
 - Run pipelines using looper `looper run project_config.yaml`.

More detailed instructions or creating a project configuration file and sample annotation sheet can be found in the [Looper documentation](http://looper.readthedocs.io).


## Pipeline steps

### Input files

Right now the pipeline takes only BAM files from paired-end sequencing as input.

### Merging input files

If given more than one file if given to each `input` parameter, the pipeline will merge them before begining processing. The merged, unmapped input BAM file will be output in `$sample_name/unmapped`. These files are temporary and will be removed if the pipeline finishes successfully.

### Sequencing read quality control with FASTQC

[FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) is ran on the unaligned input FASTQ files for quality control.
HTML reports and accompaning zip files will be output in the root directory of each sample.

In paired-end mode, a FastQC report for each read will be produced but aggregated statistics are also available.

### HiC-Pro pipeline

See https://github.com/nservant/HiC-Pro

The main output of HiC-Pro is a ["allValidPairs"](https://github.com/nservant/HiC-Pro/blob/master/doc/RESULTS.rst) file containing valid read pairs fragment-aware annotated.
In addition, raw and normalized matrix files are produced.

> TODO: parse HiC-Pro statistics as pipeline runs.

### Format convertion

Due to the recent explosion in file formats to store contact data, output from HiC-Pro is converted to several formats:

 - [Juicebox's ".hic" format](https://github.com/theaidenlab/juicebox/wiki)
 - [Cooler format](http://cooler.readthedocs.io/)
 - [Pairix-indexed bedGraph format](https://github.com/4dn-dcic/pairix)


### Peak calling

We use MACS2 to call peaks from paired-end tags. This can be used for quality control inspection and comparing with a parallel ChIP-seq experiment.

### Loop calling

We use [cLoops](https://github.com/YaqiangCao/cLoops) and [hichipper](http://hichipper.readthedocs.io/) to call loops from HiC data.
In the future, support for [HiCCups](https://github.com/theaidenlab/juicer/wiki/HiCCUPS) loops is planned.


## Collecting statistics from pipeline runs

You can easily collect statistics from all runs using looper: `lopper summarize project_config.yaml`


### Statistics

> TODO

## Sample quality, wet lab troubleshooting and advice

> TODO

## Contributing

Pull requests welcome. Active development should occur in the development branch.

## Contributors

* Andre Rendeiro, arendeiro@cemm.oeaw.ac.at
