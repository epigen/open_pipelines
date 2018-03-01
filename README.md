# Pipelines

These are some of the NGS pipelines developed and used in the [Bock lab](http://medical-epigenomics.org) at [CeMM](http://cemm.at).

This repository contains pipelines for processing NGS data and associated scripts used by them (in the [/pipelines/tools](pipelines/tools) subdirectory).
Pipelines here are configured to work with [`looper`](http://looper.readthedocs.io) and use [`pypiper`](http://pypiper.readthedocs.io) (see the corresponding repositories). These pipelines work with metadata formatted as a [portable encapsulated project](http://pepkit.github.io).

# Installation and usage
1. Install [`looper`](http://looper.readthedocs.io) and [`pypiper`](http://pypiper.readthedocs.io): 
  - `pip install https://github.com/epigen/looper/zipball/v0.7.2`
  - `pip install https://github.com/epigen/pypiper/zipball/v0.6`
2. Clone this repository:
  - `git clone git@github.com:epigen/open_pipelines.git`
3. Produce a configuration file for your project ([see here how to do it](http://looper.readthedocs.io/en/latest/define-your-project.html)).
4. Link the pipelines to your project by pointing to the pipeline interfaces configuration:

      ```yaml
          metadata:
              pipeline_interfaces: open_pipelines/pipeline_interface.yaml
      ```
5. Run all jobs using looper:
  - `looper run project/metadata/project_config.yaml`


If you are just _using a pipeline_ in a project, and you are not _developing the pipeline_, you should treat this cloned repo as read-only, frozen code, which should reside in a shared project workspace. There should be only one clone for the project, to avoid running data under changing pipeline versions (you should not pull any pipeline updates unless you plan to re-run the whole thing).


## Test data for pipelines

Small example data for several pipeline types is available in the [microtest repository](https://github.com/epigen/microtest/)


# Pipeline documentation
See the dedicated pages detailing the steps and output of each of the pipelines:
 - [ATAC-seq and DNAse-seq](pipelines/atacseq.md)
 - [ChIP-seq and ChIPmentation](pipelines/chipseq.md)
 - [RNA-seq](pipelines/rnaseq.md)

The remaining pipelines will get their dedicated documentation soon.


# Contributing

We appreciate and encourage contributions to existing pipelines or submitions of new ones.

Simply clone the repository, make your changes and create a pull request with the changes in the `dev` branch.
