# Pipelines

These are some of the NGS pipelines developed and used in the [Bock lab](http://medical-epigenomics.org) at [CeMM](http://cemm.at).

This repository contains pipelines for processing NGS data and associated scripts used by them.

These pipelines are meant to work with [`looper`](https://github.com/epigen/looper/) and use [`pypiper`](https://github.com/epigen/pypiper/) (see the corresponding repositories).

## Running pipelines with Looper

If you have set up `looper` and `pypiper` and [set up your project configuration](http://looper.readthedocs.io/en/latest/define-your-project.html), you simply have to [link these pipelines to your project.](http://looper.readthedocs.io/en/latest/pipeline-interface.html) (step 4 in the following list).

1. Install [`looper`](https://github.com/epigen/looper/) and [`pypiper`](https://github.com/epigen/pypiper/) (you might need to add the `--user` option in those pip commands): 
  - `pip install https://github.com/epigen/looper/zipball/master`
  - `pip install https://github.com/epigen/pypiper/zipball/master`
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


## Pipeline steps and outputs

More thorough documentation of each pipeline and its outputs will follow soon.


## Contributing

We appreciate and encourage contributions to existing pipelines or submitions of new ones.

Simply clone the repository, make your changes and create a pull request with the changes in the `dev` branch.
