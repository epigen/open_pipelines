# Pipelines

[More thorough documentation of each pipeline and its outputs will follow soon.]

This repository contains pipelines for processing NGS data and associated scripts used by them (in the [/pipelines/tools](pipelines/tools) subdirectory).
Pipelines here are configured to work with [`looper`](http://looper.readthedocs.io) and use [`pypiper`](http://pypiper.readthedocs.io) (see the corresponding repositories). These pipelines work with metadata formatted as a [portable encapsulated project](http://pepkit.github.io).

# Installing
1. Install [`looper`](http://looper.readthedocs.io) and [`pypiper`](http://pypiper.readthedocs.io): 
  - `pip install https://github.com/epigen/looper/zipball/master`
  - `pip install https://github.com/epigen/pypiper/zipball/master`
2. Clone this repository:
  - `git clone git@github.com:epigen/open_pipelines.git`
3. Produce a config file (it just has a bunch of paths).
4. Go!

If you are just _using a pipeline_ in a project, and you are not _developing the pipeline_, you should treat this cloned repo as read-only, frozen code, which should reside in a shared project workspace. There should be only one clone for the project, to avoid running data under changing pipeline versions (you should not pull any pipeline updates unless you plan to re-run the whole thing).


# Running pipelines

We use `Looper` to run pipelines. This just requires a yaml format config file passed as an argument, which contains all the settings required.

This can, for example, submit each job to SLURM (or SGE, or run them locally).

```bash
looper run metadata/config.yaml
```

### Running on test data

Small example data for several pipeline types is available in the [microtest repository](https://github.com/epigen/microtest/)


# Post-pipeline processing

Once a pipeline has been run (or is running), you can do some post-processing on the results. 
`Looper` has a command to do this: `looper summarize`, which collects statistics produced by the pipelines for all submitted samples.


# Developing pipelines
If you plan to create a new pipeline or develop existing pipelines, consider cloning this repo to your personal space, where you do the development. Push changes from there. Use this personal repo to run any tests or whatever, but consider making sure a project is run from a different (frozen) clone, to ensure uniform results.
