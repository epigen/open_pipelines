# Pipelines

Note! Documentation is under heavy work still.

# Usage

Prerequisites: These pipelines use [pypiper](https://github.com/epigen/pypiper) to run a pipeline for a single sample, and [looper](https://github.com/epigen/looper) to handle multi-samples job submission for a project. You can do a user-specific install of both like this:

```
pip install --user https://github.com/epigen/pypiper/zipball/master
pip install --user https://github.com/epigen/looper/zipball/master
```


### Option 1 (install this package)

```bash
pip install https://github.com/epigen/open_pipelines/zipball/master

```

### Option 2 (clone the repository)

```bash
git clone git@github.com:epigen/pipelines.git
```

Add the location of the cloned repository to your looper [project configuration](http://looper.readthedocs.io/en/latest/inputs.html#project-config-file) file.

# Developing pipelines

If you plan to create a new pipeline or develop existing pipelines, clone this repo to your personal space, where you do the development.
You can tell looper about which pipelines are available through the [pipeline interface configuration](http://looper.readthedocs.io/en/latest/inputs.html#pipeline-interface-yaml) file.
