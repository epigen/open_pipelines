#! /usr/bin/env python

import sys
import os


# take care of extra required modules depending on Python version
extra = {}


try:
	from setuptools import setup
	if sys.version_info < (2, 7):
		extra['install_requires'] = ['argparse']
	if sys.version_info >= (3,):
		extra['use_2to3'] = True
except ImportError:
	from distutils.core import setup
	if sys.version_info < (2, 7):
		extra['dependencies'] = ['argparse']


def get_static(name, keep=None):
	""" Determine additional files to include with the package. """
	pkg_path = os.path.dirname(os.path.realpath(__file__))
	static = [os.path.join(name, f) for f in
			  os.listdir(os.path.join(pkg_path, name))]
	return static if keep is None else list(filter(keep, static))

# Additional files to include with package
pipeline_configs = get_static("pipelines",
	keep=lambda fpath: os.path.splitext(fpath)[1] in [".yaml", ".yml"])
scripts = get_static(
	os.path.join("pipelines", "tools"), keep=lambda fpath: '.' in fpath)


with open("VERSION", 'r') as versionfile:
	version = versionfile.read().strip()


with open("requirements-pypi.txt", 'r') as reqs_file:
	reqs = [l.rstrip() for l in reqs_file.readlines() if not l.startswith("#")]


setup(
	name="pipelines",
	packages=["pipelines"],
	version=version,
	description="NGS pipelines in Python.",
	long_description=open('README.md').read(),
	classifiers=[
		"Development Status :: 3 - Alpha",
		"License :: OSI Approved :: GNU General Public License v2 or later (GPLv2+)",
		"Programming Language :: Python :: 2.7",
		"Topic :: Scientific/Engineering :: Bio-Informatics"
	],
	keywords="bioinformatics, sequencing, ngs, ATAC-Seq, ChIP-seq, RNA-seq, RRBS, WGBS",
	url="https://github.com/epigen/pipelines",
	author=u"Nathan Sheffield, Johanna Klughammer, Andre Rendeiro, Charles Dietz",
	license="GPL2",
	install_requires=reqs,
	entry_points={
		"console_scripts": [
			'atacseq_pipeline = pipelines.atacseq:main',
			'chipseq_pipeline = pipelines.chipseq:main',
			'quantseq_pipeline = pipelines.quantseq:main',
			'starrseq_pipeline = pipelines.starrseq:main'
		],
	},
	scripts=scripts,
	data_files=[("configs", pipeline_configs)],
	include_package_data=True,
	**extra
)
