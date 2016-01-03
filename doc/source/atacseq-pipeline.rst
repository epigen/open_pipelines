
ATAC-seq pipeline
=========================

Command-line arguments:

.. code-block:: text

	usage: atacseq-pipeline [-h] [-y SAMPLE_CONFIG] [-R] [-F] [-D]
	                        [-C CONFIG_FILE] -O PARENT_OUTPUT_FOLDER
	                        [-P NUMBER_OF_CORES]
	                        [-I INPUT_FILES [INPUT_FILES ...]] [-S SAMPLE_NAME]
	                        [-G GENOME_ASSEMBLY] [-Q SINGLE_OR_PAIRED]

	ATAC-seq pipeline.

	optional arguments:
	  -h, --help            show this help message and exit
	  -y SAMPLE_CONFIG, --sample-yaml SAMPLE_CONFIG
	                        Yaml config file with sample attributes.
	  -R, --recover         Recover mode, overwrite locks
	  -F, --fresh-start     Fresh start mode, overwrite all
	  -D, --dirty           Make all cleanups manual
	  -C CONFIG_FILE, --config CONFIG_FILE
	                        pipeline config file in YAML format; relative paths
	                        are considered relative to the pipeline script.
	                        defaults to atacseq.yaml
	  -O PARENT_OUTPUT_FOLDER, --output-parent PARENT_OUTPUT_FOLDER
	                        parent output directory of the project (required). The
	                        sample_name argument will be appended to this folder
	                        for output
	  -P NUMBER_OF_CORES, --cores NUMBER_OF_CORES
	                        number of cores to use for parallel processes
	  -I INPUT_FILES [INPUT_FILES ...], --input INPUT_FILES [INPUT_FILES ...]
	                        one or more input files (required)
	  -S SAMPLE_NAME, --sample-name SAMPLE_NAME
	                        unique name for output subfolder and files (required)
	  -G GENOME_ASSEMBLY, --genome GENOME_ASSEMBLY
	                        identifier for genome assempbly (required)
	  -Q SINGLE_OR_PAIRED, --single-or-paired SINGLE_OR_PAIRED
	                        single or paired end? default: single

Output files:

- unmapped:
	- ``sample``.fastq [1]_ :  unmapped reads in fastq format.
	- ``sample``.trimmed.fastq [1]_ : unmapped reads with trimmed adapters in fastq format.
- mapped:
	- ``sample``.trimmed.bowtie2.bam : mapped reads in bam format. Sorted by read position. 
	- ``sample``.trimmed.bowtie2.bai : index of above.
	- ``sample``.trimmed.bowtie2.filtered.bam : mapped reads filtered for high quality [2]_ in bam format. Sorted by read position.
	- ``sample``.trimmed.bowtie2.filtered.bai : index of above.
- peaks:
	- ``sample`` _peaks.narrowPeak :  MACS2 output of called peaks in NarrowPeak format.
	- ``sample`` _peaks.xls : MACS2 output of called peaks in MS Excel format.
- ``sample``.alnStats.txt : Summary file of the flags in the mapped reads.
- ``sample``.alnMetrics.txt : Log file of the mapping tool.
- ``sample``.dups.txt : Log file of the duplicate removal tool.

.. rubric:: Footnotes

.. [1] File is removed if run is completed successfully.
.. [2] High quality reads non-duplicate, uniquely mapping reads with mapping quality >=30 and for paired-end reads only proper pairs.
