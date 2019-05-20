# CircRNA SnakeMake Workflow

This workflow performs circRNA detection and Quantification using three public
available tools:

	- CIRI v2.0
	- DCC
	- CircExplorer2

## Authors
	- André M. Ribeiro dos Santos (@andremrsantos)

## Requirements
	- snakemake
	- singularity

The current version requires singularity to run CIRI and DCC pipeline, since they
either lack a conda package or their package fails on execution. To avoid this issue,
I have provided a docker container to run both applications.

## Usage
	1. If you simplily want to use this workflow, download and extract
	[latest release](https://github.com/andremrsantos/circrna_sml/releases) or
	copy the current master version.

	git clone https://github.com/andremrsantos/circrna_smk.git

	2. On the directory you intend to run the pipeline (or the workflow directory)
	make sure you have a `cluster.yaml` and `config.yaml`, according to the
	templates provided. You can copy those files and edit for your environment
	specifications.

	cp <workflow-dir>/cluster.yaml <workflow-dir>/config.yaml ./
	## When setting the `seq_pe` or `seq_se`, if you don't intend to use either
	## paied or single-end files leave the path pointing to an empty or non-existent
	## directory or snakemake may try to look at all subfolders.

	3. Run the pipeline

	snakemake --snakefile <workflow-dir>/Snakefile --use-singularity --use-conda
	## Or you can use HPC profile provided to run on a torque environment
	snakemake --snakefile <workflow-dir>/Snakefile --profile <workflow-dir>/profile/hpc

## Results

the results will be distributed in the directories:
	- seq (filtered and trimmed sequences)
	- align (alignment results)
	- circrna/circexplorer (circexplorer results)
	- circrna/ciri (ciri results)
	- circrna/dcc (dcc results)
