# ONT Sequencing Summary Statistics

This pipeline reads a colletion of ONT FASTQ files and generates summary statistics from them, such as the number
of reads, number of bases, and average read lengths.

## Run design

The pipeline is designed to be installed in a separate directory from where it is run. This documentation will
refer to `PIPELINE_DIRECTORY` as the directory where the pipeline is installed (location where `git pull` was run) and
`RUN_DIRECTORY` as the directory where the pipeline is run. `RUN_DIRECTORY` should be initially an empty directory, and
configuration files will be placed here.

The pipeline is run with Snakemake. Examples use `-j 12` for Snakemake commands to run 12 concurrent jobs. You may
adjust this number for your needs.

## Installation

### Prerequisites

The pipeline requires Snakemake and a Python environment to run.

Snakemake should be at least version 5.4, but version 7 or later is recommended (tested on 7.25.0). The pipeline
uses checkpoints and may have unexpected behavior on earlier Snakemake versions.

Python 3 and the following Python libraries are required:
* pysam (0.20.0)
* pandas (1.5.2)
* numpy (1.22.4)
* matplotlib (3.6.2)
* BioPython (1.80)
* Excel backend such as openpyxl (3.0.10)

The pipeline was tested with the versions noted above. It should work with earlier versions of these libraries
since it is not likely using features or bug fixes new in these versions.

Additional libraries may need to be present to generate figures (PNG and PDF). To explicitly set a backend, use
set `MPLBACKEND` in the shell environment before launching Snakemake (e.g. `export MPLBACKEND=Agg`).

### Install

Change to the `PIPELINE_DIRECTORY` and pull the pipeline from GitHub:
`git pull https://github.com/EichlerLab/ont_stats.git`

Other than prerequisites, there is nothing else to install.


## Configuration

There are two configuration files that the pipeline will need, the configuration JSON and sample table.

### Configuration JSON
The configuration JSON file, `config.json`, is optional and may be omitted if there are no parameters to set.
Parameters can also be set on the command line after `--config`, but this is only recommended for some parameters
(most should go in `config.json`).

Configuration options:
* sample_table [str]: Name of the sample table to be read. Only needed if the sample table does not use a default
  file name (see "Sample table" below) or is not in the `RUN_DIRECTORY`.
* cell_replace [str]: Each cell is given an ID based on its file name. This option reformats cell names based on a
  regular expression. The value of this option should be in the format "s#FIND#REPLACE#" where "FIND" and "REPLACE" are
  strings to replace in cell names. For example, "s#_guppy-5.0.11-sup-prom_fastq_pass##" will strip extra information
  off each file name (FIND="_guppy-5.0.11-sup-prom_fastq_pass", REPLACE="").
* fig [bool, off]: Generate figures for samples and cells. Do not generate figures if not set.
* force_update [bool, off]: Force the pipeline to re-generate it's cache of samples and cells (see "Cell table tache"
  below). Use if there are changes to FOFN files, but not the input cell table itself. Should be set on the
  command-line (e.g. `snakemake -j 12 --config force_update=on`) and used only when there are known or suspected
  changes, it is not recommended to add this parameter to `config.json`.

For boolean (bool, true/false) parameters, the following string are recognized (others will generate an error):
* True (on): "true", "1", "yes", "t", "y", "on"
* False (off): "false", "0", "no", "f", "n", "off"

Example:
```
{
    "cell_replace": "s#_guppy-5.0.11-sup-prom_fastq_pass##",
    "fig": "true"
}
```

### Sample table

A sample table is read containing a list of sample names and input files. The table has two columns:

1) SAMPLE: Name of the sample.
1) DATA: Path to a FASTQ or an FOFN pointing to one or more FASTQ files.

A sample may have more than one entry, and all data will be processed for each sample.

The optional `sample_table` configuration option points to the sample table file name. If absent, the files
"samples.xlsx", "samples.tsv", "samples.tsv.gz", "samples.csv", "samples.csv.gz" are searched in that order using
the first one found. Any table filename set by the `sample_table` parameter must end in ".xlsx", ".tsv",
".tsv.gz", ".csv", or ".csv.gz" (not case sensitive) or this pipeline does not know how to read it. For Excel
files, the first sheet is read and other sheets (if present) are ignored.

#### DATA field and FOFN files

The `DATA` field points to either a single FASTQ file or an FOFN ("file of file names") containing one or more paths to
input FASTQ files. Paths in an FOFN should be absolute. The pipeline assumes that each FASTQ file is one ONT sequencing
cell, and the cell name is extracted from the FASTQ file name (modifiable by the `cell_replace` configuration
parameter). All paths in one FOFN belong to a single sample, multiple samples should be mixed in an FOFN (create
one entry in the sample table for each sample).

For example, if sample "Sample1" has three ONT cells, "Sample1.fofn" might look like:
```
/path/to/data/ont/Sample1/Cell1.fastq.gz
/path/to/data/ont/Sample1/Cell2.fastq.gz
/path/to/data/ont/Sample1/Cell3.fastq.gz
```

In this example, three cells, "Cell1", "Cell2", and "Cell3" are processed.

## Cell table cache

The *sample table* in the table of input samples and paths to FASTQ files or FOFN files. The pipeline first
generates a *cell table* with one line per cell by processing each entry in the sample table, recursing into
FOFN files, and writing one record for every cell found. The cell table is stored in `data/cell_table.tsv.gz` and
an Excel copy in `data/cell_table.xlsx`. The pipeline only uses the TSV, and the Excel is provided for easy access
to the paths the pipeline ultimately resolves for each cell.

The pipeline maintains MD5 checksums of input sample tables. This allows it to detect changes to the sample table
and force steps to re-run. For example, if more samples are added to a project, the pipeline will detect changes
to the sample table, update the cell table, and process new sample on its next run.

This cache will ***not*** detect changes to FOFN files pointed to by in the sample table since that would not
change the sample table's MD5 checksum. In this case, it's necessary to run the pipeline with the `force_update`
parameter one time after FOFN changes are made (e.g. `snakemake -j 12 --config force_update=on`).

The MD5 checksum cache is stored in `data/cell_table_md5_cache`. The MD5 checksum of the current cell table is
at the top, and if it changed since the first time the pipeline was run in this `RUN_DIRECTORY`, a history of
MD5 checksums will follow (one per line, newest first). Each time the pipeline runs, it checks the top checksum
with the current input sample table and regenerates the cell table if it differs.

## Running the pipeline

### Setup

Before the first run, change to the `RUN_DIRECTORY` and link `Snakefile` from the `PIPELINE_DIRECTORY`.

For example:
```
cd /path/to/run_directory
ln -s /path/to/pipeline_directory/Snakefile
```

The pipeline uses the link back to the pipeline directory to find its libraries.

### Run profiles

The pipeline comes with runtime profiles to set Snakemake options, and optionally, distribute jobs over a Slurm cluster.

Two default profiles are provided, `default` for running jobs on the current machine and `slurm` for distributing
over a Slurm cluster.

To use these profiles without modification, link the "profiles" directory from the `PIPELINE_DIRECTORY` to the
`RUN_DIRECTORY`:

```
cd /path/to/run_directory
ln -s /path/to/pipeline_directory/profiles
```

If you need to make modifications to the profile, copy the profiles directory instead and edit the profile YAML
files. To setup for other clusters, use the `slurm` profile as a guide. Snakemake has some guidance for cluster
configuration (https://snakemake.readthedocs.io/en/stable/executing/cluster.html).

The profile can be selected on the command line. For example, to run the Slurm profile: 
`snakemake --profile profiles/slurm -j 12`

### Run

Before running the pipeline, setup your current environment so the `snakemake` command and Python 3 environment
with all dependencies are accessible.

Then execute Snakemake on slurm:
```
`snakemake --profile profiles/slurm -j 12`
```

To execute Snakemake on the current machine:
```
`snakemake --profile profiles/default -j 12`
```

The pipeline will read the sample table and process all samples. If the `fig` parameter is set, it will also
generate a default set of figures.

Note these examples use 12 concurrent jobs (`-j 12`), adjust for your environment and needs.

## Acknowledgments and citing

This pipeline was developed in Dr. Evan Eichler's lab at the University of Washington and is now maintained
through Dr. Christine Beck's lab at The Jackson Laboratory.

The 2021 HGSVC marker paper can be cited for the pipeline since it's current form was developed for this work:
> Ebert et al., “Haplotype-Resolved Diverse Human Genomes and Integrated Analysis of Structural Variation,”
Science, February 25, 2021, eabf7117, https://doi.org/10.1126/science.abf7117 (PMID: 33632895).

### Contact

Please open a case on the Github page for problems.

You may also contact Peter Audano directly (e-mail omitted to reduce SPAM).
