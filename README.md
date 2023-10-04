# ONT read stats

Process ONT reads and generate summary statistics.


## Run modes

The pipeline can be run in two modes:
1) Multi-sample mode reads from a table of samples and input files.
1) Single-sample mode by specifying "fofn" and "sample" (optional) in the config from the command-line.


### Multi-sample mode

A sample table is read containing a list of sample names and input files. The table has two columns:

1) SAMPLE: Name of the sample.
1) DATA: Path to a FASTQ or an FOFN pointing to one or more FASTQ files.

A sample may have more than one entry, and all data will be processed for each sample.
