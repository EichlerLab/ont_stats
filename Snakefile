"""
Get stats for raw sequencing data in each sample.
"""

import collections
import os
import pysam
import re

import pandas as pd
import numpy as np

import ontstatlib.plot
import ontstatlib.stats

import matplotlib.pyplot as plt

localrules: ont_stats_get_tables

# Read config
if os.path.isfile('config.json'):
    configfile: config.json


#
# Get cell table
#

# Cached cell table. Used to avoid re-reading the table.
df_cell_cache = None

def get_cell_table(sample=None):
    """
    Get the cell table.

    If the cell table was already read, the existing cell table is returned. If the cell table was not read, it is read,
    processed, and returned.

    :param sample: Return the sample table subset for this sample. If the cell table cache is read, the full table will
        still be read and cached, this only affects the subset returned by this function.

    :return: Cell table.
    """

    # Return sample table if it was already read
    global df_cell_cache

    # Get table name
    if df_cell_cache is None:
        sample_table_name = config.get('sample_table', None)

        if sample_table_name is None:
            for table_name in ('samples.tsv', 'samples.tsv.gz', 'samples.xlsx'):
                if os.path.isfile(table_name):
                    sample_table_name = table_name
                    break

        if sample_table_name is None:
            raise RuntimeError('Sample table not found in "sample_table" config option or default table names')

        # Read
        if sample_table_name.lower().endswith('.xlsx'):
            df_sample = pd.read_excel(sample_table_name)
        else:
            df_sample = pd.read_csv(sample_table_name, sep='\t')

        missing_cols = [col for col in ['SAMPLE', 'DATA'] if col not in df_sample.columns]

        if missing_cols:
            raise RuntimeError(f'Sample table is missing columns: {",".join(missing_cols)}')

        if np.any(pd.isnull(df['SAMPLE'])) or np.any(df['SAMPLE'].apply(lambda val: val.strip() == '')):
            raise RuntimeError('Blank sample names in sample table')

        if np.any(pd.isnull(df['DATA'])) or np.any(df['DATA'].apply(lambda val: val.strip() == '')):
            raise RuntimeError('Blank data entries in sample table')

        df_sample = df_sample[['SAMPLE', 'DATA']]

        # Cell re-write parameters (alters the cell name)
        if 'cell_replace' in config:
            cell_replace = re.search('s#([^#]+)#([^#]*)#', config['cell_replace'])

            if cell_replace is not None:
                cell_replace = (cell_replace[1], cell_replace[2])
            else:
                raise RuntimeError(f'Parameter "cell_replace" should be in the format "s#FIND#REPLACE#" where "FIND" and "REPLACE" are strings to replace in cell names: {config["cell_replace"]}')

        else:
            cell_replace=None

        # Find all cells
        df_cell_list = list()

        for index, row in df_sample.iterrows():
            sample = row['SAMPLE']

            file_list_in = [row['DATA']]

            while file_list_in:
                file_name = file_list_in[0]
                file_list_in = file_list_in[1:]

                if file_name.lower().endswith('.fofn'):
                    # FOFN file

                    new_file_list = list()

                    with open(file_name, 'rt') as in_file:
                        for line in in_file:
                            line = line.strip()

                            if not line or line.startswith('#'):
                                continue

                            new_file_list.append(line)

                    file_list_in = new_file_list + file_list_in

                else:
                    # File name entry

                    # Check format
                    ext_match = re.match(r'.+(\.fastq|\.fastq\.gz)$', file_name, re.IGNORECASE)

                    if ext_match is None:
                        raise RuntimeError(f'Unrecognized file type in input: Must end with ".fastq" or ".fastq.gz": {file_name}')

                    # Get cell name
                    cell_name = os.path.basename(line)[:-len(ext_match[1])]

                    if cell_replace is not None:
                        cell_name = re.sub(cell_replace[0], cell_replace[1], cell_name).strip()

                        if not cell_name:
                            raise RuntimeError(f'Cell name empty after config[cell_replace] was applied: {file_name}')

                    # Add entry
                    df_cell_list.append(
                        pd.Series(
                            [sample, cell_name, file_name],
                            index=['SAMPLE', 'CELL', 'DATA']
                        )
                    )

        # Merge records
        df_cell_cache = pd.concat(df_cell_list, axis=1).T

    # Subset
    df_cell = df_cell_cache.copy()

    if sample is not None:
        df_cell = df_cell.loc[df_cell['SAMPLE'] == sample]

    return df_cell


#
# Definitions
#

def get_all_input(wildcards):

    df_cell = get_cell_table()

    file_list = list()

    for sample in set(df_cell['SAMPLE']):
        file_list += [
            f'{sample}/cell_summary.tsv.gz',
            f'{sample}/plot/cdf_cell.png',
            f'{sample}/plot/hist-line_cell_qv.png',
            f'{sample}/plot/density_cell_len.png'
        ]

    for index, row in df_cell:

        file_list += [
            f'{row["SAMPLE"]}/cells/{row["CELL"]}/cell_summary.tsv.gz',
            f'{row["SAMPLE"]}/cells/{row["CELL"]}/zmw_summary.tsv.gz'
        ]

    return file_list

def get_file_per_cell(sample, file_pattern):

    df_cell = get_cell_table()

    df_cell = df_cell.loc[df_cell['SAMPLE'] == sample]

    return [
        file_pattern.format(sample=sample, cell=cell)
            for cell in df_cell['CELL']
    ]

#
# Rules
#

localrules: ont_stats_all

# Get all tables.
rule ont_stats_all:
    input:
        files=get_all_input

rule ont_stats_cell_table:
    output:
        tsv='cell_table.tsv.gz'
    run:
        get_cell_table().to_csv(output.tsv, sep='\t', index=False, compression='gzip')

include: 'rules/ont_stats.snakefile'
include: 'rules/plot.snakefile'
