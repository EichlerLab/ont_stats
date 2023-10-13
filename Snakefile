"""
Get stats for raw sequencing data in each sample.
"""

import collections
import os
import pysam
import re

import pandas as pd
import numpy as np

import sys

import matplotlib.pyplot as plt

PIPELINE_DIR = os.path.dirname(os.path.realpath(workflow.snakefile))

sys.path.append(PIPELINE_DIR)

import ontstatlib

# Read config
if os.path.isfile('config.json'):
    configfile: 'config.json'

# Constants
DATA_CELL_FILE_NAME = 'data/cell_table.tsv.gz'
DATA_CELL_FILE_NAME_XLSX = 'data/cell_table.xlsx'

INPUT_TABLE_MD5_CACHE = 'data/cell_table_md5_cache'

# Check the input table and generate a cache record for each cell table change (most recent at the top). Forces the
# pipeline to re-run the ont_stats_cell_table checkpoint for each input table change.
ontstatlib.util.check_and_update_cell_table_md5_cache(config, INPUT_TABLE_MD5_CACHE)

#
# Rules
#

# Get input rules
include: os.path.join(PIPELINE_DIR, 'rules/input.sm')


# Get all tables.
localrules: ont_stats_all

rule ont_stats_all:
    input:
        files=get_all_input

# Write cell table
checkpoint ont_stats_cell_table:
    input:
        md5_cache=INPUT_TABLE_MD5_CACHE
    output:
        tsv=DATA_CELL_FILE_NAME,
        xlsx=DATA_CELL_FILE_NAME_XLSX
    run:

        cell_table = ontstatlib.util.make_cell_table(config)

        cell_table.to_csv(output.tsv, sep='\t', index=False, compression='gzip')

        cell_table.to_excel(output.xlsx, index=False)


include: os.path.join(PIPELINE_DIR, 'rules/ont_stats.sm')
include: os.path.join(PIPELINE_DIR, 'rules/fig.sm')
include: os.path.join(PIPELINE_DIR, 'rules/merged.sm')
