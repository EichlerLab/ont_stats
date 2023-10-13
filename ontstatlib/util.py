# Utility functions

import collections
import hashlib
import numpy as np
import os
import pandas as pd
import pathlib
import re

def get_cell_table_file_name(config):
    sample_table_name = config.get('sample_table', None)

    if sample_table_name is None:
        for table_name in ('samples.xlsx', 'samples.tsv', 'samples.tsv.gz', 'samples.csv', 'samples.csv.gz'):
            if os.path.isfile(table_name):
                sample_table_name = table_name
                break

    if sample_table_name is None:
        raise RuntimeError('Sample table not found in "sample_table" config option or default table names')

    return sample_table_name


def check_and_update_cell_table_md5_cache(config, cache_file_name):
    """
    Check the cell table for changes and keep a record of MD5 checksums for all versions of the cell table processed.

    :param config: Pipeline config. Reads cell table location and the "force_update" parameter from the config.
    :param cache_file_name: Name of the cache file to check and update.
    """

    md5_pattern = re.compile(r'^[0-9a-f]{32}$')

    # Get the force_update parameter
    try:
        force_update = as_bool(config.get('force_update', False))
    except ValueError as e:
        raise RuntimeError(f'Error translating "force_update" configuration parameter as a boolean value: {config["force_update"]}')

    # Read the MD5 of the cell table the last pipeline run used.
    if os.path.isfile(cache_file_name):
        with open(cache_file_name, 'rt') as in_file:

            # Read existing cache and get top md5
            existing_cache = in_file.read().strip()
            existing_md5 = existing_cache.split('\n', 1)[0]

            # Check the MD5 pattern to make sure it looks like a real MD5
            if re.search(md5_pattern, existing_md5) is None:
                # Not the expected MD5 pattern, something else written to this file?
                existing_md5 = None

    else:
        existing_md5 = None
        existing_cache = None

    # MD5 of the current input table
    cell_table_file_name = get_cell_table_file_name(config)

    with open(cell_table_file_name, 'rt') as in_file:
        input_md5 = hashlib.md5(in_file.read().encode()).hexdigest()

    if re.search(md5_pattern, input_md5) is None:
        raise RuntimeError(f'BUG: Computed MD5 for input table {cell_table_file_name}, but MD5 string does not match the expected pattern of 32 characters: {input_md5}')

    # Update the MD5 cache
    if existing_md5 is None or existing_md5 != input_md5:
        parent_directory = os.path.dirname(cache_file_name).strip()

        if parent_directory not in {'.', '..', ''}:
            os.makedirs(parent_directory, exist_ok=True)

        with open(cache_file_name, 'wt') as out_file:

            # Write top MD5
            out_file.write(input_md5)
            out_file.write('\n')

            # If previous MD5s were found (changes to the input table), record them below the current MD5.
            if existing_cache is not None:
                out_file.write(existing_cache)
                out_file.write('\n')

    elif force_update:
        pathlib.Path(cache_file_name).touch()


def make_cell_table(config):
    """
    Read the sample table and generate the cell table from it.

    :param config: Pipeline configuration parameters.

    :return: Cell table (one entry per sample/cell) for all samples.
    """

    # Get sample table name
    sample_table_name = get_cell_table_file_name(config)
    print(f'Reading sample table: {sample_table_name}')

    # Read sample table
    if sample_table_name.lower().endswith('.xlsx'):
        df_sample = pd.read_excel(sample_table_name)
    elif sample_table_name.lower().endswith('.tsv') or sample_table_name.lower().endswith('.tsv.gz'):
        df_sample = pd.read_csv(sample_table_name, sep='\t')
    elif sample_table_name.lower().endswith('.csv') or sample_table_name.lower().endswith('.csv.gz'):
        df_sample = pd.read_csv(sample_table_name)
    else:
        raise RuntimeError(f'Sample table {sample_table_name} does not have a recognized file extension (.xlsx, .tsv, .tsv.gz, .csv, .csv.gz)')

    # Check table columns and data integrity
    missing_cols = [col for col in ['SAMPLE', 'DATA'] if col not in df_sample.columns]

    if missing_cols:
        raise RuntimeError(f'Sample table is missing columns: {",".join(missing_cols)}')

    if np.any(pd.isnull(df_sample['SAMPLE'])) or np.any(df_sample['SAMPLE'].apply(lambda val: val.strip() == '')):
        raise RuntimeError('Blank sample names in sample table')

    if np.any(pd.isnull(df_sample['DATA'])) or np.any(df_sample['DATA'].apply(lambda val: val.strip() == '')):
        raise RuntimeError('Blank data entries in sample table')

    dup_sample = [sample for sample, count in collections.Counter(df_sample['SAMPLE']).items() if count > 1]

    if dup_sample:
        n_dup = len(dup_sample)
        dup_sample = dup_sample[:3]

        raise RuntimeError(f'Sample table {sample_table_name} has {n_dup} duplicate sample names: {", ".join(dup_sample)}{"..." if n_dup > 3 else ""}')

    df_sample = df_sample[['SAMPLE', 'DATA']]

    # Cell re-write parameters (alters the cell name)
    if 'cell_replace' in config:
        cell_replace = re.search('s#([^#]+)#([^#]*)#', config['cell_replace'])

        if cell_replace is not None:
            cell_replace = (re.compile(cell_replace[1]), cell_replace[2])
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
                cell_name = os.path.basename(file_name)[:-len(ext_match[1])]

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

    df_cell = pd.concat(df_cell_list, axis=1).T

    # Warn about duplicate cell names for each sample
    dup_cell_list = [
        (sample, cell) for (sample, cell), count in collections.Counter(
            df_cell[['SAMPLE', 'CELL']].apply(tuple, axis=1)
        ).items() if count > 1
    ]

    if dup_cell_list:
        n_cell = len(dup_cell_list)
        dup_cell_list = ', '.join([f'{sample}:{cell}' for sample, cell in dup_cell_list[:3]])

        raise RuntimeError(f'Found {n_cell} duplicate sample/cell entries (each sample must not have multiple of the same cell name): {dup_cell_list}{"..." if n_cell > 3 else ""}')

    # Merge records
    return df_cell


def as_bool(val):
    """
    Translate value as a boolean. If `val` is boolean, return `val`. If val is a string, `True` if lower-case string
    is "true", "1", "yes", "t", or "y". `False` if lower-case string is "false", "0", "no", "f", or "n". All other
    values throw a RuntimeError. If `val` is not bool or string, it is converted to a string and the rules above are
    applied.

    :param val: Value to interpret as a boolean.

    :return: `True` or `False` (see above).
    """

    if issubclass(val.__class__, bool):
        return val

    val = str(val).lower()

    if val in {'true', '1', 'yes', 't', 'y', 'on'}:
        return True

    if val in {'false', '0', 'no', 'f', 'n', 'off'}:
        return False

    raise ValueError('Cannot interpret as boolean value: {}'.format(val))

def get_n50(vals):
    """
    Get N50 from a list of lengths.

    :param vals: List of lengths.

    :return: N50.
    """

    vals = vals.sort_values(ascending=False)
    vals_csum = np.cumsum(vals)

    return vals.iloc[np.sum(vals_csum <= (vals_csum.iloc[-1] // 2)) + 1]
