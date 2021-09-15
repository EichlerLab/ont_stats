"""
Get raw stats and tables.
"""

import Bio.SeqIO

import gzip
import numpy as np
import pandas as pd

MAX_QV = 60


def stats_table_fastq(seq_file):
    """
    Get a table of stats from a CCS FASTQ file.

    :param seq_file: FASTQ file name.

    :return: A dataframe of read stats.
    """

    # Read each line in BAM
    record_list = list()

    with PlainOrGzReader(seq_file, 'rt') as in_file:
        for record in Bio.SeqIO.parse(in_file, 'fastq'):

            # Add record
            record_list.append(
                pd.Series(
                    [
                        record.id,
                        len(record.seq),
                        calc_exp_acc(record.letter_annotations['phred_quality'])
                    ],
                    index=['ID', 'LEN', 'QV']
                )
            )

    # Merge
    df = pd.concat(record_list, axis=1).T

    df['QV'] = df['QV'].apply(lambda qv: MAX_QV if (qv == np.inf or qv > MAX_QV) else qv)

    return df


class PlainOrGzReader:
    """
    Read a plain or a gzipped file using context guard.
    Code taken from SV-Pop (Audano/Eichler lab).

    Example:
        with PlainOrGzReader('path/to/file.gz'): ...
    """

    def __init__(self, file_name, mode='rt'):

        self.file_name = file_name

        self.is_gz = file_name.lower().endswith('.gz')

        self.mode = mode

        self.file_handle = None

    def __enter__(self):

        if self.is_gz:
            self.file_handle = gzip.open(self.file_name, self.mode)
        else:
            self.file_handle = open(self.file_name, self.mode)

        return self.file_handle

    def __exit__(self, exc_type, exc_value, traceback):

        if self.file_handle is not None:
            self.file_handle.__exit__()
            self.file_handle = None


def calc_exp_acc(qv_list, qv_trim_5=0, qv_trim_3=0):
    """
    Estimate read QV from per-base QV values. Thanks to Katy Munson for the template this code was developed from.

    :param qv_list: List of QV values.
    :param qv_trim_5: 5' trimming (remove this many bases from the left end)
    :param qv_trim_3: 3' trimming (remove this many bases from the right end)
    """

    assert qv_trim_5 >= 0 and qv_trim_3 >= 0

    # Convert PHRED QV score to error probability
    base_err = 10 ** -(
            np.array(
                qv_list[qv_trim_5:-qv_trim_3] if qv_trim_3 > 0 else qv_list[qv_trim_5:]
            ) / 10.0
    )

    # Average error probability and convert back to PHRED scale
    mean_err = np.sum(base_err) / base_err.shape[0]

    return -10 * np.log10(mean_err) if mean_err > 0.0 else np.inf
