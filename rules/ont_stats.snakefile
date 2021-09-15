"""
Generate subread stat tables.
"""

#
# Definitions
#

def get_n50(vals):
    """
    Get N50 from a list of lengths.

    :param vals: List of lengths.

    :return: N50.
    """

    vals = vals.sort_values(ascending=False)
    vals_csum = np.cumsum(vals)

    return vals.iloc[np.sum(vals_csum <= (vals_csum.iloc[-1] // 2)) + 1]


#
# Rules
#

# ont_stats_merge_cell
#
# Merge cell stats
rule ont_stats_merge_cell:
    input:
        tsv=lambda wildcards: [
            '{0}/cells/{1}/cell_summary.tsv.gz'.format(SAMPLE_NAME, cell)
            for cell in get_cell_dict(True)
        ]
    output:
        tsv='{0}/cell_summary.tsv.gz'.format(SAMPLE_NAME),
        xlsx='{0}/cell_summary.xlsx'.format(SAMPLE_NAME)
    run:

        # Merge and write
        df = pd.concat(
            [pd.read_csv(tsv_file, sep='\t', header=0) for tsv_file in input.tsv]
        )

        df.to_csv(output.tsv, sep='\t', index=False, compression='gzip')
        df.to_excel(output.xlsx, index=False)

# ont_stats_merge_sample
#
# Calculate stats for the whole sample
rule ont_stats_merge_sample:
    input:
        tsv=lambda wildcards: [
            '{0}/cells/{1}/zmw_summary.tsv.gz'.format(SAMPLE_NAME, cell)
            for cell in get_cell_dict(True)
        ]
    output:
        tsv='{0}/sample_summary.tsv.gz'.format(SAMPLE_NAME),
        xlsx='{0}/sample_summary.xlsx'.format(SAMPLE_NAME)
    run:

        # Read
        df_list = list()

        for file_name in input.tsv:
            df_list.append(
                pd.read_csv(file_name, sep='\t', usecols=('LEN', 'QV'))
            )

        df = pd.concat(df_list, axis=0)

        # Summarize
        df_summary = pd.DataFrame(pd.Series(
            [
                SAMPLE_NAME,
                len(df_list),
                df.shape[0],
                np.sum(df['LEN']),

                np.mean(df['LEN']),
                np.median(df['LEN']),
                get_n50(df['LEN']),

                np.min(df['LEN']),
                np.max(df['LEN']),
                np.std(df['LEN']),

                np.sum(df['LEN'] >= 10000),  # 10 kbp
                np.sum(df.loc[df['LEN'] >= 10000, 'LEN']),

                np.sum(df['LEN'] >= 100000), # 100 kbp
                np.sum(df.loc[df['LEN'] >= 100000, 'LEN']),

                np.sum(df['LEN'] >= 1000000), # 1 Mbp
                np.sum(df.loc[df['LEN'] >= 1000000, 'LEN']),

                np.mean(df['QV']),
                np.median(df['QV']),
                np.max(df['QV'])
            ],
            index=[
                'SAMPLE',
                'N_CELL', 'N', 'SUM',
                'MEAN', 'MED', 'N50',
                'MIN', 'MAX', 'SD',
                'N_10K', 'SUM_10K',
                'N_100K', 'SUM_100K',
                'N_1M', 'SUM_1M',
                'QV_MEAN', 'QV_MED', 'QV_MAX'
            ]
        )).T

        # Write
        df_summary.to_csv(output.tsv, sep='\t', index=False, compression='gzip')
        df_summary.to_excel(output.xlsx, index=False)



# ont_stats
#
# Get stats per subread.
rule ont_stats:
    output:
        tsv_summary=protected('{0}/cells/{{cell}}/cell_summary.tsv.gz'.format(SAMPLE_NAME)),
        tsv_zmw=protected('{0}/cells/{{cell}}/zmw_summary.tsv.gz'.format(SAMPLE_NAME))
    run:

        # Get subread file
        seq_file = get_cell_dict().get(wildcards.cell, None)

        if seq_file is None:
            raise RuntimeError('No sequence data file for cell {}'.format(wildcards.cell))

        # Get stats table
        if seq_file.lower().endswith('.bam'):
            df = ontstatlib.stats.stats_table_bam(seq_file)

        elif seq_file.lower().endswith('.fastq') or seq_file.lower().endswith('.fastq.gz'):
            df = ontstatlib.stats.stats_table_fastq(seq_file)

        else:
            raise RuntimeError(f'Sequence file does not end with ".bam", ".fastq", or ".fastq.gz": {seq_file}')

        # Write ZMW table
        df.to_csv(output.tsv_zmw, sep='\t', index=False, compression='gzip')

        # Summarize by cell
        df_summary = pd.DataFrame(pd.Series(
            [
                wildcards.cell,
                df.shape[0],
                np.sum(df['LEN']),

                np.mean(df['LEN']),
                np.median(df['LEN']),
                get_n50(df['LEN']),

                np.min(df['LEN']),
                np.max(df['LEN']),
                np.std(df['LEN']),

                np.sum(df['LEN'] >= 10000),  # 10 kbp
                np.sum(df.loc[df['LEN'] >= 10000, 'LEN']),

                np.sum(df['LEN'] >= 100000), # 100 kbp
                np.sum(df.loc[df['LEN'] >= 100000, 'LEN']),

                np.sum(df['LEN'] >= 1000000), # 1 Mbp
                np.sum(df.loc[df['LEN'] >= 1000000, 'LEN']),

                np.mean(df['QV']),
                np.median(df['QV']),
                np.max(df['QV'])
            ],
            index=[
                'CELL', 'N', 'SUM',
                'MEAN', 'MED', 'N50',
                'MIN', 'MAX', 'SD',
                'N_10K', 'SUM_10K',
                'N_100K', 'SUM_100K',
                'N_1M', 'SUM_1M',
                'QV_MEAN', 'QV_MED', 'QV_MAX'
            ]
        )).T

        # Write
        df_summary.to_csv(output.tsv_summary, sep='\t', index=False, compression='gzip')
