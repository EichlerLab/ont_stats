"""
Make cumulative read size distribution.
"""

# ont_stats_plot_cdf
#
# Make CDF plot.
rule ont_stats_plot_cdf:
    input:
        tsv=lambda wildcards: [
            '{0}/cells/{1}/zmw_summary.tsv.gz'.format(SAMPLE_NAME, cell)
            for cell in get_cell_dict().keys()
        ]
    output:
        pdf='{sample}/plot/cdf_{per_cell}.pdf',
        png='{sample}/plot/cdf_{per_cell}.png'
    wildcard_constraints:
        per_cell='cell|sample'
    run:

        # Get stat list
        size_list = ontstatlib.plot.get_metric_list(
            wildcards.sample, metric='len', per_cell=wildcards.per_cell == 'cell'
        )

        # Get figure
        fig = ontstatlib.plot.get_cdf_plot(size_list)

        # Write
        fig.savefig(output.pdf, bbox_inches='tight')
        fig.savefig(output.png, bbox_inches='tight')

        # Close
        plt.close(fig)

# ont_stats_plot_line_hist
#
# Make a line histogram
rule ont_stats_plot_line_hist:
    input:
        tsv=lambda wildcards: [
            '{0}/cells/{1}/zmw_summary.tsv.gz'.format(SAMPLE_NAME, cell)
            for cell in get_cell_dict().keys()
        ]
    output:
        pdf='{sample}/plot/hist-line_{per_cell}_{metric}.pdf',
        png='{sample}/plot/hist-line_{per_cell}_{metric}.png'
    wildcard_constraints:
        per_cell='cell|sample',
        metric='qv|passes'
    run:

        # Get list of values
        stat_list = ontstatlib.plot.get_metric_list(
            wildcards.sample, metric=wildcards.metric, per_cell=wildcards.per_cell == 'cell'
        )

        # Get figure
        fig = ontstatlib.plot.get_hist_plot(stat_list, wildcards.metric)

        # Write
        fig.savefig(output.pdf, bbox_inches='tight')
        fig.savefig(output.png, bbox_inches='tight')

        # Close
        plt.close(fig)

# ont_stats_plot_read_len_density
#
# Read length density plot.
rule ont_stats_plot_read_len_density:
    input:
        tsv=lambda wildcards: [
            '{0}/cells/{1}/zmw_summary.tsv.gz'.format(SAMPLE_NAME, cell)
            for cell in get_cell_dict().keys()
        ]
    output:
        pdf='{sample}/plot/density_{per_cell}_len.pdf',
        png='{sample}/plot/density_{per_cell}_len.png'
    wildcard_constraints:
        per_cell='cell|sample',
        metric='len'
    run:

        # Get list of values
        stat_list = ontstatlib.plot.get_metric_list(
            wildcards.sample, metric='len', per_cell=wildcards.per_cell == 'cell'
        )

        fig = ontstatlib.plot.get_density_plot(
            stat_list, z_cut=2.5, xlab='Read length (Kbp)', scale_val=1000
        )

        # Write
        fig.savefig(output.pdf, bbox_inches='tight')
        fig.savefig(output.png, bbox_inches='tight')

        # Close
        plt.close(fig)
