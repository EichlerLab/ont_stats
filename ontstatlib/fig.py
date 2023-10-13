"""
Compute read-length CDF distributions and create plots.
"""

import collections
import numpy as np
import os
import pandas as pd

import matplotlib as mpl
#mpl.use('Agg')
#If needed: export MPLBACKEND=Agg

import matplotlib.pyplot as plt

from scipy.stats import gaussian_kde
from scipy.stats import zscore



#
# Definitions
#

def z_min_max(vals, z_cut):
    """
    Subset values by z-score cutoff. Takes absolute value of z_scores and removes values greater than z_cut.

    :param vals: Value list.
    :param z_cut: Z-score cutoff.

    :return: Z-score subset values.
    """

    vals_pass = vals[np.abs(zscore(vals)) <= z_cut]

    return np.min(vals_pass), np.max(vals_pass)


#
# Make Plot
#

def get_metric_list(read_table_list, cells=None, sample=None, metric='len', dtype=None):
    """
    Return a list of tuples where the first element is a label and the second
    is an array of read sizes associated with that label.

    If `per_cell` is `True`, then the list has one entry per cell, and the label
    is the cell name. If `per_cell` is `False`, then the list has one entry for
    all reads in the sample, and the label is the sample name.

    :param read_table_list: List of "read_table" (one row per read) for each cell to be read.
    :param cells: List of cells to process. If `None`, get all cells for the sample.
    :param sample: Sample name. If None, concatenate metric into one array and label with the sample name.
    :param metric: Metric to list ("len", "qv", "passes").
    :param dtype: Type of value to expect from the `upper(metric)` column of the read tables. Defaults to `int` if the
        parameter is `None`.
    """

    # Check arguments
    if sample is not None:
        sample = sample.strip()
        if sample == '':
            sample = None

    if sample is None:
        if cells is None:
            raise RuntimeError('Must set "cells" (list of cells for each element in read_table_list) or "sample" (concatenate per sample)')

        if cells.__class__ not in {list, tuple}:
            raise RuntimeError(f'Argument "cells" must be a list or tuple: {cells.__class__}')

        if len(cells) != len(read_table_list):
            raise RuntimeError(f'Length of cells ({len(cells)}) does not match the length of read_table_list ({len(read_table_list)})')

    else:
        if cells is not None:
            raise RuntimeError('Both "sample" and "cells" arguments are set. Set one to return a list labeled by cells or one list labeled wit the sample.')

    # Check file existence (takes a long time to read, fail early)
    missing_list = [file_name for file_name in read_table_list if not os.path.isfile(file_name)]

    if missing_list:
        raise RuntimeError(
            f'Missing {len(missing_list)} cell file(s): {",".join(file_name[:3])}{"..." if len(missing_list) > 3 else ""}'
        )

    # Get size array list
    stat_array_list = get_cell_metric(read_table_list, metric, dtype)

    if sample is None:
        return list(zip(cells, stat_array_list))
    else:
        return [(sample, np.concatenate(stat_array_list))]


def get_cdf_plot(
        stat_list, width=7, height=7, dpi=300, z_cut=None, legend=True,
        grid_spec = {'color': '#b8b8c8', 'which': 'both', 'lw': 0.5},
        xlim=(None, None)
    ):
    """
    Get a CDF plot.

    :param stat_list: A list of 2-element tuples, [0] sample name, [1] a numpy array of sizes. Sizes do not need to
        be sorted.
    :param width: Figure width.
    :param height: Figure height.
    :param dpi: Figure resolution.
    :param z_cut: Cut high values at this z-score cutoff. This prevents large outliers from expanding the x-axis
        and compressing the useful visualizing to the left side of x. `None` disables `z_cut`.
    :param legend: Display a legend with line colors and sample names.
    :param grid_spec: A dictionary of grid keywords (for axes.grid) to specify how the grid should be constructed. Set
        to `None` to disable the grid and show a blank background.
    :param xlim: X-axis limits.
    """

    # Make figure
    fig, ax = plt.subplots(1, 1, figsize=(width, height), dpi=dpi)

    for label, stat in stat_list:

        stat = np.sort(stat) / 1e3  # Sort sizes and convert to kbp
        cdist = np.flip(np.cumsum(np.flip(stat))) / 1e6  # Cumulative distribution in Gbp (sizes already in kbp)

        # Trim distribution (remove skew from ultra-longreads)
        if z_cut is not None:
            max_index = np.sum(((stat - np.mean(stat)) / np.std(stat)) > z_cut)

            stat = stat[:-max_index]
            cdist = cdist[:-max_index]

        # Add to plot
        ax.plot(stat, cdist, '-', label=label)

    # Add labels
    ax.set_ylabel('Cumulative throughput (Gbp)', fontsize=12)
    ax.set_xlabel('Read size (kbp)', fontsize=12)

    ax.get_xaxis().set_major_formatter(
        mpl.ticker.FuncFormatter(lambda x, p: format(int(x), ','))
    )

    ax.get_yaxis().set_major_formatter(
        mpl.ticker.FuncFormatter(lambda x, p: format(int(x), ',.2f'))
    )

    # Set x limits
    if xlim is not None:
        if xlim[0] is not None:
            ax.set_xlim(left=xlim[0])

        if xlim[1] is not None:
            ax.set_xlim(right=xlim[1])

    # Add legend
    if legend:
        ax.legend(prop={'size': 10})

    # Add grid
    if grid_spec is not None:
        ax.grid(**grid_spec)

    # Return plot
    return fig


def get_hist_plot(
        stat_list, metric, width=7, height=7, dpi=300, legend=True,
        grid_spec={'color': '#b8b8c8', 'which': 'both', 'lw': 0.5},
        x_lim=None
    ):
    """
    :param stat_list: A list of 2-element tuples, [0] sample name, [1] a numpy array of sizes. Sizes do not need to
        be sorted.
    :param metric: Metric plotted. May be "qv" or "passes".
    :param width: Figure width.
    :param height: Figure height.
    :param dpi: Figure resolution.
    :param legend: Display a legend with line colors and sample names.
    :param grid_spec: A dictionary of grid keywords (for axes.grid) to specify how the grid should be constructed. Set
        to `None` to disable the grid and show a blank background.
    :param xlim: X-axis limits.

    :return: Tuple of plot figure and axes. Axes is not needed unless the plot is modified.
    """

    if metric == 'qv':
        x_label = 'Estimated QV'

        if x_lim is None:
            x_lim = (20, 61)

    elif metric == 'passes':
        x_label = 'CCS Passes'

    else:
        raise RuntimeError('Unknown metric: {}'.format(metric))

    # Get figure
    fig, ax = plt.subplots(1, 1, figsize=(width, height), dpi=dpi)

    for label, stats in stat_list:
        counter = collections.Counter(stats)
        x_vals = sorted(counter.keys())
        ax.plot(x_vals, [counter[xval] / 1000 for xval in x_vals], '-', label=label)

    ax.set_xlabel(x_label, fontsize=12)
    ax.set_ylabel('Read count (k)', fontsize=12)

    # Aestetics
    ax.get_yaxis().set_major_formatter(
        mpl.ticker.FuncFormatter(lambda x, p: format(int(x), ','))
    )

    ax.grid(grid_spec)

    if x_lim is not None:
        ax.set_xlim(x_lim[0], x_lim[1])

    if legend:
        ax.legend(prop={'size': 10})

    return fig


def get_density_plot(
        stat_list, scale_val=1, width=7, height=7, dpi=300, z_cut=None, legend=True,
        grid_spec = {'color': '#b8b8c8', 'which': 'both', 'lw': 0.5},
        bins=500, covariance_factor=0.25, xlim=(None, None), xlab=None, ylab=None,
        color=dict(), line_style=dict()
    ):
    """
    Get a density plot.

    :param stat_list: A list of 2-element tuples, [0] sample name, [1] a numpy array of sizes. Sizes do not need to
        be sorted.
    :param width: Figure width.
    :param height: Figure height.
    :param dpi: Figure resolution.
    :param z_cut: Cut high values at this z-score cutoff. This prevents large outliers from expanding the x-axis
        and compressing the useful visualizing to the left side of x. `None` disables `z_cut`.
    :param legend: Display a legend with line colors and sample names.
    :param grid_spec: A dictionary of grid keywords (for axes.grid) to specify how the grid should be constructed. Set
        to `None` to disable the grid and show a blank background.
    :params bins: Number of points along x to calculate density.
    :params covariance_factor: Density covariance factor.
    :param xlim: X-axis limits.
    :param xlab: X-axis label.
    :param ylab: Y-axis label.
    :param color: Color dictionary keyed by first tuple element for each entry in stat_list.
    :param line_style: Dictionary of line styles keyed by first tuple element for each entry in stat_list.
    """

    # Make figure
    #fig, ax = plt.figure(1, figsize=(width, height), dpi=dpi)
    # ax = fig.add_subplot(1, 1, 1)

    fig, ax = plt.subplots(1, 1, figsize=(width, height), dpi=dpi)

    try:

        # Set x-axis limits
        if xlim is None:
            xlim = (None, None)

        if z_cut is not None and (xlim[0] is None or xlim[1] is None):

            min_vals, max_vals = list(zip(*[z_min_max(stats, z_cut) for label, stats in stat_list]))

            if xlim[0] is None:
                x_min = np.min(min_vals)
            else:
                x_min = x_lim[0]

            if xlim[1] is None:
                x_max = np.max(max_vals)
            else:
                x_max = x_lim[1]

        else:
            x_min = np.min([np.min(stats) for label, stats in stat_list]) if xlim[0] is not None else xlim[0]
            x_max = np.max([np.max(stats) for label, stats in stat_list]) if xlim[1] is not None else xlim[1]

        # Get x-axis limits
        x = np.linspace(x_min, x_max, bins)

        # Get densities
        for label, stats in stat_list:

            stats = np.sort(stats)  # Sort sizes and convert to kbp

            # Get density
            density = gaussian_kde(stats)

            density.covariance_factor = lambda : covariance_factor
            density._compute_covariance()

            ax.plot(
                x, density(x) * np.sum(stats),
                line_style.get(label, '-'),
                color=color.get(label, None),
                label=label
            )

        # Add labels
        if xlab is not None:
            ax.set_xlabel(xlab, fontsize=12)

        if ylab is not None:
            ax.set_ylabel(ylab, fontsize=12)

        # Aestetics

        ax.get_yaxis().set_visible(False)

        ax.get_xaxis().set_major_formatter(
            mpl.ticker.FuncFormatter(lambda x, p: format(int(x / scale_val), ','))
        )

        # Add legend
        if legend:
            ax.legend(prop={'size': 10})

        # Add grid
        if grid_spec is not None:
            ax.grid(**grid_spec)

    except Exception as ex:
        plt.close(fig)
        raise ex

    # Return plot
    return fig


def get_cell_metric(read_table_list, metric='len', dtype=None):
    """
    Get a list of Numpy arrays where each array is the read sizes from one cell.

    :param read_table_list: List of "read_table" (one row per read) for each cell to be read.
    :param metric: Metric to list ("len", "qv", "passes").
    :param dtype: Type of value to expect from the `upper(metric)` column of the read tables. Defaults to `int` if the
        parameter is `None`.
    """

    if dtype is None:
        dtype = int

    # Check metric
    if metric not in {'len', 'qv', 'passes'}:
        raise RuntimeError('Unknown metric to retrieve: {}: Expected "len", "qv", or "passes"'.format(metric))

    # Check file existence (takes a long time to read, fail early)
    missing_list = [file_name for file_name in read_table_list if not os.path.isfile(file_name)]

    if missing_list:
        raise RuntimeError(
            f'Missing {len(missing_list)} cell file(s): {",".join(file_name[:3])}{"..." if len(missing_list) > 3 else ""}'
        )

    # Read
    metric_array_list = list()

    for file_name in read_table_list:
        df = pd.read_csv(file_name, sep='\t', header=0)

        # Get array
        metric_array_list.append(np.asarray(df[metric.upper()], dtype))

    # Return list
    return metric_array_list
