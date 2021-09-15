"""
Compute read-length CDF distributions and create plots.
"""

import collections
import numpy as np
import os
import pandas as pd

import matplotlib as mpl
mpl.use('Agg')

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

def get_metric_list(sample, cells=None, per_cell=False, metric='len', rel_path='.'):
    """
    Return a list of tuples where the first element is a label and the second
    is an array of read sizes associated with that label.

    If `per_cell` is `True`, then the list has one entry per cell, and the label
    is the cell name. If `per_cell` is `False`, then the list has one entry for
    all reads in the sample, and the label is the sample name.

    :param sample: Sample name.
    :param cells: List of cells to process. If `None`, get all cells for the sample.
    :param per_cell: Plot one distribution per cell.
    :param metric: Metric to list ("len", "qv", "passes").
    :param rel_path: Search for samples from this relative path.
    """

    # Get cell list
    if cells is None:
        cells = os.listdir(os.path.join(rel_path, sample, 'cells'))

    # Check all cell files before reading (takes a long time to read, fail early)
    cells = [cell.strip() if cell is not None else '' for cell in cells]

    for cell in cells:
        if not cell:
            raise RuntimeError('Empty cell name in cell list')

        file_name = os.path.join(rel_path, sample, 'cells', cell, 'zmw_summary.tsv.gz')

        if not os.path.isfile(file_name):
            raise RuntimeError('Missing cell file: {}'.format(file_name))

    # Get size array list
    stat_array_list = get_cell_metric(sample, cells, metric, rel_path)

    if per_cell:
        stat_list = list(zip(cells, stat_array_list))
    else:
        stat_list = [(sample, np.concatenate(stat_array_list))]

    # Return CDF list
    return stat_list


def get_cdf_plot(
        size_list, width=7, height=7, dpi=300, z_cut=None, legend=True,
        grid_spec = {'color': '#b8b8c8', 'which': 'both', 'lw': 0.5},
        xlim=(None, None)
    ):
    """
    Get a CDF plot.

    :param size_list: A list of 2-element tuples, [0] sample name, [1] a numpy array of sizes. Sizes do not need to
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

    for label, sizes in size_list:

        sizes = np.sort(sizes) / 1e3  # Sort sizes and convert to kbp
        cdist = np.flip(np.cumsum(np.flip(sizes))) / 1e6  # Cumulative distribution in Gbp (sizes already in kbp)

        # Trim distribution (remove skew from ultra-longreads)
        if z_cut is not None:
            max_index = np.sum(((sizes - np.mean(sizes)) / np.std(sizes)) > z_cut)

            sizes = sizes[:-max_index]
            cdist = cdist[:-max_index]

        # Add to plot
        ax.plot(sizes, cdist, '-', label=label)

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


def get_cell_metric(sample, cells, metric='len', rel_path='.'):
    """
    Get a list of Numpy arrays where each array is the read sizes from one cell.

    :param sample: Sample name.
    :param cells: List of cells to process. If `None`, get all cells for the sample.
    :param metric: Metric to list ("len", "qv", "passes").
    :param rel_path: Search for samples from this relative path.
    """

    # Check sample
    if sample is None:
        raise RuntimeError('Sample is None')

    sample = sample.strip()

    if not sample:
        raise RuntimeError('Sample name is empty')

    if metric not in {'len', 'qv', 'passes'}:
        raise RuntimeError('Unknown metric to retrieve: {}: Expected "len", "qv", or "passes"'.format(metric))

    # Get cells
    cells = os.listdir(os.path.join(rel_path, sample, 'cells'))

    # Read
    metric_array_list = list()

    for cell in cells:
        file_name = os.path.join(rel_path, sample, 'cells', cell, 'zmw_summary.tsv.gz')

        if not os.path.isfile(file_name):
            raise RuntimeError('Missing cell file in get_cell_sizes(): {}'.format(file_name))

        df = pd.read_csv(file_name, sep='\t', header=0)

        # Get array
        metric_array_list.append(np.asarray(df[metric.upper()], np.int32))

    # Return list
    return metric_array_list
