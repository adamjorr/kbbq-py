"""
Utilities for plotting calibration and other relevant info.
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns

def plot_benchmark(fhin, outfile, plottype = 'calibration'):
    """
    Plot predicted vs. true probability or size of each quality bin.

    This plots the predicted quality score (x) against the true quality score (y)
    for each quality bin. This is done by calculating the mean probability of error
    for each bin and converting that probability back to a quality score.

    With plottype = 'sample-size', instead plots the size of each predicted quality
    score bin.

    :param python.IOBase fhin: file containing benchmark data
    :param str outfile: name of file to save plot to
    """
    if plottype == 'calibration':
        cols = [0,1,2]
        ylab = 'Actual Quality Score'
        title = 'Substitution Error Calibration'
    elif plottype == 'sample-size':
        cols = [0,3,2]
        ylab = 'Number of Bases'
        title = 'Sample Size'
    else:
        raise ValueError("Invalid plot type specified. Type should be calibration or sample-size")
    xlab = "Predicted Quality Score"
    sns.set()
    x, y, labels = zip(*[[line.rstrip().split()[c] for c in cols] for line in fhin])
    x = np.array([int(i) for i in x], dtype = np.int)
    y = np.array([int(i) for i in y], dtype = np.int)
    labels = np.array([str(i) for i in labels], dtype = np.unicode)
    plot, ax = plt.subplots(figsize = (9,9))
    plt.title(title)
    if plottype == 'calibration':
        ax.plot(np.arange(max(x) + 1), 'k:', label = "Perfect")
    for l in np.unique(labels):
        ax.plot(x[labels == l], y[labels == l], 'o-', alpha = .8, label = l)
    plt.xlabel(xlab)
    plt.ylabel(ylab)
    ax.legend(loc = 'upper left')
    plot.savefig(outfile, bbox_inches = 'tight')
