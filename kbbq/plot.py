"""
Utilities for plotting calibration and other relevant info.
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
import sklearn
import kbbq.compare_reads as cr

def plot_calibration(data, truth, labels, plotname, plottitle = None, maxscore = 42):
    """
    Plot true probability of each quality bin.

    This plots the predicted quality score (x) against the true quality score (y)
    for each quality bin. This is done by calculating the mean probability of error
    for each bin and converting that probability back to a quality score.

    :param list(np.array) data: list of 1d arrays of quality scores
    :param np.array truth: boolean array where errors are True.
    :param list(str) labels: list of labels to add to plot corresponding
     to each data array
    :param str plotname: name of file to save plot to
    :param str plottitle: title of plot. If not provided, plotname is used instead.
    :param int maxscore: max score to plot.

    """
    print(tstamp(), "Making Quality Score Plot . . .", file = sys.stderr)
    assert np.ndim(truth) == 1
    plottitle = (plottitle if plottitle is not None else plotname)
    sns.set()
    qualplot, ax = plt.subplots(figsize = (9,9))
    plt.title(plottitle)
    ax.plot(np.arange(maxscore + 1), 'k:', label = "Perfect")
    assert np.all([np.array_equal(data[0].shape, data[i].shape) for i in range(len(data))])
    for d, label in zip(data, labels):
        print(tstamp(), "Plotting %s . . ." % (label), file = sys.stderr)
        assert np.ndim(d) == 1
        est_p = cr.q_to_p(d)
        bscore = sklearn.metrics.brier_score_loss(truth.reshape(-1), est_p.reshape(-1))
        numtotal = np.bincount(d.reshape(-1), minlength = (maxscore+1))
        numerrs = np.bincount(d[truth].reshape(-1), minlength = len(numtotal)) #errors
        p = np.true_divide(numerrs[numtotal != 0],numtotal[numtotal != 0])
        q = cr.p_to_q(p)
        ax.plot(np.arange(len(numtotal))[numtotal != 0], q, 'o-', alpha = .8, label ="%s, %1.5f" % (label, bscore))
    plt.xlabel("Predicted Quality Score")
    plt.ylabel("Actual Quality Score")
    ax.legend(loc = "upper left")
    qualplot.savefig(plotname, bbox_inches='tight')

def plot_samplesize(data, labels, plotname, plottitle = None, maxscore = 42):
    """
    Plot the size of each quality bin.

    This plots each quality bin (x) against the number
    of bases assigned to that quality bin (y) as a 
    line plot.

    :param list(np.array) data: list of 1d arrays of quality scores
    :param list(str) labels: list of labels to add to plot corresponding
     to each data array
    :param str plotname: name of file to save plot to
    :param str plottitle: title of plot. If not provided, plotname is used instead.
    :param int maxscore: max score to plot.

    """
    print(tstamp(), "Making Sample Size Plot . . .", file = sys.stderr)
    plottitle = (plottitle if plottitle is not None else plotname)
    sns.set()
    ssplot, ax = plt.subplots(figsize = (9,9))
    ssplot.suptitle(plottitle)
    for d, label in zip(data, labels):
        numtotal = np.bincount(d.reshape(-1), minlength = (maxscore+1))
        ax.plot(np.arange(len(numtotal))[numtotal != 0], numtotal[numtotal != 0], 'o-', alpha = .6, label = f"{label}")
    plt.xlabel("Predicted Quality Score")
    plt.ylabel("Sample Size")
    ax.legend(loc = "upper left")
    ssplot.savefig(plotname)

def plot_benchmark(fhin, outfile, plottype = 'calibration'):
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
