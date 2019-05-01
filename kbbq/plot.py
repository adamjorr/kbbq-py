#!/usr/bin/env python3
import numpy as np
import matplotlib
matplotlib.use('Agg')
import seaborn as sns

def plot_calibration(data, truth, labels, plotname, plottitle = None, maxscore = 42):
    print(tstamp(), "Making Quality Score Plot . . .", file = sys.stderr)
    assert np.ndim(truth) == 1
    plottitle = (plottitle if plottitle is not None else plotname)
    sns.set()
    qualplot, ax = plt.subplots(figsize = (9,9))
    plt.title(plottitle)
    ax.plot(np.arange(maxscore + 1), 'k:', label = "Perfect")
    assert np.all([np.array_equal(data[0].shape, data[i].shape) for i in range(len(data))])
    for i in range(len(data)):
        label = labels[i]
        print(tstamp(), "Plotting %s . . ." % (label), file = sys.stderr)
        assert np.ndim(data[i]) == 1
        est_p = q_to_p(data[i])
        bscore = sklearn.metrics.brier_score_loss(truth.reshape(-1), est_p.reshape(-1))
        numtotal = np.bincount(data[i].reshape(-1), minlength = (maxscore+1))
        numerrs = np.bincount(data[i][truth].reshape(-1), minlength = len(numtotal)) #errors
        p = np.true_divide(numerrs[numtotal != 0],numtotal[numtotal != 0])
        q = p_to_q(p)
        ax.plot(np.arange(len(numtotal))[numtotal != 0], q, 'o-', alpha = .8, label ="%s, %1.5f" % (label, bscore))
    plt.xlabel("Predicted Quality Score")
    plt.ylabel("Actual Quality Score")
    ax.legend(loc = "upper left")
    qualplot.savefig(plotname, bbox_inches='tight')

def plot_samplesize(data, labels, plotname, plottitle = None, maxscore = 42):
    print(tstamp(), "Making Sample Size Plot . . .", file = sys.stderr)
    plottitle = (plottitle if plottitle is not None else plotname)
    sns.set()
    ssplot, ax = plt.subplots(figsize = (9,9))
    ssplot.suptitle(plottitle)
    for i in range(len(data)):
        numtotal = np.bincount(data[i].reshape(-1), minlength = (maxscore+1))
        ax.plot(np.arange(len(numtotal))[numtotal != 0], numtotal[numtotal != 0], 'o-', alpha = .6, label = f"{labels[i]}")
    plt.xlabel("Predicted Quality Score")
    plt.ylabel("Sample Size")
    ax.legend(loc = "upper left")
    ssplot.savefig(plotname)

