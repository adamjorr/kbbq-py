"""
Utilities for emulating GATK's ApplyBQSR tool

BQSR model construction hard clips soft clips and
trims adaptors. ApplyBQSR does not. So we need
different functions for each.
"""

import numpy as np
import pandas as pd
from .. import compare_reads as utils
from .. import recaltable

def table_to_vectors(table, rg_order, maxscore = 42):
    #the recal table uses the PU of the read group as the read group entry in the table
    #see vectors_to_report for more info
    dinuc_order = utils.Dinucleotide.dinuc_to_int.keys()
    rgtable = table.tables[2].data.reindex(rg_order)
    meanq = rgtable['EstimatedQReported'].values.astype(np.float64)
    global_errs = rgtable['Errors'].values.astype(np.int64)
    global_total = rgtable['Observations'].values

    qtable = table.tables[3].data.reindex(pd.MultiIndex.from_product([rg_order, np.arange(maxscore + 1)]))
    q_shape = (len(rg_order), maxscore + 1)
    q_errs = qtable['Errors'].fillna(0, downcast = 'infer').values.reshape(q_shape)
    q_total = qtable['Observations'].fillna(0, downcast = 'infer').values.reshape(q_shape)

    postable = table.tables[4].data.loc[rg_order, np.arange(maxscore + 1), 'Cycle']
    postable = postable.reset_index(level = 'CovariateValue').astype({'CovariateValue' : np.int_}).set_index('CovariateValue', append = True)
    seqlen = postable.index.get_level_values('CovariateValue').max()
    postable = postable.reindex(pd.MultiIndex.from_product(
        [rg_order, np.arange(maxscore + 1), ['Cycle'],
        np.concatenate([np.arange(seqlen)+1, np.flip(-(np.arange(seqlen)+1),axis = 0)])
        ]))
    pos_shape = (len(rg_order), maxscore + 1, 2 * seqlen)
    pos_errs = postable['Errors'].fillna(0, downcast = 'infer').values.reshape(pos_shape)
    pos_total = postable['Observations'].fillna(0, downcast = 'infer').values.reshape(pos_shape)

    dinuctable = table.tables[4].data.reindex(pd.MultiIndex.from_product([rg_order, np.arange(maxscore + 1), ['Context'], dinuc_order]))
    dinuc_shape = (len(rg_order), maxscore + 1, len(dinuc_order))
    dinuc_errs = dinuctable['Errors'].fillna(0, downcast = 'infer').values.reshape(dinuc_shape)
    dinuc_total = dinuctable['Observations'].fillna(0, downcast = 'infer').values.reshape(dinuc_shape)

    return meanq, global_errs, global_total, q_errs, q_total, pos_errs, pos_total, dinuc_errs, dinuc_total

def bamread_cycle_covariates(read):
    cycle = utils.generic_cycle_covariate(read.query_length, read.is_read2)
    if read.is_reverse:
        cycle = np.flip(cycle)
    return cycle

def bamread_dinuc_covariates(read, use_oq = True, minscore = 6):
    seq = read.query_sequence
    quals = (utils.bamread_get_oq(read) if use_oq else seq.query_qualities)
    if read.is_reverse:
        seq = ''.join([utils.Dinucleotide.complement.get(x,'N') for x in reversed(seq)])
        quals = np.flip(quals)
    fulldinuc = np.zeros(read.query_length, dtype = np.int)
    dinuccov = utils.generic_dinuc_covariate(np.array(list(seq), dtype = 'U1'), quals, minscore)
    if read.is_reverse:
        # flip back to fwd coordinates
        dinuccov = np.flip(dinuccov)
    return dinuccov

def recalibrate_bamread(read, meanq, globaldeltaq, qscoredeltaq, positiondeltaq, dinucdeltaq, rg_to_int, use_oq = True, minscore = 6, maxscore = 42):
    complement = utils.Dinucleotide.complement
    
    original_quals = (utils.bamread_get_oq(read) if use_oq else seq.query_qualities)
    recalibrated_quals = np.array(original_quals, dtype = np.int)
    rg = rg_to_int[read.get_tag('RG')]

    valid_positions = (original_quals >= minscore)
    qcov = original_quals[valid_positions]
    cycle = bamread_cycle_covariates(read)[valid_positions]
    dinuccov = bamread_dinuc_covariates(read)[valid_positions]

    recalibrated_quals[valid_positions] = (meanq[rg] + globaldeltaq[rg] + qscoredeltaq[rg, qcov] + dinucdeltaq[rg, qcov, dinuccov] + positiondeltaq[rg, qcov, cycle]).astype(np.int)
    return recalibrated_quals

def get_delta_qs(meanq, rg_errs, rg_total, q_errs, q_total, pos_errs, pos_total, dinuc_errs, dinuc_total, maxscore = 42):
    # shapes are:
    #   [rg]
    #   [rg, q]
    #   [rg, q, covariate]
    #seqlen = pos_total.shape[2] / 2
    nrgs = meanq.shape[0]
    rgdeltaq = utils.gatk_delta_q(meanq, rg_errs, rg_total)
    # the qscoredeltaq is 1d, with the index being the quality score
    prior1 = np.broadcast_to((meanq + rgdeltaq)[:,np.newaxis], q_total.shape).copy()
    qscoredeltaq = utils.gatk_delta_q( prior1 , q_errs, q_total)
    ## positiondeltaq is 2d, first dimension is quality score and second is position
    ## dinucdeltaq is 2d, first dimension is quality score and second is nucleotide context
    prior2 = np.broadcast_to((prior1 + qscoredeltaq)[...,np.newaxis], pos_total.shape).copy()
    positiondeltaq = utils.gatk_delta_q(prior2, pos_errs, pos_total)
    prior3 = np.broadcast_to((prior1 + qscoredeltaq)[...,np.newaxis], dinuc_total.shape).copy()
    dinucdeltaq = utils.gatk_delta_q(prior3, dinuc_errs, dinuc_total)

    #need to add another value of dinuc, for invalid dinuc
    pad = np.zeros((len(dinucdeltaq.shape),2), dtype = np.int_)
    pad[-1,1] = 1 #add a 0 to the last axis
    dinucdq = np.pad(dinucdeltaq, pad_width = pad, mode = 'constant', constant_values = 0)

    return rgdeltaq.copy(), qscoredeltaq.copy(), positiondeltaq.copy(), dinucdq.copy()
