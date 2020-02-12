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
import collections

def table_to_vectors(table, rg_order, maxscore = 42):
    """
    Get a tuple of recalibration data vectors from a :class:`RecalibrationReport`.

    The recal table uses the PU of the read group as the read group entry in the table.
    See :func:`vectors_to_report` for more info

    The rg_order must be given to define the output order of read groups, since there
    is no natural ordering for them. If the PUs listed in the table are 'apple' and 
    'orange', and rg_order is ['apple', 'orange'], data from the 'apple' read group
    will be in index 0 and data from the 'orange' read group will be in index 1.

    :param table: Report to get data from
    :type table: :class:`RecalibrationReport`
    :param rg_order: List of PU's controlling the order of read groups.
    :type rg_order: list(str)
    :param maxscore: largest quality score possible
    :type maxscore: int
    :return: data vectors (meanq, global_errs, global_total, q_errs, q_total, pos_errs, pos_total, dinuc_errs, dinuc_total)
    :rtype: tuple of :class:`numpy.ndarray` of int
    """
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
    """
    Get the cycle covariate from a bam read.

    :param read: read to get covariate values from
    :type read: :class:`pysam.AlignedSegment`
    :return: cycle covariate
    :rtype: :class:`numpy.ndarray` of int
    """
    cycle = utils.generic_cycle_covariate(read.query_length, read.is_read2)
    if read.is_reverse:
        cycle = np.flip(cycle)
    return cycle

def bamread_dinuc_covariates(read, use_oq = True, minscore = 6):
    """
    Get the dinuc covariates from a bam read.

    These range from 0 to 15 for each possible dinucleotide, or -1 for an invalid
    value. The first base is always invalid (since there is no base before it) and
    any base overlapping an N will be invalid. Any base with a quality score below
    minscore is also invalid.

    :param read: read to get covariate values from
    :type read: :class:`pysam.AlignedSegment`
    :return: dinuc covariate
    :rtype: :class:`numpy.ndarray` of int
    """
    seq = read.query_sequence
    quals = (utils.bamread_get_oq(read) if use_oq else np.array(read.query_qualities, dtype = np.int))
    if read.is_reverse:
        seq = ''.join([utils.Dinucleotide.complement.get(x,'N') for x in reversed(seq)])
        quals = np.flip(quals)
    fulldinuc = np.zeros(read.query_length, dtype = np.int)
    dinuccov = utils.generic_dinuc_covariate(np.array(list(seq), dtype = 'U1'), quals, minscore)
    if read.is_reverse:
        # flip back to fwd coordinates
        dinuccov = np.flip(dinuccov)
    return dinuccov

def recalibrate_bamread(read, meanq, globaldeltaq, qscoredeltaq, positiondeltaq, dinucdeltaq, rg_to_int, use_oq = True, minscore = 6):
    """
    Recalibrate a bamread and return the new quality scores.

    Note that this function doesn't alter the given read.

    :param read: read to recalibrate
    :type read: :class:`pysam.AlignedSegment`
    :param others: delta q values where the indices represent the rg, q, and dinuc or cycle covariates.
    :type others: :class:`numpy.ndarray` of int
    :param rg_to_int: dictionary to convert rg IDs to ints
    :type rg_to_int: dict(str, int)
    :return: recalibrated quality scores
    :rtype: :class:`numpy.ndarray` of int
    """
    complement = utils.Dinucleotide.complement
    
    original_quals = (utils.bamread_get_oq(read) if use_oq else np.array(read.query_qualities, dtype = np.int))
    recalibrated_quals = np.array(original_quals, dtype = np.int)
    rg = rg_to_int[read.get_tag('RG')]

    valid_positions = (original_quals >= minscore)
    qcov = original_quals[valid_positions]
    cycle = bamread_cycle_covariates(read)[valid_positions]
    dinuccov = bamread_dinuc_covariates(read)[valid_positions]

    recalibrated_quals[valid_positions] = (meanq[rg] + globaldeltaq[rg] + qscoredeltaq[rg, qcov] + dinucdeltaq[rg, qcov, dinuccov] + positiondeltaq[rg, qcov, cycle]).astype(np.int)
    return recalibrated_quals

def get_delta_qs(meanq, rg_errs, rg_total, q_errs, q_total, pos_errs, pos_total, dinuc_errs, dinuc_total):
    """
    Find the covariate delta Q's as arrays given the data arrays.

    The input vector shapes are: [rg], [rg, q], [rg, q, covariate]. Returns a tuple
    of delta Q arrays. The indices specify the covariate values, the values specify
    how much the Q score should change based on its covariate values.

    The preferred method of getting the covariate delta Q's is the 
    :func:`get_modeldqs_from_covariates` function.

    :param all: data array
    :type all: :class:`numpy.ndarray` of int
    :return: Delta Q arrays (rgdeltaq, qscoredeltaq, positiondeltaq, dinucdeltaq)
    :rtype: tuple of :class:`numpy.ndarray` of int
    """
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

ModelDQs = collections.namedtuple('ModelDQs', ['mean','rg','q','cycle','dinuc'],
    module = __name__)
ModelDQs.__doc__ += """
    A class to hold the model delta Q's and enforce a consistent ordering.
    """
ModelDQs.mean.__doc__ += """The mean Q for each read group.
    1-dimensional :class:`np.ndarray` with length equal to the number of read groups."""
ModelDQs.rg.__doc__ += """Shift from the mean for each read group. 1-dimensional
    :class:`np.ndarray` with length equal to the number of read groups."""
ModelDQs.q.__doc__ += """Shift from the mean for each quality score. 2-dimensional
    :class:`np.ndarray` with shape (number of read groups, number of quality scores)"""
ModelDQs.cycle.__doc__ += """Shift from the mean for each quality score. 3-dimensional
    :class:`np.ndarray` with shape (number of read groups, number of quality scores,
    number of cycles)"""
ModelDQs.dinuc.__doc__ += """Shift from the mean for each quality score. 3-dimensional
    :class:`np.ndarray` with shape (number of read groups, number of quality scores,
    number of cycles)"""

def get_modelvecs_from_covariates(covariates):
    """
    Find the model vectors given covariate data.

    The model vectors are: meanq, global_errs, global_total, q_errs, q_total,
    pos_errs, pos_total, dinuc_errs, dinuc_total

    :param covariates: covariate data
    :type covariates: :class:`CovariateData`
    :returns: model vectors
    :rtype: tuple(:class:`numpy.ndarray`)
    """
    expected_errs = np.sum(utils.q_to_p(np.arange(covariates.qcov.num_qs()))[np.newaxis,:] * covariates.qcov.total, axis = 1)
    meanq = utils.p_to_q(expected_errs / covariates.qcov.rgcov.total)
    global_errs, global_total = covariates.qcov.rgcov.errors, covariates.qcov.rgcov.total
    q_errs, q_total = covariates.qcov.errors, covariates.qcov.total
    pos_errs, pos_total = covariates.cyclecov.errors, covariates.cyclecov.total
    dinuc_errs, dinuc_total = covariates.dinuccov.errors, covariates.dinuccov.total
    return meanq, global_errs, global_total, q_errs, q_total, pos_errs, pos_total, dinuc_errs, dinuc_total

def get_modeldqs_from_covariates(covariates):
    """
    Find the covariate delta Q's given covariate data.

    :param covariates: covariate data
    :type covariates: :class:`CovariateData`
    :returns: delta Q's
    :rtype: :class:`.ModelDQs`
    """
    nrgs = covariates.qcov.rgcov.num_rgs()
    expected_errs = np.sum(utils.q_to_p(np.arange(covariates.qcov.num_qs()))[np.newaxis,:] * covariates.qcov.total, axis = 1)
    meanq = utils.p_to_q(expected_errs / covariates.qcov.rgcov.total)
    rgdeltaq = utils.gatk_delta_q(meanq, *covariates.qcov.rgcov[...])
    prior1 = np.broadcast_to((meanq + rgdeltaq)[:,np.newaxis], covariates.qcov.shape).copy()
    qscoredeltaq = utils.gatk_delta_q(prior1, *covariates.qcov[...])
    prior2 = np.broadcast_to((prior1 + qscoredeltaq)[...,np.newaxis], covariates.cyclecov.shape).copy()
    cycledeltaq = utils.gatk_delta_q(prior2, *covariates.cyclecov[...])
    prior3 = np.broadcast_to((prior1 + qscoredeltaq)[...,np.newaxis], covariates.dinuccov.shape).copy()
    dinucdeltaq = utils.gatk_delta_q(prior3, *covariates.dinuccov[...])
    #need to add another value of dinuc, for invalid dinuc
    pad = np.zeros((len(dinucdeltaq.shape),2), dtype = np.int_)
    pad[-1,1] = 1 #add a 0 to the last axis
    dinucdq = np.pad(dinucdeltaq, pad_width = pad, mode = 'constant', constant_values = 0)
    return ModelDQs(meanq, rgdeltaq, qscoredeltaq, cycledeltaq, dinucdq)
