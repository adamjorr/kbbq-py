"""
Utilities for emulating GATK's BQSR tool.

BQSR model construction hard clips soft clips and
trims adaptors. ApplyBQSR does not. So we need
different functions for each.
"""

import pysam
import numpy as np
import pandas as pd
import scipy.stats
from .. import compare_reads as utils
from .. import recaltable

###############################

# Covariate Functions

###############################


def bamread_bqsr_cycle(read):
    fullcycle = np.zeros(read.query_length, dtype = np.int) #full length
    cycle = utils.generic_cycle_covariate(read.query_alignment_length, read.is_read2) #excludes soft-clipped bases!
    #soft-clipped bases will be skipped in other code
    #so it's no problem that some cycles will stay at 0
    if read.is_reverse:
        cycle = np.flip(cycle)
    fullcycle[read.query_alignment_start:read.query_alignment_end] = cycle
    return fullcycle

def bamread_bqsr_dinuc(read, use_oq = True, minscore = 6):
    unclipped_start = read.query_alignment_start
    unclipped_end = read.query_alignment_end
    seq = read.query_sequence[unclipped_start:unclipped_end]
    quals = (utils.bamread_get_oq(read) if use_oq else np.array(read.query_qualities, dtype = np.int))
    quals = quals[unclipped_start:unclipped_end]
    if read.is_reverse:
        seq = ''.join([utils.Dinucleotide.complement.get(x,'N') for x in reversed(seq)])
        quals = np.flip(quals)
    fulldinuc = np.zeros(read.query_length, dtype = np.int)
    dinuccov = utils.generic_dinuc_covariate(
        np.array(list(seq), dtype = 'U1'),
        quals, minscore).copy()
    if read.is_reverse:
        # flip back to fwd coordinates
        dinuccov = np.flip(dinuccov)
    fulldinuc[unclipped_start:unclipped_end] = dinuccov
    return fulldinuc

def bam_to_bqsr_covariates(bamfileobj, fastafilename, var_pos, minscore = 6, maxscore = 42):
    """
    Given a BAM file object, FASTA reference file name and var_pos dict,
    get the standard covariate arrays.
    """
    rg_to_pu = utils.get_rg_to_pu(bamfileobj)
    nrgs = len(rg_to_pu.keys())
    rg_to_int = dict(zip(rg_to_pu, range(len(rg_to_pu))))
    fasta = pysam.FastaFile(fastafilename)
    #the below can probably be spun out to a function, i think we only use fullskips
    ref = {chrom : np.array(list(fasta.fetch(reference = chrom)), dtype = np.unicode) for chrom in fasta.references}
    varsites = {chrom : np.array(var_pos[chrom], dtype = np.int) for chrom in var_pos.keys()}
    fullskips = {chrom : np.zeros(len(ref[chrom]), dtype = np.bool) for chrom in ref.keys()}
    for chrom in fullskips.keys():
        variable_positions = varsites[chrom]
        fullskips[chrom][variable_positions] = True

    nreads = np.sum([s.total for s in bamfileobj.get_index_statistics()])
    read = next(bamfileobj)
    seqlen = len(read.query_qualities)

    rgs = np.zeros(seqlen, dtype = np.int_)
    meanq = np.zeros(nrgs, dtype = np.int_)
    expected_errs = np.zeros(nrgs, dtype = np.longdouble)
    rg_errs = np.zeros(nrgs, dtype = np.int_)
    rg_total = np.zeros(nrgs, dtype = np.int_)
    q_errs = np.zeros((nrgs, maxscore + 1), dtype = np.int_)
    q_total = np.zeros((nrgs, maxscore + 1), dtype = np.int_)
    pos_errs = np.zeros((nrgs, maxscore + 1, 2 * seqlen), dtype = np.int_)
    pos_total = np.zeros((nrgs, maxscore + 1, 2 * seqlen), dtype = np.int_)
    dinuc_errs = np.zeros((nrgs, maxscore + 1, 16), dtype = np.int_)
    dinuc_total = np.zeros((nrgs, maxscore + 1, 16), dtype = np.int_)
    
    try:
        while True:
            seq = read.query_sequence
            rgs[:] = rg_to_int[read.get_tag('RG')]
            errors, skips = utils.find_read_errors(read, ref, fullskips)
            q = utils.bamread_get_oq(read)
            pos = bamread_bqsr_cycle(read)
            dinucleotide = bamread_bqsr_dinuc(read)
            seq = np.array(list(seq), dtype = 'U1')
            trimmed = trim_bamread(read)
            
            skips[q < minscore] = True
            skips[trimmed] = True
            skips[seq == 'N'] = True

            valid = np.logical_not(skips)
            dinuc_valid = np.logical_and(dinucleotide != -1, valid)
            e_and_valid = np.logical_and(errors, valid)
            e_and_dvalid = np.logical_and(errors, dinuc_valid)

            rge = rgs[e_and_valid]
            rgv = rgs[valid]
            qe = q[e_and_valid]
            qv = q[valid]

            np.add.at(expected_errs, rgv, utils.q_to_p(qv))
            np.add.at(rg_errs, rge, 1)
            np.add.at(rg_total, rgv, 1)
            np.add.at(q_errs, (rge, qe), 1)
            np.add.at(q_total, (rgv, qv), 1)
            np.add.at(pos_errs, (rge, qe, pos[e_and_valid]), 1)
            np.add.at(pos_total, (rgv, qv, pos[valid]), 1)
            np.add.at(dinuc_errs, (rgs[e_and_dvalid], q[e_and_dvalid], dinucleotide[e_and_dvalid]), 1)
            np.add.at(dinuc_total, (rgs[dinuc_valid], q[dinuc_valid], dinucleotide[dinuc_valid]), 1)
            read = next(bamfileobj)
    except StopIteration:
        pass
    meanq = utils.p_to_q(expected_errs / rg_total)
    return meanq, rg_errs, rg_total, q_errs, q_total, pos_errs, pos_total, dinuc_errs, dinuc_total

#############################################

# Trimming Functions

#############################################

def bamread_adaptor_boundary(read):
    #https://github.com/broadinstitute/gatk/blob/43b2b3bd4e723552414b32b8b2a7341b81f1f688/src/main/java/org/broadinstitute/hellbender/utils/read/ReadUtils.java#L534
    #0-based in ref coordinates
    if ( read.tlen == 0 or
        not read.is_paired or
        read.is_unmapped or
        read.mate_is_unmapped or
        read.is_reverse == read.mate_is_reverse):
            return None
    if read.is_reverse:
        #next_reference_start is 1-based
        #reference_start is 0-based
        #reference_end is 0-based but points to 1 past the last base
        #   (so essentially it's 1-based)
        if (read.reference_end - 1) > (read.next_reference_start):
            #"well-formed" insert len
            return read.next_reference_start - 1
        else:
            return None
    else:
        if read.reference_start <= (read.next_reference_start + read.tlen):
            #"well-formed" insert len
            return read.reference_start + abs(read.tlen)
        else:
            return None


def trim_bamread(read):
    #https://github.com/broadinstitute/gatk/blob/b11abd12b7305767ed505a8ff644a63659abf2cd/src/main/java/org/broadinstitute/hellbender/utils/clipping/ReadClipper.java#L388
    #return an array of seqlen which includes bases to skip
    #next_reference_start is 0-based
    #reference_start is 0-based
    #reference_end is 0-based but points to 1 past the last base
    #   (so essentially it's 1-based)
    adaptor_boundary = bamread_adaptor_boundary(read)
    skips = np.zeros(len(read.query_qualities), dtype = np.bool)
    if adaptor_boundary is None:
        return skips
    else:
        if read.is_reverse:
            if adaptor_boundary >= read.reference_start:
                #clip from start (left)
                #we need to get the boundary in read coordinates rather than ref
                boundary_reached = False
                for readidx, refidx in reversed(read.get_aligned_pairs()):
                    if refidx is not None and refidx <= adaptor_boundary:
                        boundary_reached = True
                    if boundary_reached and readidx is not None:
                        adaptoridx = readidx + 1 #slice syntax
                        break
                else:
                    #couldn't find boundary
                    #I think this can only happen if the boundary lies in a
                    #deletion that covers the rest of the read.
                    adaptoridx = 0
                skips[:adaptoridx] = True #skip first x bases
            return skips
        else:
            if adaptor_boundary <= (read.reference_end - 1):
                #clip from end (right)
                #reference_end is 1 past the end
                #reference_end - 1 - adaptor_boundary + 1
                boundary_reached = False
                for readidx, refidx in read.get_aligned_pairs():
                    if refidx is not None and refidx >= adaptor_boundary:
                        boundary_reached = True
                    if boundary_reached and readidx is not None:
                        adaptoridx = readidx
                        break
                else:
                    #couldn't find boundary
                    #I think this can only happen if the boundary lies in a
                    #deletion that covers the rest of the read.
                    adaptoridx = len(skips)
                skips[adaptoridx:] = True #skip last x bases
            return skips

##################################################

# RecalibrationReport Creation Functions

##################################################

def quantize(q_errs, q_total, maxscore = 93):
    """
    This function doesn't match the GATK version, but
    it's not used so it's not a priority.
    """
    qe = np.sum(q_errs, axis = 0)
    qt = np.sum(q_total, axis = 0)
    unobserved = (qt == 0)
    quantizer = np.arange(maxscore + 1)
    quantizer[0:qt.shape[0]][unobserved] = maxscore
    quantizer[qt.shape[0]:] = maxscore
    return quantizer

def vectors_to_report(meanq, global_errs, global_total, q_errs, q_total,
    pos_errs, pos_total, dinuc_errs, dinuc_total, rg_order, maxscore = 42):
    """
    Turn the set of recalibration vectors into a
    :class:`kbbq.recaltable.RecalibrationReport` object.

    For the recalibration vectors, each dimension corresponds to a covariate.
    The first index is always the read group, and the second (if it exists)
    represents the raw quality score, the final index is either the cycle or
    dinucleotide covariate.

    :param meanq: Mean q for each read group
    :type meanq: :class:`numpy.ndarray` [:]
    :param global_errs: Number of errors for each read group
    :type global_errs: :class:`numpy.ndarray` [:]
    :param global_total: Number of observations for each read group
    :type global_total: :class:`numpy.ndarray` [:]
    :param q_errs: Number of errors for each read group and q score subset.
    :type q_errs: :class:`numpy.ndarray` [:,:]
    :param q_total: Number of observations for each read group and q score subset.
    :type q_total: :class:`numpy.ndarray` [:,:]
    :param pos_errs: Number of errors for each read group, q, and cycle subset.
    :type pos_errs: :class:`numpy.ndarray` [:,:,:]
    :param pos_total: Number of observations for each read group, q, and cycle subset.
    :type pos_total: :class:`numpy.ndarray` [:,:,:]
    :param dinuc_errs: Number of errors for each read group, q, and dinucleotide subset.
    :type dinuc_errs: :class:`numpy.ndarray` [:,:,:]
    :param dinuc_total: Number of observations for each read group, q, and dinucleotide
        subset.
    :type dinuc_total: :class:`numpy.ndarray` [:,:,:]
    :param rg_order: The order of read groups
    :type rg_order: list(str)
    :param int maxscore: The maximum possible quality score
    :return: the recalibration table
    :rtype: :class:`kbbq.recaltable.RecalibrationReport`
    """

    #these will be mostly default values, except quantization
    #which I don't attempt to implement.
    #I'm afraid bad things will happen if I don't include at least null values
    #for all the args so I'll just include them all.
    #This may need to be cleaned up later.

    args = {
        'binary_tag_name' : 'null',
        'covariate' : 'ReadGroupCovariate,QualityScoreCovariate,ContextCovariate,CycleCovariate',
        'default_platform' : 'null',
        'deletions_default_quality' : '45',
        'force_platform' : 'null',
        'indels_context_size' : '3',
        'insertions_default_quality' : '45',
        'low_quality_tail' : '2',
        'maximum_cycle_value' : '500',
        'mismatches_context_size' : '2',
        'mismatches_default_quality' : '-1',
        'no_standard_covs' : 'false',
        'quantizing_levels' : '16',
        'recalibration_report' : 'null',
        'run_without_dbsnp' : 'false',
        'solid_nocall_strategy' : 'THROW_EXCEPTION',
        'solid_recal_mode' : 'SET_Q_ZERO'
        }
    argdata = {'Argument' : list(args.keys()),
    'Value' : list(args.values())
    }
    argtable = pd.DataFrame(data = argdata)

    rg_est_q = -10.0 * np.log10(np.sum(utils.q_to_p(np.arange(q_total.shape[1])) * q_total, axis = 1) / global_total).round(decimals = 5).astype(np.float)
    rg_est_q[np.isnan(rg_est_q)] = 0
    rgdata = {'ReadGroup' : rg_order,
        'EventType' : 'M',
        'EmpiricalQuality' : (utils.gatk_delta_q(rg_est_q, global_errs.copy(), global_total.copy()) + rg_est_q).astype(np.float),
        'EstimatedQReported' : rg_est_q,
        'Observations' : global_total,
        'Errors' : global_errs.astype(np.float)
        }
    rgtable = pd.DataFrame(data = rgdata)
    rgtable = rgtable[rgtable.Observations != 0]

    qualscore = np.broadcast_to(np.arange(q_total.shape[1]), (q_total.shape)).copy()
    qualdata = {'ReadGroup' : np.repeat(rg_order, q_total.shape[1]),
        'QualityScore' : qualscore.flatten(),
        'EventType' : np.broadcast_to('M', (q_total.shape)).flatten(),
        'EmpiricalQuality' : (utils.gatk_delta_q(qualscore.flatten(), q_errs.flatten(), q_total.flatten()) + qualscore.flatten()).astype(np.float),
        'Observations' : q_total.flatten(),
        'Errors' : q_errs.flatten().astype(np.float)
        }
    qualtable = pd.DataFrame(data = qualdata)
    qualtable = qualtable[qualtable.Observations != 0]

    #no quantization, but still have to make the quantization table
    #TODO: actual quant algo
    quantscores = np.arange(94)
    qcount = np.zeros(quantscores.shape)
    qcount[qualscore[0,]] = np.sum(q_total, axis = 0)
    quantized = quantize(q_errs, q_total) #TODO: actually quantize
    quantdata = {'QualityScore' : quantscores,
        'Count' : qcount,
        'QuantizedScore' : quantized
        }
    quanttable = pd.DataFrame(data = quantdata)

    dinuc_q = np.repeat(np.broadcast_to(np.arange(dinuc_total.shape[1]), (dinuc_total.shape[0:2])), dinuc_total.shape[2])
    dinuc_to_int = utils.Dinucleotide.dinuc_to_int
    covtable_colorder = ['ReadGroup','QualityScore','CovariateName','CovariateValue']
    dinucdata = {'ReadGroup' : np.repeat(rg_order, np.prod(dinuc_total.shape[1:])),
        'QualityScore' : dinuc_q.flatten(),
        'CovariateValue' : np.broadcast_to(np.array(utils.Dinucleotide.dinucs), dinuc_total.shape).flatten(),
        'CovariateName' : np.broadcast_to('Context', dinuc_total.shape).flatten(),
        'EventType' : np.broadcast_to('M',dinuc_total.shape).flatten(),
        'EmpiricalQuality' : (utils.gatk_delta_q(dinuc_q.flatten(), dinuc_errs.flatten(), dinuc_total.flatten()) + dinuc_q.flatten()).astype(np.float),
        'Observations' : dinuc_total.flatten(),
        'Errors' : dinuc_errs.flatten().astype(np.float)
        }
    dinuctable = pd.DataFrame(data = dinucdata)
    dinuctable = dinuctable[dinuctable.Observations != 0]

    cycle_q = np.repeat(np.broadcast_to(np.arange(pos_total.shape[1]), (pos_total.shape[0:2])), pos_total.shape[2])
    ncycles = pos_total.shape[2] / 2
    cycle_values = np.concatenate([np.arange(ncycles) + 1, np.flip(-(np.arange(ncycles)+1),axis=0)]).astype(np.int)
    cycledata = {'ReadGroup' : np.repeat(rg_order, np.prod(pos_total.shape[1:])).flatten(),
        'QualityScore' : cycle_q.flatten(),
        'CovariateValue' : np.broadcast_to(cycle_values, pos_total.shape).astype(np.unicode).flatten(),
        'CovariateName' : np.broadcast_to('Cycle',pos_total.shape).flatten(),
        'EventType' : np.broadcast_to('M',pos_total.shape).flatten(),
        'EmpiricalQuality' : (utils.gatk_delta_q(cycle_q.flatten(), pos_errs.flatten(), pos_total.flatten()) + cycle_q.flatten()).astype(np.float),
        'Observations' : pos_total.flatten(),
        'Errors' : pos_errs.flatten().astype(np.float)
        }
    cycletable = pd.DataFrame(data = cycledata)
    covariatetable = dinuctable.append(cycletable)
    covariatetable = covariatetable.set_index(covtable_colorder)
    covariatetable = covariatetable[covariatetable.Observations != 0]
    covariatetable = covariatetable.swaplevel('CovariateValue','CovariateName')
    covariatetable = covariatetable.sort_index(level = 0, sort_remaining = True)
    covariatetable = covariatetable.reset_index()
    #we do this to fix ordering because concatenating the tables ruins it

    titles = ['Arguments','Quantized','RecalTable0','RecalTable1','RecalTable2']
    descriptions = ['Recalibration argument collection values used in this run',
        'Quality quantization map', '' , '' , '']
    gatktables = [recaltable.GATKTable(title, desc, table) for title, desc, table in \
        zip(titles, descriptions, [argtable, quanttable, rgtable, qualtable, covariatetable])]

    return recaltable.RecalibrationReport(gatktables)

def bam_to_report(bamfileobj, fastafilename, var_pos):
    rgs = list(utils.get_rg_to_pu(bamfileobj).values())
    *vectors, = bam_to_bqsr_covariates(bamfileobj, fastafilename, var_pos)
    return vectors_to_report(*vectors, rgs)
