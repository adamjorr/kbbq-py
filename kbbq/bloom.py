"""
Methods to create and manipulate khmer bloom filters.

The base khmer bloom filter object is a Nodegraph.
Nodegraph.get(kmer) returns 0 or 1 depending on whether the kmer is present
Nodegraph.count(kmer) sets the value to 1 (present)
Nodegraph.add(khmer) sets the value to 1 (present) and returns true if it's new.
Nodegraph.consume(str) does add on each kmer in the string and returns the number of kmers added
Nodegraph.get_kmer_counts(str) does get on each kmer in the string and returns a list of counts
Nodegraph.save(filename) saves to filename
Nodegraph.load(filename) loads from filename
Nodegraph.hash(kmer) returns the hash of the k-mer
Nodegraph.get_kmer_hashes(kmer) hashes each k
khmer.calc_expected_collisions(nodegraph, force = False, max_false_pos = .15) will return false positive rate and issue a warning if it's too high.
when force is true, don't exit if false positive is higher than max.
#may also want to try a Nodetable at some point. It should have the same fns.

ksize = 32
maxmem = khmer.khmer_args.memory_setting('8G')
fake_args = namedtuple(max_memory_usage = maxmem, n_tables = 4)
tablesize = calculate_graphsize(fake_args, 'nodegraph', multiplier = 1.0)
nodegraph = khmer.Nodegraph(ksize, tablesize, fake_args.n_tables) 

"""

import khmer
import khmer.khmer_args
import numpy as np
import scipy.stats
import collections

overlapcache = {}

def create_empty_nodegraph(ksize = 32, max_mem = '8G'):
    """
    Create an empty nodegraph. From this point, fill it from a file or start filling
    it with count.
    """
    mem = khmer.khmer_args.memory_setting(max_mem)
    fake_args = collections.namedtuple('fake_args', ['max_memory_usage', 'n_tables'])
    fargs = fake_args(max_memory_usage = mem, n_tables = 4)
    tablesize = khmer.khmer_args.calculate_graphsize(fargs, 'nodegraph', multiplier = 1.0)
    return khmer.Nodegraph(ksize, tablesize, fargs.n_tables)

def count_read(read, graph, sampling_rate):
    """
    Put ReadData read in the countgraph. Each k-mer will be added with p = sampling_rate
    """
    # ksize = graph.ksize()
    # #https://rigtorp.se/2011/01/01/rolling-statistics-numpy.html
    # kmers = numpy.lib.stride_stricks.as_strided(read.seq,
    #     shape = (len(read.seq) - ksize + 1, ksize),
    #     strides = read.seq.strides * 2) #2D array of shape (nkmers, ksize)
    #probably will be fastest to get all hashes in C then select the ones i want
    hashes = np.array(graph.get_kmer_hashes(read.seq.tostring().decode('utf32')), dtype = np.ulonglong)
    sampled = np.random.choice([True, False],
        size = hashes.shape,
        replace = True,
        p = [sampling_rate, 1.0 - sampling_rate])
    contains_n = np.any(rolling_window(read.seq, graph.ksize()) == 'N', axis = 1)
    sampled[contains_n] = False
    for h in hashes[sampled]:
        graph.count(h.item())

def kmers_in_graph(read, graph):
    """
    Query the graph for each kmer in read and return a :class:`numpy.ndarray` of bools.

    The returned array has length len(read) - ksize + 1
    """
    ingraph = np.array(graph.get_kmer_counts(read.seq.tostring().decode('utf32')), dtype = np.bool)
    ingraph[np.any(rolling_window(read.seq, graph.ksize()) == 'N', axis = 1)] = False
    return ingraph

def overlapping_kmers_in_graph(read, graph):
    """
    Get the number of kmers overlapping each read position that are in the graph.

    The returned array has length len(read).
    """
    ksize = graph.ksize()
    kmers = kmers_in_graph(read, graph)
    num_in_graph = np.zeros(len(read), dtype = np.int)
    num_in_graph_windowed = rolling_window(num_in_graph, ksize)
    np.add.at(num_in_graph_windowed, (kmers,), 1) #this will modify num_in_graph
    return num_in_graph

def overlapping_kmers_possible(read, ksize):
    """
    Get the possible number of kmers overlapping each read position.

    For example, position 1 will always have only 1 overlapping kmer, position 2 will
    have 2, etc.

    This is cached.
    """
    if ksize > len(read):
        raise ValueError(f"ksize {ksize} too small for read with length {len(read)}. \
            ksize must be <= the read length.")
    if len(read) not in overlapcache:
        num_possible = np.zeros(len(read), dtype = np.int)
        koverlaps = rolling_window(num_possible, ksize)
        np.add.at(koverlaps, Ellipsis, 1)
        num_possible.flags.writeable = False
        overlapcache[len(read)] = num_possible
    return overlapcache[len(read)]

def p_kmer_added(sampling_rate, graph):
    """
    The probability a kmer was added to the graph, including the false positive rate.

    P(A) =  1 - ( 1 - sampling_rate ) ^ n, where n is the assumed highest multiplicity
    of a weak kmer. When sampling_rate >= .1, n = 2. Otherwise, n = .2 / sampling_rate.
    This function returns P*(A) = P(A) + B - B*P(A), where B is the false positive rate.

    The point is to perform a binomial test, with p = p_kmer_added, k = n_kmers_in_graph,
    N = n_kmers_possible. The null hypothesis is that the kmer is not an error, and this
    is p under that assumption. Since the kmer is not an error, it was observed 2 or greater
    times. The alternative hypothesis is that the kmer is an error, so the true
    probability the kmer was added to the graph was less than the number returned
    by this function.

    I think the better model would be a negative binomial. In this case, we have no idea
    how many times an overlapping kmer appeared in our dataset. We do know that there are
    r kmers that failed to be added to the hash and that each time an overlapping kmer
    appeared, there was a *sampling_rate* probability that it was added to the hash. Then
    the number of overlapping kmers in the dataset can be predicted with a negative binomial
    distribution with r = (# in hash) and p = 1 - sampling_rate. If this distribution is
    X, the read coverage at that site is X / (# k-mers that could overlap the site); in
    most cases, X / k.

    Thus with this technique, we can calculate the distribution of site-by-site coverage
    in the dataset. In marginalizing, we sum the probability vector of each X then divide
    by the number of sites. Once we're done, we must also divide by X again, since each base
    with coverage x (element of X) was counted x times.

    In our sequencing model, the number of sequenced reads at any position is Poisson.
    Because errors are correlated, the number of sequenced erroneous reads is an overdispersed
    Poisson. This can be accurately modeled with a different negative binomial.
    Thus from our inference of coverage we should be able to fit a mixture of negative
    binomials that represents the coverage expected in the dataset, and this binomial will
    tell us the probability that a read is an error given its coverage. Note that we use
    the negative binomial for two distinct purposes: one to estimate the coverage given a
    site at a read, and one to model the total coverage of the dataset.

    Notably, this approach suggests an optimal choice for the sampling rate since the
    distribution becomes less predictive for some coverages at extreme points. A strategy
    for picking the theoretically optimal sampling rate can likely be found.
    """
    fpr = khmer.calc_expected_collisions(graph, force = False, max_false_pos = .15)
    exp = .2 / sampling_rate if sampling_rate < .1 else 2
    p_a = 1 - ( 1 - sampling_rate ) ** exp
    p_added = p_a + fpr - fpr * p_a
    return p_added

def calculate_thresholds(p_added, ksize):
    """
    Calculate the thresholds. If the number of overlapping kmers is less than this number,
    the base is inferred to be erroneous.
    """
    dists = [scipy.stats.binom(n = i, p = p_added) for i in range(1, ksize + 1, 1)] 
    return np.array([0] + [int(d.ppf(.995)) for d in dists])

def infer_errors(overlapping, possible, thresholds):
    """
    Perform a binomial test to infer which bases are erroneous.

    Given a vector of the number of overlapping kmers in the hash and the number of
    possible kmers in the hash for each position, infer whether each base is erroneous.

    This is done with essentially a binomial test. We evaluate the probability the sample
    came from a binomial distribution with p = p_added. We do a right-tailed test with
    the null hypothesis that all kmers overlapping the site are erroneous. We assume
    the multiplicity of a weak kmer is less than ~2, so this is reflected in p_added.
    We only use the right tail because a lower p parameter supports the null, while a
    high value supports the alternative.
    """
    #n = 0 doesn't make sense
    return overlapping <= thresholds[possible] #if overlapping < thresholds, it's an error.

def infer_read_errors(read, graph, thresholds):
    """
    Return an array of errors given a graph and the thresholds.
    """
    overlapping = overlapping_kmers_in_graph(read, graph)
    possible = overlapping_kmers_possible(read, graph.ksize())
    errors = infer_errors(overlapping, possible, thresholds)
    assert len(errors) == len(read)
    return errors

def add_trusted_kmers(read, graph):
    """
    Add trusted kmers to graph.
    """
    ksize = graph.ksize()
    hashes = np.array(graph.get_kmer_hashes(read.seq.tostring().decode('utf32')), dtype = np.ulonglong)
    errors = rolling_window(read.errors, ksize)
    trusted_kmers = np.all(~errors, axis = 1)
    for h in hashes[trusted_kmers]:
        graph.count(h.item())

def find_longest_trusted_block(trusted_kmers):
    """
    Given a boolean array describing whether each kmer is trusted, return a pair of
    indices to the start and end of the block.

    The indices will be standard for python ranges; inclusive on the left side and
    exclusive on the right side.
    """
    transitions = np.nonzero(np.diff(trusted_kmers) != 0)[0] + 1 #indices where changes occur
    segments = np.concatenate([np.array([0]), transitions, np.array([len(trusted_kmers)])]) #indices including beginning and end
    segment_pairs = rolling_window(segments, 2) #numsegments, 2
    segment_lens = np.diff(segments) #lengths of each segment
    trusted_segments = trusted_kmers[segment_pairs[:,0]] #hack
    trusted_segment_lens = np.array(segment_lens, copy = True)
    trusted_segment_lens[~trusted_segments] = 0
    longest_trusted = np.argmax(trusted_segment_lens)
    #argmax will pick the first in case of tie
    return segment_pairs[longest_trusted,0], segment_pairs[longest_trusted,1]

def infer_errors_from_trusted_kmers(read, graph):
    """
    Return an array of errors and a bool describing whether multiple corrections were
    made.
    """
    trusted_kmers = kmers_in_graph(read, graph)
    errors = np.zeros(len(read), dtype = np.bool)
    multiple = False
    if np.all(trusted_kmers) or np.all(~trusted_kmers): #do nothing
        return errors, multiple
    else:
        longest_trusted = find_longest_trusted_block(trusted_kmers)
        ksize = graph.ksize()
        #right side
        #kmer trusted_kmers[longest_trusted[1]] is an error
        #for kmer k in range(len(seq) + ksize - 1)
            #base to look at is base k + ksize - 1
        # longest_trusted_len = longest_trusted[1] - longest_trusted[0]
        # if longest_trusted_len < ksize:
        #     k = longest_trusted[1]
        # else:
        #     k = longest_trusted[1] - 1
        #     trusted_kmers[k] = False
        k = longest_trusted[1]
        while k < len(trusted_kmers):
            if trusted_kmers[k]:
                k = k + 1
            else:
                cor_len, base, m = correction_len(read.seq[k:], graph, right = True)
                multiple = multiple or m
                if cor_len is not None: #correction found
                    # if read.seq[k + ksize - 1] != base:  #not equal to current base
                    errors[k + ksize - 1] = True
                    read.seq[k + ksize - 1] = base
                    k = k + cor_len[0]
                else:
                    #could not find a fix; try chopping up the read and trying again
                    #need to make sure this doesn't loop forever somehow; in lighter this
                    #happens once and only once
                    #it should be OK since we will eventually run out of trusted kmers
                    # print('k+ksize:',k + ksize - 1)
                    subread = read[(k+ksize-1):]
                    if len(subread) > (len(read) / 2) or (len(subread) > ksize * 2):
                        # print('trying again')
                        errors[k+ksize-1:], m = infer_errors_from_trusted_kmers(subread, graph)
                        multiple = multiple or m
                    break
                
        #left side
        k = longest_trusted[0] - 1 # = -1 if the trusted block is at the start
        while k >= 0:
            if trusted_kmers[k]:
                k = k - 1
            else:
                cor_len, base, m = correction_len(read.seq[:(k+ksize)], graph, right = False)
                multiple = multiple or m
                if cor_len is not None: #correction found
                    errors[k] = True
                    read.seq[k] = base
                    k = k - cor_len[0]
                else: #could not find a fix; try chopping up the read and trying again
                    #need to make sure this doesn't loop forever somehow; in lighter this
                    #happens once and only once
                    #it should be OK since we will eventually run out of trusted kmers
                    subread = read[:(k+ksize)]
                    if len(subread) > (len(read) / 2) or (len(subread) > ksize * 2):
                        errors[:(k+ksize)], m = infer_errors_from_trusted_kmers(subread, graph)
                        multiple = multiple or m
                    break
        return errors, multiple

def correction_len(seq, graph, right = True):
    """
    Get the number of corrections given an ndarray of sequence.

    This value is an array of values between between 1 and ksize, or len(seq) if no
    correction can be made. The length of the array represents the number of results if
    there is a tie.

    The last element returned says whether there were multiple corrections that led
    to a trusted kmer.

    This function is in desperate need of a refactor.
    """
    ksize = graph.ksize()
    kmers = rolling_window(seq.copy(), ksize)
        #note the memory is preserved across windows,
        #so changing one base in one kmer will change every kmer!
        #this is exactly the behavior we want so we can exploit this for efficiency
    largest_possible_fix = min(ksize, len(kmers))
    bases = list("ACGT")
    counts = np.zeros(len(bases), dtype = np.int)
    if right:
        idx = (0,-1)
        possible_fixes = list(range(largest_possible_fix))
    else:
        idx = (-1,0)
        possible_fixes = list(range(-1,-largest_possible_fix-1, -1))
    for b, base in enumerate(bases):
        kmers[idx] = base #kmers[0,-1] or kmers[-1,0]
        for i in possible_fixes:
            # print(np.str.join('',kmers[i]), i)
            if not graph.get(kmers[i].tostring().decode('utf32')):
                if right:
                    counts[b] = i
                else:
                    counts[b] = (-i - 1)
                break
        else: #we made it through every possible fix
            if largest_possible_fix == len(kmers): #we ran out of kmers, try to extend
                if right:
                    last_kmer = np.str.join('',kmers[-1])
                    for j in range(ksize - len(kmers)): #ksize - len(kmers) is the number of kmers we didn't see at the end
                        for extra in bases:
                            if graph.get(last_kmer[1:] + extra):
                                last_kmer = last_kmer[1:] + extra
                                break #stop looking at more bases because we found one
                        else: #we didn't find an appropriate base; we're done extending
                            counts[b] = largest_possible_fix + j
                            break
                    else: #we extended and got to see all ksize kmers
                        counts[b] = ksize
                else: #left side
                    last_kmer = np.str.join('',kmers[0])
                    for j in range(ksize - len(kmers)): #ksize - len(kmers) is the number of kmers we didn't see at the end
                        for extra in bases:
                            if graph.get(extra + last_kmer[:-1]):
                                last_kmer = extra + last_kmer[:-1]
                                break #stop looking at more bases because we found one
                        else: #we didn't find an appropriate base; we're done extending
                            counts[b] = largest_possible_fix + j
                            break
                    else: #we extended and got to see all ksize kmers
                        counts[b] = ksize
            else: #if every kmer is corrected and we saw ksize kmers, we move forward k
                counts[b] = largest_possible_fix
    if np.all(counts == 0):
        #end correction if we cannot find any
        #https://github.com/mourisl/Lighter/blob/df39031f8254f8351852f9f8b51b643475226ea0/ErrorCorrection.cpp#L574
        return None, None, False
    else:
        m = np.amax(counts)
        #we may also want to test if there are multiple maxima
        largest = counts[counts == m] #if there are multiple this will be an array
        largest[largest > largest_possible_fix] = largest_possible_fix #if we extended we only want to continue k
        if len(largest) > 1 and largest[0] < largest_possible_fix: #if there's a tie and the tie can't fix all k kmers
            largest = None
        return largest, bases[counts.argmax()], len(counts[counts != 0]) > 1

def fix_one(seq, graph):
    """
    If we don't start with an anchor we just fix one base that produces the largest
    number of trusted kmers.
    """
    original_seq = seq.copy()
    ksize = graph.ksize()
    bases = list("ACGT")
    best_fix_len = 0
    best_fix_base = None
    best_fix_pos = None
    for i in range(len(seq)):
        modified_seq = original_seq.copy()
        kmers = rolling_window(modified_seq, ksize)
        for b in bases:
            modified_seq[i] = b
            trusted_kmers = np.array(graph.get_kmer_counts(modified_seq.tostring().decode('utf32')), dtype = np.bool)
            start_pos = int(min(max(i - ksize/2 + 1, 0), len(seq) - ksize))
            if trusted_kmers[start_pos]:
                num_in_graph = np.sum(trusted_kmers[start_pos:])
                if num_in_graph > best_fix_len:
                    best_fix_base = b
                    best_fix_pos = i
                    best_fix_len = num_in_graph
    return best_fix_len, best_fix_base, best_fix_pos


def fix_overcorrection(read, ksize, minqual = 6, window = 20, threshold = 4, adjust = False):
    """
    The threshold is adjusted when there are not multiple options for any correction.
    It is only adjusted for positions [window:-window]
    """
    corrections = read.errors.copy()
    corrections_windowed = rolling_window(corrections, window)
    correction_count = np.array(rolling_window(corrections, window), dtype = np.double)
    seq = rolling_window(read.seq, window)
    quals = rolling_window(read.qual, window)
    correction_count[np.logical_and(seq == 'N', corrections_windowed)] = 0
    correction_count[np.logical_and(quals < minqual, corrections_windowed)] = .5

    thresh_ary = np.repeat(threshold, len(correction_count))
    if adjust:
        thresh_ary[window:-(window-1)] += 1
    overcorrected = np.greater(np.sum(correction_count, axis = 1), thresh_ary)
    overcorrected_sites = np.zeros(len(corrections), dtype = np.bool)
    overcorrected_sites_windowed = rolling_window(overcorrected_sites, window)
    overcorrected_sites_windowed[overcorrected] = corrections_windowed[overcorrected] #overcorrected_sites is modified
    
    #now anything within k of an overcorrected site we call overcorrected
    num_fixed = np.sum(overcorrected_sites)
    fixed_before = 0
    corrections_windowed = rolling_window(corrections, ksize)
    overcorrected_sites_windowed = rolling_window(overcorrected_sites, ksize)
    while num_fixed > 0:
        #mark anything within k of an overcorrected site as overcorrected
        #i'm guessing most of the time this just gets rid of all corrections
        any_overcorrected = np.any(overcorrected_sites_windowed, axis = 1)
        overcorrected_sites_windowed[any_overcorrected,:] = corrections_windowed[any_overcorrected]
        fixed_now = np.sum(overcorrected_sites)
        num_fixed = fixed_now - fixed_before
        fixed_before = fixed_now
    corrections[overcorrected_sites] = False
    return corrections

def rolling_window(a, window):
    """
    Use stride tricks to reshape an array into an array of sliding windows.

    Different values in the array will point to the same place in memory, so
    copy the array before altering it if that behavior is undesired.

    From https://rigtorp.se/2011/01/01/rolling-statistics-numpy.html
    """
    if window > a.shape[-1]:
        raise ValueError(f"Window size {window} is too large for array with last \
            dimension of size {a.shape[-1]}")
    shape = a.shape[:-1] + (a.shape[-1] - window + 1, window)
    strides = a.strides + (a.strides[-1],)
    return np.lib.stride_tricks.as_strided(a, shape=shape, strides=strides)
