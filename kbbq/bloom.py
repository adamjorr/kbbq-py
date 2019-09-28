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
    hashes = np.array(graph.get_kmer_hashes(np.str.join('',read.seq)), dtype = np.ulonglong)
    sampled = np.random.choice([True, False],
        size = hashes.shape,
        replace = True,
        p = [sampling_rate, 1.0 - sampling_rate])
    for h in hashes[sampled]:
        graph.count(h.item())

def kmers_in_graph(read, graph):
    """
    Query the graph for each kmer in read and return a :class:`numpy.ndarray` of bools.

    The returned array has length len(read) - ksize + 1
    """
    return np.array(graph.get_kmer_counts(np.str.join('',read.seq)), dtype = np.bool)

def overlapping_kmers_in_graph(read, graph):
    """
    Get the number of kmers overlapping each read position that are in the graph.

    The returned array has length len(read).
    This only behaves well for len(seq) > 2ksize.
    """
    ksize = graph.ksize()
    assert ksize <= len(read)
    assert len(read) > 2 * ksize
    kmers = kmers_in_graph(read, graph)
    num_in_graph = np.zeros(len(read), dtype = np.int)
    for i in range(ksize - 1):
        #each base is overlapped by < k kmers
        num_in_graph[i] = np.sum(kmers[0:i])
        num_in_graph[-1 - i] = np.sum(kmers[-1 - i:])
    for kidx, readidx in zip(range(len(kmers) - ksize + 1), range(ksize - 1, len(read) - ksize + 1)): #kidx len: len(read) - 2 ksize + 2; readidx range: len(read) - 2 ksize + 2
        #each base is overlapped by k kmers
        num_in_graph[readidx] = np.sum(kmers[kidx:kidx + ksize])
    return num_in_graph

def overlapping_kmers_possible(read, ksize):
    """
    Get the possible number of kmers overlapping each read position.

    For example, position 1 will always have only 1 overlapping kmer, position 2 will
    have 2, etc.
    """
    assert ksize <= len(read)
    assert len(read) > 2 * ksize
    num_possible = np.zeros(len(read), dtype = np.int)
    #each base is overlapped by < k kmers
    num_possible[np.arange(ksize - 1)] = np.arange(ksize - 1) + 1
    num_possible[-np.arange(ksize - 1)-1] = np.arange(ksize - 1) + 1
    #each base is overlapped by k kmers
    num_possible[ksize - 1 : len(read) - ksize + 1] = ksize
    return num_possible

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
    # #https://rigtorp.se/2011/01/01/rolling-statistics-numpy.html
    ksize = graph.ksize()
    hashes = np.array(graph.get_kmer_hashes(np.str.join('',read.seq)), dtype = np.ulonglong)
    errors = np.lib.stride_tricks.as_strided(read.errors,
        shape = (len(read.errors) - ksize + 1, ksize),
        strides = read.errors.strides * 2) #2D array of shape (nkmers, ksize)
    trusted_kmers = np.all(~errors, axis = 1)
    for h in hashes[trusted_kmers]:
        graph.count(h.item())

def infer_errors_from_trusted_kmers(read, graph):
    trusted_kmers = kmers_in_graph(read, graph)
    errors = np.zeros(len(read), dtype = np.bool)
    if np.all(trusted_kmers) or np.all(~trusted_kmers): #do nothing
        return errors
    else:
        transitions = np.nonzero(np.diff(trusted_kmers) != 0)[0] + 1 #indices where changes occur
        segments = np.concatenate([np.array([0]), transitions, np.array([len(trusted_kmers)])]) #indices including beginning and end
        segment_pairs = np.lib.stride_tricks.as_strided(segments,
            shape = (len(segments) - 2 + 1, 2),
            strides = segments.strides * 2) #numsegments, 2
        segment_lens = np.diff(segments) #lengths of each segment
        trusted_segments = trusted_kmers[segment_pairs[:,0]] #hack
        # trusted_segments = [np.all(trusted_kmers[s[0] : s[1]]) for s in segment_pairs]
        trusted_segment_lens = np.array(segment_lens, copy = True)
        trusted_segment_lens[~trusted_segments] = 0
        longest_trusted = np.argmax(trusted_segment_lens)
        #note argmax will pick the first in case of a tie
        #right side
        #kmer segments[longest + 1] is an error
        #for kmer k in range(len(seq) + ksize - 1)
            #base to look at is base k + ksize - 1
        k = segment_pairs[longest_trusted, 1]
        while k < len(read) - graph.ksize() + 1:
            if trusted_kmers[k]:
                k = k + 1
            else:
                cor_len = correction_len(read.seq[k:], graph, right = True)
                if len(cor_len) == 1: #there was not a tie
                    errors[k + graph.ksize() - 1] = True
                k = k + cor_len[0]

        #left side
        #fix below
        k = segment_pairs[longest_trusted, 0] - 1
        while k >= 0:
            if trusted_kmers[k]:
                k = k - 1
            else:
                cor_len = correction_len(read.seq[:(k+graph.ksize())], graph, right = False)
                if len(cor_len) == 1: #no tie
                    errors[k] = True
                k = k - cor_len[0]
        return errors

def correction_len(seq, graph, right = True):
    """
    Get the number of corrections given an ndarray of sequence.

    This value is an array of values between between 1 and ksize, or len(seq) if no
    correction can be made. The length of the array represents the number of results if
    there is a tie.
    """
    ksize = graph.ksize()
    kmers = np.lib.stride_tricks.as_strided(seq.copy(),
        shape = (len(seq) - ksize + 1, ksize),
        strides = 2 * seq.strides) #note the memory is preserved across windows,
        #so changing one base in one kmer will change every kmer!
        #this is exactly the behavior we want so we can exploit this for efficiency
    largest_possible_fix = min(ksize, len(kmers))

    bases = list("ACGT")
    counts = np.zeros(len(bases), dtype = np.int)
    idx = np.array([0,-1])
    possible_fixes = range(largest_possible_fix)
    if not right:
        np.flip(idx)
        possible_fixes = reversed(possible_fixes)
    for b, base in enumerate(bases):
        kmers[idx] = base #kmers[0,-1] or kmers[-1,0]
        for i in possible_fixes:
            if not graph.get(np.str.join('',kmers[i])):
                counts[b] = i if right else largest_possible_fix - i - 1
                break
        else: #we made it through every possible fix
            if right and largest_possible_fix == len(kmers): #we ran out of kmers, try to extend
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
            else: #if every kmer is corrected and we saw ksize kmers, we move forward k
                counts[b] = largest_possible_fix
    if np.all(counts == 0):
        #end correction if we cannot find any
        #https://github.com/mourisl/Lighter/blob/df39031f8254f8351852f9f8b51b643475226ea0/ErrorCorrection.cpp#L574
        return np.array(len(seq), dtype = np.int)
    else:
        m = np.amax(counts)
        #we may also want to test if there are multiple maxima
        return counts[counts == m] #if there are multiple this will be an array

def fix_overcorrection(read, ksize, minqual = 6, window = 20, threshold = 4):
    corrections = read.errors.copy()
    corrections_windowed = rolling_window(corrections, window)
    correction_count = np.array(rolling_window(corrections, window), dtype = np.int)
    seq = rolling_window(read.seq, window)
    quals = rolling_window(read.qual, window)
    correction_count[seq == 'N'] = 0
    correction_count[quals < minqual] = .5
    overcorrected = np.sum(correction_count, axis = 1) > threshold
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
    shape = a.shape[:-1] + (a.shape[-1] - window + 1, window)
    strides = a.strides + (a.strides[-1],)
    return np.lib.stride_tricks.as_strided(a, shape=shape, strides=strides)
