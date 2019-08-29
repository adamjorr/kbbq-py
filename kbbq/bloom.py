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
khmer.calc_expected_collisions(nodegraph, force = False, max_false_pos = .15) will return false positive rate and issue a warning if it's too high.
    when force is true, don't exit if false positive is higher than max.

ksize = 32
maxmem = khmer.khmer_args.memory_setting('8G')
fake_args = namedtuple(max_memory_usage = maxmem, n_tables = 4)
tablesize = calculate_graphsize(fake_args, 'nodegraph', multiplier = 1.0)
nodegraph = khmer.Nodegraph(ksize, tablesize, fake_args.n_tables) 
"""

import khmer
import khmer.khmer_args
import numpy as np

binomial_thresholds = 

def create_empty_nodegraph(ksize = 32, max_mem = '8G'):
    """
    Create an empty nodegraph. From this point, fill it from a file or start filling
    it with count.
    """
    mem = khmer.khmer_args.memory_setting(max_mem)
    fake_args = namedtuple(max_memory_usage = mem, n_tables = 4)
    tablesize = khmer.khmer_args.calculate_graphsize
    return khmer.Nodegraph(ksize, tablesize, fake_args.n_tables)

def count_read(read, graph):
    """
    Put ReadData read in the countgraph
    """
    graph.consume(np.str.join('',read.seq))

def kmers_in_graph(read, graph):
    """
    Query the graph for each kmer in read and return a :class:`numpy.ndarray` of bools.

    The returned array has length len(read) - ksize + 1
    """
    return np.array(graph.get_kmer_counts(np.str.join('',read.seq)), dtype = np.bool)

def n_kmers_in_graph(read, graph):
    """
    Get the number of kmers overlapping each read position that are in the graph.

    The returned array has length len(read).
    """

def n_kmers_possible(read):
    """
    Get the possible number of kmers overlapping each read position.

    For example, position 1 will always have only 1 overlapping kmer, position 2 will
    have 2, etc.
    """

def p_kmer_added(sampling_rate):
    """
    The probability a kmer was added to the graph.

    This is 1 - ( 1 - sampling_rate ) ^ n, where n is the assumed highest multiplicity
    of a weak kmer. When sampling_rate >= .1, n = 2. Otherwise, n = .2 / sampling_rate.

    The point is to perform a binomial test, with p = p_kmer_added, k = n_kmers_in_graph,
    N = n_kmers_possible. The null hypothesis is that the kmer is not an error, and this
    is p under that assumption. Since the kmer is not an error, it was observed 2 or greater
    times. The alternative hypothesis is that the kmer is an error, so the true
    probability the kmer was added to the graph was less than the number returned
    by this function.

    I think the better model would be a negative binomial. In this case, we have no idea
    how many times an overlapping kmer appeared in our dataset. We do know that there are
    r kmers that failed to be added to the hash and that each time an overlapping kmer
    appeared, there was a *sampling_rate* probability that it was added to the hash.
    """




