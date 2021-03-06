""" DNA related algorithms
"""
try:
    import matplotlib.pyplot as plt
except ImportError:
    pass
import itertools
import Queue
import sys

from bioinformatics_algorithms.data_structures import dictionaries

def pattern_count(t, p, start=0, end=0, m=0):
    """ Counts number of overlapping occurrences of pattern p in text t.

    :param t: String - Text (DNA) to look for
    :param p: String - Pattern (K-mer) to find in t
    :param start: Integer - Start in position <start> in the DNA
    :param end: Integer - End in position <end> in the DNA
    :param mismatches: Integer - Allow m mismatches
    :returns: Integer - n number of occurrences of p in t
    :raises: ValueError - If start < 0 or >= len(t)
    """
    if start < 0 or start >= len(t):
        raise ValueError('The starting position should be between 0 and the size ' + \
            'of the DNA')
    k = len(p)
    count = 0
    end = len(t) - k + 1 if end == 0 else end
    for i in range(0, end):
        if hamming_distance(t[i:i+k], p) <= m:
            count += 1
    return count


def frequent_kmers(DNA, k, m=0, reverse=False):
    """ Return a list of most frequent k-mers in DNA.

    :param DNA: String - DNA
    :param k: Integer - Length of the K-mer
    :param m: Allow for m mismatches
    :returns: Set - Set of most frequent K-mers in DNA
    """
    fd = dictionaries.FrequencyDict(DNA, k, m)

    freq = 0
    frequent = set()
    for kmer, frequency in fd.iteritems():
        rev = reverse_complement(kmer, as_string=True)
        if reverse and fd.has_key(rev):
            frequency += fd[rev]
        if frequency > freq:
            freq = frequency
            frequent = set([kmer])
        elif frequency == freq:
                frequent.add(kmer)
    return frequent


def reverse_complement(DNA, as_string=False):
    """ Returns the Reverse Complement of DNA
    """
    reverse_dic = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
    reverse_comp = [reverse_dic[nuc] for nuc in DNA[::-1]]
    if as_string:
        reverse_comp = ''.join(reverse_comp)
    return reverse_comp


def find_clumps(DNA, k, L, t):
    """ Find kmers forming clumps in DNA

    For a determined K-mer we say it forms an (L, t)-Clump, if K-mer appears at
    least t times in a region of length L in DNA

    :param DNA: String - Genome
    :param k: Integer - Length of the K-mer
    :param L: Integer - Size of the clump
    :param t: Integer - Number of times kmer must appear in the DNA region

    :return: List - K-mers forming (L, t)-Clumps in DNA
    """
    assert len(DNA) >= L

    clumps = set()
    # Construct the frequency dict for the first region of size L in the DNA
    fa = dictionaries.FrequencyDict(DNA[:L], k)

    # For each kmer in the first window, check if frequency >= t and correspondingly
    # add the kmer to the clumps set
    kmers = set()
    for i in range(0, L - k + 1):
        kmer = DNA[i:i+k]
        if not kmer in kmers:
            kmers.add(kmer)
            _t = fa[kmer]
            if _t >= t:
                clumps.add(kmer)

    # Decrease the frequency of the first kmer for the next iteration
    first_kmer = DNA[0:k]
    fa[first_kmer] -= 1

    # Go through all other regions of length L in the DNA
    for i in range(1, len(DNA)-L+1):
        # If not the first iteration, increase the frequency of the recently added
        # last kmer. If that frequency >= t, add the kmer to the set of clumps
        last_kmer = DNA[i+L-k:i+L]
        fa[last_kmer] += 1
        if fa[last_kmer] >= t:
            clumps.add(last_kmer)

        # Decrese the frequency of the first kmer in the region, as the sliding
        # window will remove it
        first_kmer = DNA[i:i+k]
        fa[first_kmer] -= 1
    return clumps


def skew(DNA, chart=False):
    """ Give all values of Skew(DNA)

    We define skew(DNA) as the difference between the total number of occurrences
    of G and the total number of occurrences of C in DNA.

    :param DNA: String - DNA to calculate skew
    :param chart: Boolean - If True, will save a skew.png chart in the current directory.

    :return: Tuple - With skew array and list of positions where the skew minimizes.
    """
    res = [0]
    G_C = 0
    _min = 0
    indexes = []
    for i, n in enumerate(DNA):
        if n == 'G':
            G_C += 1
        elif n == 'C':
            G_C -= 1
        # Compute the min. NOTE: We have to sum one to the indexes because we already
        # start with an extra element in the res (a 0)
        if G_C < _min:
            _min = G_C
            indexes = [i+1]
        elif G_C == _min:
            indexes.append(i+1)
        res.append(G_C)
    if chart:
        if sys.modules.get('matplotlib', None):
            plt.plot(res)
            plt.ylabel('G - C diff')
            plt.title('Skew diagram')
            plt.savefig('skew.png')
        else:
            print "No matplotlib module found, not creating skew diagram"
    return (res, indexes)


def hamming_distance(s1, s2):
    """ Computes the Hamming Distance between s1 and s2
    """
    assert len(s1) == len(s2)
    return sum([1 if c1 != c2 else 0 for c1, c2 in zip(s1, s2)])
