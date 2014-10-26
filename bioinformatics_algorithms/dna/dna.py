""" DNA related algorithms
"""
import Queue

from bioinformatics_algorithms.data_structures import arrays

def pattern_count(t, p, start=0, end=0):
    """ Counts number of overlapping occurrences of pattern p in text t.

    :param t: String - Text (DNA) to look for
    :param p: String - Pattern (K-mer) to find in t
    :param start: Integer - Start in position <start> in the DNA
    :param end: Integer - End in position <end> in the DNA
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
        if t[i:i+k] == p:
            count += 1
    return count


def frequent_kmers(DNA, k):
    """ Return a list of most frequent k-mers in DNA

    :param DNA: String - DNA
    :param k: Integer - Length of the K-mer
    :returns: Set - Set of most frequent K-mers in DNA
    """
    kmers = set()
    most_frequent = Queue.PriorityQueue()
    # Go through all the K-mers in DNA and, if the K-mer is found, frist check
    # if it has been found previously and, if not, count its occurences.
    for i in range(0, len(DNA) - k + 1):
        kmer = ''.join(DNA[i:i+k])
        if kmer in kmers:
            continue
        else:
            kmers.add(kmer)
            # Priority queue will return the lowest first, so we can negate it if we
            # want to get the higher priorite (frequence) first
            most_frequent.put((-pattern_count(DNA, kmer), kmer))

    # Extract most frequent K-mers
    result = set()
    first = True
    freq = -1
    if not most_frequent.empty():
        while not most_frequent.empty():
            kmer = most_frequent.get()
            if not first and kmer[0] != freq:
                break
            else:
                result.add(kmer[1])
                freq = kmer[0]
                first = False
    return result


def reverse_complement(DNA, as_string=False):
    """ Returns the Reverse Complement of DNA
    """
    reverse_dic = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
    reverse_comp = [reverse_dic[nuc] for nuc in DNA[::-1]]
    if as_string:
        reverse_comp = ''.join(reverse_comp)
    return reverse_comp


def find_kmer(kmer, DNA):
    """ Find all occurences of K-mer in DNA

    :param kmer: String - K-mer to find in DNA
    :param DNA: String - DNA

    :return: List - Starting positions of kmer in DNA
    """
    res = []
    k = len(kmer)
    for i in range(len(DNA) - k + 1):
        if DNA[i:i+k] == kmer:
            res.append(i)
    return res


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
    # Construct the frequency array for the first region of size L in the DNA
    fa = arrays.FrequencyArray(DNA[:L], k)

    # For each kmer in the first window, check if frequency >= t and correspondingly
    # add the kmer to the clumps set
    kmers = set()
    for i in range(0, L - k + 1):
        kmer = DNA[i:i+k]
        if not kmer in kmers:
            kmers.add(kmer)
            _t = fa.get_frequency(kmer)
            if _t >= t:
                clumps.add(kmer)

    # Decrease the frequency of the first kmer for the next iteration
    first_kmer = DNA[0:k]
    f = fa.get_frequency(first_kmer) - 1
    fa.set_frequency(first_kmer, f)

    # Go through all other regions of length L in the DNA
    for i in range(1, len(DNA)-L+1):
        region = DNA[i:i+L]
        # If not the first iteration, increase the frequency of the recently added
        # last kmer. If that frequency >= t, add the kmer to the set of clumps
        last_kmer = DNA[i+L-k:i+L]
        f = fa.get_frequency(last_kmer) + 1
        fa.set_frequency(last_kmer, f)
        if f >= t:
            clumps.add(last_kmer)

        # Decrese the frequency of the first kmer in the region, as the sliding
        # window will remove it
        first_kmer = DNA[i:i+k]
        f = fa.get_frequency(first_kmer) - 1
        fa.set_frequency(first_kmer, f)
        first_window = False
    return clumps


def skew(DNA, chart=False):
    """ Give all values of Skew(DNA)

    We define skew(DNA) as the difference between the total number of occurrences
    of G and the total number of occurrences of C in DNA.
    """
    res = [0]
    G_C = 0
    for i, n in enumerate(DNA):
        if n == 'G':
            G_C += 1
        elif n == 'C':
            G_C -= 1
        res.append(G_C)
    return res