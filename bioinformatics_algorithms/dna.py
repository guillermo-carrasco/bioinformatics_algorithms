""" DNA related algorithms
"""
import Queue

def pattern_count(t, p, start=0):
    """ Counts number of overlapping occurrences of pattern p in text t.

    :param t: String - Text (DNA) to look for
    :param p: String - Pattern (K-mer) to find in t
    :param start: Integer - Start in position <start> in the DNA
    :returns: Integer - n number of occurrences of p in t
    :raises: ValueError - If start < 0 or >= len(t)
    """
    if start < 0 or start >= len(t):
        raise ValueError('The starting position should be between 0 and the size ' + \
            'of the DNA')
    k = len(p)
    count = 0
    for i in range(0, len(t) - k + 1):
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

