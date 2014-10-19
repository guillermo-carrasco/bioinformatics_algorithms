""" DNA related algorithms
"""

def pattern_count(t, p):
    """ Counts number of overlapping occurrences of pattern p in text t.

    :param t: String - Text (DNA) to look for
    :param p: String - Pattern (K-mer) to find in t
    :returns: Integer - n number of occurrences of p in t
    """
    k = len(p)
    count = 0
    for i in range(0, len(t) - k + 1):
        if t[i:i+k] == p:
            count += 1
    return count
