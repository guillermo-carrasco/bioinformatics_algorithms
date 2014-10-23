""" Simple tests for the DNA related methods
"""
import unittest

from bioinformatics_algorithms.dna import dna

class DNATests(unittest.TestCase):
    """ Test all methods related with DNA manipulation
    """
    def test_1_pattern_count(self):
        """ Tests for pattern_count method...
        """
        DNA = 'GCGCG'
        PATTERN = 'GCG'
        with self.assertRaises(ValueError):
            dna.pattern_count(DNA, PATTERN, -1)
        with self.assertRaises(ValueError):
            dna.pattern_count(DNA, PATTERN, 10)
        self.assertEqual(dna.pattern_count(DNA, PATTERN), 2)

    def test_2_frequent_kmers(self):
        """ Tests for the frequen words (K-mers) method...
        """
        DNA = 'ACGTTGCATGTCGCATGATGCATGAGAGCT'
        K = 4
        self.assertEqual(dna.frequent_kmers(DNA, K), set(['CATG', 'GCAT']))

    def test_3_reverse_complement(self):
        """ Testing reverse complement method...
        """
        DNA = 'AAAACCCGGT'
        reverse = 'ACCGGGTTTT'
        self.assertEqual(dna.reverse_complement(DNA), list(reverse))
        self.assertEqual(dna.reverse_complement(DNA, as_string=True), reverse)

    def test_4_find_kmer(self):
        """ Testing finding k-mer in DNA method...
        """
        self.assertEqual(dna.find_kmer('ATAT', 'GATATATGCATATACTT'), [1, 3, 9])

    def test_5_find_clumps(self):
        """ Testing find clumps method...
        """
        DNA = 'CGGACTCGACAGATGTGAAGAACGACAATGTGAAGACTCGACACGACAGAGTGAAGAGAAGAGGAAACATTGTAA'
        k = 5
        L = 50
        t = 4
        res = {'CGACA', 'GAAGA'}
        self.assertEqual(dna.clumps(DNA, k, L, t), res)
