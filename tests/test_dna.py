""" Simple tests for the DNA related methods
"""
import unittest

from bioinformatics_algorithms import dna

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

