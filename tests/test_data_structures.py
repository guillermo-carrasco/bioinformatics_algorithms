""" Tests for Data Structures
"""
import unittest

from bioinformatics_algorithms.data_structures import arrays


class TestArrayDataStructures(unittest.TestCase):
    """ Test array-based data structures
    """
    def test_1_frequency_array(self):
        """ Testing frequency array...
        """
        DNA = 'ACGCGGCTCTGAAA'
        K = 2
        fa = arrays.FrequencyArray(DNA, K)
        self.assertEqual(fa.get_frequency_array(), [2, 1, 0, 0, 0, 0, 2, 2, 1, 2, 1, 0, 0, 1, 1, 0])
