""" Tests for Data Structures
"""
import unittest

from bioinformatics_algorithms.data_structures import arrays


class TestArrayDataStructures(unittest.TestCase):
    """ Test array-based data structures
    """
    def test_1_frequency_array(self):
        """ Testing frequency array process...
        """
        fa = arrays.FrequencyArray('ACGTTGCATGTCGCATGATGCATGAGAGCT', 4)
        with open('test_frequency_array.txt', 'r') as f:
            freq_array = [int(x) for x in f.readline().split()]

        self.assertEqual(fa.get_frequency_array(), freq_array)
        self.assertEqual(fa.get_most_frequent(), set(['CATG', 'GCAT']))
