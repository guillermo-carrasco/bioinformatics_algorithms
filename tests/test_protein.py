import unittest

from bioinformatics_algorithms.protein import protein

class TestProteinMethods(unittest.TestCase):
    """ Tests for the protein submodule
    """
    def test_1_rna_to_amino(self):
        """ Testing rna_to_amino method...
        """
        rna = 'AUGGCCAUGGCGCCCAGAACUGAGAUCAAUAGUACCCGUAUUAACGGGUGA'
        expected = 'MAMAPRTEINSTRING'
        res = []
        for a in protein.rna_to_amino(rna):
            res.append(a)
        self.assertEqual(expected, ''.join(res))
