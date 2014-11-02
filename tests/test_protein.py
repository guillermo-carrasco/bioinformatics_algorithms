import unittest

from bioinformatics_algorithms.protein import protein

class TestProteinMethods(unittest.TestCase):
    """ Tests for the protein submodule
    """
    def test_1_rna_to_amino(self):
        """ Testing rna_to_amino method...
        """
        rna = 'AUGGCCAUGGCGCCCAGAACUGAGAUCAAUAGUACCCGUAUUAACGGGUGA'
        res = 'MAMAPRTEINSTRING'
        self.assertEqual(res, ''.join(protein.rna_to_amino(rna)))


    def test_2_encoded_bt(self):
        """ Testing encoded_by method...
        """
        dna = 'ATGGCCATGGCCCCCAGAACTGAGATCAATAGTACCCGTATTAACGGGTGA'
        amino = 'MA'
        res = ['ATGGCC', 'ATGGCC', 'GGCCAT']
        self.assertEqual(res, protein.encoded_by(dna, amino))