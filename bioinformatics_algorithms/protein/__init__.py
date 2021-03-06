""" Module for protein related methods
"""

RNA_TO_AMINO={"UUU":"F", "UUC":"F", "UUA":"L", "UUG":"L",
              "UCU":"S", "UCC":"S", "UCA":"S", "UCG":"S",
              "UAU":"Y", "UAC":"Y", "UAA":"", "UAG":"",
              "UGU":"C", "UGC":"C", "UGA":"", "UGG":"W",
              "CUU":"L", "CUC":"L", "CUA":"L", "CUG":"L",
              "CCU":"P", "CCC":"P", "CCA":"P", "CCG":"P",
              "CAU":"H", "CAC":"H", "CAA":"Q", "CAG":"Q",
              "CGU":"R", "CGC":"R", "CGA":"R", "CGG":"R",
              "AUU":"I", "AUC":"I", "AUA":"I", "AUG":"M",
              "ACU":"T", "ACC":"T", "ACA":"T", "ACG":"T",
              "AAU":"N", "AAC":"N", "AAA":"K", "AAG":"K",
              "AGU":"S", "AGC":"S", "AGA":"R", "AGG":"R",
              "GUU":"V", "GUC":"V", "GUA":"V", "GUG":"V",
              "GCU":"A", "GCC":"A", "GCA":"A", "GCG":"A",
              "GAU":"D", "GAC":"D", "GAA":"E", "GAG":"E",
              "GGU":"G", "GGC":"G", "GGA":"G", "GGG":"G"}

AMINIO_TO_RNA = {'M': ['AUG'],
                 'H': ['CAU', 'CAC'],
                 'Q': ['CAA', 'CAG'],
                 'P': ['CCU', 'CCG', 'CCC', 'CCA'],
                 'R': ['CGU', 'CGG', 'CGC', 'CGA', 'AGA', 'AGG'],
                 'L': ['CUU', 'CUG', 'CUC', 'CUA', 'UUG', 'UUA'],
                 'D': ['GAU', 'GAC'],
                 'E': ['GAA', 'GAG'],
                 'A': ['GCA', 'GCC', 'GCU', 'GCG'],
                 'G': ['GGA', 'GGC', 'GGU', 'GGG'],
                 'V': ['GUA', 'GUC', 'GUU', 'GUG'],
                 'Y': ['UAC', 'UAU'],
                 '*': ['UAG', 'UAA', 'UGA'],
                 'S': ['UCA', 'UCC', 'UCU', 'UCG', 'AGU', 'AGC'],
                 'C': ['UGC', 'UGU'],
                 'W': ['UGG'],
                 'F': ['UUC', 'UUU'],
                 'N': ['AAC', 'AAU'],
                 'K': ['AAA', 'AAG'],
                 'T': ['ACA', 'ACC', 'ACG', 'ACU'],
                 'I': ['AUU', 'AUC', 'AUA']
                 }

MASS_TABLE = {'G': 57, 'A': 71, 'S': 87, 'P': 97, 'V': 99, 'T': 101, 'C': 103, 
              'I': 113, 'L': 113, 'N': 114, 'D': 115, 'K': 128, 'Q': 128, 'E': 129,
              'M': 131, 'H': 137, 'F': 147, 'R': 156, 'Y': 163, 'W': 186}