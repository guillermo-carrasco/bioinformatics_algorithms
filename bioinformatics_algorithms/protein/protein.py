from collections import Counter
from itertools import permutations
from math import factorial
from string import maketrans

from bioinformatics_algorithms.data_structures import dictionaries
from bioinformatics_algorithms.dna import dna
from bioinformatics_algorithms.protein import RNA_TO_AMINO, AMINIO_TO_RNA, MASS_TABLE
from bioinformatics_algorithms.utils import factorial


def rna_to_amino(rna):
    """ Translate an RNA string into an aminoacids string. Stop codons won't be translated.
    """
    translated = []
    for i in range(0, len(rna), 3):
        translated.append(RNA_TO_AMINO[rna[i:i+3]])
    return translated


def dna_to_rna(dna):
    """ Transcribe a DNA string into an RNA string
    """
    transcribe = maketrans('T', 'U')
    return dna.translate(transcribe)


def _encoded_by_rec(encodings, amino_seq):
    """ Generate all RNA sequences that encode a given aminoacid sequence
    """
    if not encodings:
        [encodings.add(t) for t in AMINIO_TO_RNA[amino_seq[0]]]
        encodings = _encoded_by_rec(encodings, amino_seq[1:])
    elif amino_seq:
        aux = set()
        for t in AMINIO_TO_RNA[amino_seq[0]]:
            for e in encodings:
                e += t
                aux.add(e)
        encodings = _encoded_by_rec(aux, amino_seq[1:])
    return encodings


def encoded_by(DNA, amino_seq):
    """ Find substrings of a genome DNA encoding a given amino acid sequence
    """
    encodings = _encoded_by_rec(set(), amino_seq)
    fd = dictionaries.FrequencyDict(DNA, len(amino_seq*3))
    res = []
    for e in encodings:
        e_dna = e.translate(maketrans('U', 'T'))
        e_rev = dna.reverse_complement(e_dna, as_string=True)
        f = fd.get(e_dna, 0)
        res.extend([e_dna]*f)
        f = fd.get(e_rev, 0)
        res.extend([e_rev]*f)
    return res


def count_subpeptides(n):
    """ Return the number of sub-peptides of a cyclic peptide of length n.

    As it is cyclic, it can always form n groups of length 1, 2, 3, ..., n-1. I.e
    the peptide NQEL can form: N, Q, E, L, NQ, QE, EL, LN, NQE, QEL, ELN, LNQ == 12

    So basically... n*n-1
    """
    return n*(n-1)


def spectrum(peptide, cyclic=False):
    """ Returns the mass spectrum of a peptide 

    :param bool cyclic: Indicates that the peptide is a cyclic peptide
    """
    # Calculate mass prefixes
    prefix_mass = [0]
    for index, base in enumerate(peptide):
        prefix_mass.append(MASS_TABLE[base] + prefix_mass[index])
    peptide_mass = prefix_mass[-1]
    # And now calculate the linear spectrum using the prefix mass
    spectrum = [0]
    for i in range(len(peptide)):
        for j in range(i+1, len(peptide) + 1):
            spectrum.append(prefix_mass[j] - prefix_mass[i])
            if cyclic and (i > 0 and j < len(peptide)):
                spectrum.append(peptide_mass - (prefix_mass[j] - prefix_mass[i]))
    return sorted(spectrum)


def score(peptide, theoretical_spectrum):
    """ Compute the score of a cyclic peptide against a spectrum
    """
    _spectrum = spectrum(peptide, cyclic=True)
    i = 0
    j = 0
    score = 0
    while i < len(theoretical_spectrum) and j < len(_spectrum):
        if theoretical_spectrum[i] == _spectrum[j]:
            score += 1
            i += 1
            j += 1
        elif theoretical_spectrum[i] > _spectrum[j]:
            j += 1
        else:
            i += 1
    return score


def peptides_with_mass(mass, recursive=False):
    """ Compute the number of peptides of given mass.


    The order DOES matter, so we need to save the correct results and sum to the 
    solution the total number of permutations of each result.
    """
    amin_table = ['G', 'A', 'S', 'P', 'V', 'T', 'C', 'L', 'N', 'D', 'Q', 'E', 'M', 'H', 'F', 'R', 'Y', 'W']
    mass_table = [57, 71, 87, 97, 99, 101, 103, 113, 114, 115, 128, 129, 131, 137, 147, 156, 163, 186]
    return _peptides_with_mass_rec(mass_table, amin_table, mass, '')


def _peptides_with_mass_rec(mass_table, amin_table, mass, peptide):
    """ Recursive version of peptides_with_mass
    """
    if len(mass_table) == 0 and mass > 0:
        return 0
    elif mass < 0:
        return 0
    elif mass == 0:
        # We have to return all the possible permutations of the branch. Each branch
        # will have a unique set of weights, so there is no possibility of repeated
        # set of weights across branches
        c = Counter(peptide)
        div = 1
        for v in c.values():
            div *= factorial(v)
        return factorial(len(peptide))/div
    else:
        return _peptides_with_mass_rec(mass_table, amin_table, mass - mass_table[0], peptide + amin_table[0]) + \
               _peptides_with_mass_rec(mass_table[1:], amin_table[1:], mass, peptide)

