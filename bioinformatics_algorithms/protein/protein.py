from string import maketrans

from bioinformatics_algorithms.data_structures import dictionaries
from bioinformatics_algorithms.dna import dna
from bioinformatics_algorithms.protein import RNA_TO_AMINO, AMINIO_TO_RNA


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