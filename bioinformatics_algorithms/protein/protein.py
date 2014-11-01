from bioinformatics_algorithms.protein import RNA_TO_AMINO

def rna_to_amino(rna):
    """ Translate an RNA string into an aminoacids string. Stop codons won't be translated.
    """
    for i in range(0, len(rna), 3):
        yield RNA_TO_AMINO[rna[i:i+3]]
