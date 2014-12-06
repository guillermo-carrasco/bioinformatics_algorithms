""" Methods and classes to manage inputs
"""

def read_fasta(fasta_file):
    """ Return a generator to read a FASTA file

    :param str fasta_file: Path to the FASTA file
    """
    try:
        with open(fasta_file, 'r') as f:
            l = f.readline()
            while l:
                l = l.rstrip().lstrip()
                yield l
                l = f.readline()
    except IOError as e:
        e.message = "{} not found or don't have correct permissions :-(".format(fasta_file)
        raise e


def generic_input(file_path, types):
    """ Reads a generic input file

    :param str file_path: Path to the input file
    _param list types: Type of each line to expect in the input file, i.e int, str, etc.
    """
    try:
        with open(file_path, 'r') as f:
            for t in types:
                l = f.readline().lstrip().rstrip()
                if t == 'str':
                    yield l
                elif t == 'int':
                    yield int(l)
                elif t == 'float':
                    yield float(l)
                elif t == 'int_list':
                    yield [int(element) for element in l.split()]
                elif t == 'float_list':
                    yield [float(element) for element in l.split()]
                else:
                    raise NotImplementedError('This kind of input data is not contemplated, sorry')
    except IOError as e:
        e.message = "{} not found or don't have correct permissions :-(".format(file_path)
        raise e
