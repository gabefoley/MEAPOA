from Bio.SubsMat import MatrixInfo, SeqMat
from Bio.Alphabet import NucleotideAlphabet
blosum62 = MatrixInfo.blosum62
blosum50 = MatrixInfo.blosum50


def score_match(pair, matrix):
    if pair not in matrix:
        return matrix[(tuple(reversed(pair)))]
    else:
        return matrix[pair]


# borodovsky_data = {('A', 'A') : 0.5, ('A', 'C') : 0.15, ('A', 'G') : 0.05, ('A', 'T') : 0.3,
#                    ('C', 'A'): 0.15, ('C', 'C'): 0.5, ('C', 'G'): 0.3, ('C', 'T'): 0.05,
#                    ('G', 'A'): 0.05, ('G', 'C'): 0.3, ('G', 'G'): 0.5, ('G', 'T'): 0.15,
#                    ('T', 'A') : 0.3, ('T', 'C') : 0.05, ('T', 'G') : 0.15, ('T', 'T') : 0.5,}

borodovsky_data = {('A', 'A') : 0.5, ('A', 'C') : 0.15, ('A', 'G') : 0.05, ('A', 'T') : 0.3,
                   ('C', 'C'): 0.5, ('C', 'G'): 0.3, ('C', 'T'): 0.05,
                   ('G', 'G'): 0.5, ('G', 'T'): 0.15,
                   ('T', 'T') : 0.5,}

borodovsky_4_7 = SeqMat(borodovsky_data, mat_name='borodovsky_4_7')

