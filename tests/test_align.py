import pytest
from mea_poa.align import align_seqs
import mea_poa.pair_hmm  as ph
import math
import numpy as np
from pandas import DataFrame

def compare_matrices(query, correct):
    for i in range(query.shape[0]):
        for j in range(query.shape[1]):
            # print (i,j)
            # print (query[i][j])
            # print (correct[i][j])
            # print()
            if (query[i][j] == 0.0 or query[i][j] == float("-Inf")) and (correct[i][j] == 0.0 or correct[i][j] ==
                float(\
                    "-Inf")):
                # print ('it was')
                continue

            elif math.isclose(math.log(query[i][j]), math.log(correct[i][j]) ,rel_tol = 1e-4):
                continue
            else:
                return f'False - at {i} and {j} and {math.log(query[i][j])} is not the same as {math.log(correct[i][j])}'

    return True





def test_align_simple_blosum_62_2(simple_blosum_62_2):

    aligned_profile = simple_blosum_62_2.performViterbiAlignment()

    print (aligned_profile)

    assert "".join(aligned_profile.seqs[0].sequence) == "RTAG"
    assert "".join(aligned_profile.seqs[1].sequence) == "-TA-"


def test_align_simple_blosum_50_2(simple_blosum_50_2):

    aligned_profile = simple_blosum_50_2.performViterbiAlignment()

    print (aligned_profile)

    assert "".join(aligned_profile.seqs[0].sequence) == "RTAG"
    assert "".join(aligned_profile.seqs[1].sequence) == "-TA-"


def test_align_durbin_blosum_50_2(durbin_blosum_50_2):
    aligned_profile = durbin_blosum_50_2.performViterbiAlignment()

    print (aligned_profile)

    print (DataFrame(durbin_blosum_50_2.vM))
    print (DataFrame(durbin_blosum_50_2.vX))
    print (DataFrame(durbin_blosum_50_2.vY))

    assert "".join(aligned_profile.seqs[0].sequence) == "HEAG-AWGHE-E"
    assert "".join(aligned_profile.seqs[1].sequence) == "----PAW-HEAE"

def test_align_borodovsky_2(borodovsky_borodovsky_4_7_2):

    correct_vM = [[1, 0, 0, 0, 0 , 0],
           [0, 0.25, 0, 0, 0, 0],
           [0, 0, 0.0375, 0.005, 0.0000375, 0.0000003125],
           [0, 0, 0.0015, 0.0009375, 0.00075, 0.0001]]

    correct_vX = [[0, 0, 0, 0 , 0, 0],
                  [0, 0, 0.0125, 0.0003125, 0.000007813, 0.0000001953],
                  [0, 0, 0, 0.001875, 0.00025, 0.00000625],
                  [0, 0, 0, 0.000075, 0.00004688, 0.0000375]]

    correct_vY = [[0, 0, 0, 0, 0, 0],
                  [0, 0, 0, 0, 0, 0],
                  [0, 0.0125, 0, 0, 0, 0],
                  [0, 0.0003125, 0.001875, 0.00025, 0.000001875, 0.00000001563]]

    aligned_profile = borodovsky_borodovsky_4_7_2.performViterbiAlignment()

    print ()
    print (aligned_profile)
    print (DataFrame(borodovsky_borodovsky_4_7_2.vM.transpose()))
    print()
    print (DataFrame(borodovsky_borodovsky_4_7_2.vX.transpose()))
    print()
    print (DataFrame(borodovsky_borodovsky_4_7_2.vY.transpose()))
    print()

    assert "".join(aligned_profile.seqs[0].sequence) == "TTACG"
    assert "".join(aligned_profile.seqs[1].sequence) == "T-A-G"

    # print (DataFrame(correct_vM))
    # print()
    # print (DataFrame(np.log(correct_vM)))
    # print ()
    # print(DataFrame(borodovsky_borodovsky_4_7_2.vM.transpose()))



    assert compare_matrices(borodovsky_borodovsky_4_7_2.vM.transpose(), correct_vM) == True
    assert compare_matrices(borodovsky_borodovsky_4_7_2.vX.transpose(), correct_vY) == True
    assert compare_matrices(borodovsky_borodovsky_4_7_2.vY.transpose(), correct_vX) == True

def test_borodovsky_forward(borodovsky_borodovsky_4_7_2):

    correct_vM = [[1, 0, 0, 0, 0 , 0],
           [0, 0.25, 0.02, 0.0003, 0.00000125, 0.00000009375],
           [0, 0.012, 0.0375, 0.01, 0.00018, 0.000001944],
           [0, 0.00015, 0.0024, 0.0010015, 0.00196, 0.0002639]]

    correct_vX = [[0, 0.5, 0.00125, 0.00003125 , 0.0000007813, 0.000000001953],
                  [0, 0, 0.0125, 0.0013125, 0.00004781, 0.000001258],
                  [0, 0, 0.0006, 0.00189, 0.0005473, 0.00002268],
                  [0, 0, 0.0000075, 0.0001202, 0.00005308, 0.00009919]]

    correct_vY = [[0, 0, 0, 0, 0, 0],
                  [0.05, 0, 0, 0, 0, 0],
                  [0.00125, 0.0125, 0.001, 0.000015, 0.0000000625, 0.000000004688],
                  [0.00003125, 0.0009125, 0.00025, 0.0005004, 0.000009002, 0.00000009730]]

def test_borodovsky_backward(borodovsky_borodovsky_4_7_2):

    correct_vM = [[0, 0, 0, 0, 0 , 0],
           [0, 0.1, 0.005, 0.000125, 0.000003125, 0.00000007813],
           [0, 0.005, 0.25, 0.0275, 0.0001131, 0.000003234],
           [0, 0.000125, 0.00213, 0.00195, 0.000838, 0.00007574]]

    correct_vX = [[0, 0, 0, 0, 0 , 0],
           [0, 0.1, 0.0025, 0.0000625, 0.000001563, 0.0000000391],
           [0, 0, 0.04, 0.0022, 0.00006, 0.000001875],
           [0, 0, 0.0002, 0.00301, 0.001175, 0.00005653]]

    correct_vY = [[0, 0, 0, 0, 0 , 0],
           [0, 0.1, 0, 0, 0, 0],
           [0, 0.0025, 0.04, 0.0012, 0.000005, 0.000000375],
           [0, 0.0000625, 0.0012, 0.00303, 0.0011, 0.00002716]]