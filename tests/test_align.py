import pytest
from mea_poa.align import align_seqs
import mea_poa.pair_hmm  as ph
import math
import numpy as np
from pandas import DataFrame

def compare_matrices(query, correct, check_log=True):
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



            elif check_log and math.isclose(math.log(query[i][j]), math.log(correct[i][j]) ,rel_tol = 1e-4):
                continue

            elif not check_log and math.isclose(query[i][j], correct[i][j] ,rel_tol = 1e-2):
                continue
            else:
                if check_log:
                    return f'False - at {i} and {j} and {math.log(query[i][j])} is not the same as {math.log(correct[i][j])}'
                else:
                    return f'False - at {i} and {j} and {query[i][j]} is not the same as {correct[i][j]}'

    return True





def test_align_simple_blosum_62_2(simple_blosum_62_2):

    simple_blosum_62_2.performViterbiAlignment()

    aligned_profile = simple_blosum_62_2.get_alignment('viterbi')

    print (aligned_profile)

    assert "".join(aligned_profile.seqs[0].sequence) == "RTAG"
    assert "".join(aligned_profile.seqs[1].sequence) == "-TA-"


def test_align_simple_blosum_50_2(simple_blosum_50_2):

    simple_blosum_50_2.performViterbiAlignment()

    aligned_profile = simple_blosum_50_2.get_alignment('viterbi')

    print (aligned_profile)

    assert "".join(aligned_profile.seqs[0].sequence) == "RTAG"
    assert "".join(aligned_profile.seqs[1].sequence) == "-TA-"


def test_align_durbin_blosum_50_2(durbin_blosum_50_2):
    durbin_blosum_50_2.performViterbiAlignment()
    aligned_profile = durbin_blosum_50_2.get_alignment('viterbi')

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

    borodovsky_borodovsky_4_7_2.performViterbiAlignment()

    aligned_profile = borodovsky_borodovsky_4_7_2.get_alignment('viterbi')

    print ()
    print (aligned_profile)
    print (DataFrame(borodovsky_borodovsky_4_7_2.vM))
    print()
    print (DataFrame(borodovsky_borodovsky_4_7_2.vX))
    print()
    print (DataFrame(borodovsky_borodovsky_4_7_2.vY))
    print()

    assert "".join(aligned_profile.seqs[0].sequence) == "T-A-G"
    assert "".join(aligned_profile.seqs[1].sequence) == "TTACG"

    # print (DataFrame(correct_vM))
    # print()
    # print (DataFrame(np.log(correct_vM)))
    # print ()
    # print(DataFrame(borodovsky_borodovsky_4_7_2.vM.transpose()))



    assert compare_matrices(borodovsky_borodovsky_4_7_2.vM, correct_vM, False) == True
    assert compare_matrices(borodovsky_borodovsky_4_7_2.vX, correct_vX, False) == True
    assert compare_matrices(borodovsky_borodovsky_4_7_2.vY, correct_vY, False) == True

def test_borodovsky_forward(borodovsky_borodovsky_4_7_2):

    correct_fM = [[1, 0, 0, 0, 0 , 0],
           [0, 0.25, 0.02, 0.0003, 0.00000125, 0.00000009375],
           [0, 0.012, 0.0375, 0.01, 0.00018, 0.000001944],
           [0, 0.00015, 0.0024, 0.0010015, 0.00196, 0.0002639]]

    correct_fX = [[0, 0.05, 0.00125, 0.00003125 , 0.0000007813, 0.00000001953],
                  [0, 0, 0.0125, 0.0013125, 0.00004781, 0.000001258],
                  [0, 0, 0.0006, 0.00189, 0.0005473, 0.00002268],
                  [0, 0, 0.0000075, 0.0001202, 0.00005308, 0.00009919]]

    correct_fY = [[0, 0, 0, 0, 0, 0],
                  [0.05, 0, 0, 0, 0, 0],
                  [0.00125, 0.0125, 0.001, 0.000015, 0.0000000625, 0.000000004688],
                  [0.00003125, 0.0009125, 0.0019, 0.0005004, 0.000009002, 0.00000009730]]

    aligned_profile = borodovsky_borodovsky_4_7_2.performMEAAlignment()

    print()
    print ("FINAL RESULTS")
    print ()
    print (aligned_profile)
    print (DataFrame(borodovsky_borodovsky_4_7_2.fM))
    print()
    print (DataFrame(borodovsky_borodovsky_4_7_2.fX))
    print()
    print (DataFrame(borodovsky_borodovsky_4_7_2.fY))
    print()


    assert math.isclose(borodovsky_borodovsky_4_7_2.full_prob, 0.00003632, rel_tol=1e-4) == True

    assert compare_matrices(borodovsky_borodovsky_4_7_2.fM, correct_fM, False) == True
    assert compare_matrices(borodovsky_borodovsky_4_7_2.fX, correct_fX, False) == True
    assert compare_matrices(borodovsky_borodovsky_4_7_2.fY, correct_fY, False) == True


def test_borodovsky_backward(borodovsky_borodovsky_4_7_2):

    correct_bM = [[0, 0, 0, 0, 0, 0, 0],
                  [0, 0.00007574, 0.000838, 0.00195, 0.00213, 0.000125, 0],
                  [0, 0.000003234, 0.0001131, 0.00275, 0.025, 0.005, 0],
                  [0, 0.00000007813, 0.000003125, 0.000125, 0.005, 0.1, 0],
                  [0, 0, 0, 0, 0, 0, 0]]

    correct_bX = [[0, 0, 0, 0, 0, 0, 0],
                  [0, 0.00005653, 0.001175, 0.00301, 0.0002, 0, 0],
                  [0, 0.000001875, 0.00006, 0.0022, 0.04, 0, 0],
                  [0, 0.0000000391, 0.000001563, 0.0000625, 0.0025, 0.1, 0],
                  [0, 0, 0, 0, 0, 0, 0]]

    correct_bY = [[0, 0, 0, 0, 0, 0, 0],
                  [0, 0.00002716, 0.0011, 0.00303, 0.0012, 0.0000625, 0],
                  [0, 0.000000375, 0.000005, 0.0012, 0.04, 0.0025, 0],
                  [0, 0, 0, 0, 0, 0.1, 0],
                  [0, 0, 0, 0, 0, 0, 0]]


    aligned_profile = borodovsky_borodovsky_4_7_2.performMEAAlignment()

    print()
    print ("FINAL RESULTS")
    print ()
    print (aligned_profile)
    print (DataFrame(borodovsky_borodovsky_4_7_2.bM))
    print()
    print (DataFrame(borodovsky_borodovsky_4_7_2.bX))
    print()
    print (DataFrame(borodovsky_borodovsky_4_7_2.bY))
    print()

    assert compare_matrices(borodovsky_borodovsky_4_7_2.bM, correct_bM, False) == True
    assert compare_matrices(borodovsky_borodovsky_4_7_2.bX, correct_bX, False) == True
    assert compare_matrices(borodovsky_borodovsky_4_7_2.bY, correct_bY, False) == True

def test_posterior_matrix(borodovsky_borodovsky_4_7_2):

    aligned_profile = borodovsky_borodovsky_4_7_2.calc_posterior_for_viterbi()


    assert borodovsky_borodovsky_4_7_2.aligned_positions == [(1,1), (2,3), (3,5)]

    assert math.isclose(borodovsky_borodovsky_4_7_2.pM[1][1], 0.521, rel_tol=1e-3)
    assert math.isclose(borodovsky_borodovsky_4_7_2.pM[2][3], 0.757, rel_tol=1e-3)
    assert math.isclose(borodovsky_borodovsky_4_7_2.pM[3][5], 0.727, rel_tol=1e-3)
    print()
    pM = borodovsky_borodovsky_4_7_2.pM

    print ('Lets do it')

    seq1_matches = []
    seq2_matches = []

    seq1_idx = 0
    seq2_idx = 0
    last_state = ""

    i = 3
    j = 5

    while i > 0 and j > 0:
        print (i,j)

        print (DataFrame(pM))

        print (pM[i-1][j-1])
        print (pM[i][j])
        print ('**')
        print (pM[i-1][j-1] * pM[i][j] )
        print (pM[i-1][j])
        print (pM[i][j-1])





    print (seq1_matches)
    print (seq2_matches)





    # self.aligned_positions = self.get_aligned_positions(seq1_matches, seq2_matches)
    #
    # if type == "viterbi":
    #     self.viterbi_matches1 = seq1_matches
    #     self.viterbi_matches2 = seq2_matches
    #
    # elif type == 'mea':
    #     self.mea_matches1 = seq1_matches
    #     self.mea_matches2 = seq2_matches
    #
    #
    #     # alignment1.add_gaps(seq1_matches)
    #     # alignment2.add_gaps(seq2_matches)
    #
    #     # print (self.profile1.profile)
    #     # print (self.profile2.profile)
    #
    #     # self.profile1.add_profile(self.profile2)
    #
    #     # return self.profile1







