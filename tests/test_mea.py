import math
import numpy as np
from pandas import DataFrame

def test_mea_align_borodovsky_2(borodovsky_borodovsky_4_7_2):

    borodovsky_borodovsky_4_7_2.performMEAAlignment()

    aligned_profile = borodovsky_borodovsky_4_7_2.get_alignment('mea')

    print ()
    print (aligned_profile)
    print (DataFrame(borodovsky_borodovsky_4_7_2.pM))

    def test_mea_align_borodovsky_2(borodovsky_borodovsky_4_7_2):
        borodovsky_borodovsky_4_7_2.performMEAAlignment()

        aligned_profile = borodovsky_borodovsky_4_7_2.get_alignment('mea')

        print()
        print(aligned_profile)
        print(DataFrame(borodovsky_borodovsky_4_7_2.pM))



    assert "".join(aligned_profile.seqs[0].sequence) == "T-A-G"
    assert "".join(aligned_profile.seqs[1].sequence) == "TTACG"

def test_mea_align_borodovsky_blosum_50_2(borodovsky_blosum_50_2):
    borodovsky_blosum_50_2.performMEAAlignment()

    aligned_profile = borodovsky_blosum_50_2.get_alignment('mea')

    print ()
    print (aligned_profile)
    print (DataFrame(borodovsky_blosum_50_2.pM))

    def test_mea_align_borodovsky_2(borodovsky_borodovsky_4_7_2):
        borodovsky_borodovsky_4_7_2.performMEAAlignment()

        aligned_profile = borodovsky_borodovsky_4_7_2.get_alignment('mea')

        print()
        print(aligned_profile)
        print(DataFrame(borodovsky_borodovsky_4_7_2.pM))



    assert "".join(aligned_profile.seqs[0].sequence) == "T-A-G"
    assert "".join(aligned_profile.seqs[1].sequence) == "TTACG"