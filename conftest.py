import pytest
import mea_poa.pair_hmm as ph
import mea_poa.alignment_profile as aln_profile
import mea_poa.sub_matrix as sub_matrix
import mea_poa.align as align
import mea_poa.parameters as params
import sequence

@pytest.fixture
def simple_blosum_62_2():
    seqs = sequence.readFastaFile("./files/simple_seqs/simple_2.fasta")

    profile1 = aln_profile.AlignmentProfile([seqs[0]])
    profile2 = aln_profile.AlignmentProfile([seqs[1]])

    phmm = align.load_params(params.basic_params, [profile1, profile2], sub_matrix.blosum62, log_transform=True)

    return phmm

@pytest.fixture
def probabilities_blosum_62_2():
    seqs = sequence.readFastaFile("./files/simple_seqs/simple_2.fasta")

    profile1 = aln_profile.AlignmentProfile([seqs[0]])
    profile2 = aln_profile.AlignmentProfile([seqs[1]])

    phmm = align.load_params(params.basic_params, [profile1, profile2], sub_matrix.blosum62LatestProbs, log_transform=True)

    return phmm

@pytest.fixture
def simple_blosum_50_2():
    seqs = sequence.readFastaFile("./files/simple_seqs/simple_2.fasta")

    profile1 = aln_profile.AlignmentProfile([seqs[0]])
    profile2 = aln_profile.AlignmentProfile([seqs[1]])

    phmm = align.load_params(params.basic_params, [profile1, profile2], sub_matrix.blosum50, log_transform=True)

    return phmm

@pytest.fixture
def durbin_blosum_50_2():

    seqs = sequence.readFastaFile("./files/simple_seqs/durbin_2.fasta")
    profile1 = aln_profile.AlignmentProfile([seqs[0]])
    profile2 = aln_profile.AlignmentProfile([seqs[1]])

    phmm = align.load_params(params.basic_params, [profile1, profile2], sub_matrix.blosum62, log_transform=False)

    return phmm

@pytest.fixture
def borodovsky_borodovsky_4_7_2():

    seqs = sequence.readFastaFile("./files/simple_seqs/borodovsky.fasta")
    profile1 = aln_profile.AlignmentProfile([seqs[0]])
    profile2 = aln_profile.AlignmentProfile([seqs[1]])

    phmm = align.load_params(params.borodovsky_4_7, [profile1, profile2], sub_matrix.borodovsky_4_7,
                             log_transform=False)

    return phmm


@pytest.fixture
def borodovsky_blosum_50_2():

    seqs = sequence.readFastaFile("./files/simple_seqs/borodovsky.fasta")
    profile1 = aln_profile.AlignmentProfile([seqs[0]])
    profile2 = aln_profile.AlignmentProfile([seqs[1]])

    phmm = align.load_params(params.borodovsky_4_7, [profile1, profile2], sub_matrix.blosum62LatestProbs,
                             log_transform=False)

    return phmm