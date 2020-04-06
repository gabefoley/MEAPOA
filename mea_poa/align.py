import sequence
from sym import Alphabet
import mea_poa.alignment_profile as aln_profile
import mea_poa.guide_tree as gt
import mea_poa.pair_hmm as ph
import mea_poa.sub_matrix as sub_matrix
import mea_poa.parameters as parameters


Protein_Alphabet_wB_X_Z = Alphabet('ABCDEFGHIKLMNPQRSTVWYXZ')



# Align sequences
def align_seqs(inpath, outpath, aln_type, params=parameters.basic_params, subsmat=sub_matrix.blosum62EstimatedWithX_dict,
               log_transform=True):

    # Read sequences in
    seqs = sequence.readFastaFile(inpath, alphabet=Protein_Alphabet_wB_X_Z)

    # print (len(seqs))

    # Calculate guide tree
    guide_tree = gt.get_guide_tree(seqs)
    # print (guide_tree.ascii_art())

    # Get the alignment order
    aln_order = gt.get_aln_order(guide_tree)
    # print (aln_order)

    seq_dict = {x.name : x for x in seqs}

    # Create alignment in order from guide tree
    for node in aln_order:
        curr_node = node[0]
        if type(seq_dict[node[1][0]]) == aln_profile.AlignmentProfile:
            profile1 = seq_dict[node[1][0]]
        else:
            profile1 = aln_profile.AlignmentProfile([seq_dict[node[1][0]]])

        if type(seq_dict[node[1][1]]) == aln_profile.AlignmentProfile:
            profile2 = seq_dict[node[1][1]]
        else:
            profile2 = aln_profile.AlignmentProfile([seq_dict[node[1][1]]])

        seqs = [profile1, profile2]

        pair_hmm = load_params(params, seqs, subsmat, log_transform)

        if aln_type == 'viterbi':

            pair_hmm.performViterbiAlignment()
            aligned_profile = pair_hmm.get_alignment(type_to_get='viterbi')

        elif aln_type == 'mea':

            pair_hmm.performMEAAlignment()
            aligned_profile = pair_hmm.get_alignment(type_to_get='mea')

        seq_dict[curr_node] = aligned_profile

        print (aligned_profile)




    with open(outpath, 'w') as outfile:
        outfile.write(str(aligned_profile))


    return aligned_profile

def load_params(params, seqs, subsmat, log_transform):

    pair_hmm = ph.PairHMM(seqs, params['tau'], params['epsilon'], params['delta'], params['emissionX'],
                      params['emissionY'], subsmat, log_transform)
    return pair_hmm

# alignment = align_seqs("../tests/files/simple_seqs/bananas_5.fasta", "../tests/files/simple_seqs/bananas_5.aln",
#                        aln_type='viterbi',
#                        params=parameters.basic_params, log_transform=True)

# alignment = align_seqs("../tests/files/simple_seqs/mea_test3.fasta", "../tests/files/simple_seqs/mea_test3.aln",
#                        aln_type='mea',
#                        params=parameters.test_params2, log_transform=False)
#
#
#
# print('Final alignment')
# print(alignment)

# alignment = align_seqs("../tests/files/simple_seqs/mea_test.fasta", "../tests/files/simple_seqs/mea_test.aln",
#                        aln_type='mea',
#                        params=parameters.test_params2, log_transform=True)
#
#
#
# print('Final alignment')
# print(alignment)

# alignment = align_seqs("../tests/files/simple_seqs/mea_test.fasta", "../tests/files/simple_seqs/mea_test.aln",
#                        aln_type='viterbi',
#                        params=parameters.basic_params, log_transform=True)
#
#
#
# print('Final alignment')
# print(alignment)

# alignment = align_seqs("../tests/files/simple_seqs/mea_test2.fasta", "../tests/files/simple_seqs/mea_test2.aln",
#                        aln_type='mea',
#                        params=parameters.test_params2, subsmat=sub_matrix.blosum62LatestProbs, log_transform=False)
#
#
#
# print('Final alignment')
# print(alignment)
#
#
#
# alignment = align_seqs("../tests/files/simple_seqs/mea_test2.fasta", "../tests/files/simple_seqs/mea_test2.aln",
#                        aln_type='mea',
#                        params=parameters.test_params2, subsmat=sub_matrix.blosum62LatestProbs, log_transform=True)
#
#
#
# print('Final alignment')
# print(alignment)

# alignment = align_seqs("../tests/files/simple_seqs/simple_4.fasta", "../tests/files/simple_seqs/simple_4.aln",
#                        aln_type='mea',
#                        params=parameters.changed_params, subsmat=sub_matrix.blosum62LatestProbs, log_transform=True)
#
#
#
# print('Final alignment')
# print(alignment)

# alignment = align_seqs("../tests/files/simple_seqs/simple_8.fasta", "../tests/files/simple_seqs/simple_8.aln",
#                        aln_type='mea',
#                        params=parameters.basic_params, subsmat=sub_matrix.blosum62LatestProbs, log_transform=False)
#
#
#
# print('Final alignment')
# print(alignment)
#
# print ('R and R')
# print (sub_matrix.score_match(('R', 'R'), sub_matrix.blosum62EstimatedWithX_dict))

# print ('R and A')
#
# print (sub_matrix.score_match(('R', 'A'), sub_matrix.blosum62LatestProbs))
# print ('R and S')
#
# print (sub_matrix.score_match(('R', 'S'), sub_matrix.blosum62LatestProbs))
# print ('R and G')
#
# print (sub_matrix.score_match(('R', 'G'), sub_matrix.blosum62LatestProbs))
# print ('N and S')
# print (sub_matrix.score_match(('N', 'S'), sub_matrix.blosum62LatestProbs))
#
# print ('S and S')
# print (sub_matrix.score_match(('S', 'S'), sub_matrix.blosum62LatestProbs))
#
# print ('S and D')
# print (sub_matrix.score_match(('S', 'D'), sub_matrix.blosum62LatestProbs))