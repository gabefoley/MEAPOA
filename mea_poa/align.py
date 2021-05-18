import sequence
from sym import Alphabet
import mea_poa.alignment_profile as aln_profile
import mea_poa.guide_tree as gt
import mea_poa.pair_hmm as ph
import mea_poa.sub_matrix as sub_matrix
import mea_poa.parameters as parameters
import sequence

Protein_Alphabet_wB_X_Z = Alphabet('ABCDEFGHIKLMNPQRSTVWYXZ')


# Align sequences
def align_seqs(inpath, outpath, aln_type, params=parameters.basic_params,
               subsmat=sub_matrix.blosum62EstimatedWithX_dict,
               log_transform=True):

    print("params are")
    print(params)

    # Read sequences in
    seqs = sequence.readFastaFile(inpath, alphabet=Protein_Alphabet_wB_X_Z)

    print(len(seqs))

    if len(seqs) == 2:
        aln_order = [("N0", [seqs[0].name, seqs[1].name])]

    else:

        # Calculate guide tree
        guide_tree = gt.get_guide_tree(seqs, random=False)
        print(guide_tree.ascii_art())

        # Get the alignment order
        aln_order = gt.get_aln_order(guide_tree)
        # print (aln_order)

    print(aln_order)

    seq_dict = {x.name: x for x in seqs}

    # Predecessors start off blank
    predecessors = [{}, {}]

    # Create alignment in order from guide tree
    for node in aln_order:

        # Get the current node name and list of sequences under that node
        curr_node = node[0]
        curr_seqs = node[1]

        # List to store the aligned sequences in
        aligned = []

        # While the node has sequences underneath yet to be aligned
        while curr_seqs:

            # Get a sequence
            seq = curr_seqs.pop()

            # Make it into a profile if it isn't one already
            if type(seq_dict[seq]) != aln_profile.AlignmentProfile:
                profile = aln_profile.AlignmentProfile([seq_dict[seq]])
            else:
                profile = seq_dict[seq]



            # Add sequence to the aligned list
            aligned.append(profile)

            # if len(alns) > 1:
            #     new_align = "-align-".join(alns)
            #     alns = []
            #     alns.append(new_align)

            # If we have two profiles it is time to align
            if len(aligned) > 1:

                pair_hmm = load_params(params, aligned, subsmat, log_transform, predecessors)

                if aln_type == 'viterbi':

                    pair_hmm.performViterbiAlignment(po=False)
                    aligned_profile = pair_hmm.get_alignment(type_to_get='viterbi')

                elif aln_type == 'poviterbi':

                    pair_hmm.performViterbiAlignment(po=True)
                    aligned_profile = pair_hmm.get_alignment(type_to_get='viterbi')

                elif aln_type == 'mea':

                    pair_hmm.performMEAAlignment(po=False)
                    aligned_profile = pair_hmm.get_alignment(type_to_get='mea')

                elif aln_type == 'pomea':

                    pair_hmm.performMEAAlignment(po=True)
                    aligned_profile = pair_hmm.get_alignment(type_to_get='mea')


                # Clear the previous unaligned sequences
                aligned = []

                # Add the aligned sequences
                aligned.append(aligned_profile)

        # print ('wowza')
        # print (aligned[0])
        # print(aligned[0].predecessors)

        seq_dict[curr_node] = aligned[0]

        # print('alignment is ')
        # print(aligned_profile)

    with open(outpath, 'w') as outfile:
        outfile.write(str(aligned_profile))

    return aligned_profile


# Align sequences
# def align_seqs(inpath, outpath, aln_type, params=parameters.basic_params, subsmat=sub_matrix.blosum62EstimatedWithX_dict,
#                log_transform=True):
#
#     # Read sequences in
#     seqs = sequence.readFastaFile(inpath, alphabet=Protein_Alphabet_wB_X_Z)
#
#     # print (len(seqs))
#
#     # Calculate guide tree
#     guide_tree = gt.get_guide_tree(seqs)
#     print (guide_tree.ascii_art())
#
#     # Get the alignment order
#     aln_order = gt.get_aln_dict(guide_tree)
#     print (aln_order)
#
#     seq_dict = {x.name : x for x in seqs}
#
#     curr_seqs = {x.name : x for x in seqs}
#
#     finished = False
#
#
#
#     while not finished:
#
#
#
#         smllst_dist_aln_order = gt.get_smallest_distance_aln_order(curr_seqs.values(), aln_order)
#
#         ancestor = guide_tree.lowest_common_ancestor(smllst_dist_aln_order)
#
#         print(guide_tree.ascii_art())
#
#         print ('ancestor is')
#         print (ancestor.name)
#
#
#
#         # seq_dict[node[1][0]]
#
#
#         # Create alignment in order from guide tree
#
#         curr_node = guide_tree.lowest_common_ancestor(smllst_dist_aln_order)
#
#
#         child1 = aln_order[curr_node.name][0]
#         child2 = aln_order[curr_node.name][1]
#
#         print ('children are ')
#         print (child1)
#         print (child2)
#
#         if type(seq_dict[child1]) == aln_profile.AlignmentProfile:
#             profile1 = seq_dict[child1]
#         else:
#             profile1 = aln_profile.AlignmentProfile([seq_dict[child1]])
#
#         if type(seq_dict[child2]) == aln_profile.AlignmentProfile:
#             profile2 = seq_dict[child2]
#         else:
#             profile2 = aln_profile.AlignmentProfile([seq_dict[child2]])
#
#         profiles = [profile1, profile2]
#
#         pair_hmm = load_params(params, profiles, subsmat, log_transform)
#
#         if aln_type == 'viterbi':
#
#             pair_hmm.performViterbiAlignment()
#             aligned_profile = pair_hmm.get_alignment(type_to_get='viterbi')
#
#         elif aln_type == 'mea':
#
#             pair_hmm.performMEAAlignment()
#             aligned_profile = pair_hmm.get_alignment(type_to_get='mea')
#
#         seq_dict[curr_node] = aligned_profile
#
#         curr_seqs.pop(child1)
#         curr_seqs.pop(child2)
#
#         print ('alignment is ')
#         print (aligned_profile)
#
#         if ancestor.name == "N0":
#             print ('finishing')
#             finished = True
#
#
#
#
#     with open(outpath, 'w') as outfile:
#         outfile.write(str(aligned_profile))
#
#     print (aligned_profile)
#     return aligned_profile

def load_params(params, seqs, subsmat, log_transform, predecessors):
    pair_hmm = ph.PairHMM(seqs, params['tau'], params['epsilon'], params['delta'], params['emissionX'],
                          params['emissionY'], subsmat, log_transform)
    return pair_hmm

# alignment = align_seqs("../tests/files/simple_seqs/bananas_5.fasta", "../tests/files/simple_seqs/bananas_5.aln",
#                        aln_type='viterbi',
#                        params=parameters.basic_params, log_transform=True)

# alignment = align_seqs("../tests/files/simple_seqs/mea_test5.fasta", "../tests/files/simple_seqs/mea_test5.aln",
#                        aln_type='mea',
#                        params=parameters.test_params3, log_transform=True)
#
#
#
# print('Final alignment')
# print(alignment)

# alignment = align_seqs("../tests/files/simple_seqs/mea_test3.fasta", "../tests/files/simple_seqs/mea_test3.aln",
#                        aln_type='mea',
#                        params=parameters.test_params3, log_transform=False)
#
#
#
# print('Final alignment')
# print(alignment)


# alignment = align_seqs("../tests/files/simple_seqs/mea_test5.fasta", "../tests/files/simple_seqs/mea_test5.aln",
#                        aln_type='viterbi',
#                        params=parameters.test_params3, log_transform=True)
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

# alignment = align_seqs("../tests/files/custom_seqs/monkey.fasta", "../tests/files/custom_seqs/monkey.aln",
#                        aln_type='mea',
#                        params=parameters.changed_params, subsmat=sub_matrix.blosum62LatestProbs, log_transform=True)
#
#
#
# print('Final alignment')
# print(alignment)

# alignment = align_seqs("../tests/files/custom_seqs/col_3.fasta", "../tests/files/custom_seqs/col_3.aln",
#                        aln_type='mea',
#                        params=parameters.changed_params, subsmat=sub_matrix.blosum62LatestProbs, log_transform=True)
#
#
#
# print('Final alignment')
# print(alignment)
#
#
#
# print (sub_matrix.blosum62LatestProbs[('S','S')])
#
#
# print (sub_matrix.blosum62LatestProbs[('A','S')])
# print (sub_matrix.blosum62LatestProbs[('A','A')])


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