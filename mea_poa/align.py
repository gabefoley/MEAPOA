import sequence
import mea_poa.alignment_profile as aln_profile
import mea_poa.guide_tree as gt
import mea_poa.pair_hmm as ph
import mea_poa.sub_matrix as sub_matrix
import mea_poa.parameters as parameters

# Align sequences
def align_seqs(inpath, outpath, params=parameters.basic_params, subsmat=sub_matrix.blosum62, log_transform=True):

    # Read sequences in
    seqs = sequence.readFastaFile(inpath)

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
        # print ('curr node is ' + str(curr_node))
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

        aligned_profile = pair_hmm.performViterbiAlignment()
        # print(aligned_profile)

        seq_dict[curr_node] = aligned_profile




    with open(outpath, 'w') as outfile:
        outfile.write(str(aligned_profile))


    return aligned_profile

def load_params(params, seqs, subsmat, log_transform):

    pair_hmm = ph.PairHMM(seqs, params['tau'], params['epsilon'], params['delta'], params['emissionX'],
                      params['emissionY'], subsmat, log_transform)
    return pair_hmm

# alignment = align_seqs("../tests/files/simple_seqs/bananas_5.fasta", "../tests/files/simple_seqs/bananas_5.aln",
#                        params=parameters.basic_params, log_transform=True)
#
# print('Final alignment')
# print(alignment)
#
# print (sub_matrix.score_match(('N', 'D'), sub_matrix.blosum62))
#
# print (sub_matrix.score_match(('P', 'V'), sub_matrix.blosum62))