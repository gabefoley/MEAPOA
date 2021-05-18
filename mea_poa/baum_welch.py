import mea_poa.align as align
import mea_poa.sub_matrix as sub_matrix
import mea_poa.guide_tree as gt
from sym import Alphabet
import sequence
import mea_poa.alignment_profile as aln_profile
import itertools
import copy

Protein_Alphabet_wB_X_Z = Alphabet('ABCDEFGHIKLMNPQRSTVWYXZ')


params = {'tau' : 0.0002, 'epsilon' : 0.072, 'delta' : 0.08, 'emissionX': 0.2, 'emissionY': 0.2}
change_params = {'tau' : 0.0002, 'epsilon' : 0.072, 'delta' : 0.08, 'emissionX': 0.2, 'emissionY': 0.2}

predecessors = {}
def runBaumWelch(parameters, profiles, aln_type, count=0, stop=2):

    if count == stop:
        return parameters

    pair_hmm = align.load_params(parameters, profiles, subsmat=sub_matrix.blosum62EstimatedWithX_dict,
                                 log_transform=True, predecessors=predecessors)

    print ('Delta value is ')
    print (pair_hmm.delta)

    print ('Epsilon value is')
    print (pair_hmm.epsilon)

    print ("Emission X value is ")
    print (pair_hmm.emissionX)

    print ("Emission Y value is ")
    print (pair_hmm.emissionY)

    if aln_type == 'viterbi':

        pair_hmm.performViterbiAlignment()

        copy_hmm = copy.deepcopy(pair_hmm)

        aligned_profile = copy_hmm.get_alignment(type_to_get='viterbi')

        print (aligned_profile)

        # Update parameters
        parameters['epsilon'] = float(pair_hmm.epsilon + 0.00001)
        parameters['delta'] = float(pair_hmm.delta + 0.00001)
        parameters['emissionX'] = float(pair_hmm.emissionX + 0.00001)
        parameters['emissionY'] = float(pair_hmm.emissionY + 0.00001)

        return runBaumWelch(parameters, profiles, aln_type, count +1)

    elif aln_type == 'mea':

        pair_hmm.performMEAAlignment()

        copy_hmm = copy.deepcopy(pair_hmm)

        aligned_profile = copy_hmm.get_alignment(type_to_get='mea')

        print (aligned_profile)

        # Update parameters
        parameters['epsilon'] = float(pair_hmm.epsilon + 0.00001)
        parameters['delta'] = float(pair_hmm.delta + 0.00001)
        parameters['emissionX'] = float(pair_hmm.emissionX + 0.00001)
        parameters['emissionY'] = float(pair_hmm.emissionY + 0.00001)

        return runBaumWelch(parameters, profiles, aln_type, count +1)


            # seqs = sequence.readFastaFile(sequences, alphabet=Protein_Alphabet_wB_X_Z)
            #
            # aln_order = [("N0", [seqs[0].name, seqs[1].name])]
            #
            # seqs = get_profiles(seqs, aln_order)
            #
            # # Reload parameters
            #
            # pair_hmm = align.load_params(parameters, seqs, subsmat=sub_matrix.blosum62EstimatedWithX_dict,
            #                              log_transform=True)

# seq_file = '../tests/files/custom_seqs/2_col.fasta'
# output_file = '../tests/files/custom_seqs/col_3.aln'
# aln_type = 'mea'
# seqs = sequence.readFastaFile(seq_file, alphabet=Protein_Alphabet_wB_X_Z)
#
#
# for seq_order in list(itertools.combinations(seqs, 2)):
#
#     profiles = [aln_profile.AlignmentProfile([x]) for x in seq_order]
#     print (seq_order)
#
#     change_params = runBaumWelch(change_params, profiles, aln_type)
#
#
# print (params)
# print (change_params)

# alignment = align.align_seqs(seq_file,
#                              output_file, aln_type,
#                        params=params, subsmat=sub_matrix.blosum62EstimatedWithX_dict,log_transform=True)
#
# change_alignment = align.align_seqs(seq_file,
#                              output_file, aln_type,
#                        params=change_params, subsmat=sub_matrix.blosum62EstimatedWithX_dict,log_transform=True)