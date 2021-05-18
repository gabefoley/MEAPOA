import mea_poa.align as align
import mea_poa.parameters as parameters
import mea_poa.sub_matrix as sub_matrix
import itertools
import mea_poa.alignment_profile as aln_profile
import mea_poa.baum_welch as bw
import sequence
from sym import Alphabet

Protein_Alphabet_wB_X_Z = Alphabet('ABCDEFGHIKLMNPQRSTVWYXZ')
alignment = 'Not calculated'
po_alignment = 'Not calculated'
# seq = "../../tests/files/simple_seqs/borodovsky.fasta"

#
seq = "../../tests/files/custom_seqs/tree_check.fasta"
seqs = sequence.readFastaFile(seq, alphabet=Protein_Alphabet_wB_X_Z)

change_params = {'tau': 0.02, 'epsilon': 0.05, 'delta': 0.02, 'emissionX': 0.92, 'emissionY':
    0.2}

change_params = {'tau': 0.002, 'epsilon': 0.05, 'delta': 0.02, 'emissionX':
    0.5,
                 'emissionY':
                     0.5}




for seq_order in list(itertools.combinations(seqs, 2)):
    profiles = [aln_profile.AlignmentProfile([x]) for x in seq_order]


    # change_params = bw.runBaumWelch(change_params, profiles, "viterbi")


alignment = align.align_seqs(seq,
                             "../../tests/files/custom_seqs/3_col.aln",
                       aln_type='mea',
                       params=change_params, subsmat=sub_matrix.blosum62EstimatedWithX_dict,
                             log_transform=True)

po_alignment = align.align_seqs(seq,
                             "../../tests/files/custom_seqs/3_col.aln",
                       aln_type= 'pomea',
                       params=change_params, subsmat=sub_matrix.blosum62EstimatedWithX_dict,
                             log_transform=True)

# alignment = align.align_seqs("../../tests/files/custom_seqs/emission.fasta",
#                              "../../tests/files/custom_seqs/emission.aln",
#                        aln_type='viterbi',


#                        params=parameters.changed_params_again, subsmat=sub_matrix.blosum62EstimatedWithX_dict,
#                              log_transform=True)

# alignment = align.align_seqs("../../tests/files/custom_seqs/start_gap.fasta",
#                              "../../tests/files/custom_seqs/start_gap.aln",
#                        aln_type='mea',
#                        params=parameters.changed_params_again, subsmat=sub_matrix.blosum62EstimatedWithX_dict,
#                              log_transform=True)

# alignment = align.align_seqs("../../tests/files/custom_seqs/col_3.fasta",
#                              "../../tests/files/custom_seqs/col_3.aln",
#                        aln_type='viterbi',
#                        params=parameters.changed_params_again, subsmat=sub_matrix.blosum62LatestProbs,
#                              log_transform=True)

# alignment = align.align_seqs("../../tests/files/custom_seqs/monkey.fasta",
#                              "../../tests/files/custom_seqs/monkey.aln",
#                        aln_type='mea',
#                        params=parameters.changed_params_again, subsmat=sub_matrix.blosum62EstimatedWithX_dict,
#                              log_transform=True)



print('Final alignment')
print(alignment)
print ()
print (po_alignment)
