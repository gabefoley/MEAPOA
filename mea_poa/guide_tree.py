from skbio import DistanceMatrix
from skbio.tree import nj
from alfpy.utils import distmatrix
from alfpy.utils import seqrecords
from alfpy import ncd
import numpy as np

def calc_random_distances(seqs):
    """
    Calculate the distances between a set of sequences randomly
    :param seqs: Sequences to calculate for
    :return: Distance matrix of all distances between sequences
    """

    # Lets just create with a random distance matrix
    a = np.random.random_integers(0, 40, size=(len(seqs), len(seqs)))

    # Make it symmetric and make the diagonal equal zero
    a_symm = (a + a.T)/2
    np.fill_diagonal(a_symm, 0)

    return a_symm

def calc_distances(seqs):
    seq_records = seqrecords.SeqRecords()

    for seq in seqs:
        seq_records.add(seq.name, "".join(seq.sequence))


    dist = ncd.Distance(seq_records)
    matrix = distmatrix.create(seq_records.id_list, dist)

    return matrix.data


def get_aln_order(tree):
    """
    Given a guide tree, get the order that we should be aligning the sequences
    :param tree:
    :return: Dictionary mapping tree node : children_names in the order they should be aligned
    """
    aln_order = []
    for node in tree.postorder():
        if not node.is_tip():
            aln_order.append((node.name, [child.name.replace(" ", "_") for child in node.children]))
            # aln_order.append({node.name : [child.name for child in node.children]})
    return aln_order

def get_guide_tree(seqs):
    """
    Get a guide tree representing distances between sequences
    :param seqs: Sequences to create a tree for
    :return: Guide tree
    """

    # Get distances and ids
    distances = calc_distances(seqs)
    ids = [x.name for x in seqs]

    # distances = [[ 0,  16,  22,  26.5],
    #              [16,   0,  25.5, 24.5],
    #              [22,  25.5,  0,  22.5],
    #              [26.5, 24.5, 22.5,  0. ]]
    #

    # Make a distance matrix and Neighbour-Joining tree
    dm = DistanceMatrix(distances, ids)
    tree = nj(dm)


    # Mid-point root and then label the internal nodes
    tree = tree.root_at_midpoint()
    label_internal_nodes(tree)


    return tree

def label_internal_nodes(tree):
    """
    Label the internal nodes of a tree
    :param tree: Tree to label
    """
    count = 0
    for node in tree.preorder():
        if not node.is_tip():
            node.name = "N" + str(count)
            count +=1





