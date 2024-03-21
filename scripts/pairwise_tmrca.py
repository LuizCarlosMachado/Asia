import sys
import ete3
from ete3 import Tree
from random import shuffle
from scipy.stats import pearsonr
from scipy.stats import spearmanr

def pairwise_tmrca(t1, t2, n_comps):
    """
    Calculate the Time to Most Recent Common Ancestor (TMRCA) for pairs of leaves in two phylogenetic trees.

    Args:
        t1 (str): The file path of the first phylogenetic tree.
        t2 (str): The file path of the second phylogenetic tree.
        n_comps (int): The number of leaf pairs to calculate TMRCA for.

    Returns:
        tuple: A tuple containing the Pearson correlation coefficient and p-value of the TMRCA values.

    """
    t1 = Tree(t1)
    t2 = Tree(t2)

    all_leaves = []
    for leaf in t1:
        all_leaves.append(leaf.name)
    all_pairs = [(a, b) for idx, a in enumerate(all_leaves) for b in all_leaves[idx + 1:]]
    shuffle(all_pairs)
    all_pairs = all_pairs[0:n_comps]

    a_list = []
    b_list = []
    for pair in all_pairs:
        leaf1 = pair[0]
        leaf2 = pair[1]
        A1 = t1&leaf1
        B1 = t1&leaf2
        tmrca1 = A1.get_distance(B1) / 2
        a_list.append(tmrca1)
        A2 = t2&leaf1
        B2 = t2&leaf2
        tmrca2 = A2.get_distance(B2) / 2
        b_list.append(tmrca2)
        print("pair", sys.argv[1], sys.argv[2], leaf1, leaf2, tmrca1, tmrca2)

    correlation = pearsonr(a_list, b_list)
    return correlation


# Example usage:
correlation = pairwise_tmrca(sys.argv[1], sys.argv[2], int(sys.argv[3]))
print("stats", sys.argv[1], sys.argv[2], correlation[0], correlation[1])
