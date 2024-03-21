import numpy as np
import gzip
import tskit
import arg_needle_lib
import sys

def add_mutation_to_tree_sequence():
    """
    Adds a mutation to a tree sequence based on the given inputs.

    Args:
        sys.argv[1] (str): Path to the needle tree sequence file.
        sys.argv[2] (str): Path to the true tree sequence file.
        sys.argv[3] (str): Path to the haplotype file.
        sys.argv[4] (int): Position of the mutation.
        sys.argv[5] (str): Path to save the mutated tree sequence.

    Returns:
        None
    """
    needle_ts = tskit.load(sys.argv[1])  # .trees
    true_ts = tskit.load(sys.argv[2])  # .ts
    hap = sys.argv[3] # .haps.gz
    POS = int(sys.argv[4]) #500000

    mut_time = true_ts.mutations_time[true_ts.mutations_site == true_ts.site(position=POS).id][0]

    with gzip.open(hap, "rt") as hapF:
        for line in hapF:
            fields = line.strip().split()
            pos = int(fields[2])
            if pos == POS:
                gt = list(map(int, fields[5:]))
                break

    der_leaves = np.nonzero(gt)[0]
    center_tree = needle_ts.at(POS)
    mut_node = center_tree.mrca(*der_leaves)
    max_mut_time = center_tree.time(center_tree.parent(mut_node))-1
    min_mut_time = center_tree.time(mut_node)
    tree_time = np.maximum(np.minimum(mut_time, max_mut_time), min_mut_time)
    tables = needle_ts.dump_tables()
    site_id = tables.sites.add_row(position=POS, ancestral_state="0")
    tables.mutations.add_row(site=site_id, node=mut_node, derived_state="1", time=tree_time)

    mut_ts = tables.tree_sequence()
    mut_ts.dump(sys.argv[5])
