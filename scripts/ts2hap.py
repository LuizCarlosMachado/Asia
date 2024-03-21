import tskit
import numpy as np
import sys

def convert_ts_to_hap(true_ts, hap_file):
    """
    Convert a Tree Sequence (tskit object) to a haplotype file.

    Args:
        true_ts (tskit.TreeSequence): The input Tree Sequence object.
        hap_file (str): The path to the output haplotype file.

    Returns:
        None
    """
    gtm = true_ts.genotype_matrix()
    pos = true_ts.sites_position

    with open(hap_file, 'w') as hapF:
        for i in range(pos.shape[0]):
            if np.all(np.isin(gtm[i], [0, 1])):
                print(str(1), "SNP"+str(i+1), int(pos[i]), "A", "T", *gtm[i], sep=" ", file=hapF)

