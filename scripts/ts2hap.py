import tskit
import numpy as np
import sys

true_ts = tskit.load(sys.argv[1])
gtm = true_ts.genotype_matrix()
pos = true_ts.sites_position
hap_file = sys.argv[2]

with open(hap_file, 'w') as hapF:
    for i in range(pos.shape[0]):
        if np.all(np.isin(gtm[i], [0, 1])):
            print(str(1), "SNP"+str(i+1), int(pos[i]), "A", "T", *gtm[i], sep = " ", file = hapF)

