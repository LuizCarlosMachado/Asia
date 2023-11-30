import tskit
import sys
import os

ts = tskit.load(sys.argv[1])
snps = ts.num_mutations
#print(snps)

if snps == 0:
	os.remove(sys.argv[1])



