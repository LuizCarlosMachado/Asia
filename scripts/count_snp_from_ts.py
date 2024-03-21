import tskit
import sys
import os

def count_snps_from_ts(ts_file):
	"""
	Counts the number of SNPs (Single Nucleotide Polymorphisms) in a Tree Sequence file.

	Args:
		ts_file (str): The path to the Tree Sequence file.

	Returns:
		None
	"""
	ts = tskit.load(ts_file)
	snps = ts.num_mutations
	#print(snps)

	if snps == 0:
		os.remove(ts_file)



