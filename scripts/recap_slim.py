import sys
import tskit
import numpy as np
from numpy import random
import pyslim
import msprime
import warnings

warnings.simplefilter('ignore', msprime.TimeUnitsMismatchWarning)

x = random.randint(1e5)
r=1e-8
chrlen = 1e6
N = 1e4
mu = 1e-8
taxa = int(sys.argv[2])

ts = tskit.load(sys.argv[1])
ts = pyslim.update(ts)

#outfile = sys.argv[1] + ".slim_recap_id" + str(x) + "_samples" + str(taxa) + "_r" + str(r) + "_len" + str(int(chrlen)) + "_N" + str(int(N)) + "_mu" + str(mu) + ".ts"
outfile = sys.argv[1][:-3]+ "_samples" + str(taxa) + ".ts"

s = pyslim.recapitate(ts,recombination_rate=r,ancestral_Ne=N)

# Subsample
rng = np.random.default_rng()
alive_inds = pyslim.individuals_alive_at(s, 0)
keep_indivs = rng.choice(alive_inds, taxa, replace=False)
keep_nodes = []
for i in keep_indivs:
	    keep_nodes.extend(s.individual(i).nodes)

# Compute mutation times if these are missing
tables = s.dump_tables()
tables.compute_mutation_times()
tables.sort()
timed_s = tables.tree_sequence()

# Simplify to reduced sample set
sts = timed_s.simplify(keep_nodes, keep_input_roots=True)
# Add neutral mutations
rts = msprime.sim_mutations(sts, rate=mu, keep=True,model="jc69")
rts.dump(outfile)
with open(outfile[:-3] + ".vcf", "w") as vcf_file:
    rts.write_vcf(vcf_file)







