import sys
import msprime
import tskit
from numpy import random

def simulate_msprime_neutral(taxa):
    """
    Simulates neutral genetic variation using msprime.

    Parameters:
    - taxa (int): Number of samples to simulate.

    Returns:
    - None

    Saves the simulated genetic variation to a file with a name based on the simulation parameters.
    """

    x = random.randint(1e5)
    r = 1e-8
    chrlen = 1e7
    N = 1e4
    mu = 1e-8
    outfile = "msprime_id" + str(x) + "_samples" + str(taxa) + "_r" + str(r) + "_len" + str(int(chrlen)) + "_N" + str(int(N)) + "_mu" + str(mu) + ".ts"

    ts = msprime.sim_ancestry(samples=taxa, recombination_rate=r, ploidy=2, sequence_length=chrlen, population_size=N)
    mts = msprime.sim_mutations(ts, rate=mu, model="jc69")
    mts.dump(outfile)
    #print(mts)
