import sys
import os
import subprocess
import tskit
import tsinfer
import tsdate

ts = tskit.load(sys.argv[1])
stem = os.path.splitext(sys.argv[1])[0]

#out_file = stem + "_tsinfer.samples"
#tsinfer.SampleData.from_tree_sequence(ts, path=out_file, num_flush_threads=2)
#with tsinfer.SampleData(path=out_file, sequence_length=ts.sequence_length, num_flush_threads=4) as sample_data:
#    for var in ts.variants():
#        sample_data.add_site(var.site.position, var.genotypes, var.alleles)

out_file = stem + "_tsinfer.samples"
with tsinfer.SampleData(path=out_file, sequence_length=ts.sequence_length, num_flush_threads=4) as sample_data:
    pop = sample_data.add_population()
    for ind in range(ts.num_individuals):
        sample_data.add_individual(ploidy=2, population=pop)
    for var in ts.variants():
        sample_data.add_site(var.site.position, var.genotypes, var.alleles)

subprocess.run(["tsinfer", "infer", out_file, "-p","-t", "4"])
ts = tskit.load(stem + "_tsinfer.trees")
ts = ts.simplify(keep_unary=False)
dated_ts = tsdate.date(ts, Ne=10000, mutation_rate=1e-8)
dated_ts.dump(stem + "_tsinfer_tsdate.trees")


os.remove(stem + "_tsinfer.samples")
