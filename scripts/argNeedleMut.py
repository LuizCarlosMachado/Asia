import numpy as np
import gzip
import tskit
import arg_needle_lib
import sys

#needle_ts = tskit.load(sys.argv[1])
#true_ts = tskit.load(sys.argv[1][:-11] + ".ts")
#hap = sys.argv[1][:-11] + ".haps.gz"
#POS = int(sys.argv[2])

needle_ts = tskit.load(sys.argv[1])  # .trees
true_ts = tskit.load(sys.argv[2])  # .ts
hap = sys.argv[3] # .haps.gz
POS = int(sys.argv[4]) #500000




mut_time = true_ts.mutations_time[true_ts.mutations_site == true_ts.site(position=POS).id][0]

#with open() as hapF:
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
#tree_time = np.minimum(mut_time, max_mut_time)
tables = needle_ts.dump_tables()
site_id = tables.sites.add_row(position=POS, ancestral_state="0")
tables.mutations.add_row(site=site_id, node=mut_node, derived_state="1", time=tree_time)
#tables.mutations.add_row(site=site_id, node=mut_node, derived_state="1", time=mut_time)


mut_ts = tables.tree_sequence()
mut_ts.dump(sys.argv[5])
#mut_ts.dump(sys.argv[1] + "_mut.ts")
