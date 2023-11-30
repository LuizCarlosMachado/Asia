import sys
import tskit

ts = tskit.load(sys.argv[1])
tree = ts.at(int(sys.argv[2]))
print(tree.as_newick())
