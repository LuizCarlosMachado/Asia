import arg_needle_lib
import sys

arg = arg_needle_lib.deserialize_arg(sys.argv[1])
arg.populate_children_and_roots()
newick = arg_needle_lib.arg_to_newick(arg)
print(newick)

