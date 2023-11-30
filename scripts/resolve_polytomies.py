import sys
from ete3 import Tree

def remove_outer(tree):
    tree = tree[1:-2]
    tree = tree + ";"
    return tree

t = Tree(sys.argv[1])
t.resolve_polytomy(recursive=True)
#t_out = remove_outer(t.write())
t_out = t.write()
print(t_out)
