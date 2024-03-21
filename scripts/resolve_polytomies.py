import sys
from ete3 import Tree

def remove_outer(tree):
    """
    Removes the outer brackets from a Newick tree string.

    Args:
        tree (str): The Newick tree string.

    Returns:
        str: The Newick tree string without the outer brackets.
    """
    tree = tree[1:-2]
    tree = tree + ";"
    return tree

def resolve_polytomies(tree_file):
    """
    Resolves polytomies in a Newick tree.

    Args:
        tree_file (str): The path to the file containing the Newick tree.

    Returns:
        str: The resolved Newick tree string.
    """
    t = Tree(tree_file)
    t.resolve_polytomy(recursive=True)
    t_out = t.write()
    return t_out

if __name__ == "__main__":
    tree_file = sys.argv[1]
    resolved_tree = resolve_polytomies(tree_file)
    print(resolved_tree)
