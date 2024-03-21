import sys
import tskit

def get_tree(filepath: str, index: int) -> str:
    """
    Load a tree sequence from a file and return the Newick representation of a specific tree.

    Args:
        filepath (str): The path to the tree sequence file.
        index (int): The index of the tree to retrieve.

    Returns:
        str: The Newick representation of the specified tree.
    """
    ts = tskit.load(filepath)
    tree = ts.at(index)
    return tree.as_newick()

# Example usage:
#tree_filepath = sys.argv[1]
#tree_index = int(sys.argv[2])
#print(get_tree(tree_filepath, tree_index))

