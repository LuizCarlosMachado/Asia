import arg_needle_lib
import sys

def needleARG2newick():
    """
    Converts a Needle ARG object to Newick format and prints the result.

    Usage: python needleARG2newick.py <serialized_arg>

    Args:
        <serialized_arg> (str): The serialized Needle ARG object.

    Returns:
        None
    """
    arg = arg_needle_lib.deserialize_arg(sys.argv[1])
    arg.populate_children_and_roots()
    newick = arg_needle_lib.arg_to_newick(arg)
    print(newick)

needleARG2newick()

