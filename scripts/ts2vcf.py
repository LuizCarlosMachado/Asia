import sys
import tskit

def convert_ts_to_vcf(input_file):
    """
    Converts a TreeSequence file to VCF format.

    Args:
        input_file (str): Path to the input TreeSequence file.

    Returns:
        None
    """
    ts = tskit.load(input_file)
    outfile = input_file

    ts.dump(outfile)
    with open(outfile[:-3] + ".vcf", "w") as vcf_file:
        ts.write_vcf(vcf_file)

