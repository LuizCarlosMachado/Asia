import sys
import tskit

ts = tskit.load(sys.argv[1])
outfile = sys.argv[1]

ts.dump(outfile)
with open(outfile[:-3] + ".vcf", "w") as vcf_file:
    ts.write_vcf(vcf_file)

