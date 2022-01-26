import os, sys

import orthoani
import argparse
from Bio.SeqIO import read





def main(args):

    outfile = args.output
    anchor_gfile = args.anchor_genome
    compare_genomes_pathlist = args.compare_genome_pathlist
    compare_genomes = [x.strip() for x in open(compare_genomes_pathlist, "r")]
    path_name = args.path

    with open(outfile, 'w') as out:
        header = ["comparison_name", "path", "anchor_name", "compare_name", "orthoani"]
        out.write(",".join(header) + "\n")
        anchor_g = read(anchor_gfile, "fasta")
        import pdb;pdb.set_trace()
        anchor_name = os.path.basename(anchor_gfile).rsplit('_genomic', 1)[0]
        for cg in compare_genomes:
            compare_g = read(cg, "fasta")
            compare_name = os.path.basename(cg).rsplit('_genomic', 1)[0]
            comparison_name = f"{anchor_name}_x_{compare_name}"

            ani = orthoani.orthoani(anchor_g, compare_g)
            ani = str(ani)

            out.write(f"{comparison_name},{path_name},{anchor_name},{compare_name},{ani}\n")


def cmdline(sys_args):
    "Command line entry point w/argparse action."
    p = argparse.ArgumentParser()
    p.add_argument("--anchor-genome", required=True)
    p.add_argument("--compare-genome-pathlist", required=True)
    p.add_argument("--path")
    p.add_argument("--output", required=True)
    args = p.parse_args()
    return main(args)

if __name__ == '__main__':
    returncode = cmdline(sys.argv[1:])
    sys.exit(returncode)

