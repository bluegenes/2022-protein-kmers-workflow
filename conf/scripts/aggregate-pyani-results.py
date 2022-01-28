import os
import sys
import argparse
import glob
import pprint

import numpy as np
import pandas as pd

from collections import defaultdict, namedtuple

anchorpyani = namedtuple('anchorpyani',
                           'comparison_name, anchor_name, compare_name, path, lowest_common_rank, pyani_ident, pyani_coverage, pyani_aln_length, pyani_sim_errors, pyani_hadamard')

compare_ranklist = ["genus", "family", "order", "class", "order", "class", "phylum", "superkingdom"]

def main(args):
    # get basename for these sequences
    anchor_results = []
    pathInfo = pd.read_csv(args.path_info, dtype=str, sep="\t", header=0)
    if args.path_name:
        select_path = args.path_name
        pathInfo = pathInfo[pathInfo["paths"] == select_path] # subset
    allpaths = pathInfo['paths'].unique().tolist()

    #argh, passing in diff classes/labels isn't working. just convert here
    acc2label = {}
    if args.labels:
        with open(args.labels, 'r') as lb:
            for line in lb:
                md5,fn,label=line
                acc = fn.rsplit('_')[0]
                acc2label[acc] = label
    else:
        acc2label[acc] = acc

    print(acc2label)

    for path in allpaths:
        anchor_acc = pathInfo.loc[(pathInfo["path"] == path) & (pathInfo["rank"] == "anchor")]["accession"].values[0]
        anchor_label = acc2label[anchor_acc]
        compare_accs = pathInfo.loc[(pathInfo["path"] == path) & (pathInfo["rank"] != "anchor")]["accession"].values.tolist()

        path_results_dir = os.path.join(args.pyani_results_dir, path, results)
        lenF = os.path.join(path_results_dir, "matrix_aln_lengths_2.tab")
        covF = os.path.join(path_results_dir,"matrix_coverage_2.tab")
        hadamardF = os.path.join(path_results_dir,"matrix_hadamard_2.tab")
        identF = os.path.join(path_results_dir,"matrix_identity_2.tab")
        simerrF = os.path.join(path_results_dir,"matrix_sim_errors_2.tab")

        #pathfiles = {"pyani_ident": identF, "pyani_coverage": covF, "pyani_aln_length": lenF, "pyani_sim_errors": simerrF, "pyani_hadamard": hadamardF}
        # read in all matrices
        lenD = pd.read_csv(lenF, sep="\t", header=0, index_col=0)
        covD = pd.read_csv(covF, sep="\t", header=0, index_col=0)
        hadD = pd.read_csv(hadamardF, sep="\t", header=0, index_col=0)
        idD = pd.read_csv(identF, sep="\t", header=0, index_col=0)
        seD = pd.read_csv(simerrF, sep="\t", header=0, index_col=0)
        # get info for each comparison
        for rank in compare_ranklist:
            compare_acc = pathInfo.loc[(pathInfo["path"] == path) & (pathInfo["rank"] == rank)]["accession"].values[0]
            compare_label = acc2label[compare_acc]
            comparison_name = f"{anchor_acc}_x_{compare_acc}"
            pyani_ident, pyani_coverage, pyani_aln_length, pyani_sim_errors, pyani_hadamard = np.nan, np.nan, np.nan, np.nan, np.nan
            # pyani SHOULD produce symmetrical matrix, so all accs should exist regardless of value, right?
            pyani_ident = idD.at[anchor_label, compare_label]
            pyani_coverage = covD.at[anchor_label, compare_label]
            pyani_aln_length = lenD.at[anchor_label, compare_label]
            pyani_sim_errors = seD.at[anchor_label, compare_label]
            pyani_hadamard = hadD.at[anchor_label, compare_label]

            this_info = anchorpyani(comparison_name, anchor_acc, compare_acc, path, rank, pyani_ident, pyani_coverage, pyani_aln_length, pyani_sim_errors, pyani_hadamard)
            anchor_results.append(this_info)

    # convert path ANI comparison info to pandas dataframe
    anchor_aniDF = pd.DataFrame.from_records(anchor_results, columns = anchorpyani._fields)

    # print to csv
    anchor_aniDF.to_csv(args.output_csv, index=False)
    print(f"done! path comparison info written to {args.output_csv}")

def cmdline(sys_args):
    "Command line entry point w/argparse action."
    p = argparse.ArgumentParser()
    p.add_argument("pyani_results_dir")
    p.add_argument("--path-info", default="gtdb-r95-reps.pathinfo.tsv")
    p.add_argument("--path-name") # just do the one path
    p.add_argument("--labels")
    p.add_argument("--output-csv", required=True)
    args = p.parse_args()
    return main(args)

if __name__ == '__main__':
    returncode = cmdline(sys.argv[1:])
    sys.exit(returncode)

