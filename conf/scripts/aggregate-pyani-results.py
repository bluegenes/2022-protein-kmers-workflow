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

compare_ranklist = ["genus", "family", "order", "class", "phylum", "superkingdom"]

def main(args):
    # get basename for these sequences
    anchor_results = []
    pathInfo = pd.read_csv(args.path_info, dtype=str, sep="\t", header=0)

    if "3" in str(args.pyani_version):
        len_fn= "matrix_aln_lengths_1.tab"
        cov_fn = "matrix_coverage_1.tab"
        had_fn = "matrix_hadamard_1.tab"
        id_fn = "matrix_identity_1.tab"
        se_fn = "matrix_sim_errors_1.tab"
    elif "2" in str(args.pyani_version):
        len_fn= "ANIb_alignment_lengths.tab"
        cov_fn = "ANIb_alignment_coverage.tab"
        had_fn = "ANIb_hadamard.tab"
        id_fn = "ANIb_percentage_identity.tab"
        se_fn = "ANIb_similarity_errors.tab"
    if args.path_name:
        select_path = args.path_name
        pathInfo = pathInfo[pathInfo["path"] == select_path] # subset
    allpaths = pathInfo['path'].unique().tolist()

    for path in allpaths:
        anchor_acc = pathInfo.loc[(pathInfo["path"] == path) & (pathInfo["rank"] == "anchor")]["accession"].values[0]
        compare_accs = pathInfo.loc[(pathInfo["path"] == path) & (pathInfo["rank"] != "anchor")]["accession"].values.tolist()

        path_results_dir = args.pyani_results_dir
        lenF = os.path.join(path_results_dir, len_fn)
        covF = os.path.join(path_results_dir, cov_fn)
        hadamardF = os.path.join(path_results_dir, had_fn)
        identF = os.path.join(path_results_dir, id_fn)
        simerrF = os.path.join(path_results_dir, se_fn)

        # read in all matrices
        lenD = pd.read_csv(lenF, sep="\t", header=0, index_col=0)
        covD = pd.read_csv(covF, sep="\t", header=0, index_col=0)
        hadD = pd.read_csv(hadamardF, sep="\t", header=0, index_col=0)
        idD = pd.read_csv(identF, sep="\t", header=0, index_col=0)
        seD = pd.read_csv(simerrF, sep="\t", header=0, index_col=0)

        # use headers on one file to get full column names:
        names = lenD.columns.tolist()
        anchor_label = [x for x in names if x.startswith(anchor_acc)][0]

        # get info for each comparison
        for rank in compare_ranklist:
            compare_acc = pathInfo.loc[(pathInfo["path"] == path) & (pathInfo["rank"] == rank)]["accession"].values[0]
            compare_label = [x for x in names if x.startswith(compare_acc)][0]
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
    p.add_argument("--pyani-version", default = 'v0.3')
    p.add_argument("--path-info", default="gtdb-rs202.taxonomy.paths.tsv")
    p.add_argument("--path-name") # just do the one path
    p.add_argument("--labels")
    p.add_argument("--output-csv", required=True)
    args = p.parse_args()
    return main(args)

if __name__ == '__main__':
    returncode = cmdline(sys.argv[1:])
    sys.exit(returncode)

