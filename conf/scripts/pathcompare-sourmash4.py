import os
import sys
import screed
import argparse
from collections import namedtuple
import csv
import pandas as pd

import sourmash
from sourmash import load_one_signature

# catch warnings
import warnings


CompareResult = namedtuple('CompareResult','comparison_name, anchor_name, compare_name, path, lowest_common_rank, alphabet, ksize, scaled, jaccard, max_containment, anchor_containment, compare_containment, anchor_hashes, compare_hashes, num_common, jaccard_ani, containA_ani, containAani_low, containAani_high, containBani, containBani_low, containBani_high, mcANI, acANI') #, jaccard_warning')

def compare_mh(mhA, mhB, A_name, B_name, lowest_common_rank, path_name, alpha, ksize, scaled):
    comparison_name = f"{A_name}_x_{B_name}"
    A_hashes = len(mhA)
    B_hashes = len(mhB)
    intersect_numhashes = mhA.count_common(mhB)
    jaccard = mhA.jaccard(mhB)
    containA = mhA.contained_by(mhB)
    containB = mhB.contained_by(mhA)
    max_contain = max(containA,containB)
    jaccard_warning_encountered=False
    # get ANI values
    query_cANI_result  = mhA.containment_ani(mhA, containment=containA, estimate_ci=True)
    # go around any zeroing out
    containAani = 1- query_cANI_result.dist
    containAani_high = 1- query_cANI_result.dist_low
    containAani_low = 1- query_cANI_result.dist_high
    #containAani,containAani_low,containAani_high = mhA.containment_ani(mhA, containment=containA, estimate_ci=True)
    match_cANI_result = mhB.containment_ani(mhA, containment=containB, estimate_ci=True)
    # go around any zeroing out
    containBani = 1- match_cANI_result.dist
    containBani_high = 1- match_cANI_result.dist_low
    containBani_low = 1- match_cANI_result.dist_high
    #containBani,containBani_low,containBani_high= mhB.containment_ani(mhA, containment=containB)
    mcANI = mhA.max_containment_ani(mhB, max_containment=max_contain).ani
    #mcANI,mcANI_low,mcANI_high = mhA.max_containment_ani(mhB, max_containment=max_contain)
    # figure out when jaccard RuntimeWarnings happen by catching warnings
    #with warnings.catch_warnings(record=True) as w:
    jaccard_ani = mhA.jaccard_ani(mhB, jaccard=jaccard).ani
    acANI = mhA.avg_containment_ani(mhB)
        #if len(w) > 0:
        #    jaccard_warning_encountered=True
    return CompareResult(comparison_name, A_name, B_name, path_name, lowest_common_rank, alpha, ksize,
                         scaled, jaccard, max_contain, containA, containB, A_hashes, B_hashes, intersect_numhashes,
                         jaccard_ani, containAani, containAani_low, containAani_high,
                         containBani, containBani_low, containBani_high, mcANI, acANI)#, jaccard_warning_encountered)


def find_sigfile(accession, sigdir, alphabet):
    #sigdir = os.path.abspath(sigdir)
    if alphabet in ["nucleotide", "dna", "rna"]:
        ext = "nucleotide"
    else:
        ext = "protein"
    sigfile = os.path.join(sigdir, f"{accession}.{ext}.sig")
    if not os.path.isfile(sigfile):
        print(f"ERROR: Can't find sigfile! {sigfile}. Exiting.")
        sys.exit(-1)
    return sigfile


def main(args):
    # handle alpha
    is_dna=False
    is_protein=True
    is_dayhoff=False
    is_hp=False

    alphabet = args.alphabet
    if alphabet in ['dna', 'rna', 'nucleotide']:
        is_dna=True
        is_protein=False
    elif alphabet == "dayhoff":
        is_dayhoff = True
    elif alphabet == "hp":
        is_hp = True

    scaled_vals = args.scaled
    scaled_vals.sort() # ascending order
    ksize = args.ksize
    sigdir = args.sigdir
    translated_sigdir = args.translated_sigdir

    select_ksize = ksize
    if is_protein:
        select_ksize = ksize *3


    # load path info
    pathinfo = pd.read_csv(args.paths_csv, dtype=str, sep="\t", header=0)

    # if a single path is provided, subset dataframe to just that path
    if args.path:
        pathinfo = pathinfo[pathinfo["path"] == args.path]

    paths = list(pathinfo["path"].drop_duplicates())
    groupbyPath = pathinfo.groupby('path')

    rank_order = ["genus", "family", "order", "class", "phylum", "superkingdom"]
    path_comparisons = []

    output_fieldnames = CompareResult._fields
    with open(args.output, 'w') as outF:
        w = csv.DictWriter(outF, fieldnames=output_fieldnames)
        w.writeheader()

        for n, path in enumerate(paths):
            #if n !=0 and n % 50 == 0:
            print(f"... assessing {n}th path, {path}\n")
            groupInfo = groupbyPath.get_group(path).set_index("rank")
            anchor_acc = groupInfo.at["anchor", "accession"]
            anchor_sigfile = find_sigfile(anchor_acc, sigdir, alphabet)
            anchor_sig = load_one_signature(anchor_sigfile, ksize=select_ksize)
            # compare vs rest of path
            path_comparisons = []
            for lowest_common_rank in rank_order:
                compare_acc = groupInfo.at[lowest_common_rank, "accession"]
                if args.compare_translated:
                    compare_sigfile = find_sigfile(compare_acc, translated_sigdir, alphabet)
                else:
                    compare_sigfile = find_sigfile(compare_acc, sigdir, alphabet)
                # select and load comparison sig
                compare_sig = load_one_signature(compare_sigfile, ksize=select_ksize)
                compare_sig.minhash = compare_sig.minhash.to_mutable()
                for sc in scaled_vals:
                    sc = int(sc)
                    # downsample to scaled
                    anc_mh = anchor_sig.minhash.downsample(scaled=sc)
                    cmp_mh = compare_sig.minhash.downsample(scaled=sc)
                    comparison = compare_mh(anc_mh, cmp_mh, anchor_acc, compare_acc, lowest_common_rank, path, alphabet, ksize, sc)
                    path_comparisons.append(comparison)

            for c in path_comparisons:
                w.writerow(c._asdict())



def cmdline(sys_args):
    "Command line entry point w/argparse action."
    p = argparse.ArgumentParser()
    p.add_argument("--paths-csv", required=True)
    p.add_argument('--sigdir', default = "")
    p.add_argument('--translated-sigdir', default = "")
    p.add_argument("-k", "--ksize",  type=int)
    p.add_argument("--alphabet")
    p.add_argument("--scaled", nargs='*', default=[1])
    p.add_argument("--path")
    p.add_argument("--compare-translated", action="store_true")
    p.add_argument("--output")
    args = p.parse_args()
    return main(args)

if __name__ == '__main__':
    returncode = cmdline(sys.argv[1:])
    sys.exit(returncode)

