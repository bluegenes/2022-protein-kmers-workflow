import os
import sys
import argparse
from collections import namedtuple
import csv
import pandas as pd

import sourmash
from sourmash import load_one_signature
from sourmash.logging import notify
from sourmash.distance_utils import set_size_chernoff


REL_ERR = [0.05, 0.1, 0.15, 0.2, 0.25, 0.3]
SizeACResult = namedtuple('SizeACResult', "name, alphabet, ksize, scaled, rel_err, probability")

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

def assess_size_accuracy(name, alphabet, ksize, mh, scaled_list, rel_err_list = REL_ERR):
    res_list = []
    for rel_err in rel_err_list:
        for sc in scaled_list:
            sc = int(sc)
            if sc < mh.scaled:
                notify(f"Scaled {sc} is lower than original mh scaled ({mh.scaled}). Skipping.")
                continue
            mh_sc = mh.downsample(scaled=sc)
            probability = set_size_chernoff(mh_sc.unique_dataset_hashes, sc, relative_error=rel_err)
            res = SizeACResult(name, alphabet, ksize, sc, rel_err, probability)
            res_list.append(res)
    #return probability >= confidence
    return res_list

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

    # just get the accessions
    accs = list(pathinfo["accession"].drop_duplicates())
    anchor_accs = pathinfo[pathinfo["rank"] == "anchor"]["accession"].tolist()
    n_accs = float(len(accs))

    output_fieldnames = SizeACResult._fields
    with open(args.output, 'w') as outF:
        w = csv.DictWriter(outF, fieldnames=output_fieldnames)
        w.writeheader()

        for n, acc in enumerate(accs):
            #if n % 100 == 0:
            perc = n/n_accs
            notify(f"progress: {perc:.1f}%")
            sigfile = find_sigfile(acc, sigdir, alphabet)
            sig = load_one_signature(sigfile, ksize=select_ksize)

            # assess size accuracy
            acc_res = assess_size_accuracy(sig.name, alphabet, ksize, sig.minhash, scaled_vals, args.relative_error)
            for res in acc_res:
                w.writerow(res._asdict())
            # for proteins, also do translated sig
            if is_protein and acc not in anchor_accs:
                tr_sigfile = find_sigfile(acc, translated_sigdir, alphabet)
                tr_sig = load_one_signature(tr_sigfile, ksize=select_ksize)
                # assess accuracy of translated
                tr_acc_res = assess_size_accuracy(sig.name, "tr-" + alphabet, ksize, sig.minhash, scaled_vals, args.relative_error)
                for res in tr_acc_res:
                    w.writerow(res._asdict())



def cmdline(sys_args):
    "Command line entry point w/argparse action."
    p = argparse.ArgumentParser()
    p.add_argument("--paths-csv", required=True)
    p.add_argument('--sigdir', default = "")
    p.add_argument('--translated-sigdir', default = "")
    p.add_argument("-k", "--ksize",  type=int)
    p.add_argument("--alphabet")
    p.add_argument("--scaled", nargs='*', default=[1])
    p.add_argument("--relative-error", nargs='*', default=REL_ERR)
    p.add_argument("--path")
    p.add_argument("--compare-translated", action="store_true")
    p.add_argument("--output")
    args = p.parse_args()
    return main(args)

if __name__ == '__main__':
    returncode = cmdline(sys.argv[1:])
    sys.exit(returncode)

