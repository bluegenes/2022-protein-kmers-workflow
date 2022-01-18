#! /usr/bin/env python3
import argparse
import sourmash.lca
from sourmash.lca.command_index import load_taxonomy_assignments
from sourmash.lca.lca_utils import LineagePair, pop_to_rank, is_lineage_match

import sys
import collections
import csv
import random

def get_name_at_rank(lineage, rank):
    for pair in lineage:
        if pair.rank == rank:
            return pair.name

    raise Exception


def main():
    p = argparse.ArgumentParser()
    p.add_argument('lineage_csv')
    p.add_argument("pathinfo")
    p.add_argument("--all-lineages")
    p.add_argument('-N', '--num_paths', type=int, default=0)
    args = p.parse_args()


    # this should just be representatives
    #representative_assignments, n_rows = load_taxonomy_assignments(args.lineage_csv,
    #                                                start_column=2)

    #print(f'loaded {len(representative_assignments)} representative assignments.')
    alltax, n_rows = load_taxonomy_assignments(args.all_lineages,
                                                    start_column=2)

    print(f'loaded {len(alltax)} assignments total.')

    representative_assignments = alltax

    # here, 'representative_assignments' is a dictionary where the keys are the
    # identifiers in column 1, and the values are tuples of LineagePairs.

    rep_paths_to_idents = collections.defaultdict(set)
    # build paths to idents
    # connect every lineage in lineage_paths to their identifiers.
    for ident, lineage in representative_assignments.items():
    #for ident, lineage in alltax.items():
        # bring back to each rank, skipping species
        for rank in sourmash.lca.taxlist(include_strain=False):
            if rank == 'species': continue

            tup = pop_to_rank(lineage, rank)
            assert tup
            assert tup[-1].rank == rank

            # record identifer for that lineage.
            rep_paths_to_idents[tup].add(ident)

    # now, collect all genus level paths
    genus_paths = set()
    for k, v in rep_paths_to_idents.items():
        if k[-1].rank == 'genus':
            genus_paths.add(k)

    # for each genus level path, pull out actual accessions for the
    # entire path.
    genus_paths = list(sorted(genus_paths))

    species_to_idents = collections.defaultdict(set)
    for ident, lineage in alltax.items():
        # get all species-level assignments
        tup = pop_to_rank(lineage, "species")
        assert tup
        assert tup[-1].rank == "species"
        # record identifer for that lineage.
        species_to_idents[tup].add(ident)

    # now find paths and write
    with open(args.pathinfo, "w") as outP:
        header=["accession", "path", "rank", "lineage"]
        outP.write("\t".join(header) + "\n")

        extract_idents = set()
        paths_found = 0
        paths_with_species = 0
        for chosen_genus in genus_paths:
            if args.num_paths and paths_found == args.num_paths:
                break

            # for this chosen genus, find identifiers for lineages that differ at
            # each step. lineages may not have this, note - CTB explain later :)
            d = {}
            lineage = list(chosen_genus)
            assert lineage[-1].rank == 'genus'
            last_rank = 'species'

            # pick a specific identifier and grab full lineage for it:
            chosen_ident = next(iter(rep_paths_to_idents[chosen_genus]))
            chosen_lineage = representative_assignments[chosen_ident]

            # now find differences at every level.
            n_found = 0
            while lineage:
                rank = lineage[-1].rank
                idents = rep_paths_to_idents[tuple(lineage)]

                for ident in idents:
                    this_lineage = representative_assignments[ident]

                    # find ones that match at this level, and not previous level
                    if is_lineage_match(this_lineage, chosen_lineage, rank) and \
                       not is_lineage_match(this_lineage, chosen_lineage, last_rank):
                         d[rank] = ident
                         n_found += 1
                         break

                last_rank = rank
                lineage.pop()

            if n_found < 6:
                continue                      # skip on to next genus-level path

            assert n_found == 6

            paths_found += 1
            pathname = f"path{str(paths_found)}" # use path number as a unique identifier for each path
            ### now, start tracking down identifiers to extract!

            extract_idents.add(chosen_ident)

            # print some stuff out --
            print('------------- path:', paths_found)
            for rank in sourmash.lca.taxlist():
                if rank == 'species':
                    # find a random species comparison
                    species_idents = species_to_idents[chosen_lineage]
                    species_idents = species_idents.remove(chosen_ident)
                    # if we have additional members of this species:
                    if species_idents:
                        import pdb;pdb.set_trace()
                        paths_with_species +=1
                        # pick random one with random.choice()
                        species_compare_ident = next(iter(species_to_idents[chosen_genus]))
                        print(' ', "species",species_compare_ident, sourmash.lca.display_lineage(chosen_lineage))
                        outP.write("\t".join([species_compare_ident, pathname, "species", sourmash.lca.display_lineage(chosen_lineage)]) + "\n")

                    # write anchor info
                    print(' ', "anchor", chosen_ident, sourmash.lca.display_lineage(chosen_lineage))
                    outP.write("\t".join([chosen_ident, pathname, "anchor", sourmash.lca.display_lineage(chosen_lineage)]) + "\n")
                    break

                assert rank in d
                ident = d[rank]
                this_lineage = representative_assignments[ident]
                print(' ', rank, ident, sourmash.lca.display_lineage(this_lineage))
                outP.write("\t".join([ident, pathname, rank, sourmash.lca.display_lineage(this_lineage)]) + "\n")

                extract_idents.add(ident)

        print(f'Total Paths:{paths_found}')
        print(f'Paths with species:{paths_with_species}')

if __name__ == '__main__':
    main()
