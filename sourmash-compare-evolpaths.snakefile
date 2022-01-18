"""
Author: Tessa Pierce Ward, UC Davis Lab for Data Intensive Biology
Run: snakemake -j1 -n
"""

import csv
import pandas as pd
from collections import defaultdict


configfile: "conf/gtdb-rs202-evolpaths.yaml"
basename = config['basename']
out_dir = config.get('output_dir', 'outputs')
logs_dir = os.path.join(out_dir, "logs")

# read in evolpaths --> get ACCS, build path2acc; acc: [path accs] dicts
pathinfo_file = config["paths_csv"]
paths = pd.read_csv(pathinfo_file, dtype=str, sep="\t", header=0)
ACCS = paths['accession'].unique().tolist()
# path: accessions dict
path2acc = paths.groupby('path')['accession'].apply(list).to_dict()
paths.set_index("accession", inplace=True)
path_names = path2acc.keys()

#anchor_acc = pathinfo[(pathinfo["path"] == w.path) & (pathinfo["rank"] == "species")].index[0]

# Grab taxonomic info (needed for signame, fastafiles) 
tax_info = pd.read_csv(config['taxonomy_csv'], header=0)
tax_info = tax_info[tax_info['ident'].isin(ACCS)]
tax_info['signame'] = tax_info['ident'] + ' ' + tax_info['species']
# add default fastapaths
tax_info['genome_fastafile'] = '/home/ntpierce/2021-rank-compare/genbank/genomes/'+ tax_info['ident'] + "_genomic.fna.gz"
tax_info['protein_fastafile'] = '/home/ntpierce/2021-rank-compare/genbank/proteomes/'+ tax_info['ident'] + "_protein.faa.gz"
# set index
tax_info.set_index("ident", inplace=True)

# change protein fastafile name for prodigal proteomes
prodigal_info = pd.read_csv("gtdb-rs202.prodigal-filenames.csv", index_col=0)
#update with prodigal filenames
tax_info.update(prodigal_info)

# check params are in the right format, build alpha-ksize combos
alphabet_info = config['alphabet_info']
alpha_ksizes=[]
prot_alpha_ksizes =[]
for alpha, info in alphabet_info.items():
    scaled = info["scaled"]
    ksize = info["ksize"]
    if not isinstance(scaled, list):
        scaled = [scaled]
        config["alphabet_info"][alpha]["scaled"] = scaled
    if not isinstance(ksize, list):
        ksize=[ksize]
        config["alphabet_info"][alpha]["ksize"] = ksize
    # build a parameter for the right combinations
    alpha_ksizes += expand(f"{alpha}-k{{ksize}}", ksize = ksize)
    # need separate for translated sigs
    if alpha == 'protein':
        prot_alpha_ksizes += expand(f"{alpha}-k{{ksize}}", ksize = ksize)


rule all:
    input:
        os.path.join(out_dir, "path-compare", f"{basename}.pathcompare.csv.gz"),
#        os.path.join(out_dir, "path-compare-translate", f"{basename}.pathcompare.csv.gz")

def make_param_str(ksizes, scaled):
    ks = [ f'k={k}' for k in ksizes ]
    ks = ",".join(ks)
    scaled = min(scaled) #take minimum value of scaled list
    return f"{ks},scaled={scaled},abund"


rule sketch_nucl:
    input: 
        lambda w: tax_info.at[w.acc, "genome_fastafile"]
    output: os.path.join(out_dir, "sketch", "{acc}.nucleotide.sig")
    params:
        sketch_params=make_param_str(config["alphabet_info"]["nucleotide"]["ksize"], config["alphabet_info"]["nucleotide"]["scaled"]),
        signame = lambda w: tax_info.at[w.acc, "signame"],
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt *6000,
        runtime=1200,
    log: os.path.join(logs_dir, "sourmash_sketch_genomic", "{acc}.sketch.log")
    benchmark: os.path.join(logs_dir, "sourmash_sketch_genomic", "{acc}.sketch.benchmark")
    conda: "conf/env/sourmash.yml"
    shell:
        """
        sourmash sketch dna {input} -p {params.sketch_params} -o {output} --name {params.signame:q} 2> {log}
        """    

rule sketch_prot:
    input: 
        lambda w: tax_info.at[w.acc, "protein_fastafile"]
    output: os.path.join(out_dir, "sketch", "{acc}.protein.sig")
    params:
        sketch_params=make_param_str(config["alphabet_info"]["protein"]["ksize"], config["alphabet_info"]["protein"]["scaled"]),
        signame = lambda w: tax_info.at[w.acc, "signame"],
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt *3000,
        runtime=1200,
    log: os.path.join(logs_dir, "sourmash_sketch_protein", "{acc}.sketch.log")
    benchmark: os.path.join(logs_dir, "sourmash_sketch_protein", "{acc}.sketch.benchmark")
    conda: "conf/env/sourmash.yml"
    shell:
        """
        sourmash sketch protein {input} -p {params.sketch_params} -o {output} --name {params.signame:q} 2> {log}
        """    

rule sketch_translate:
    input: 
        lambda w: tax_info.at[w.acc, "genome_fastafile"]
    output: os.path.join(out_dir, "sketch_translate", "{acc}.protein.sig")
    params:
        sketch_params=make_param_str(config["alphabet_info"]["protein"]["ksize"], config["alphabet_info"]["protein"]["scaled"]),
        signame = lambda w: tax_info.at[w.acc, "signame"],
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt *3000,
        runtime=1200,
    log: os.path.join(logs_dir, "sourmash_sketch_translate", "{acc}.sketch.log")
    benchmark: os.path.join(logs_dir, "sourmash_sketch_translate", "{acc}.sketch.benchmark")
    conda: "conf/env/sourmash.yml"
    shell:
        """
        sourmash sketch translate {input} -p {params.sketch_params} -o {output} --name {params.signame:q} 2> {log}
        """    


localrules: signames_to_file
rule signames_to_file:
    input:
        expand(os.path.join(out_dir, "sketch", "{acc}.{{alphabet}}.sig"), acc=ACCS)
    output: os.path.join(out_dir, f"{basename}.{{alphabet}}.siglist.txt")
    run:
        with open(str(output), "w") as outF:
            for inF in input:
                full_filename = os.path.abspath(str(inF))
                outF.write(f"{full_filename}\n")


localrules: translate_signames_to_file
rule translate_signames_to_file:
    input:
        expand(os.path.join(out_dir, "sketch_translate", "{acc}.protein.sig"), acc=ACCS)
    output: os.path.join(out_dir, f"{basename}.translate.protein.siglist.txt")
    run:
        with open(str(output), "w") as outF:
            for inF in input:
                full_filename = os.path.abspath(str(inF))
                outF.write(f"{full_filename}\n")


rule compare_paths:
    input: 
        siglist= os.path.join(out_dir, f"{basename}.{{alphabet}}.siglist.txt"),
        paths = config["paths_csv"], 
    output: 
        f"{out_dir}/path-compare/{{path}}.{{alphabet}}-k{{ksize}}.pathcompare.csv"
    params:
        sigdir = os.path.join(out_dir, "sketch"),
        scaled = lambda w: " ".join(str(x) for x in alphabet_info[w.alphabet]["scaled"])
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt *3000,
        runtime=1200,
    conda: "conf/env/sourmash-dist.yml"
    log: f"{logs_dir}/path-compare/{{path}}.{{alphabet}}-k{{ksize}}.pathcompare.log"
    benchmark: f"{logs_dir}/path-compare/{{path}}.{{alphabet}}-k{{ksize}}.pathcompare.benchmark",
    shell:
        """
        python conf/scripts/pathcompare.py --paths-csv {input.paths} -k {wildcards.ksize} --alphabet {wildcards.alphabet} \
                                         --path {wildcards.path} --output {output} --sigdir {params.sigdir} \
                                         --scaled {params.scaled} 2> {log}
        """

localrules: aggregate_pathcompare
rule aggregate_pathcompare:
    input:
        expand(os.path.join(out_dir, "path-compare", "{path}.{alphak}.pathcompare.csv"), path=path_names, alphak=alpha_ksizes)
    output:
        os.path.join(out_dir, "path-compare", "{basename}.pathcompare.csv.gz")
    run:
        # aggregate all csvs --> single csv
        aggDF = pd.concat([pd.read_csv(str(csv), sep=",") for csv in input])
        aggDF["alpha-ksize"] = aggDF["alphabet"] + "-" + aggDF["ksize"].astype(str)
        aggDF.to_csv(str(output), index=False)


rule compare_paths_translate:
    input: 
        siglist= os.path.join(out_dir, f"{basename}.translate.protein.siglist.txt"),
        paths = config["paths_csv"], 
    output: 
        f"{out_dir}/path-compare-translate/{{path}}.{{alphabet}}-k{{ksize}}.pathcompare.csv"
    params:
        sigdir = os.path.join(out_dir, "sketch_translate"),
        scaled = lambda w: " ".join(str(x) for x in alphabet_info[w.alphabet]["scaled"])
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt *3000,
        runtime=1200,
    conda: "conf/env/sourmash-dist.yml"
    log: f"{logs_dir}/path-compare-translate/{{path}}.{{alphabet}}-k{{ksize}}.pathcompare.log"
    benchmark: f"{logs_dir}/path-compare-translate/{{path}}.{{alphabet}}-k{{ksize}}.pathcompare.benchmark",
    shell:
        """
        python conf/scripts/pathcompare.py --paths-csv {input.paths} -k {wildcards.ksize} --alphabet {wildcards.alphabet} \
                                         --path {wildcards.path} --output {output} --sigdir {params.sigdir} \
                                         --scaled {params.scaled} 2> {log}
        """

localrules: aggregate_pathcompare_translate
rule aggregate_pathcompare_translate:
    input:
        expand(os.path.join(out_dir, "path-compare-translate", "{basename}.{alphak}.pathcompare.csv"), basename=basename, alphak=prot_alpha_ksizes)
    output:
        os.path.join(out_dir, "path-compare-translate", "{basename}.pathcompare.csv.gz")
    run:
        # aggregate all csvs --> single csv
        aggDF = pd.concat([pd.read_csv(str(csv), sep=",") for csv in input])
        aggDF["alpha-ksize"] = aggDF["alphabet"] + "-" + aggDF["ksize"].astype(str)
        aggDF.to_csv(str(output), index=False)

