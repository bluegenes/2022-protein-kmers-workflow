import os, sys
import pandas as pd
import glob

from snakemake.workflow import srcdir

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
# subset to path accessions
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

onstart:
    print("------------------------------")
    print("Estimate similarity for 'evolutionary paths' genomes, proteomes")
    print("------------------------------")

onsuccess:
    print("\n--- Workflow executed successfully! ---\n")

onerror:
    print("Alas!\n")


rule all:
    input: 
        expand( os.path.join(out_dir, "compareM", "{basename}.path-compareM.{input_type}.csv.gz"), basename=basename, input_type = ["proteome"]),
        #expand( os.path.join(out_dir, "compareM", "{basename}.path-compareM.{input_type}.csv.gz"), basename=basename, input_type = ["genomic", "proteome"]),
        
        # fastani
        expand(os.path.join(out_dir, "fastani", "{basename}.path-fastani.csv.gz"),basename=basename),
        #expand(os.path.join(out_dir, "fastani-compare", "{basename}.fastani.tsv"), basename=basename),

        # compareM
        #os.path.join(out_dir, "compareM", "aai/aai_summary.tsv")


#### Protein CompareM Rules ###
localrules: build_protein_filepaths_for_compareM
rule build_protein_filepaths_for_compareM:
    output:
        os.path.join(out_dir, "compareM/proteome", "{path}", "{path}.protein-filepaths.txt")
    run:
        with open(str(output), "w") as out:
            acc_list = path2acc[wildcards.path]
            for acc in acc_list:
                fn = tax_info.at[acc, 'protein_fastafile']
                out.write(f"{fn}\n")


rule protein_AAI_via_compareM:
    input: 
        os.path.join(out_dir, "compareM/proteome", "{path}/{path}.protein-filepaths.txt")
    output:
        os.path.join(out_dir, "compareM/proteome/paths", "{path}/aai/aai_summary.tsv"),
    params:
        #proteins_cmd = "--proteins" if input_type == "protein" else "",
        #file_ext = ".faa.gz" if input_type == "protein" else ".fna",
        proteins_cmd = "--proteins",
        file_ext = ".faa.gz",
        outdir = lambda w: os.path.join(out_dir, "compareM/proteome/paths", w.path),
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt *3000,
        runtime=100,
    log: os.path.join(logs_dir, "compareM/proteome/paths", "{path}.compareM.log")
    benchmark: os.path.join(logs_dir, "compareM/proteome/paths", "{path}.compareM.benchmark")
    conda: "conf/env/compareM.yml"
    shell:
        """
        comparem aai_wf --cpus {threads} {params.proteins_cmd} --file_ext {params.file_ext:q}  --sensitive {input} {params.outdir} > {log} 2>&1
        """


localrules: write_protein_compareM_result_csv
rule write_protein_compareM_result_csv:
    input:
        expand(os.path.join(out_dir, "compareM/proteome/paths", "{path}/aai/aai_summary.tsv"), path=path_names)
    output:
        os.path.join(out_dir, "compareM/proteome", "{basename}.path-compareM.filecsv"),
    run:
        with open(str(output), "w") as out:
            for path in path_names:
                out.write(f"{path},{out_dir}/compareM/proteome/paths/{path}/aai/aai_summary.tsv\n")


localrules: aggregate_proteome_compareM_results
rule aggregate_proteome_compareM_results:
    input:
        compareM_proteome=os.path.join(out_dir, "compareM/proteome", "{basename}.path-compareM.filecsv"),
        paths=config["paths_csv"],
    output: 
        proteome=os.path.join(out_dir, "compareM", "{basename}.path-compareM.proteome.csv.gz"),
    log: os.path.join(logs_dir, "aggregate-compareM/proteome", "{basename}.path-compareM.aggregate.log")
    benchmark: os.path.join(logs_dir, "aggregate-compareM/proteome", "{basename}.path-compareM.aggregate.benchmark")
    shell:
        """
        python conf/scripts/aggregate-compareM-results.py --comparem-tsv-filecsv {input.compareM_proteome} \
                                             --path-info {input.paths} \
                                             --output-csv {output} > {log} 2>&1
        """

## nucleotide compareM ##
rule build_genome_filepaths_for_compareM:
    output:
        filepaths = temp(os.path.join(out_dir, "compareM/genomic/paths", "{path}", "{path}.genomic-filepaths.txt")),
        #fna_files= temp(expand(os.path.join(out_dir, "compareM/genomic", "{path}", "{acc}_genomic.fna"),acc=path2acc[{path}])),
    params: 
        outdir =  lambda w: os.path.join(out_dir, "compareM/genomic/paths", w.path)
    run:
        with open(str(output.filepaths), "w") as out:
            acc_list = path2acc[wildcards.path]
            for acc in acc_list:
                fn = tax_info.at[acc, 'genome_fastafile']
                fn_gunzip = os.path.join(params.outdir, os.path.basename(fn).rsplit(".gz")[0])
                shell("gunzip -c {fn} > {fn_gunzip}")
                out.write(f"{fn_gunzip}\n")


# CompareM from nucl uses Prodigal under the hood, but this prodigal cant use gzipped nucl files!!(???). Workaround by temp unzipping fna files.
rule nucl_AAI_via_compareM:
    input:
        filepaths=os.path.join(out_dir, "compareM/genomic/paths", "{path}/{path}.genomic-filepaths.txt"),
        fna_files= lambda w: expand(os.path.join(out_dir, "compareM/genomic/paths", "{path}", "{acc}_genomic.fna"),acc=path2acc[w.path])
    output:
        os.path.join(out_dir, "compareM/genomic", "{path}/aai/aai_summary.tsv"),
    params:
        proteins_cmd = "",
        file_ext = ".fna",
        outdir = lambda w: os.path.join(out_dir, "compareM/genomic/paths", w.path),
        fna_filepath = lambda w: os.path.join(out_dir, "compareM/genomic/paths", w.path, "*fna"),
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt *3000,
        runtime=100,
    log: os.path.join(logs_dir, "compareM/genomic/paths", "{path}.compareM.log")
    benchmark: os.path.join(logs_dir, "compareM/genomic/paths", "{path}.compareM.benchmark")
    conda: "conf/env/compareM.yml"
    shell:
        """
        comparem aai_wf --cpus {threads} {params.proteins_cmd} --file_ext {params.file_ext:q}  --sensitive {input.filepaths} {params.outdir} > {log} 2>&1
        # would be better to make these temp files
        rm -rf {params.fna_filepath}
        """

localrules: write_genomic_compareM_result_csv
rule write_genomic_compareM_result_csv:
    input:
        expand(os.path.join(out_dir, "compareM/genomic/paths", "{path}/aai/aai_summary.tsv"), path=path_names)
    output:
        os.path.join(out_dir, "compareM/genomic", "{basename}.path-compareM.filecsv"),
    run:
        with open(str(output), "w") as out:
            for path in path_names:
                out.write(f"{path},{out_dir}/compareM/genomic/{path}/aai/aai_summary.tsv\n")


localrules: aggregate_genomic_compareM_results
rule aggregate_genomic_compareM_results:
    input:
        compareM_genomic=os.path.join(out_dir, "compareM/genomic", "{basename}.path-compareM.filecsv"),
        paths=config["paths_csv"],
    output: 
        genomic=os.path.join(out_dir, "compareM", "{basename}.path-compareM.genomic.csv.gz"),
    log: os.path.join(logs_dir, "aggregate-compareM/genomic", "{basename}.path-compareM.aggregate.log")
    benchmark: os.path.join(logs_dir, "aggregate-compareM/genomic", "{basename}.path-compareM.aggregate.benchmark")
    shell:
        """
        python conf/scripts/aggregate-compareM-results.py --comparem-tsv-filecsv {input.compareM_genomic} \
                                             --path-info {input.paths} \
                                             --output-csv {output} > {log} 2>&1
        """


### fastANI rules ###
localrules: build_filepaths_for_fastani
rule build_filepaths_for_fastani:
    output: os.path.join(out_dir, "fastani", "{path}", "{path}.filepaths.txt")
    run:
        with open(str(output), "w") as out:
            acc_list = path2acc[wildcards.path]
            for acc in acc_list:
                fn = tax_info.at[acc, 'genome_fastafile']
                out.write(f"{fn}\n")


def get_genome_info(w):
    anchor_acc = paths[(paths["path"] == w.path) & (paths["rank"] == "anchor")].index[0]
    anchor_g = tax_info.at[anchor_acc, 'genome_fastafile']
    path_glist =  os.path.join(out_dir, "fastani", f"{w.path}/{w.path}.filepaths.txt")
    return {"anchor_genome": anchor_g, "path_genomes": path_glist}

rule compare_via_fastANI:
    input:  
        unpack(get_genome_info)
    output: os.path.join(out_dir, "fastani", "{path}/{path}.fastani.tsv"),
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt *3000,
        runtime=1200,
    log: os.path.join(logs_dir, "fastani", "{path}/{path}.fastani.log")
    benchmark: os.path.join(logs_dir, "fastani", "{path}/{path}.fastani.benchmark")
    conda: "conf/env/fastani.yml"
    shell:
        """
        fastANI -q {input.anchor_genome:q} --rl {input.path_genomes:q} -o {output} > {log} 2>&1
        """

localrules: write_fastani_result_csv
rule write_fastani_result_csv:
    input: expand(os.path.join(out_dir, "fastani", "{path}/{path}.fastani.tsv"), path=path_names)
    output: os.path.join(out_dir, "fastani", "{basename}.path-fastani.filecsv"),
    run:
        with open(str(output), "w") as out:
            for path in path_names:
                out.write(f"{path},{out_dir}/fastani/{path}/{path}.fastani.tsv\n")


localrules: aggregate_fastani_results
rule aggregate_fastani_results:
    input:
        fastani=os.path.join(out_dir, "fastani", "{basename}.path-fastani.filecsv"),
        paths=config["paths_csv"],
    output: os.path.join(out_dir, "fastani", "{basename}.path-fastani.csv.gz"),
    log: os.path.join(logs_dir, "fastani", "{basename}.path-fastani.aggregate.log")
    benchmark: os.path.join(logs_dir, "fastani", "{basename}.path-fastani.aggregate.benchmark")
    shell:
        """
        python conf/scripts/aggregate-fastani-results.py --fastani-filecsv {input.fastani} \
                                             --path-info {input.paths} \
                                             --output-csv {output} > {log} 2>&1
        """

