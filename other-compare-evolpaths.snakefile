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
        # pyani
        #expand(os.path.join(out_dir, "pyani/paths/{path}/results/matrix_identity.tab"), path=path_names)
        #os.path.join(out_dir, "pyani", f"{basename}.pyani.csv.gz"),
        os.path.join(out_dir, "pyani", f"{basename}.pyani-ANIb.csv.gz")
        #os.path.join(out_dir, "pyani", f"{basename}.pyani.csv")
        # orthoani
        #os.path.join(out_dir, "orthoani", f"{basename}.orthoani.csv")
       
        
        # fastani
        #expand(os.path.join(out_dir, "fastani", "{basename}.path-fastani.csv.gz"),basename=basename),
        #expand(os.path.join(out_dir, "fastani-compare", "{basename}.fastani.tsv"), basename=basename),


        # compareM
        #os.path.join(out_dir, "compareM", "aai/aai_summary.tsv")
       # expand( os.path.join(out_dir, "compareM", "{basename}.path-compareM.{input_type}.csv.gz"), basename=basename, input_type = ["proteome"]),
        #expand( os.path.join(out_dir, "compareM", "{basename}.path-compareM.{input_type}.csv.gz"), basename=basename, input_type = ["genomic", "proteome"]),


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


def get_genomes_for_orthoani(w):
    anchor_acc = paths[(paths["path"] == w.path) & (paths["rank"] == "anchor")].index[0]
    anchor_g = tax_info.at[anchor_acc, 'genome_fastafile']
    path_accs = path2acc[w.path]
    compare_genomes = []
    for acc in path_accs:
        if acc != anchor_acc:
            compare_genomes.append(tax_info.at[acc, 'genome_fastafile'])
    return {"anchor_genome": anchor_g, "path_genomes": compare_genomes}


## orthoANI ##
rule compare_via_orthoANI:
    input:  
        unpack(get_genome_info) # use fastani compare pathlist
        #unpack(get_genomes_for_orthoani)
    output: os.path.join(out_dir, "orthoani", "paths", "{path}.orthoani.csv"),
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt *3000,
        runtime=1200,
    log: os.path.join(logs_dir, "orthoani", "{path}.orthoani.log")
    benchmark: os.path.join(logs_dir, "orthoani", "{path}.orthoani.benchmark")
    conda: "conf/env/orthoani.yml"
        #orthoani -q {input.anchor_genome:q} --r {input.path_genomes:q} -o {output} > {log} 2>&1 # commandline
    shell:
        """
        python conf/scripts/orthoani.py --path {wildcards.path} --anchor-genome {input.anchor_genome} \
               --compare-genome-pathlist {input.path_genomes} --output {output} > {log}
        """
                
rule aggregate_orthoani:
    input: expand(os.path.join(out_dir, "orthoani", "paths", "{path}.orthoani.csv"), path=path_names)
    output: os.path.join(out_dir, "orthoani", "{basename}.orthoani.csv")
    run:
        # aggregate all csvs --> single csv
        aggDF = pd.concat([pd.read_csv(str(csv), sep=",") for csv in input])
        aggDF.to_csv(str(output), index=False)



def get_genomes_for_pyani(w):
    path_accs = path2acc[w.path]
    path_genomes = []
    for acc in path_accs:
        path_genomes.append(tax_info.at[acc, 'genome_fastafile'])
    return path_genomes 



### split into 6-genome folders + unzip fna files and generate classes/labels, then run pyANI index and compare 
# make a folder with just the genomes in a single path
localrules: make_pyani_path_folder
rule make_pyani_path_folder:
    input:
        get_genomes_for_pyani
    output: 
        classes=os.path.join(out_dir, "pyani", "paths", "{path}", "py.classes.txt"),
        labels=os.path.join(out_dir, "pyani", "paths", "{path}", "py.labels.txt")
    params:
        acc_list = lambda w: path2acc[w.path],
        pathdir = lambda w: os.path.join(out_dir, 'pyani', 'paths', w.path),
    run:
        import hashlib
        os.makedirs(params.pathdir, exist_ok=True)
        with open(output.classes, 'w') as out_classes:
            with open(output.labels, 'w') as out_labels:
                for acc in params.acc_list:
                    fn = tax_info.at[acc, 'genome_fastafile']
                    dest = os.path.join(params.pathdir, f"{acc}_genomic.fna.gz")
                    dest_unz = os.path.join(params.pathdir, f"{acc}_genomic.fna")
                    md5_file = os.path.join(params.pathdir, f"{acc}_genomic.md5")
                    shell("cp {fn} {dest}")
                    shell("gunzip -c {fn} > {dest_unz}")
                    # get md5 of unzipped fna
                    with open(dest_unz, "rb") as f:
                        bytes = f.read()
                        md5 = hashlib.md5(bytes).hexdigest()
                    # write to md5 file
                    with open(md5_file, 'w') as mfile:
                        mfile.write(f"{md5}\t{dest_unz}\n")
                    fna_base = os.path.basename(dest_unz).rsplit('.fna')[0]
                    out_classes.write(f"{md5}\t{fna_base}\t{acc}\n")
                    out_labels.write(f"{md5}\t{fna_base}\t{acc}\n")


rule pyani_index_and_createdb:
    input: 
        os.path.join(out_dir, "pyani", "paths", "{path}", "py.classes.txt"),
        os.path.join(out_dir, "pyani", "paths", "{path}", "py.labels.txt")
    output:
        classes=os.path.join(out_dir, "pyani/paths", "{path}/classes.txt"),
        labels=os.path.join(out_dir, "pyani/paths", "{path}/labels.txt"),
        db=os.path.join(out_dir, "pyani/paths/{path}/.pyani-{path}/pyanidb")
    params:
        pathdir = lambda w: os.path.join(out_dir, 'pyani/paths', w.path),
        pyanidb = lambda w: os.path.join(out_dir, 'pyani/paths', w.path, f".pyani-{w.path}/pyanidb"),
        classes_basename = "classes.txt",
        labels_basename = "labels.txt"
    log: os.path.join(logs_dir, "pyani", "{path}.index-and-createdb.log")
    benchmark: os.path.join(logs_dir, "pyani", "{path}.index-and-createdb.benchmark")
    conda: "conf/env/pyani0.3.yml"
    shell:
        """
        pyani index -i {params.pathdir} --classes {params.classes_basename} --labels {params.labels_basename} 
        pyani createdb --dbpath {params.pyanidb} -v -l {log}
        """

    
rule pyANI_ANIm:
    input:  
        classes=os.path.join(out_dir, "pyani/paths", "{path}/py.classes.txt"),
        labels=os.path.join(out_dir, "pyani/paths", "{path}/py.labels.txt"),
        idx_classes=os.path.join(out_dir, "pyani/paths", "{path}/classes.txt"),
        idx_labels=os.path.join(out_dir, "pyani/paths", "{path}/labels.txt")
    output: 
        directory(os.path.join(out_dir, "pyani/paths", "{path}/results/nucmer_output")),
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt *5000,
        runtime=1200,
    params:
        pyanidb = lambda w: os.path.join(out_dir, 'pyani/paths', w.path, f".pyani-{w.path}/pyanidb"),
        genome_dir = lambda w: os.path.join(out_dir, 'pyani/paths', w.path),
        output_dir = lambda w: os.path.join(out_dir, 'pyani/paths', w.path, "results"),
    log: os.path.join(logs_dir, "pyani", "paths/{path}/{path}.pyANI-aniM.log")
    benchmark: os.path.join(logs_dir, "pyani", "paths/{path}/{path}.pyANI-aniM.benchmark")
    conda: "conf/env/pyani0.3.yml"
    shell:
        """
        pyani anim --dbpath {params.pyanidb} --name {wildcards.path} \
             --classes {input.classes} --labels {input.labels} \
            -i {params.genome_dir} -o {params.output_dir} \
            -l {log} -v
        """

rule pyANI_ANIb:
    input:  
        classes=os.path.join(out_dir, "pyani/paths", "{path}/py.classes.txt"),
        labels=os.path.join(out_dir, "pyani/paths", "{path}/py.labels.txt"),
        idx_classes=os.path.join(out_dir, "pyani/paths", "{path}/classes.txt"),
        idx_labels=os.path.join(out_dir, "pyani/paths", "{path}/labels.txt")
    output: 
        covF= os.path.join(out_dir, "pyani/paths/{path}/ANIb_results", "ANIb_alignment_coverage.tab"),
        lenF= os.path.join(out_dir, "pyani/paths/{path}/ANIb_results","ANIb_alignment_lengths.tab"),
        hadF= os.path.join(out_dir, "pyani/paths/{path}/ANIb_results","ANIb_hadamard.tab"),
        idF=  os.path.join(out_dir, "pyani/paths/{path}/ANIb_results","ANIb_percentage_identity.tab"),
        seF=  os.path.join(out_dir, "pyani/paths/{path}/ANIb_results","ANIb_similarity_errors.tab"),
        bn =  os.path.join(out_dir, "pyani/paths/{path}/ANIb_results","blastn_output.tar.gz"),
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt *5000,
        runtime=1200,
    params:
        #pyanidb = lambda w: os.path.join(out_dir, 'pyani/paths', w.path, f".pyani-{w.path}/pyanidb"),
        genome_dir = lambda w: os.path.join(out_dir, 'pyani/paths', w.path),
        output_dir = lambda w: os.path.join(out_dir, 'pyani/paths', w.path, "ANIb_results"),
    log: os.path.join(logs_dir, "pyani", "paths/{path}/{path}.pyANI-anib.log")
    benchmark: os.path.join(logs_dir, "pyani", "paths/{path}/{path}.pyANI-anib.benchmark")
    conda: "conf/env/pyani0.2.yml"
    shell:
        """
        average_nucleotide_identity.py -i {params.genome_dir} \
             -o {params.output_dir} -m ANIb -v \
             --labels {input.labels} --classes {input.classes} \
             --force > {log}
        """

rule pyANI_report_ANIm:
    input:  
        os.path.join(out_dir, "pyani/paths", "{path}/results/nucmer_output"),
    output: 
        idF = os.path.join(out_dir, "pyani/paths/{path}/results/matrix_identity_1.tab"),
        lenF = os.path.join(out_dir, "pyani/paths/{path}/results/matrix_aln_lengths_1.tab"),
        covF = os.path.join(out_dir, "pyani/paths/{path}/results/matrix_coverage_1.tab"),
        seF = os.path.join(out_dir, "pyani/paths/{path}/results/matrix_sim_errors_1.tab"),
        hadF = os.path.join(out_dir, "pyani/paths/{path}/results/matrix_hadamard_1.tab"),
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt *5000,
        runtime=1200,
    params:
        ANIm_dir = lambda w: os.path.join(out_dir, 'pyani/paths', w.path, "results"),
        pyanidb = lambda w: os.path.join(out_dir, 'pyani/paths', w.path, f".pyani-{w.path}/pyanidb"),
    log: os.path.join(logs_dir, "pyani", "paths/{path}/{path}.pyANI-aniM.log")
    benchmark: os.path.join(logs_dir, "pyani", "paths/{path}/{path}.pyANI-aniM.benchmark")
    conda: "conf/env/pyani0.3.yml"
    shell:
        """
        pyani report -o {params.ANIm_dir} --dbpath {params.pyanidb} --runs --run_results 1 --formats stdout -l {log}
        pyani report -v -o {params.ANIm_dir} --formats stdout --run_matrices 1  \
        --dbpath {params.pyanidb} -l {log}
        """

localrules: aggregate_path_anim
rule aggregate_path_anim:
    input:
        idF = os.path.join(out_dir, "pyani/paths/{path}/results/matrix_identity_1.tab"),
        lenF = os.path.join(out_dir, "pyani/paths/{path}/results/matrix_aln_lengths_1.tab"),
        covF = os.path.join(out_dir, "pyani/paths/{path}/results/matrix_coverage_1.tab"),
        seF = os.path.join(out_dir, "pyani/paths/{path}/results/matrix_sim_errors_1.tab"),
        hadF = os.path.join(out_dir, "pyani/paths/{path}/results/matrix_hadamard_1.tab"),
    output:
        os.path.join(out_dir, "pyani/paths", "{path}/results/{path}.pyani.csv"),
    params:
        results_dir =  os.path.join(out_dir, "pyani/paths/{path}/results")
    shell:
        """
        python conf/scripts/aggregate-pyani-results.py {params.results_dir} --path-name {wildcards.path} --output-csv {output}
        """

localrules: aggregate_path_anib
rule aggregate_path_anib:
    input:
        covF= os.path.join(out_dir, "pyani/paths/{path}/ANIb_results", "ANIb_alignment_coverage.tab"),
        lenF= os.path.join(out_dir, "pyani/paths/{path}/ANIb_results","ANIb_alignment_lengths.tab"),
        hadF= os.path.join(out_dir, "pyani/paths/{path}/ANIb_results","ANIb_hadamard.tab"),
        idF=  os.path.join(out_dir, "pyani/paths/{path}/ANIb_results","ANIb_percentage_identity.tab"),
        seF=  os.path.join(out_dir, "pyani/paths/{path}/ANIb_results","ANIb_similarity_errors.tab"),
    output:
        os.path.join(out_dir, "pyani/paths", "{path}/ANIb_results/{path}.pyani.csv"),
    params:
        results_dir =  os.path.join(out_dir, "pyani/paths/{path}/ANIb_results")
    shell:
        """
        python conf/scripts/aggregate-pyani-results.py {params.results_dir} --path-name {wildcards.path} --output-csv {output} --pyani-version v0.2
        """

localrules: aggregate_all_anim
rule aggregate_all_anim:
    input:
        expand(os.path.join(out_dir, "pyani/paths", "{path}/results/{path}.pyani.csv"), path=path_names)
    output:
        os.path.join(out_dir, "pyani", "{basename}.pyani.csv.gz")
    run:
        # aggreate all csv.gzs --> single csv
        aggDF = pd.concat([pd.read_csv(str(csv), sep=",") for csv in input])
        aggDF.to_csv(str(output), index=False)

localrules: aggregate_all_anib
rule aggregate_all_anib:
    input:
        expand(os.path.join(out_dir, "pyani/paths", "{path}/ANIb_results/{path}.pyani.csv"), path=path_names)
    output:
        os.path.join(out_dir, "pyani", "{basename}.pyani-ANIb.csv.gz")
    run:
        # aggreate all csv.gzs --> single csv
        aggDF = pd.concat([pd.read_csv(str(csv), sep=",") for csv in input])
        aggDF.to_csv(str(output), index=False)


#localrules: make_pyani_path_folder
#rule make_pyani_path_folder_maybe_use_for_orthoani:
#    input:
#        get_genomes_for_pyani
#        #lambda w: expand(tax_info.at[{acc}, 'genome_fastafile'], acc=path2acc[w.path])
#    output:
#        classes=os.path.join(out_dir, "pyani", "paths", "{path}", "py.classes.txt"),
#        labels=os.path.join(out_dir, "pyani", "paths", "{path}", "py.labels.txt")
#    params:
#        acc_list = lambda w: path2acc[w.path],
#        pathdir = lambda w: os.path.join(out_dir, 'pyani', 'paths', w.path),
#    run:
#        import hashlib
#        import screed
#        os.makedirs(params.pathdir, exist_ok=True)
#        with open(output.classes, 'w') as out_classes:
#            with open(output.labels, 'w') as out_labels:
#                for acc in params.acc_list:
#                    fn = tax_info.at[acc, 'genome_fastafile']
#                    dest = os.path.join(params.pathdir, f"{acc}_genomic.fna.gz")
#                    #dest_unz = os.path.join(params.pathdir, f"{acc}_genomic.fna")
#                    unzipped_single_fasta = os.path.join(params.pathdir, f"{acc}_genomic.fna")
#                    md5_file = os.path.join(params.pathdir, f"{acc}_genomic.md5")
#                    shell("cp {fn} {dest}")
#                    #shell("gunzip -c {fn} > {dest_unz}")
#                    # get md5 of unzipped fna
#                    #with open(dest_unz, "rb") as f:
#                    #    bytes = f.read()
#                    #    md5 = hashlib.md5(bytes).hexdigest()
#                    # write to md5 file
#                    #with open(md5_file, 'w') as mfile:
#                    #    mfile.write(f"{md5}\t{dest_unz}\n")
#                    #fna_base = os.path.basename(dest_unz).rsplit('.fna')[0]
#                    # now handle single fasta aggregation
#                    with open(unzipped_single_fasta, 'w') as sfa:
#                        sequences=[]
#                        with screed.open(fn) as ff:
#                            for read in ff:
#                                sequences.append(read.sequence)
#                        sfa.write(f'>{acc}\n')
#                        sfa.write("".join(sequences) + "\n")
#                    # get md5 of this file
#                    with open(unzipped_single_fasta, "rb") as f:
#                        bytes = f.read()
#                        md5 = hashlib.md5(bytes).hexdigest()
#                    # write to md5 file
#                    with open(single_md5_file, 'w') as mfile:
#                        mfile.write(f"{md5}\t{unzipped_single_fasta}\n")
#                    fna_base = os.path.basename(unzipped_single_fasta).rsplit('.fna')[0]
#                    # write labels and classes for these files
#                    out_classes.write(f"{md5}\t{fna_base}\t{acc}\n")
#                    out_labels.write(f"{md5}\t{fna_base}\t{acc}\n")
#
