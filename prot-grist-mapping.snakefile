## simplest approach:
# sketch translate metgenome
# sourmash gather --> prot gtdb
# using pre-downloaded /generated gtdb proteomes,
# index proteomes
# merge reads using bbmerge/pear
# map  reads --> proteomes
# generate summary stats
# generate grist-style plots
import pandas as pd
#configfile: prot-mapping.yml
configfile: "conf/protein-grist.yml"
SAMPLES = config["samples"]
out_dir = config["outdir"]
logs_dir = os.path.join(out_dir, "logs")

ABUNDTRIM_MEMORY = float(config['metagenome_trim_memory'])

# make dictionary of proteomes
proteome_inf = pd.read_csv("/home/ntpierce/2021-rank-compare/output.rank-compare/gtdb-rs207.fromfile.csv")
proteome_inf.set_index('ident', inplace=True)
# to do-- read /get these from gather results! 
proteomes_to_map = ['GCF_014131715.1']

rule all:
    input: 
        # grist outputs
        #expand(out_dir + 'reports/report-{sample}.html', sample=SAMPLES)
        # read mapping outputs
        expand(out_dir + '/merge_reads/{sample}.merged.fq.gz', sample=SAMPLES),
        expand(out_dir + "/proteomes/{ident}_protein.faa.gz.bwt",ident=proteomes_to_map),
        expand(out_dir + "/protein_mapping/{sample}.x.{ident}.mapped.fq.gz", sample=SAMPLES, ident=proteomes_to_map),
        expand(out_dir + "/protein_mapping/{sample}.x.{ident}.regions.bed", sample=SAMPLES, ident=proteomes_to_map),
        expand(out_dir + "/protein_mapping/{sample}.x.{ident}.depth.txt", sample=SAMPLES, ident=proteomes_to_map),
        expand(out_dir + "/protein_mapping/{sample}.x.{ident}.vcf.gz", sample=SAMPLES, ident=proteomes_to_map),
        expand(out_dir + "/protein_mapping/{sample}.summary.csv", sample=SAMPLES),

rule run_genome_grist_til_gather:
    input: "conf/protein-grist.yml",
        #config=config["grist-config"] #os.path.join(out_dir, f"config.grist.{basename}.yml"),
    output: expand(os.path.join(out_dir, "reports/report-{sample}.html"), sample=SAMPLES)
    log: os.path.join(logs_dir, "genome-grist.log")
    benchmark: os.path.join(logs_dir, "genome-grist.benchmark")
    threads: 32
    resources:
        mem_mb=150000,
        #mem_mb=lambda wildcards, attempt: attempt*145000
    conda: "conf/env/grist.yml"
    shell:
        """
        genome-grist run {input} --resources mem_mb={resources.mem_mb} \
                          -j {threads} summarize_gather --nolock \
                          --profile slurm -n 2> {log}
        """

#rule protein_mapping:
# bbmerge reads for use with prodigal protein mapper
# does this work with gzipped reads?
rule bbmerge_paired_reads:
    input:
        interleaved = ancient(out_dir + '/trim/{sample}.trim.fq.gz'),
    output:
        merged = protected(out_dir + '/merge_reads/{sample}.merged.fq.gz'),
        unmerged = protected(out_dir + '/merge_reads/{sample}.unmerged.fq.gz'),
        insert_size_hist = protected(out_dir + '/merge_reads/{sample}.insert-size-histogram.txt'),
    conda: 'conf/env/bbmap.yml'
    threads: 6
    resources:
        mem_mb = int(ABUNDTRIM_MEMORY / 1e6),
    params:
        mem = ABUNDTRIM_MEMORY,
    shell: 
        """
        bbmerge.sh -Xmx25G in={input} \
            out={output.merged} outu={output.unmerged} \
            ihist={output.insert_size_hist}
        """

## TO DO: 
#1. use GTDB proteomes available in 2021-rank-compare/genbank/proteomes and/or 2021-rank-compare/prodigal/
#2. checkpoint rule to read in gather results and load all idents for mapping. Does this need to be in order??
#3. test index and alignment w/paladin


rule paladin_index_wc:
    input:
        #Checkpoint_ProteomeFiles(f"{out_dir}/proteomes/{{ident}}_protein.faa.gz"),
        proteome = lambda w: "/home/ntpierce/2021-rank-compare/" + proteome_inf.at[f"{w.ident}", "protein_filename"]
    output:
        idx = out_dir + "/proteomes/{ident}_protein.faa.gz.bwt",
    params:
        prot_dir = out_dir + '/proteomes',
        proteome = lambda w: f"{out_dir}/proteomes/{w.ident}_protein.faa.gz",
    conda: 'conf/env/paladin.yml'
    shell:
        """
        mkdir -p {params.prot_dir}
        cp {input.proteome} {params.prot_dir}
        paladin index -r3 {params.proteome}
        """
        #paladin index -r3 {input.proteome}

rule paladin_align_wc:
    input:
        merged_reads= out_dir + '/merge_reads/{sample}.merged.fq.gz',
        index= out_dir + "/proteomes/{ident}_protein.faa.gz.bwt",
    output:
        bam = out_dir + "/protein_mapping/{sample}.x.{ident}.bam",
    params:
        index_base=lambda w: f"{out_dir}/proteomes/{w.ident}_protein.faa.gz"
    log:
        out_dir + "logs/{sample}.x.{ident}.paladin-align.log"
    benchmark:
        out_dir + "logs/{sample}.x.{ident}.paladin-align.benchmark"
    threads: 4
    conda: "conf/env/paladin.yml"
    shell:
        """
        paladin align -t {threads} -T 20 {params.index_base} {input.merged_reads} | samtools view -Sb - > {output}
        """

# sort bam 
rule sort_bam:
    input:
        bam = out_dir + "/{dir}/{bam}.bam",
    output:
        sort = out_dir + "/{dir}/sorted/{bam}.sorted.bam",
    conda: "conf/env/minimap2.yml"
    shell: 
        """
        samtools sort {input.bam} -o {output.sort} -O bam
        """

# extract FASTQ from BAM
rule bam_to_fastq_wc:
    input:
        sorted_bam = out_dir + "/protein_mapping/sorted/{bam}.sorted.bam",
    output:
        mapped = out_dir + "/protein_mapping/{bam}.mapped.fq.gz",
    conda: "conf/env/minimap2.yml"
    shell: 
        """
        samtools bam2fq {input.sorted_bam} | gzip > {output.mapped}
        """


# get per-base depth information from BAM
rule bam_to_depth_wc:
    input:
        sorted_bam = out_dir + "/{dir}/sorted/{bam}.sorted.bam",
    output:
        depth = out_dir + "/{dir}/{bam}.depth.txt",
    conda: "conf/env/minimap2.yml"
    shell: 
        """
        samtools depth -aa {input.sorted_bam} > {output.depth}
        """

# wild card rule for getting _covered_ regions from BAM
rule bam_covered_regions_wc:
    input:
        sorted_bam = out_dir + "/{dir}/sorted/{bam}.sorted.bam",
    output:
        regions = out_dir + "/{dir}/{bam}.regions.bed",
    conda: "conf/env/covtobed.yml"
    shell: 
        """
        covtobed {input.sorted_bam} -l 100 -m 1 | \
            bedtools merge -d 5 -c 4 -o mean > {output.regions}
        """

# calculating SNPs/etc.
rule mpileup_wc:
    input:
        query = out_dir + "/proteomes/{ident}_protein.faa.gz",
        sorted_bam = out_dir + "/{dir}/sorted/{sample}.x.{ident}.sorted.bam",
    output:
        bcf = out_dir + "/{dir}/{sample}.x.{ident}.bcf",
        vcf = out_dir + "/{dir}/{sample}.x.{ident}.vcf.gz",
        vcfi = out_dir + "/{dir}/{sample}.x.{ident}.vcf.gz.csi",
    conda: "conf/env/bcftools.yml"
    shell: 
        """
        proteomefile=$(mktemp -t g.proteome.XXXXXXX)
        gunzip -c {input.query} > $proteomefile
        bcftools mpileup -Ou -f $proteomefile {input.sorted_bam} | bcftools call -mv -Ob -o {output.bcf}
        rm $proteomefile
        bcftools view {output.bcf} | bgzip > {output.vcf}
        bcftools index {output.vcf}
        """

# summarize depth into a CSV
rule summarize_samtools_depth_wc:
    input:
        depth = expand(out_dir + "/{{dir}}/{{sample}}.x.{ident}.depth.txt", ident = proteomes_to_map),
        vcf = expand(out_dir + "/{{dir}}/{{sample}}.x.{ident}.vcf.gz", ident = proteomes_to_map)
        #depth = Checkpoint_GatherResults(outdir + f"/{{dir}}/{{sample}}.x.{{ident}}.depth.txt"),
        #vcf = Checkpoint_GatherResults(outdir + f"/{{dir}}/{{sample}}.x.{{ident}}.vcf.gz")
    output:
        csv = f"{out_dir}/{{dir}}/{{sample}}.summary.csv"
    shell: 
        """
        python summarize_mapping.py {wildcards.sample} \
             {input.depth} -o {output.csv}
        """

#rule write_grist_config:
#    input:
#        database=os.path.join(out_dir, "databases", "{basename}.zip"),
#        metagenomes=config["metagenome_list"],
#    output: os.path.join(out_dir, "config.grist.{basename}.yml")
#    params:
#        metagenome_trim_memory=config.get("metagenome_trim_memory", "1e9"),
#        ksize=ksize,
#        search_ksize=config["search_ksize"]
#    run:
#        with open(str(output), 'w') as out:
#            out.write(f"out_dir: {out_dir}\n")
#            out.write(f"metagenome_trim_memory: {params.metagenome_trim_memory}\n")
#            out.write(f"sourmash_database_glob_pattern: {input.database}\n") # can this have filepath, or does it need to be the basename only?
#            out.write(f"sample:\n")
#            for mg in metagenomes:
#                out.write(f"  - {mg}\n")
#            out.write(f"sourmash_database_ksize: {params.search_ksize}\n")
#            out.write(f"sourmash_compute_ksizes:\n")
#            for k in params.ksize:
#                out.write(f"  - {k}\n")
#            out.write("tempdir:")
#            out.write("  - /scratch")



