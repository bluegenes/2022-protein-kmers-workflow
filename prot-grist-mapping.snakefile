## simplest approach:
# sketch translate metgenome
# sourmash gather --> prot gtdb
# using pre-downloaded /generated gtdb proteomes,
# index proteomes
# merge reads using bbmerge/pear
# map  reads --> proteomes
# generate summary stats
# generate grist-style plots

#configfile: prot-mapping.yml
configfile: "conf/protein-grist.yml"
SAMPLES = config["samples"]
out_dir = config["outdir"]
logs_dir = os.path.join(out_dir, "logs")

rule all:
    input: expand(out_dir + 'reports/report-{sample}.html', sample=SAMPLES)

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
#rule bbmerge_paired_reads_wc:
#    input:
#        interleaved = ancient(outdir + '/trim/{sample}.trim.fq.gz'),
#    output:
#        merged = protected(outdir + '/merge_reads/{sample}.merged.fq.gz'),
#        unmerged = protected(outdir + '/merge_reads/{sample}.unmerged.fq.gz'),
#        insert_size_hist = protected(outdir + '/merge_reads/{sample}.insert-size-histogram.txt'),
#    conda: 'env/bbmap.yml'
#    threads: 6
#    resources:
#        mem_mb = int(ABUNDTRIM_MEMORY / 1e6),
#    params:
#        mem = ABUNDTRIM_MEMORY,
#    shell: """
#            bbmerge.sh -t {threads} -Xmx {params.mem} in={input} \
#            out={output.merged} outu={output.unmerged} \
#            ihist={output.ihist}
#    """
#
#rule paladin_index_wc:
#    input:
#        proteome = Checkpoint_ProteomeFiles(f"{outdir}/proteomes/{{ident}}_protein.faa.gz"),
#    output:
#        idx = outdir + "/proteomes/{ident}_protein.faa.bwt",
#    conda: 'env/paladin.yml'
#    shell:
#        """
#        paladin index -r3 {input.proteome}
#        """

#rule paladin_align_wc:
#    input:
##        reads='{sample}.assembled.fastq',
#        merged_reads= outdir + '/merge_reads/{sample}.merged.fq.gz',
#        index= outdir + "/proteomes/{ident}_protein.faa.bwt",
#    output:
#        bam = outdir + "/protein_mapping/{sample}.x.{ident}.bam",
#    params:
#        index_base=lambda w: f"{outdir}/proteomes/{w.ident}_protein.faa"
##    log:
##        "logs/{sample}.paladin-align.log"
##    threads: 4
#    conda: "paladin.yml"
#    shell:
#        """
#        paladin align -t {threads} -T 20 {params.index_base} {input.reads} | samtools view -Sb - > {output}
#        """

# extract FASTQ from BAM
#rule bam_to_fastq_wc:
#    input:
#        bam = outdir + "/protein_mapping/{bam}.bam",
#    output:
#        mapped = outdir + "/protein_mapping/{bam}.mapped.fq.gz",
#    conda: "env/minimap2.yml"
#    shell: """
#        samtools bam2fq {input.bam} | gzip > {output.mapped}
#    """

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
#            out.write(f"outdir: {out_dir}\n")
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



