## simplest approach:
# sketch translate metgenome
# sourmash gather --> prot gtdb
# using pre-downloaded /generated gtdb proteomes,
# index proteomes
# merge reads using bbmerge/pear
# map  reads --> proteomes
# generate summary stats
# generate grist-style plots
import csv
import pandas as pd
import glob
#configfile: prot-mapping.yml
configfile: "conf/protein-grist-reps.yml"
SAMPLES = config["samples"]

out_dir = config["outdir"]
logs_dir = os.path.join(out_dir, "logs")

ABUNDTRIM_MEMORY = float(config['metagenome_trim_memory'])

# make dictionary of proteomes
proteome_inf = pd.read_csv("/home/ntpierce/2021-rank-compare/output.rank-compare/gtdb-rs207.fromfile.csv")
proteome_inf.set_index('ident', inplace=True)


rule all:
    input: 
        # grist outputs
        ancient(expand(out_dir + '/reports/report-{sample}.html', sample=SAMPLES)),
        # read mapping outputs
        #expand(out_dir + '/merge_reads/{sample}.merged.fq.gz', sample=SAMPLES),
        ancient(expand(out_dir + "/protein_mapping/{sample}.summary.csv", sample=SAMPLES)),
        expand(f"{out_dir}/protein_leftover/{{sample}}.summary.csv", sample=SAMPLES),

# use existed rs202 gather results for now
#checkpoint gather_reads_wc:
#    input:
#        gather_csv = f'output.protein-grist-rs202/gather/{{sample}}.gather.csv'
#        #gather_csv = f'{out_dir}/gather/{{sample}}.gather.csv'
#    output:
#        touch(f"output.protein-grist-rs202/gather/.gather.{{sample}}")   # checkpoints need an output ;)
#       #touch(f"{out_dir}/gather/.gather.{{sample}}")   # checkpoints need an output ;)

class Checkpoint_GatherResults:
    """Given a pattern containing {ident} and {sample}, this class
    will generate the list of {ident} for that {sample}'s gather results.
    Alternatively, you can omit {sample} from the pattern and include the
    list of sample names in the second argument to the constructor.
    """
    def __init__(self, pattern, samples=None):
        self.pattern = pattern
        self.samples = samples

    def get_genome_idents(self, sample):
        # NTP: use rs202 gather results for now
        gather_csv = f'{out_dir}/gather/{sample}.gather.csv'
        #gather_csv = f'output.protein-grist-rs202/gather/{sample}.gather.csv'
        assert os.path.exists(gather_csv), "gather output does not exist!?"

        genome_idents = []
        with open(gather_csv, 'rt') as fp:
           r = csv.DictReader(fp)
           for row in r:
               ident = row['name'].split()[0]
               genome_idents.append(ident)
        print(f'loaded {len(genome_idents)} identifiers from {gather_csv}.',
              file=sys.stderr)

        return genome_idents

    def __call__(self, w):
        # get 'sample' from wildcards?
        if self.samples is None:
            return self.do_sample(w)
        else:
            assert not hasattr(w, 'sample'), "if 'samples' provided to constructor, cannot also be in rule inputs"

            ret = []
            for sample in self.samples:
                d = dict(sample=sample)
                w = snakemake.io.Wildcards(fromdict=d)

                x = self.do_sample(w)
                ret.extend(x)

            return ret

    def do_sample(self, w):
        # wait for the results of 'gather_reads_wc'; this will trigger
        # exception until that rule has been run.
        #checkpoints.gather_reads_wc.get(**w)
        
        # wait for the results of 'run_genome_grist_til_gather'; this will trigger
        checkpoints.run_genome_grist_til_gather.get()

        # parse hitlist_genomes,
        genome_idents = self.get_genome_idents(w.sample)

        p = expand(self.pattern, ident=genome_idents, **w)

        return p



# NTP use this later to re-run gather with rs207
#rule run_genome_grist_til_gather:
checkpoint run_genome_grist_til_gather:
    input: "conf/protein-grist.yml",
        #config=config["grist-config"] #os.path.join(out_dir, f"config.grist.{basename}.yml"),
    output: 
        expand(os.path.join(out_dir, "reports/report-{sample}.html"), sample=SAMPLES),
        expand(f'{out_dir}/gather/{{sample}}.gather.csv', sample=SAMPLES),
        #touc:h(f"output.protein-grist/gather/.gather.{{sample}}")   # checkpoints need an output ;)
    log: os.path.join(logs_dir, "genome-grist.log")
    benchmark: os.path.join(logs_dir, "genome-grist.benchmark")
    threads: 32
    resources:
        #mem_mb=lambda wildcards, attempt: attempt*145000
        mem_mb=150000,
        time=400, #240,
        partition="med2",
    conda: "conf/env/grist.yml"
    shell:
        """
        genome-grist run {input} --resources mem_mb={resources.mem_mb} \
                          -j {threads} summarize_gather --nolock \
                          2> {log}
        """

# bbmerge reads for use with prodigal protein mapper
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
        time=400, #240,
        partition="med2",
    params:
        mem = ABUNDTRIM_MEMORY,
    wildcard_constraints:
        sample="/w+"
    log: os.path.join(logs_dir, "bbmerge", "{sample}.log")
    benchmark: os.path.join(logs_dir, "bbmerge", "{sample}.benchmark")
    shell: 
        """
        bbmerge.sh -Xmx25G in={input} \
            out={output.merged} outu={output.unmerged} \
            ihist={output.insert_size_hist} 2> {log}
        """

## TO DO: 
#1. DONE use GTDB proteomes available in 2021-rank-compare/genbank/proteomes and/or 2021-rank-compare/prodigal/
#2. DONE checkpoint rule to read in gather results and load all idents for mapping.
#3. DONE test index and alignment w/paladin
#4. integrate protein mapping and smash gather results (notebook)


def find_proteome(w):
    dl_protdir = "/home/ntpierce/2021-rank-compare/genbank/proteomes"
    prodigal_protdir = "/home/ntpierce/2021-rank-compare/genbank/prodigal"
    filename = f"{w.ident}_protein.faa.gz"
    for pd in [prodigal_protdir, dl_protdir]: # must check prodigal dir first for this to work properly
        fn = glob.glob(os.path.join(pd, filename))
        if len(fn) == 1:
            break
    return fn[0]
   


rule paladin_index_wc:
    input:
        #Checkpoint_ProteomeFiles(f"{out_dir}/proteomes/{{ident}}_protein.faa.gz"),
        #lambda w: "/home/ntpierce/2021-rank-compare/" + proteome_inf.at[f"{w.ident}", "protein_filename"]
        find_proteome
    output:
        idx = out_dir + "/proteomes/{ident}_protein.faa.gz.bwt",
        proteome = out_dir + "/proteomes/{ident}_protein.faa.gz",
    params:
        prot_dir = out_dir + '/proteomes',
        proteome = lambda w: f"{out_dir}/proteomes/{w.ident}_protein.faa.gz",
    resources:
        mem_mb = lambda wildcards, attempt: attempt*3000,
        time=400, #240,
        partition="med2",
    log: os.path.join(logs_dir, "paladin_index", "{ident}.log")
    benchmark: os.path.join(logs_dir, "paladin_index", "{ident}.benchmark")
    conda: 'conf/env/paladin.yml'
    shell:
        """
        mkdir -p {params.prot_dir}
        cp {input} {output.proteome}
        paladin index -r3 {output.proteome} 2> {log}
        """

rule paladin_align_wc:
    input:
        merged_reads= out_dir + '/merge_reads/{sample}.merged.fq.gz',
        index= out_dir + "/proteomes/{ident}_protein.faa.gz.bwt",
    output:
        bam = out_dir + "/protein_mapping/{sample}.x.{ident}.bam",
    params:
        index_base=lambda w: f"{out_dir}/proteomes/{w.ident}_protein.faa.gz"
    resources:
        mem_mb = lambda wildcards, attempt: attempt*15000,
        time=6000,
        partition="med2",
    log: os.path.join(logs_dir, "paladin_align", "{sample}.x.{ident}.paladin-align.log")
    benchmark: os.path.join(logs_dir, "paladin_align", "{sample}.x.{ident}.paladin-align.benchmark")
    threads: 20
    conda: "conf/env/paladin.yml"
    shell:
        """
        paladin align -t {threads} -T 20 {params.index_base} {input.merged_reads} | samtools view -Sb - > {output} 2> {log}
        """

# sort bam 
rule sort_bam:
    input:
        bam = out_dir + "/{dir}/{bam}.bam",
    output:
        sort = out_dir + "/{dir}/sorted/{bam}.bam",
    resources:
        mem_mb = lambda wildcards, attempt: attempt*3000,
        time=90,
        partition="med2",
    log: os.path.join(logs_dir, "sort_bam", "{dir}/{bam}.log")
    benchmark: os.path.join(logs_dir, "sort_bam", "{dir}/{bam}.benchmark")
    conda: "conf/env/minimap2.yml"
    shell: 
        """
        samtools sort {input.bam} -o {output.sort} -O bam 2> {log}
        """

# extract FASTQ from BAM
rule bam_to_fastq_wc:
    input:
        sorted_bam = out_dir + "/protein_mapping/sorted/{bam}.bam",
    output:
        mapped = out_dir + "/protein_mapping/sorted/{bam}.mapped.fq.gz",
    resources:
        mem_mb = lambda wildcards, attempt: attempt*3000,
        time=240,
        partition="med2",
    log: os.path.join(logs_dir, "bam_to_fastq", "{bam}.log")
    benchmark: os.path.join(logs_dir, "bam_to_fastq", "{bam}.benchmark")
    conda: "conf/env/minimap2.yml"
    shell: 
        """
        samtools bam2fq {input.sorted_bam} | gzip > {output.mapped} 2> {log}
        """


# get per-base depth information from BAM
rule bam_to_depth_wc:
    input:
        sorted_bam = out_dir + "/{dir}/sorted/{bam}.bam",
    output:
        depth = out_dir + "/{dir}/{bam}.depth.txt",
    resources:
        mem_mb = lambda wildcards, attempt: attempt*3000,
        time=240,
        partition="med2",
    log: os.path.join(logs_dir, "bam_to_depth", "{dir}/{bam}.log")
    benchmark: os.path.join(logs_dir, "bam_to_depth", "{dir}/{bam}.benchmark")
    conda: "conf/env/minimap2.yml"
    shell: 
        """
        samtools depth -aa {input.sorted_bam} > {output.depth} 2> {log}
        """

# wild card rule for getting _covered_ regions from BAM
rule bam_covered_regions_wc:
    input:
        sorted_bam = out_dir + "/{dir}/sorted/{bam}.bam",
    output:
        regions = out_dir + "/{dir}/{bam}.regions.bed",
    resources:
        mem_mb = lambda wildcards, attempt: attempt*3000,
        time=240,
        partition="med2",
    log: os.path.join(logs_dir, "bam_covered_regions", "{dir}/{bam}.log")
    benchmark: os.path.join(logs_dir, "bam_covered_regions", "{dir}/{bam}.benchmark")
    conda: "conf/env/covtobed.yml"
    shell: 
        """
        covtobed {input.sorted_bam} -l 100 -m 1 | \
            bedtools merge -d 5 -c 4 -o mean > {output.regions} 2> {log}
        """

# calculating SNPs/etc.
rule mpileup_wc:
    input:
        query = out_dir + "/proteomes/{ident}_protein.faa.gz",
        sorted_bam = out_dir + "/{dir}/sorted/{sample}.x.{ident}.bam",
    output:
        bcf = out_dir + "/{dir}/{sample}.x.{ident}.bcf",
        vcf = out_dir + "/{dir}/{sample}.x.{ident}.vcf.gz",
        vcfi = out_dir + "/{dir}/{sample}.x.{ident}.vcf.gz.csi",
    conda: "conf/env/bcftools.yml"
    resources:
        mem_mb = lambda wildcards, attempt: attempt*3000,
        time=240,
        partition="med2",
    log: os.path.join(logs_dir, "mpileup", "{dir}/{sample}.x.{ident}.log")
    benchmark: os.path.join(logs_dir, "mpileup", "{dir}/{sample}.x.{ident}.benchmark")
    shell: 
        """
        proteomefile=$(mktemp -t g.proteome.XXXXXXX)
        gunzip -c {input.query} > $proteomefile
        bcftools mpileup -Ou -f $proteomefile {input.sorted_bam} | bcftools call -mv -Ob -o {output.bcf} 2> {log}
        rm $proteomefile
        bcftools view {output.bcf} | bgzip > {output.vcf} 2>> {log}
        bcftools index {output.vcf} 2>> {log}
        """

# summarize depth into a CSV
rule summarize_samtools_depth_wc:
    input:
        depth = Checkpoint_GatherResults(out_dir + f"/{{dir}}/{{sample}}.x.{{ident}}.depth.txt"),
        vcf = Checkpoint_GatherResults(out_dir + f"/{{dir}}/{{sample}}.x.{{ident}}.vcf.gz"),
    output:
        csv = f"{out_dir}/{{dir}}/{{sample}}.summary.csv"
    resources:
        mem_mb = lambda wildcards, attempt: attempt*3000,
        time=240,
        partition="med2",
    log: os.path.join(logs_dir, "summarize_samtools_depth", "{dir}/{sample}.log")
    benchmark: os.path.join(logs_dir, "summarize_samtools_depth", "{dir}/{sample}.benchmark")
    shell: 
        """
        python summarize_mapping.py {wildcards.sample} \
             {input.depth} -o {output.csv} 2> {log}
        """

# convert mapped reads to protein_leftover reads
rule extract_protein_leftover_reads_wc:
# THIS WRITES {outdir}/mapping/{sample}.x.{ident}.protein_leftover.fq.gz
    input:
        csv = ancient(f'{out_dir}/gather/{{sample}}.gather.csv'),
        mapped = Checkpoint_GatherResults(f"{out_dir}/protein-mapping/sorted/{{sample}}.x.{{ident}}.mapped.fq.gz"),
    output:
        touch(f"{out_dir}/protein_leftover/.protein_leftover.{{sample}}")
    #conda: "conf/env/sourmash.yml" # i think i have sourmash in main env
    params:
        outdir = out_dir,
    shell:
        """
        python -Werror -m genome_grist.subtract_gather \
            {wildcards.sample:q} {input.csv} --outdir={params.out_dir:q}
        """

# rule for mapping protein_leftover reads to genomes -> BAM
rule map_protein_leftover_reads_wc:
    input:
        all_csv = f"{out_dir}/mapping/{{sample}}.summary.csv",
        #query = Checkpoint_GenomeFiles(f"{out_dir}/genomes/{{ident}}_genomic.fna.gz"),
        protein_leftover_reads_flag = f"{out_dir}/protein_leftover/.protein_leftover.{{sample}}",
        index= out_dir + "/proteomes/{ident}_protein.faa.gz.bwt",
    output:
        bam=out_dir + "/protein_leftover/{sample}.x.{ident}.bam",
    params:
        index_base=lambda w: f"{out_dir}/proteomes/{w.ident}_protein.faa.gz"
    conda: "conf/env/paladin.yml"
    threads: 4
    shell: 
       """
       paladin align -t {threads} -T 20 {params.index_base} \
         {out_dir}/protein_mapping/{wildcards.sample}.x.{wildcards.ident}.protein_leftover.fq.gz \
         | samtools view -b -F 4 - | samtools sort - > {output.bam} 2> {log}
       """
 #   shell:
        #minimap2 -ax sr -t {threads} {input.query} \
     #{outdir}/mapping/{wildcards.sample}.x.{wildcards.ident}.protein_leftover.fq.gz | \
     #       samtools view -b -F 4 - | samtools sort - > {output.bam}
 #       """
 #       """

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



