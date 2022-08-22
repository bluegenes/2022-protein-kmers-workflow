## simplest approach:
# using sourmash gather results and
# using pre-downloaded /generated gtdb proteomes,
# index proteomes
# merge reads using bbmerge/pear
# map reads --> proteomes
# generate summary stats
# generate grist-style plots

import csv
import pandas as pd
import glob
configfile: "conf/prot-mapping.yml"
SAMPLES = config["samples"]
out_dir = config.get('output_dir', 'output.pgather-vs-pmapping')
logs_dir = os.path.join(out_dir, "logs")
benchmarks_dir = os.path.join(out_dir, "benchmarks")

basename = config.get("basename", 'phhz')

sample_info = pd.read_csv(config['sample_info'])
SAMPLES = sample_info["name"].to_list()
# hacky, but useful: drop trimmed read names into sample_info dict for easier finding
#sample_info['trim'] = f'{out_dir}/trim/' + sample_info['name'] + '.trim.fq.gz'
#sample_info['abundtrim'] = f'{out_dir}/abundtrim/' + sample_info['name'] + '.abundtrim.fq.gz'
# set name as index for easy access
sample_info.set_index('name', inplace=True)

#search_databases = config['search_databases'] # must be dictionary
gather_dir = config["gather_dir"]
gather_alpha = config.get("alphabet", 'protein')
gather_scaled = config.get("scaled", 200)
gather_ksize = config.get("ksize", 10)
out_dir = out_dir + f".{gather_alpha}-k{gather_ksize}-sc{gather_scaled}"

onstart:
    print("------------------------------------------------------------------------------------------")
    print("  compare sourmash protein gather to protein mapping for metagenomic taxonomic profiling")
    print("------------------------------------------------------------------------------------------")

onsuccess:
    print("--------------------")
    print("workflow completed!")
    print("--------------------")


# make dictionary of gtdb proteomes
proteome_inf = pd.read_csv("/home/ntpierce/2021-rank-compare/output.rank-compare/gtdb-rs207.fromfile.csv")
proteome_inf.set_index('ident', inplace=True)
phylodb_inf = pd.read_csv("/group/ctbrowngrp/sourmash-db/phylodb/phylodb_1.076.fromfile.wident.csv")
phylodb_inf.set_index('ident', inplace=True)
# join these dicts
proteomes= pd.concat([proteome_inf, phylodb_inf])

rule all:
    input: 
        # read mapping outputs
        #expand(out_dir + '/merge_reads/{sample}.merged.fq.gz', sample=SAMPLES),
        expand(out_dir + "/{dir}/{sample}.summary.csv", sample=SAMPLES, dir=['protein_mapping']), #, 'protein_leftover']),
        #expand(f"{out_dir}/protein_leftover/{{sample}}.summary.csv", sample=SAMPLES),

# use gather results from separate workflow - rule instead of checkpoint? 
#rule copy_reads_gather:
checkpoint copy_reads_gather:
    input:
        #gather_csv = lambda w: sample_info.at[w.sample, 'gather_csv']
        gather_csv= lambda w: os.path.join(gather_dir, 'abund-gather', f"{w.sample}.{gather_alpha}-k{gather_ksize}-sc{gather_scaled}.gather.csv")
        #gather_csv = f'output.protein-grist-rs202/gather/{{sample}}.gather.csv'
        #gather_csv = f'{out_dir}/gather/{{sample}}.gather.csv'
    output: 
        os.path.join(out_dir, "gather", "{sample}.gather.csv"),
        #touch(f"{out_dir}/gather/.gather.{{sample}}")   # checkpoints need an output ;)
    shell:
        """
        cp {input} {output}
        """
   

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
        gather_csv=os.path.join(out_dir, 'gather', f'{sample}.gather.csv')
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
        # wait for the results of 'copy_reads_gather'; this will trigger
        # exception until that rule has been run.
        checkpoints.copy_reads_gather.get(**w)

        # parse hitlist_genomes,
        genome_idents = self.get_genome_idents(w.sample)

        p = expand(self.pattern, ident=genome_idents, **w)

        return p


# bbmerge reads for use with prodigal protein mapper
rule bbmerge_paired_reads:
    input:
        interleaved = ancient(gather_dir + '/trim/{sample}.trim.fq.gz'),
    output:
        merged = protected(out_dir + '/merge_reads/{sample}.merged.fq.gz'),
        unmerged = protected(out_dir + '/merge_reads/{sample}.unmerged.fq.gz'),
        insert_size_hist = protected(out_dir + '/merge_reads/{sample}.insert-size-histogram.txt'),
    conda: 'conf/env/bbmap.yml'
    threads: 6
    resources:
        mem_mb = int(20e9 / 1e6),
        time=240,
        partition="low2",
    params:
        mem = 20e9,
    #wildcard_constraints:
        #sample="/w+"
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
    phylodb_protdir = "/group/ctbrowngrp/sourmash-db/phylodb/phylodb_1.076_taxsplit.ff"
    filename = f"{w.ident}_protein.faa.gz"
    phydb_filename = f"{w.ident}.fa"
    for pd in [prodigal_protdir, dl_protdir]: # must check prodigal dir first for this to work properly
        fn = glob.glob(os.path.join(pd, filename))
        if len(fn) == 1:
            return fn[0]
    if not fn:
        fn = glob.glob(os.path.join(phylodb_protdir, phydb_filename))
        if len(fn) == 1:
            return fn[0]
   
rule paladin_index:
    input:
        find_proteome
    output:
        idx = out_dir + "/proteomes/{ident}_protein.faa.gz.bwt",
        proteome = out_dir + "/proteomes/{ident}_protein.faa.gz",
    params:
        prot_dir = out_dir + '/proteomes',
        proteome = lambda w: f"{out_dir}/proteomes/{w.ident}_protein.faa.gz",
    resources:
        mem_mb = lambda wildcards, attempt: attempt*3000,
        time=240,
        partition="low2",
    log: os.path.join(logs_dir, "paladin_index", "{ident}.log")
    benchmark: os.path.join(logs_dir, "paladin_index", "{ident}.benchmark")
    conda: 'conf/env/paladin.yml'
    shell:
        """
        mkdir -p {params.prot_dir}
        cp {input} {output.proteome}
        paladin index -r3 {output.proteome} 2> {log}
        """

rule paladin_align:
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
        paladin align -t {threads} -T 20 {params.index_base} {input.merged_reads} \
                      | samtools view -b -F 4 - | samtools sort  - > {output} 2> {log}
        """

# extract FASTQ from BAM
rule bam_to_fastq:
    input:
        bam = out_dir + "/protein_mapping/{bam}.bam",
    output:
        mapped = out_dir + "/protein_mapping/{bam}.mapped.fq.gz",
    resources:
        mem_mb = lambda wildcards, attempt: attempt*3000,
        time=240,
        partition="med2",
    log: os.path.join(logs_dir, "bam_to_fastq", "{bam}.log")
    benchmark: os.path.join(logs_dir, "bam_to_fastq", "{bam}.benchmark")
    conda: "conf/env/minimap2.yml"
    shell: 
        """
        samtools bam2fq {input.bam} | gzip > {output.mapped} 2> {log}
        """


# get per-base depth information from BAM
rule bam_to_depth:
    input:
        sorted_bam = out_dir + "/{dir}/{bam}.bam",
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
rule bam_covered_regions:
    input:
        sorted_bam = out_dir + "/{dir}/{bam}.bam",
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
rule mpileup:
    input:
        query = out_dir + "/proteomes/{ident}_protein.faa.gz",
        sorted_bam = out_dir + "/{dir}/{sample}.x.{ident}.bam",
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
rule summarize_samtools_depth:
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
#rule extract_protein_leftover_reads:
# THIS WRITES {outdir}/mapping/{sample}.x.{ident}.protein_leftover.fq.gz
#    input:
#        csv = ancient(f'{out_dir}/gather/{{sample}}.gather.csv'),
#        mapped = Checkpoint_GatherResults(f"{out_dir}/protein-mapping/{{sample}}.x.{{ident}}.mapped.fq.gz"),
#        # out_dir + "/protein_mapping/{bam}.mapped.fq.gz"
#    output:
#        touch(f"{out_dir}/protein_leftover/.protein_leftover.{{sample}}")
#    conda: "conf/env/sourmash.yml"
#    params:
#        outdir = out_dir,
#    shell:
#        """
#        python -Werror -m genome_grist.subtract_gather \
#            {wildcards.sample:q} {input.csv} --outdir={params.outdir:q}
#        """
#
## rule for mapping protein_leftover reads to genomes -> BAM
#rule map_protein_leftover_reads:
#    input:
#        all_csv = f"{out_dir}/protein_mapping/{{sample}}.summary.csv",
#        #query = Checkpoint_GenomeFiles(f"{out_dir}/genomes/{{ident}}_genomic.fna.gz"),
#        protein_leftover_reads_flag = f"{out_dir}/protein_leftover/.protein_leftover.{{sample}}",
#        index= out_dir + "/proteomes/{ident}_protein.faa.gz.bwt",
#    output:
#        bam=out_dir + "/protein_leftover/{sample}.x.{ident}.bam",
#    params:
#        index_base=lambda w: f"{out_dir}/proteomes/{w.ident}_protein.faa.gz"
#    conda: "conf/env/paladin.yml"
#    threads: 4
#    shell: 
#        """
#        paladin align -t {threads} -T 20 {params.index_base} \
#         {out_dir}/protein_mapping/{wildcards.sample}.x.{wildcards.ident}.protein_leftover.fq.gz \
#         | samtools view -b -F 4 - | samtools sort - > {output.bam} 2> {log}
#        """
 #   shell:
        #minimap2 -ax sr -t {threads} {input.query} \
     #{outdir}/mapping/{wildcards.sample}.x.{wildcards.ident}.protein_leftover.fq.gz | \
     #       samtools view -b -F 4 - | samtools sort - > {output.bam}
 #       """
 #       """
