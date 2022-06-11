import pandas as pd

# read in samples
configfile: "/home/ntpierce/2022-np-mgx/conf.yml"
sample_info_csv = config['samples']
sample_info = pd.read_csv(sample_info_csv)
essential_cols = ["name", "read1", "read2"]
#assert essential_cols in sample_info.columns.to_list()

sample_info.set_index("name", inplace=True)
print(sample_info)
SAMPLES= sample_info.index.to_list()
print("samples:", SAMPLES)
out_dir  = config['outdir']
data_dir  = config['datadir']
logs_dir = os.path.join(out_dir, "logs")

ABUNDTRIM_MEMORY = float(config['metagenome_trim_memory'])

SOURMASH_DB_KSIZE = config.get('sourmash_database_ksize', ['10'])
SOURMASH_COMPUTE_SCALED = config.get('sourmash_scaled', '200')
SOURMASH_COMPUTE_TYPE = config.get('sourmash_sigtype', 'protein')
SOURMASH_COMPUTE_KSIZES = config.get('sourmash_compute_ksizes', ['7',' 10'])

rule all:
    input: 
        expand(os.path.join(out_dir, "trim", "{sample}.trim.fq.gz"), sample=SAMPLES),
        expand(os.path.join(out_dir, "abundtrim", "{sample}.abundtrim.fq.gz"), sample=SAMPLES),
        expand(os.path.join(out_dir, "abundtrim", "{sample}.abundtrim.fq.gz.kmer-report.txt"), sample=SAMPLES),
        expand(os.path.join(out_dir, "abundtrim", "{sample}.abundtrim.fq.gz.reads-report.txt"), sample=SAMPLES),
        expand(os.path.join(out_dir, "sigs", "{sample}.abundtrim.sig.zip"), sample=SAMPLES)



# fastp trim
rule trim_adapters_wc:
    input:
        r1 = lambda w: os.path.join(data_dir, sample_info.at[w.sample, "read1"]), #ancient(out_dir + "/raw/{sample}_1.fastq.gz"),
        r2 = lambda w: os.path.join(data_dir, sample_info.at[w.sample, "read2"]),
    output:
        interleaved = protected(out_dir + '/trim/{sample}.trim.fq.gz'),
        json=out_dir + "/trim/{sample}.trim.json",
        html=out_dir + "/trim/{sample}.trim.html",
    conda: 'conf/env/trim.yml'
    threads: 4
    resources:
        mem_mb=5000,
        runtime_min=600,
    shadow: "shallow"
    log: os.path.join(logs_dir, "trim", "{sample}.fastp.log")
    benchmark: os.path.join(logs_dir, "trim", "{sample}.fastp.benchmark")
    shell: 
        """
        fastp --in1 {input.r1} --in2 {input.r2} \
             --detect_adapter_for_pe  --qualified_quality_phred 4 \
             --length_required 25 --correction --thread {threads} \
             --json {output.json} --html {output.html} \
             --low_complexity_filter --stdout | gzip -9 > {output.interleaved} 2> {log}
        """

# abundtrim with khmer
# k-mer abundance trimming
rule kmer_trim_reads_wc:
    input: 
        interleaved = ancient(out_dir + '/trim/{sample}.trim.fq.gz'),
    output:
        protected(out_dir + "/abundtrim/{sample}.abundtrim.fq.gz")
    conda: 'conf/env/trim.yml'
    resources:
        mem_mb = int(ABUNDTRIM_MEMORY / 1e6),
    params:
        mem = ABUNDTRIM_MEMORY,
    log: os.path.join(logs_dir, "abundtrim", "{sample}.kmer_trim.log")
    benchmark: os.path.join(logs_dir, "abundtrim", "{sample}.kmer_trim.benchmark")
    shell: 
        """
            trim-low-abund.py -C 3 -Z 18 -M {params.mem} -V \
            {input.interleaved} -o {output} --gzip 2> {log}
        """


# count k-mers
rule estimate_distinct_kmers_wc:
    message:
        """
        Count distinct k-mers for {wildcards.sample} using 'unique-kmers.py' from the khmer package.
        """
    conda: 'conf/env/trim.yml'
    input:
        out_dir + "/abundtrim/{sample}.abundtrim.fq.gz"
    output:
        report = out_dir + "/abundtrim/{sample}.abundtrim.fq.gz.kmer-report.txt",
    params:
        ksize = SOURMASH_DB_KSIZE,
    log: os.path.join(logs_dir, "abundtrim", "{sample}.estimate_distinct_kmers.log")
    benchmark: os.path.join(logs_dir, "abundtrim", "{sample}.estimate_distinct_kmers.benchmark")
    shell: 
        """
        unique-kmers.py {input} -k {params.ksize} -R {output.report} 2> {log}
        """

# count reads and bases
rule count_trimmed_reads_wc:
    message: 
        """
        Count reads & bp in trimmed file for {wildcards.sample}.
        """
    input:
        out_dir + "/abundtrim/{sample}.abundtrim.fq.gz"
    output:
        report = out_dir + "/abundtrim/{sample}.abundtrim.fq.gz.reads-report.txt",
    # from Matt Bashton, in https://bioinformatics.stackexchange.com/questions/935/fast-way-to-count-number-of-reads-and-number-of-bases-in-a-fastq-file
    #log: os.path.join(logs_dir, "abundtrim", "{sample}.count_trimmed_reads.log")
    benchmark: os.path.join(logs_dir, "abundtrim", "{sample}.count_trimmed_reads.benchmark")
    shell:
        """
        gzip -dc {input} |
             awk 'NR%4==2{{c++; l+=length($0)}}
                  END{{
                        print "n_reads,n_bases"
                        print c","l
                      }}' > {output.report}
        """

# make a `sourmash sketch` -p param string.
def make_param_str(ksizes, scaled):
    ks = [ f'k={k}' for k in ksizes ]
    ks = ",".join(ks)
    return f"{ks},scaled={scaled},abund"

# compute sourmash signature from abundtrim reads
rule smash_abundtrim_wc:
    input:
        metagenome = ancient(out_dir + "/abundtrim/{sample}.abundtrim.fq.gz"),
    output:
        sig = out_dir + "/sigs/{sample}.abundtrim.sig.zip"
    conda: "conf/env/sourmash.yml"
    params:
        param_str = make_param_str(ksizes=SOURMASH_COMPUTE_KSIZES,
                                   scaled=SOURMASH_COMPUTE_SCALED),
        action = "translate" if SOURMASH_COMPUTE_TYPE == "protein" else "dna"
    shell: 
        """
        sourmash sketch {params.action} -p {params.param_str} \
           {input} -o {output} \
           --name {wildcards.sample}
        """
