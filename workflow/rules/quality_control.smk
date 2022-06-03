rule fastqc:
    input: "data/from_sequencer/{sample}.fastq.gz"
    output:
        html="results/qc/fastqc/{sample}.html",
        zip="results/qc/fastqc/{sample}_fastqc.zip"
    log: "results/logs/fastqc/{sample}.log"
    threads: 1
    params: '--quiet'
    wrapper: "v1.5.0/bio/fastqc"

def multiqc_input(wildcards):
    samples = glob_wildcards('data/from_sequencer/{sample}.fastq.gz').sample
    return expand("results/qc/fastqc/{sample}.html", sample=samples)

rule multiqc:
    input: multiqc_input
    output: "results/qc/multiqc.html"
    params: ""
    log: "results/logs/multiqc.log"
    wrapper: "v1.5.0/bio/multiqc"
