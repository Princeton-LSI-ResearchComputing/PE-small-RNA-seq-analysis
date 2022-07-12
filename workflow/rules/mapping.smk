rule fastp_pe:
    input:
        sample=[
            "results/fastq/{sample}_{unit}.1.fastq.gz",
            "results/fastq/{sample}_{unit}.2.fastq.gz",
        ],
    output:
        trimmed=[
            "results/trimmed/{sample}_{unit}.1.fastq.gz",
            "results/trimmed/{sample}_{unit}.2.fastq.gz",
        ],
        unpaired="results/trimmed/{sample}_{unit}.singletons.fastq.gz",
        failed="results/trimmed/{sample}_{unit}.failed.fastq.gz",
        html="results/report/fastp/{sample}_{unit}.html",
        json="results/report/fastp/{sample}_{unit}.fastp.json",
    log:
        "results/logs/fastp/{sample}_{unit}.log",
    params:
        adapters="--detect_adapter_for_pe",
    threads: 2
    wrapper:
        "v1.7.0/bio/fastp"


rule join_references:
    input:
        targets="data/references/targets.fa",
        reference=config["reference"],
    output:
        "data/references/grch38_p12_targets.fa",
    log:
        "results/logs/join_references.log",
    conda:
        "../envs/fastx_toolkit.yaml"
    shell:
        "cat {input.reference:q} {input.targets:q} | fasta_formatter > {output:q} &2> {log:q}"


rule samtools_fasta_index:
    input:
        "data/references/grch38_p12_targets.fa",
    output:
        "data/references/grch38_p12_targets.fa.fai",
    log:
        "results/logs/faidx/grch38_p12_targets.log",
    params:
        extra="",  # optional params string
    wrapper:
        "v1.7.0/bio/samtools/faidx"


rule bowtie2_build_large:
    input:
        ref="data/references/grch38_p12_targets.fa",
    output:
        multiext(
            "data/references/grch38_p12_targets",
            ".1.bt2l",
            ".2.bt2l",
            ".3.bt2l",
            ".4.bt2l",
            ".rev.1.bt2l",
            ".rev.2.bt2l",
        ),
    log:
        "results/logs/bowtie2_build/build.log",
    params:
        extra="--large-index",  # optional parameters
    threads: 12
    resources:
        time="06:00:00",
    wrapper:
        "v1.7.0/bio/bowtie2/build"


rule bowtie2:
    input:
        sample=[
            "results/trimmed/{sample}_{unit}.1.fastq.gz",
            "results/trimmed/{sample}_{unit}.2.fastq.gz",
        ],
        idx=multiext(
            "data/references/grch38_p12_targets",
            ".1.bt2l",
            ".2.bt2l",
            ".3.bt2l",
            ".4.bt2l",
            ".rev.1.bt2l",
            ".rev.2.bt2l",
        ),
    output:
        "results/mapped/{sample}_{unit}.bam",
    log:
        "results/logs/bowtie2/{sample}_{unit}.log",
    params:
        extra="",  # optional parameters
    threads: 8  # Use at least two threads
    wrapper:
        "v1.7.0/bio/bowtie2/align"
