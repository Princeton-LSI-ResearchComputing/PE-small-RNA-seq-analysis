def input_fastp_pe(wildcards):
    sequence = units.loc[wildcards.sample, wildcards.unit]
    return [
        sequence.fq1,
        sequence.fq2,
    ]


rule fastp_pe:
    input:
        sample=input_fastp_pe,
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
        "v1.5.0/bio/fastp"


rule join_references:
    input:
        targets="data/references/targets.fa",
        reference=config["reference"],
    output:
        "data/references/grch38_p12_targets.fa",
    shell:
        "cat {input} > {output}"


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
    threads: 8
    wrapper:
        "v1.5.0/bio/bowtie2/build"


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
        "v1.5.0/bio/bowtie2/align"
