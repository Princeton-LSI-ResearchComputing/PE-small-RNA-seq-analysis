rule bowtie2_build_large:
    input:
        ref="data/references/grch38_pjy103_pjy300.fa",
    output:
        multiext(
            "data/references/grch38_pjy103_pjy300",
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
        mem_mb=16000,
    wrapper:
        "v1.18.3/bio/bowtie2/build"


rule bowtie2:
    input:
        sample=[
            "results/trimmed/{sample}_{unit}.1.fastq.gz",
            "results/trimmed/{sample}_{unit}.2.fastq.gz",
        ],
        idx=multiext(
            "data/references/grch38_pjy103_pjy300",
            ".1.bt2l",
            ".2.bt2l",
            ".3.bt2l",
            ".4.bt2l",
            ".rev.1.bt2l",
            ".rev.2.bt2l",
        ),
    output:
        temp("results/mapped/{sample}_{unit}.bam"),
    log:
        "results/logs/bowtie2/{sample}_{unit}.log",
    params:
        extra="",  # optional parameters
    threads: 12  # Use at least two threads
    wrapper:
        "v1.18.3/bio/bowtie2/align"


rule samtools_sort:
    input:
        "results/mapped/{sample}_{unit}.bam",
    output:
        "results/mapped_sorted/{sample}_{unit}.bam",
    log:
        "results/samtools_sort/{sample}_{unit}.log",
    params:
        extra="",
    threads: 8
    resources:
        mem_mb=16000,
    wrapper:
        "v1.18.3/bio/samtools/sort"
