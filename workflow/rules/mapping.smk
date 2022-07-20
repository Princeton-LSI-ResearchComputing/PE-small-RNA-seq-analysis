rule bowtie2_build_large:
    input:
        ref="data/references/grch38_pjy142_pjy151_pjy209.fa",
    output:
        multiext(
            "data/references/grch38_pjy142_pjy151_pjy209",
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
            "data/references/grch38_pjy142_pjy151_pjy209",
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
        "v1.7.0/bio/bowtie2/align"
