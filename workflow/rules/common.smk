def input_fastq(wildcards):
    sequence = units.loc[wildcards.sample, wildcards.unit]
    return {
        "fq1": sequence.fq1,
        "fq2": sequence.fq2,
    }


def input_exogenous_mapping(wildcards):
    unit_info = units.loc[wildcards.sample, wildcards.unit]

    return {
        "sample": [unit_info.fq1, unit_info.fq2],
        "idx": multiext(
            f"data/references/{wildcards.genome}",
            ".1.bt2l",
            ".2.bt2l",
            ".3.bt2l",
            ".4.bt2l",
            ".rev.1.bt2l",
            ".rev.2.bt2l",
        ),
    }


def exogenous_unmapped_fastq(wildcards):
    sample_info = samples.loc[wildcards.sample]

    return {
        "sample": [
            f"results/alignments/{sample_info.pegRNA}/unmapped/{wildcards.sample}_{wildcards.unit}.1.fq.gz",
            f"results/alignments/{sample_info.pegRNA}/unmapped/{wildcards.sample}_{wildcards.unit}.2.fq.gz",
        ],
    }


rule unzip:
    input:
        "data/references/{genome}.fa.gz",
    output:
        "data/references/{genome}.fa",
    log:
        "results/logs/unzip/{genome}.txt",
    conda:
        "../envs/gzip.yaml"
    shell:
        "gzip -dc {input:q} > {output:q} 2> {log:q}"
