def input_fastq(wildcards):
    sequence = units.loc[wildcards.sample, wildcards.unit]
    return {
        "fq1": sequence.fq1,
        "fq2": sequence.fq2,
    }


def input_exogenous_mapping(wildcards):
    unit_info = units.loc[wildcards.sample, wildcards.unit]
    sample_info = samples.loc[wildcards.sample]

    return {
        "sample": [
            f"results/trimmed/{wildcards.sample}_{wildcards.unit}.1.fastq.gz",
            f"results/trimmed/{wildcards.sample}_{wildcards.unit}.2.fastq.gz",
        ],
        "idx": multiext(
            f"data/references/{sample_info.exogenous_rna}",
            ".1.bt2l",
            ".2.bt2l",
            ".3.bt2l",
            ".4.bt2l",
            ".rev.1.bt2l",
            ".rev.2.bt2l",
        ),
    }


def input_sample_reference_fasta(wildcards):
    sample_info = samples.loc[wildcards.sample]

    return f"data/references/{sample_info.exogenous_rna}.fa"


def exogenous_unmapped_fastq(wildcards):
    sample_info = samples.loc[wildcards.sample]

    return {
        "sample": [
            f"results/alignments/{sample_info.exogenous_rna}/unmapped/{wildcards.sample}_{wildcards.unit}.1.fq.gz",
            f"results/alignments/{sample_info.exogenous_rna}/unmapped/{wildcards.sample}_{wildcards.unit}.2.fq.gz",
        ],
    }
