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
        "sample": [unit_info.fq1, unit_info.fq2],
        "idx": multiext(
            f"data/references/{sample_info.pegRNA}",
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

    # TODO Can use standardized names here, no need for function
    return {
        "sample": [
            f"results/alignments/{sample_info.pegRNA}/unmapped/{wildcards.sample}_{wildcards.unit}.1.fq.gz",
            f"results/alignments/{sample_info.pegRNA}/unmapped/{wildcards.sample}_{wildcards.unit}.2.fq.gz",
        ],
    }
