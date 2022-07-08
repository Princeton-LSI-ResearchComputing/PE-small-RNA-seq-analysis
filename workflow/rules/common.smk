def input_fastq(wildcards):
    sequence = units.loc[wildcards.sample, wildcards.unit]
    return {
        "fq1": sequence.fq1,
        "fq2": sequence.fq2,
    }
