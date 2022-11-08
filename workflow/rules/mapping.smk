rule bowtie2_build_large:
    input:
        ref="data/references/{genome}.fa",
    output:
        multiext(
            "data/references/{genome}",
            ".1.bt2l",
            ".2.bt2l",
            ".3.bt2l",
            ".4.bt2l",
            ".rev.1.bt2l",
            ".rev.2.bt2l",
        ),
    log:
        "results/logs/bowtie2_build/{genome}.log",
    params:
        extra="--large-index",  # optional parameters
    threads: 12
    resources:
        time="06:00:00",
        mem_mb=16000,
    wrapper:
        "v1.18.3/bio/bowtie2/build"


rule bowtie2_exogenous_rna:
    input:
        unpack(input_exogenous_mapping),
    output:
        bam=temp("results/alignments/{genome}/unsorted/{sample}_{unit}.bam"),
        unmapped_fq1="results/alignments/{genome}/unmapped/{sample}_{unit}.1.fq.gz",
        unmapped_fq2="results/alignments/{genome}/unmapped/{sample}_{unit}.2.fq.gz",
    log:
        "results/logs/bowtie2/{genome}/{sample}_{unit}.log",
    params:
        extra=lambda wildcards, output: f"--un-conc-gz {output['unmapped_fq1'][:-8]}.%.fq.gz",  # optional parameters
    wildcard_constraints:
        genome="PJY\d+",
    threads: 12  # Use at least two threads
    wrapper:
        "v1.18.3/bio/bowtie2/align"


rule bowtie2_hg38:
    input:
        unpack(exogenous_unmapped_fastq),
        idx=multiext(
            "data/references/Homo_sapiens.GRCh38.dna.primary_assembly",
            ".1.bt2l",
            ".2.bt2l",
            ".3.bt2l",
            ".4.bt2l",
            ".rev.1.bt2l",
            ".rev.2.bt2l",
        ),
    output:
        bam=temp(
            "results/alignments/Homo_sapiens.GRCh38.dna.primary_assembly/unsorted/{sample}_{unit}.bam"
        ),
        unmapped_fq1="results/alignments/Homo_sapiens.GRCh38.dna.primary_assembly/unmapped/{sample}_{unit}.1.fq.gz",
        unmapped_fq2="results/alignments/Homo_sapiens.GRCh38.dna.primary_assembly/unmapped/{sample}_{unit}.2.fq.gz",
    log:
        "results/logs/bowtie2/Homo_sapiens.GRCh38.dna.primary_assembly/{sample}_{unit}.log",
    params:
        extra=lambda wildcards, output: f"--un-conc-gz {output['unmapped_fq1'][:-8]}.%.fq.gz",  # optional parameters
    threads: 12  # Use at least two threads
    wrapper:
        "v1.18.3/bio/bowtie2/align"


rule samtools_sort:
    input:
        "results/alignments/{genome}/unsorted/{sample}_{unit}.bam",
    output:
        "results/alignments/{genome}/sorted/{sample}_{unit}.bam",
    log:
        "results/samtools-sort/{genome}-{sample}_{unit}.log",
    params:
        extra="",
    threads: 8
    resources:
        mem_mb=16000,
    wrapper:
        "v1.18.3/bio/samtools/sort"
