rule riboswitch_bam:
    input:
        "results/alignments/exogenous_rna/sorted/{sample}_{unit}.bam",
    output:
        bam="results/alignments/exogenous_rna/riboswitch/{sample}_{unit}.bam",
        idx="results/alignments/exogenous_rna/riboswitch/{sample}_{unit}.bam.bai",
    log:
        "results/logs/riboswitch/samtools_view/{sample}_{unit}.log",
    params:
        extra="",  # optional params string
        region="PJY300:454-490",  # optional params string
    threads: 2
    wrapper:
        # "v1.19.2/bio/samtools/view"
        # TODO Update once PR is merged and released
        # https://github.com/snakemake/snakemake-wrappers/pull/940
        "file:../../snakemake-wrappers/bio/samtools/view"


rule riboswitch_snp_freq:
    input:
        bam="results/alignments/exogenous_rna/sorted/{sample}_{unit}.bam",
        ref=input_sample_reference_fasta,
    output:
        "results/alignments/exogenous_rna/riboswitch/{sample}_{unit}.brc.tsv",
    log:
        "results/logs/riboswitch/bam_readcount/{sample}_{unit}.log",
    params:
        region="",
    conda:
        "../envs/bam_readcount.yaml"
    shell:
        "bam-readcount -w1 -f {input.ref:q} {input.bam:q} {params.region} > {output:q} 2> {log:q}"


rule riboswitch_reads:
    input:
        "results/alignments/exogenous_rna/riboswitch/{sample}_{unit}.bam",
    output:
        "results/alignments/exogenous_rna/riboswitch/{sample}_{unit}.1.fq",
        "results/alignments/exogenous_rna/riboswitch/{sample}_{unit}.2.fq",
    log:
        "results/logs/riboswitch/samtools_view/{sample}_{unit}.log",
    params:
        sort="",
        fastq="",
    threads: 3
    wrapper:
        "v1.19.2/bio/samtools/fastq/separate"


rule riboswtich_reads_collapsed:
    input:
        "results/alignments/exogenous_rna/riboswitch/{sample}_{unit}.{read}.fq",
    output:
        "results/alignments/exogenous_rna/riboswitch/{sample}_{unit}.{read}.collapsed.fq",
    log:
        "results/logs/riboswitch/fastx_collapser/{sample}_{unit}.{read}.log",
    conda:
        "../envs/fastx_toolkit.yaml"
    shell:
        """
        if [ -s {input:q} ]
        then
            fastx_collapser -i {input:q} -o {output:q} 2> {log:q}
        else
            echo "Input file \"{input}\" empty" > {log:q}
            touch {output:q}
        fi
        """


rule cross_alignment_check:
    input:
        sample=[
            "results/trimmed/{sample}_{unit}.1.fastq.gz",
            "results/trimmed/{sample}_{unit}.2.fastq.gz",
        ],
        idx=multiext(
            "data/references/{genome}",
            ".1.bt2l",
            ".2.bt2l",
            ".3.bt2l",
            ".4.bt2l",
            ".rev.1.bt2l",
            ".rev.2.bt2l",
        ),
    output:
        bam=temp(
            "results/alignments/exogenous_rna/cross_alignment_check/unsorted/{sample}_{unit}_to_{genome}.bam"
        ),
    log:
        "results/logs/bowtie2/exogenous_rna/cross_alignment_check/{sample}_{unit}_to_{genome}.log",
    params:
        extra="--no-unal",  # optional parameters
    threads: 12  # Use at least two threads
    wrapper:
        "v1.18.3/bio/bowtie2/align"
