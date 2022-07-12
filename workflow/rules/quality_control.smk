rule rename_fastq:
    input:
        unpack(input_fastq),
    output:
        fq1="results/fastq/{sample}_{unit}.1.fastq.gz",
        fq2="results/fastq/{sample}_{unit}.2.fastq.gz",
    log:
        fq1="results/logs/fastq/{sample}_{unit}.1.log",
        fq2="results/logs/fastq/{sample}_{unit}.2.log",
    conda:
        "../envs/coreutils.yaml"
    threads: 1
    shell:
        """
        ln -sr {input.fq1:q} {output.fq1:q} &> {log.fq1:q}
        ln -sr {input.fq2:q} {output.fq2:q} &> {log.fq2:q}
        """


rule fastqc:
    input:
        "results/fastq/{sample}_{unit}.{readnum}.fastq.gz",
    output:
        html="results/qc/fastqc/{sample}_{unit}.{readnum}.html",
        zip="results/qc/fastqc/{sample}_{unit}.{readnum}_fastqc.zip",
    log:
        "results/logs/fastqc/{sample}_{unit}.{readnum}.log",
    threads: 1
    params:
        "--quiet",
    wrapper:
        "v1.7.0/bio/fastqc"


rule samtools_stats:
    input:
        bam="results/mapped/{sample}_{unit}.bam",
    output:
        "results/samtools_stats/{sample}_{unit}.txt",
    params:
        extra="",  # Optional: extra arguments.
        region="",  # Optional: region string.
    log:
        "results/logs/samtools_stats/{sample}_{unit}.log",
    wrapper:
        "v1.7.0/bio/samtools/stats"


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
        "v1.7.0/bio/samtools/sort"


rule samtools_index:
    input:
        "results/mapped_sorted/{sample}_{unit}.bam",
    output:
        "results/mapped_sorted/{sample}_{unit}.bam.bai",
    log:
        "results/logs/samtools_index/{sample}_{unit}.log",
    params:
        extra="",  # optional params string
    threads: 4  # This value - 1 will be sent to -@
    wrapper:
        "v1.7.0/bio/samtools/index"


rule samtools_idxstats:
    input:
        bam="results/mapped_sorted/{sample}_{unit}.bam",
        bai="results/mapped_sorted/{sample}_{unit}.bam.bai",
    output:
        "results/samtools_idxstats/{sample}_{unit}.bam.idxstats",
    log:
        "results/logs/samtools_idxstats/{sample}_{unit}.log",
    params:
        extra="",  # optional params string
    wrapper:
        "v1.7.0/bio/samtools/idxstats"


rule multiqc:
    input:
        expand(
            "results/qc/fastqc/{unit.sample_name}_{unit.unit_name}.{readnum}.html",
            unit=units.itertuples(),
            readnum=(1, 2),
        ),
        expand(
            "results/trimmed/{unit.sample_name}_{unit.unit_name}.qc",
            unit=units.itertuples(),
        ),
        expand(
            "results/logs/bowtie2/{unit.sample_name}_{unit.unit_name}.log",
            unit=units.itertuples(),
        ),
        expand(
            "results/samtools_stats/{unit.sample_name}_{unit.unit_name}.txt",
            unit=units.itertuples(),
        ),
        expand(
            "results/samtools_idxstats/{unit.sample_name}_{unit.unit_name}.bam.idxstats",
            unit=units.itertuples(),
        ),
    output:
        "results/qc/multiqc.html",
    params:
        "",
    log:
        "results/logs/multiqc.log",
    wrapper:
        "v1.7.0/bio/multiqc"
