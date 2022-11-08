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
        "v1.18.3/bio/fastqc"


rule samtools_stats:
    input:
        bam="results/alignments/{genome}/sorted/{sample}_{unit}.bam",
    output:
        "results/samtools_stats/{genome}/{sample}_{unit}.txt",
    params:
        extra="",  # Optional: extra arguments.
        region="",  # Optional: region string.
    log:
        "results/logs/samtools_stats/{genome}/{sample}_{unit}.log",
    wrapper:
        "v1.18.3/bio/samtools/stats"


rule samtools_index:
    input:
        "results/alignments/{genome}/sorted/{sample}_{unit}.bam",
    output:
        "results/alignments/{genome}/sorted/{sample}_{unit}.bam.bai",
    log:
        "results/logs/samtools_index/{genome}/{sample}_{unit}.log",
    params:
        extra="",  # optional params string
    threads: 4  # This value - 1 will be sent to -@
    wrapper:
        "v1.18.3/bio/samtools/index"


rule samtools_idxstats:
    input:
        bam="results/alignments/{genome}/sorted/{sample}_{unit}.bam",
        bai="results/alignments/{genome}/sorted/{sample}_{unit}.bam.bai",
    output:
        "results/samtools_idxstats/{genome}/{sample}_{unit}.bam.idxstats",
    log:
        "results/logs/samtools_idxstats/{genome}/{sample}_{unit}.log",
    params:
        extra="",  # optional params string
    wrapper:
        "v1.18.3/bio/samtools/idxstats"


rule deeptools_bampefragmentsize:
    input:
        bam=expand(
            "results/alignments/{{genome}}/sorted/{unit.sample_name}_{unit.unit_name}.bam",
            unit=units.itertuples(),
        ),
        bai=expand(
            "results/alignments/{{genome}}/sorted/{unit.sample_name}_{unit.unit_name}.bam.bai",
            unit=units.itertuples(),
        ),
    output:
        table="results/deeptools_bampefragmentsize/{genome}/fragment_size_table.tsv",
        raw="results/deeptools_bampefragmentsize/{genome}/fragment_size_raw_lengths.tsv",
        summary="results/deeptools_bampefragmentsize/{genome}/fragment_size_summary.txt",
    log:
        "results/logs/deeptools_bampefragmentsize/{genome}/deeptools_bampefragmentsize.log",
    threads: 12
    conda:
        "../envs/deeptools.yaml"
    shell:
        """
        bamPEFragmentSize \
            -T "Fragment size of PE RNA-Seq data" \
            --table {output.table:q} \
            --outRawFragmentLengths {output.raw:q} \
            --numberOfProcessors {threads} \
            -b {input.bam:q} \
            > {output.summary:q} \
            2> {log:q}
        """


rule multiqc_library:
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
    output:
        "results/qc/multiqc_library.html",
    params:
        "",
    log:
        "results/logs/multiqc_library.log",
    wrapper:
        "v1.18.3/bio/multiqc"


rule multiqc_mapped:
    input:
        expand(
            "results/logs/bowtie2/{{genome}}/{unit.sample_name}_{unit.unit_name}.log",
            unit=units.itertuples(),
        ),
        expand(
            "results/samtools_stats/{{genome}}/{unit.sample_name}_{unit.unit_name}.txt",
            unit=units.itertuples(),
        ),
        expand(
            "results/samtools_idxstats/{{genome}}/{unit.sample_name}_{unit.unit_name}.bam.idxstats",
            unit=units.itertuples(),
        ),
        "results/deeptools_bampefragmentsize/{genome}/fragment_size_table.tsv",
        "results/deeptools_bampefragmentsize/{genome}/fragment_size_raw_lengths.tsv",
    output:
        "results/qc/multiqc_{genome}.html",
    params:
        "",
    log:
        "results/logs/multiqc_{genome}.log",
    wrapper:
        "v1.18.3/bio/multiqc"
