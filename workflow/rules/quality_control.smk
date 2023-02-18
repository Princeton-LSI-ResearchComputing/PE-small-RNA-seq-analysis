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


rule fastqc_unmapped:
    input:
        "results/alignments/Homo_sapiens.GRCh38.dna.primary_assembly/unmapped/{sample}_{unit}.{readnum}.fq.gz",
    output:
        html="results/qc/fastqc_unmapped/{sample}_{unit}.{readnum}.html",
        zip="results/qc/fastqc_unmapped/{sample}_{unit}.{readnum}_fastqc.zip",
    log:
        "results/logs/fastqc_unmapped/{sample}_{unit}.{readnum}.log",
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
        "results/alignments/{genome}/{format}/{sample}_{unit}.bam",
    output:
        "results/alignments/{genome}/{format}/{sample}_{unit}.bam.bai",
    log:
        "results/logs/samtools_index/{genome}/{format}/{sample}_{unit}.log",
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


rule insert_size:
    input:
        "results/alignments/{genome}/sorted/{sample}_{unit}.bam",
    output:
        txt="results/picard_insert_size/{genome}/{sample}_{unit}.isize.txt",
        pdf="results/picard_insert_size/{genome}/{sample}_{unit}.isize.pdf",
    log:
        "results/logs/picard_insert_size/{genome}/{sample}_{unit}.log",
    params:
        # optional parameters (e.g. relax checks as below)
        extra="--VALIDATION_STRINGENCY LENIENT --METRIC_ACCUMULATION_LEVEL null --METRIC_ACCUMULATION_LEVEL SAMPLE",
    # optional specification of memory usage of the JVM that snakemake will respect with global
    # resource restrictions (https://snakemake.readthedocs.io/en/latest/snakefiles/rules.html#resources)
    # and which can be used to request RAM during cluster job submission as `{resources.mem_mb}`:
    # https://snakemake.readthedocs.io/en/latest/executing/cluster.html#job-properties
    resources:
        mem_mb=1024,
    wrapper:
        "v1.23.1/bio/picard/collectinsertsizemetrics"


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
        expand(
            "results/picard_insert_size/{{genome}}/{unit.sample_name}_{unit.unit_name}.isize.txt",
            unit=units.itertuples(),
        ),
    output:
        "results/qc/multiqc_{genome}.html",
    params:
        "",
    log:
        "results/logs/multiqc_{genome}.log",
    wrapper:
        "v1.18.3/bio/multiqc"


rule multiqc_featurecounts:
    input:
        "results/featureCounts/noncoding_exon_type.featureCounts.summary",
    output:
        "results/qc/multiqc_featurecounts.html",
    params:
        "",
    log:
        "results/logs/multiqc_featurecounts.log",
    wrapper:
        "v1.18.3/bio/multiqc"
