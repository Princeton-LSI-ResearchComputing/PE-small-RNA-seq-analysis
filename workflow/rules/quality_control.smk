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


rule multiqc:
    input:
        expand(
            "results/qc/fastqc/{unit.sample_name}_{unit.unit_name}.{readnum}.html",
            unit=units.itertuples(),
            readnum=(1, 2),
        ),
        expand(
            "results/report/fastp/{unit.sample_name}_{unit.unit_name}.fastp.json",
            unit=units.itertuples(),
        ),
        expand(
            "results/logs/bowtie2/{unit.sample_name}_{unit.unit_name}.log",
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
