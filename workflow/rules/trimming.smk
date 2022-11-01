rule cutadapt:
    input:
        [
            "results/fastq/{sample}_{unit}.1.fastq.gz",
            "results/fastq/{sample}_{unit}.2.fastq.gz",
        ],
    output:
        fastq1="results/trimmed/{sample}_{unit}.1.fastq.gz",
        fastq2="results/trimmed/{sample}_{unit}.2.fastq.gz",
        qc="results/trimmed/{sample}_{unit}.qc",
    params:
        # https://cutadapt.readthedocs.io/en/stable/guide.html#adapter-types
        adapters="-a " + config["adapter-read-1"] + " -A " + config["adapter-read-2"],
        # https://cutadapt.readthedocs.io/en/stable/guide.html#
        extra="--minimum-length 10 -q 20",
    log:
        "results/logs/cutadapt/{sample}_{unit}.log",
    threads: 8  # set desired number of threads here
    resources:
        mem_mb=8000,
    wrapper:
        "v1.18.3/bio/cutadapt/pe"
