import os
from snakemake.remote.FTP import RemoteProvider as FTPRemoteProvider

FTP = FTPRemoteProvider()


rule fetch_grch38:
    input:
        # only keeping the file so we can move it out to the cwd
        FTP.remote(
            "http://ftp.ensembl.org/pub/release-107/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz",
            keep_local=True,
        ),
    output:
        "data/references/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz",
    log:
        "results/logs/fetch_grch38.log",
    run:
        shell("mv {input:q} {output:q} &> {log:q}")


rule join_references:
    input:
        pjy103="data/references/PJY103.fa",
        pjy300="data/references/PJY300.fa",
        grch38="data/references/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz",
    output:
        "data/references/grch38_pjy103_pjy300.fa",
    log:
        "results/logs/join_references.log",
    conda:
        "../envs/fastx_toolkit.yaml"
    shell:
        """
        zcat {input.grch38:q} | \
        cat - {input.pjy103:q} {input.pjy300:q} | \
        fasta_formatter -o {output:q} \
        2> {log:q}
        """


rule samtools_fasta_index:
    input:
        "data/references/grch38_pjy103_pjy300.fa",
    output:
        "data/references/grch38_pjy103_pjy300.fa.fai",
    log:
        "results/logs/faidx/grch38_pjy103_pjy300_faindex.log",
    params:
        extra="",  # optional params string
    wrapper:
        "v1.18.3/bio/samtools/faidx"
