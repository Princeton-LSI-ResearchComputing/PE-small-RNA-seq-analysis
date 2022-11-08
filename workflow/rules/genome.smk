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


rule unzip:
    input:
        "data/references/{genome}.fa.gz",
    output:
        "data/references/{genome}.fa",
    log:
        "results/logs/unzip/{genome}.txt",
    conda:
        "../envs/gzip.yaml"
    shell:
        "gzip -dc {input:q} > {output:q} 2> {log:q}"


rule samtools_fasta_index:
    input:
        "data/references/{genome}.fa",
    output:
        "data/references/{genome}.fa.fai",
    log:
        "results/logs/faidx/{genome}.faidx.log",
    params:
        extra="",  # optional params string
    wrapper:
        "v1.18.3/bio/samtools/faidx"
