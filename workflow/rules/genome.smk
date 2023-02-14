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
        "results/logs/references/fetch_grch38.log",
    run:
        shell("mv {input:q} {output:q} &> {log:q}")


rule fetch_grch38_gtf:
    input:
        # only keeping the file so we can move it out to the cwd
        FTP.remote(
            "https://ftp.ensembl.org/pub/release-108/gtf/homo_sapiens/Homo_sapiens.GRCh38.108.gtf.gz",
            keep_local=True,
        ),
    output:
        "data/references/Homo_sapiens.GRCh38.108.gtf.gz",
    log:
        "results/logs/reference/fetch_grch38_gff.log",
    run:
        shell("mv {input:q} {output:q} &> {log:q}")


rule fetch_rna_central_annotations:
    input:
        # only keeping the file so we can move it out to the cwd
        FTP.remote(
            "https://ftp.ebi.ac.uk/pub/databases/RNAcentral/releases/21.0/genome_coordinates/gff3/homo_sapiens.GRCh38.gff3.gz",
            keep_local=True,
        ),
    output:
        "data/references/rna_central/genome_coordinates/homo_sapiens.GRCh38.gff3.gz",
    log:
        "results/logs/references/fetch_rna_central_annotations_gff.log",
    run:
        shell("mv {input:q} {output:q} &> {log:q}")


rule rna_central_typeName_attribute:
    input:
        "data/references/rna_central/genome_coordinates/homo_sapiens.GRCh38.gff3.gz",
    output:
        "data/references/rna_central/genome_coordinates/homo_sapiens.GRCh38.typeName.gff3",
    log:
        "results/logs/references/add_typename_rna_central_annotations_gff.log",
    shell:
        "gzip -dc {input:q} | awk 'match($9, /Name=([^;]*);.*type=([^;]*);/, a) "
        '\'{{print $0 ";typeName=" a[2] "|" a[1]}}\' '
        "> {output:q} "
        "2> {log:q}"


rule unzip:
    input:
        "data/references/{genome_file}.gz",
    output:
        "data/references/{genome_file}",
    log:
        "results/logs/unzip/{genome_file}.txt",
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
