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


rule fetch_gencode_annotations_gff:
    input:
        # only keeping the file so we can move it out to the cwd
        FTP.remote(
            "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_43/gencode.v43.primary_assembly.annotation.gff3.gz",
            keep_local=True,
        ),
    output:
        "data/references/gencode.v43.primary_assembly.annotation.gff3.gz",
    log:
        "results/logs/reference/fetch_grch38_gff.log",
    run:
        shell("mv {input:q} {output:q} &> {log:q}")


rule gff_to_tsv:
    input:
        gff="data/references/gencode.v43.primary_assembly.annotation.gff3",
    output:
        tsv="data/references/gencode.v43.primary_assembly.annotation.tsv",
    log:
        "results/logs/reference/gffread-gff-to-tsv.log",
    conda:
        "../envs/gffread.yaml"
    shell:
        'gffread --table "@chr, @start, @end, gene_name, transcript_name, transcript_type" {input.gff:q} -o {output.tsv:q} 2> {log}'


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
