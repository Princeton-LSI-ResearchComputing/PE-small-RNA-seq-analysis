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
        pjy142="data/references/pjy142.fa",
        pjy151="data/references/pjy151.fa",
        pjy209="data/references/pjy209.fa",
        grch38="data/references/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz",
    output:
        "data/references/grch38_pjy142_pjy151_pjy209.fa",
    log:
        "results/logs/join_references.log",
    conda:
        "../envs/fastx_toolkit.yaml"
    shell:
        """
        zcat {input.grch38:q} | \
        cat - {input.pjy142:q} {input.pjy151:q} {input.pjy209:q} | \
        fasta_formatter -o {output:q} \
        2> {log:q}
        """


rule samtools_fasta_index:
    input:
        "data/references/grch38_pjy142_pjy151_pjy209.fa",
    output:
        "data/references/grch38_pjy142_pjy151_pjy209.fa.fai",
    log:
        "results/logs/faidx/grch38_p12_targets.log",
    params:
        extra="",  # optional params string
    wrapper:
        "v1.7.0/bio/samtools/faidx"
