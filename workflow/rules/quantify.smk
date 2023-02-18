import os


rule bam_to_sam:
    input:
        "results/alignments/Homo_sapiens.GRCh38.dna.primary_assembly/sorted/{sample}_{unit}.bam",
    output:
        temp(
            "results/alignments/Homo_sapiens.GRCh38.dna.primary_assembly/sorted/{sample}_{unit}.sam"
        ),
    log:
        "results/logs/samtools_view/{sample}_{unit}_sam_to_bam.log",
    params:
        extra="-f 66",
        region="",
    threads: 2
    wrapper:
        "v1.23.4/bio/samtools/view"


rule count_smrna_annotations:
    input:
        script="workflow/scripts/smRNAseqv3_annotation_read_processed_gff3_primary_assembly.py",
        sam="results/alignments/Homo_sapiens.GRCh38.dna.primary_assembly/sorted/{sample}_{unit}.sam",
        annotation="data/references/gencode.v43.primary_assembly.annotation.tsv",
    output:
        multiext(
            "results/smrna_count/{sample}_{unit}",
            "_biotype_count.txt",
            "_gene_count.txt",
            "_log.txt",
            "_transcript_count.txt",
            "_unannotated_chr_count.txt",
        ),
    params:
        outdir=lambda wildcards, output: Path(output[0]).parent,
    log:
        "results/logs/smrna_count/{sample}_{unit}.log",
    conda:
        "../envs/smrnaseq_annotation_count.yaml"
    shell:
        (
            "{input.script:q} "
            "--sam_file {input.sam:q} "
            "--annotation_file {input.annotation:q} "
            "--outdir {params.outdir:q} "
            "&> {log:q}"
        )
