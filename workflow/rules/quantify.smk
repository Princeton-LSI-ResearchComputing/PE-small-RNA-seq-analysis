import os


rule filtered_alignment:
    input:
        bam="results/alignments/{genome}/sorted/{sample}_{unit}.bam",
        idx="results/alignments/{genome}/sorted/{sample}_{unit}.bam.bai",
    output:
        temp(
            "results/alignments/{genome}/filtered/{sample}_{unit}_first_proper_pair.{format,sam|bam}"
        ),
    log:
        "results/logs/samtools_view_filter/{genome}/{sample}_{unit}_{format}_filter.log",
    params:
        extra="-f 66",
        region="",
    threads: 2
    wrapper:
        "v1.23.4/bio/samtools/view"


rule count_smrna_annotations:
    input:
        script="workflow/scripts/smRNAseqv3_annotation_read_processed_gff3_primary_assembly.py",
        sam="results/alignments/Homo_sapiens.GRCh38.dna.primary_assembly/filtered/{sample}_{unit}_first_proper_pair.sam",
        annotation="data/references/gencode.v43.primary_assembly.annotation.tsv",
    output:
        multiext(
            "results/smrna_count/{sample}_{unit}_first_proper_pair",
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


rule count_exogenous_rna_alignments:
    input:
        bam="results/alignments/exogenous_rna/filtered/{sample}_{unit}_first_proper_pair.bam",
        idx="results/alignments/exogenous_rna/filtered/{sample}_{unit}_first_proper_pair.bam.bai",
    output:
        "results/exogenous_rna_count/{sample}_{unit}_idxstats.txt",
    log:
        "results/logs/exogenous_rna_count/{sample}_{unit}_idxstats.log",
    params:
        extra="",
    wrapper:
        "v1.23.4/bio/samtools/idxstats"


rule feature_counts:
    input:
        # list of sam or bam files
        samples="results/alignments/Homo_sapiens.GRCh38.dna.primary_assembly/filtered/{sample}_{unit}_first_proper_pair.bam",
        annotation="data/references/gencode.v43.primary_assembly.annotation.gff3.gz",
        # optional input
        #chr_names="",           # implicitly sets the -A flag
        #fasta="genome.fasta"    # implicitly sets the -G flag
    output:
        multiext(
            "results/smrna_featurecounts/{sample}_{unit}_first_proper_pair",
            ".featureCounts",
            ".featureCounts.summary",
        ),
    threads: 4
    resources:
        mem_mb=30720,
    params:
        strand=0,  # optional; strandness of the library (0: unstranded [default], 1: stranded, and 2: reversely stranded)
        r_path="",  # implicitly sets the --Rpath flag
        extra="-p --countReadPairs --extraAttributes gene_type",
    log:
        "results/logs/smrna_featurecounts/{sample}_{unit}.log",
    wrapper:
        "v3.3.3/bio/subread/featurecounts"


rule feature_counts_multioverlap:
    input:
        # list of sam or bam files
        samples="results/alignments/Homo_sapiens.GRCh38.dna.primary_assembly/filtered/{sample}_{unit}_first_proper_pair.bam",
        annotation="data/references/gencode.v43.primary_assembly.annotation.gff3.gz",
        # optional input
        #chr_names="",           # implicitly sets the -A flag
        #fasta="genome.fasta"    # implicitly sets the -G flag
    output:
        multiext(
            "results/smrna_featurecounts_multioverlap/{sample}_{unit}_first_proper_pair",
            ".featureCounts",
            ".featureCounts.summary",
        ),
    threads: 4
    resources:
        mem_mb=30720,
    params:
        strand=0,  # optional; strandness of the library (0: unstranded [default], 1: stranded, and 2: reversely stranded)
        r_path="",  # implicitly sets the --Rpath flag
        extra="-O -p --countReadPairs --extraAttributes gene_type",
    log:
        "results/logs/smrna_featurecounts/{sample}_{unit}.log",
    wrapper:
        "v3.3.3/bio/subread/featurecounts"
