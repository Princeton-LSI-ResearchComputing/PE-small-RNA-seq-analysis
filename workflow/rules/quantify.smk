import os

# rule feature_counts_rna_central_wrapper:
#     input:
#         # list of sam or bam files
#         samples="results/alignments/Homo_sapiens.GRCh38.dna.primary_assembly/sorted/{sample}_{unit}.bam",
#         annotation="data/references/rna_central/genome_coordinates/homo_sapiens.GRCh38.gff3.gz",
#         # optional input
#         #chr_names="",           # implicitly sets the -A flag
#         #fasta="genome.fasta"    # implicitly sets the -G flag
#     output:
#         multiext(
#             "results/featureCounts/rna_central/{sample}_{unit}_noncoding_exon_name",
#             ".featureCounts",
#             ".featureCounts.summary",
#         ),
#     threads: 12
#     resources:
#         mem_mb=32000,
#     params:
#         strand=0,  # optional; strandness of the library (0: unstranded [default], 1: stranded, and 2: reversely stranded)
#         r_path="",  # implicitly sets the --Rpath flag
#         extra="-p -B -M -O --primary -t noncoding_exon --frac-overlap 1.0 -g typeName --extraAttributes type",  # count all primary alignments, even if reads are multimapped
#     log:
#         "results/logs/featureCounts/rna_central/{sample}_{unit}_noncoding_exon_name.log",
#     wrapper:
#         "v1.21.2/bio/subread/featurecounts"


rule feature_counts_rna_central_with_core:
    input:
        # sam or bam
        sample="results/alignments/Homo_sapiens.GRCh38.dna.primary_assembly/sorted/{sample}_{unit}.bam",
        annotation="data/references/rna_central/genome_coordinates/homo_sapiens.GRCh38.typeName.gff3",
    output:
        featureCounts="results/featureCounts/rna_central/{sample}_{unit}.featureCounts",
        summary="results/featureCounts/rna_central/{sample}_{unit}.featureCounts.summary",
        core="results/featureCounts/rna_central/{sample}_{unit}.bam.featureCounts",
    threads: 12
    resources:
        mem_mb=32000,
    params:
        strand=0,  # optional; strandness of the library (0: unstranded [default], 1: stranded, and 2: reversely stranded)
        extra="-p -B -M -O --primary -t noncoding_exon --fracOverlap 1.0 -R CORE -g typeName",  # count all primary alignments, even if reads are multimapped
    log:
        "results/logs/featureCounts/rna_central/{sample}_{unit}_noncoding_exon_name.log",
    conda:
        "../envs/subread.yaml"
    shell:
        (
            "featureCounts "
            "-T {threads} "
            "-s {params.strand} "
            "-a {input.annotation:q} "
            "{params.extra} "
            "-o {output[0]:q} "
            "{input.sample:q} "
            "&> {log:q} "
            "&& "
            'mv "$(basename{input.sample}).featureCounts" {output.core:q}'
        )
