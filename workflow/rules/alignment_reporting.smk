rule exogenous_rna_bedpe:
    input:
        bam="results/alignments/exogenous_rna/sorted/{sample}_{unit}.bam",
        bai="results/alignments/exogenous_rna/sorted/{sample}_{unit}.bam.bai",
    output:
        bedpe="results/alignments/exogenous_rna/bedpe/{sample}_{unit}.bedpe",
    log:
        "results/logs/exogenous_rna_bedpe/{sample}_{unit}.log",
    params:
        samtools_options="-bf 0x2",  # only proper pairs
        regions="",
        bedtools_options="-bedpe",
    threads: 2
    group:
        "exogenous_targets_bam"
    conda:
        "../envs/bedtools_samtools.yaml"
    shell:
        """
        samtools view {params.samtools_options} --bam {input.bam:q} --threads {threads} {params.regions} | \
        samtools sort -n | \
        bedtools bamtobed {params.bedtools_options} -i stdin | \
        sort -k 1,1 -k6,6n > {output.bedpe:q} 2> {log:q}
        """


# rule feature_counts:
#     input:
#         # list of sam or bam files
#         samples="results/mapped_sorted/{sample}_{unit}.bam",
#         annotation="data/references/ncrna.gtf",
#         # optional input
#         #chr_names="",           # implicitly sets the -A flag
#         #fasta="genome.fasta"    # implicitly sets the -G flag
#     output:
#         multiext(
#             "results/{sample}",
#             ".featureCounts",
#             ".featureCounts.summary",
#             ".featureCounts.jcounts",
#         ),
#     threads: 2
#     params:
#         strand=0,  # optional; strandness of the library (0: unstranded [default], 1: stranded, and 2: reversely stranded)
#         r_path="",  # implicitly sets the --Rpath flag
#         extra="-O --fracOverlap 0.2 -J",
#     log:
#         "logs/{sample}.log",
#     wrapper:
#         "v1.7.0/bio/subread/featurecounts"
# rule exogenous_targets_bam:
#     input:
#         bam="results/mapped_sorted/{sample}_{unit}.bam",
#         bai="results/mapped_sorted/{sample}_{unit}.bam.bai",
#     output:
#         bam="results/exogenous_targets_bam/{sample}_{unit}_exogenous_targets.bam",
#         #idx="results/exogenous_targets_bam/{sample}_{unit}_exogenous_targets.bam.bai",
#     log:
#         "results/logs/exogenous_targets_bam/{sample}_{unit}.log",
#     params:
#         regions="PJY142 PJY151 PJY209",
#     threads: 2
#     group:
#         "exogenous_targets_bam"
#     conda:
#         "../envs/samtools.yaml"
#     shell:
#         """
#         samtools view --bam {input.bam:q} --output {output.bam} --thread {threads} {params.regions} &> {log:q}
#         """
#
#
# rule exogenous_targets_bam_index:
#     input:
#         "results/exogenous_targets_bam/{sample}_{unit}_exogenous_targets.bam",
#     output:
#         "results/exogenous_targets_bam/{sample}_{unit}_exogenous_targets.bam.bai",
#     log:
#         "results/logs/samtools_index/{sample}_{unit}.log",
#         "results/logs/exogenous_targets_bam_index/{sample}_{unit}.log",
#     params:
#         extra="",  # optional params string
#     threads: 2  # This value - 1 will be sent to -@
#     group:
#         "exogenous_targets_bam"
#     wrapper:
#         "v1.7.0/bio/samtools/index"
#
#
# rule exogenous_targets_genomecov:
#     input:
#         bam="results/exogenous_targets_bam/{sample}_{unit}_exogenous_targets.bam",
#         bai="results/exogenous_targets_bam/{sample}_{unit}_exogenous_targets.bam.bai",
#     output:
#         "results/exogenous_targets_genomecov/{sample}_{unit}_coverage.tsv",
#     log:
#         "results/logs/exogenous_targets_genomecov/{sample}_{unit}_coverage.log",
#     params:
#         "-pc -dz",  # optional parameters
#     resources:
#         time="06:00:00",
#     wrapper:
#         "v1.7.1/bio/bedtools/genomecov"
