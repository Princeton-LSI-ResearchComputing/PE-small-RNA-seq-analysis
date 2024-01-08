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


rule exogenous_rna_alignment_stats:
    input:
        bam=expand(
            "results/alignments/exogenous_rna/sorted/{unit.sample_name}_{unit.unit_name}.bam",
            unit=units.itertuples(),
        ),
        bai=expand(
            "results/alignments/exogenous_rna/sorted/{unit.sample_name}_{unit.unit_name}.bam.bai",
            unit=units.itertuples(),
        ),
        human_small_rna=expand(
            "results/smrna_count/{unit.sample_name}_{unit.unit_name}_first_proper_pair_gene_count.txt",
            unit=units.itertuples(),
        ),
    output:
        "results/alignments/exogenous_rna_alignment_stats/all_exogenous_rna_alignments_stats.tsv",
    log:
        "results/logs/exogenous_rna_alignment_stats/all_exogenous_rna_alignments_stats.ipynb",
    conda:
        "../envs/alignment_stats.yaml"
    notebook:
        "notebooks/alignment_stats.py.ipynb"
