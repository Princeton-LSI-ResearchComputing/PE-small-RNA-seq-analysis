rule feature_counts:
    input:
        # list of sam or bam files
        samples=expand(
            "results/alignments/Homo_sapiens.GRCh38.dna.primary_assembly/sorted/{unit.sample_name}_{unit.unit_name}.bam",
            unit=units.itertuples(),
        ),
        annotation="data/references/Homo_sapiens.GRCh38.108.gtf.gz",
        # optional input
        #chr_names="",           # implicitly sets the -A flag
        #fasta="genome.fasta"    # implicitly sets the -G flag
    output:
        multiext(
            "results/alignments/Homo_sapiens.GRCh38.dna.primary_assembly/featureCounts/all_counts",
            ".featureCounts",
            ".featureCounts.summary",
        ),
    threads: 12
    params:
        strand=0,  # optional; strandness of the library (0: unstranded [default], 1: stranded, and 2: reversely stranded)
        r_path="",  # implicitly sets the --Rpath flag
        # extra="-O --fracOverlap 0.2 -J",
    log:
        "results/logs/featureCounts/Homo_sapiens.GRCh38.dna.primary_assembly/all_counts.log",
    wrapper:
        "v1.21.2/bio/subread/featurecounts"
