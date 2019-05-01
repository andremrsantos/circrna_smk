configfile: "config.yaml"

pe_samples, = glob_wildcards(config["seq_pe"] + "/{sample}_R1.fastq.gz")
se_samples, = glob_wildcards(config["seq_se"] + "/{sample,[\w\d-]+}.fastq.gz")
samples = se_samples + pe_samples

rule all:
    input:
        # Trimming output
        expand("seq/{sample}.trim.fastq.gz", sample = se_samples),
        expand("seq/{sample}_R1.trim.fastq.gz", sample = pe_samples),
        expand("seq/{sample}_R2.trim.fastq.gz", sample = pe_samples),
        # STAR mapping
        expand("align/{sample}_star/Aligned.sortedByCoord.out.bam",
               sample = samples),
        expand("align/{sample}_star/Chimeric.out.junction",
               sample = samples),
        expand("align/{sample}_star/{mate}.Chimeric.out.junction",
               sample = pe_samples, mate = ["R1", "R2"]),
        # BWA mapping
        expand("align/{sample}_bwa.sam", sample = samples),
        # CircRNA quantification
        expand("circrna/ciri/{sample}_circRNA.tsv", sample = samples),
        expand("circrna/circexplorer/{sample}_circRNA.tsv", sample = samples),
        "circrna/dcc_se/CircRNACount",
        "circrna/dcc_pe/CircRNACount"


singularity: "docker://continuumio/miniconda3:4.5.12"

include: "rules/trim.smk"
include: "rules/align.smk"
include: "rules/quant.smk"

## Getting all rules
rule all:
    input:
        expand("circrna/ce/{samples}_back_spliced_junction.bed",
               samples = samples),
        expand("circrna/ciri/{samples}_circRNA.txt",
               samples = samples),
        "circrna/dcc/se/CircRNACount",
        "circrna/dcc/pe/CircRNACount",
        "circrna/ce/ce_circrna.tsv",
        "circrna/ciri/ciri_circrna.tsv",
        "circrna/dcc/dcc_circrna.tsv",
        "circrna/circrna_intersect.tsv"
