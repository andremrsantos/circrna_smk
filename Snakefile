configfile: "config.yaml"

pe_samples, = glob_wildcards(config["seq_pe"] + "/{sample}_R1.fastq.gz")
se_samples, = glob_wildcards(config["seq_se"] + "/{sample,[\w\d-]+}.fastq.gz")
samples = se_samples + pe_samples
seq_types = []
if len(pe_samples) > 0:
    seq_types.append("pe")

if len(se_samples) > 0:
    seq_types.append("se")

print(pe_samples)
print(se_samples)

rule all:
    input:
        # Trimming output
        expand("seq/{sample}.trim.fastq.gz", sample = se_samples),
        expand("seq/{sample}_{mate}.trim.fastq.gz",
               sample = pe_samples, mate = ["R1", "R2"]),
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
        expand("circrna/dcc_{seq_type}/CircRNACount", seq_type = seq_types)

singularity: "docker://continuumio/miniconda3:4.5.12"

include: "rules/trim.smk"
include: "rules/align.smk"
include: "rules/quant.smk"
