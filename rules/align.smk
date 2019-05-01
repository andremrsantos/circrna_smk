def get_trimmed(wc):
    if wc.sample in pe_samples:
        return expand("seq/{sample}_{mate}.trim.fastq.gz",
                      mate = ["R1", "R2"], **wc)
    return "seq/{sample}.trim.fastq.gz".format(**wc)

rule align_bwa:
    input: get_trimmed
    output: "align/{sample}_bwa.sam"
    log: "log/{sample}_bwa.log"
    threads: 16
    params:
        index = config["ref"]["bwa"]
    conda: "envs/map.yaml"
    shell: """
    bwa mem -t {threads} -T19 {params.index} {input} \
    > {output} 2> {log}
    """

rule align_star:
    input: get_trimmed
    output:
        "align/{sample}_star/Aligned.sortedByCoord.out.bam"
        "align/{sample}_star/Chimeric.out.junction"
        "align/{sample}_star/ReadsPerGene.out.tab"
    log: "log/{sample}_star.log"
    threads: 16
    params:
        index = config["ref"]["star"],
        annotation = config["ref"]["annotation"],
        prefix = "align/{sample}_star/"
    conda: "envs/map.yaml"
    shell: """
    STAR --runThreadN {threads} \
    --genomeDir {params.index} \
    --outSAMtype BAM SortedByCoordinate \
    --readFilesIn {input} \
    --readFilesCommand zcat \
    --sjdbGTFfile {params.annotation} \
    --outFileNamePrefix {params.prefix} \
    --quantMode TranscriptomeSAM \
    --outSAMunmapped Within \
    --outFilterType BySJout \
    --outSAMattributes NH HI AS NM MD \
    --outFilterMultimapNmax 20 \
    --outFilterMismatchNmax 999 \
    --outFilterScoreMin 1 \
    --outFilterMatchNmin 1 \
    --outFilterMismatchNoverReadLmax 0.04 \
    --outSJfilterOverhangMin 15 15 15 15 \
    --alignSJoverhangMin 15 \
    --alignSJDBoverhangMin 15 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --sjdbScore 1 \
    --chimSegmentMin 15 \
    --chimScoreMin 15 \
    --chimScoreSeparation 10 \
    --chimJunctionOverhangMin 15 \
    --outStd Log > {log}

    samtools index {output[1]} >> {log}
    """

rule align_star_mate:
    input:
        "seq/{sample}_{mate}.trim.fastq.gz"
    output:
        "align/{sample}_star/{mate}.Chimeric.out.junction"
    log: "log/{sample}_{mate}_star.log"
    threads: 16
    params:
        index = config["ref"]["star"],
        annotation = config["ref"]["annotation"],
        prefix = "align/{sample}_star/{mate}."
    conda: "envs/map.yaml"
    shell: """
    STAR --runThreadN {threads} \
    --genomeDir {params.index} \
    --outSAMtype None \
    --readFilesIn {input} \
    --readFilesCommand zcat \
    --sjdbGTFfile {params.annotation}
    --outFileNamePrefix {params.prefix} \
    --outFilterMultimapNmax 20 \
    --outFilterMismatchNmax 2 \
    --outFilterScoreMin 1 \
    --outFilterMatchNmin 1 \
    --outSJfilterOverhangMin 15 15 15 15 \
    --alignSJoverhangMin 15 \
    --alignSJDBoverhangMin 15 \
    --chimSegmentMin 15 \
    --chimScoreMin 15 \
    --chimScoreSeparation 10 \
    --chimJunctionOverhangMin 15 \
    --outStd Log > {log}
    """
