trim = ("--adapter_sequence {adapter_r1} " +
        "--adapter_sequence_r2 {adapter_r2} " +
        "-l {min_length} -5 {left} -3 {right} " +
        "-W {window} -M {window_mean} -g").format(**config["trim"])

rule trim_fastp_se:
    input: config["seq_se"] + "/{sample}.fastq.gz"
    output: "seq/{sample}.trim.fastq.gz"
    threads: 8
    log:
        html = "log/{sample}_fastp.html",
	json = "log/{sample}_fastp.json"
    conda: "envs/trim.yaml"
    shell: """
    fastp --thread {threads} {trim} \
    -i {input} -o {output} \
    -j {log.json}  -h {log.html}
    """

rule trim_fastp_pe:
    input:
        R1 = config["seq_pe"] + "/{sample}_R1.fastq.gz",
        R2 = config["seq_pe"] + "/{sample}_R2.fastq.gz"
    output:
        R1 = "seq/{sample}_R1.trim.fastq.gz",
        R2 = "seq/{sample}_R2.trim.fastq.gz"
    threads: 8
    log:
        html = "log/{sample}_fastp.html",
	json = "log/{sample}_fastp.json"
    conda: "envs/trim.yaml"
    shell: """
    fastp --thread {threads} {trim} \
    -i {input.R1}  -I {input.R2} \
    -o {output.R1} -O {output.R2} \
    -j {log.json} -h {log.html}
    """
