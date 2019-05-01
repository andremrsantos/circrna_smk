junc_pattern = "align/{sample}_star/{mate}Chimeric.out.junction"

def get_sheets(wc):
    if wc.seq_type == "se":
        return ["circrna/dcc_se/samplesheet"]
    return ["circrna/dcc_pe/samplesheet",
            "circrna/dcc_pe/mate1",
            "circrna/dcc_pe/mate2"]

def get_aligns(wc):
    return expand(
        "align/{sample}_star/Aligned.sortedByCoord.bam",
        sample = se_samples if wc.seq_type == "se" else pe_samples)

rule ciri:
    input: "align/{sample}_bwa.sam"
    output: "circrna/ciri/{sample}_circRNA.tsv"
    log: "log/{sample}_ciri.log"
    threads: 8
    params:
        ref = config["ref"]["bwa"],
        ann = config["ref"]["annotation"]
    singularity: "docker://andremrsantos/ciri2"
    shell: """
    CIRI2 \
    -I {input} -O {output} \
    -F {params.ref} -A {params.ann} \
    -T {threads} -G {log} -M "None"
    """

rule circexplorer_parse:
    input: "align/{sample}_star/Chimeric.out.junction"
    output: "circrna/circexplorer/{sample}_bsj.bed"
    log: "log/{sample}_ce_parse.log"
    threads: 4
    conda: "envs/circexplorer.yaml"
    shell: """
    CIRCexplorer2 parse -t STAR -b {output} {input} > {log}
    """

rule circexplorer_annotate:
    input: "circrna/circexplorer/{sample}_bsj.bed"
    output: "circrna/circexplorer/{sample}_circRNA.tsv"
    log: "log/{sample}_ce_annotate.log"
    threads: 4
    conda: "envs/circexplorer.yaml"
    params:
        ref = config["ref"]["bwa"],
        pred = config["ref"]["pred"]
    shell: """
    CIRCexplorer2 annotate \
    -g {params.ref} -r {params.pred} \
    -b {input} -o {output} > {log}
    """

rule dcc_se_samplesheet:
    input: expand(junc_pattern, sample = se_samples, mate = "")
    output: "circrna/dcc_se/samplesheet"
    run:
        with open(output[0], "w") as out:
            out.write("\n".join(input))

rule dcc_pe_samplesheet:
    input:
        expand(junc_pattern, sample = pe_samples, mate = ""),
        expand(junc_pattern, sample = pe_samples, mate = "R1."),
        expand(junc_pattern, sample = pe_samples, mate = "R2.")
    output:
        "circrna/dcc_pe/samplesheet",
        "circrna/dcc_pe/mate1",
        "circrna/dcc_pe/mate2"
    run:
        for i in range(3):
            with open(output[i], "w") as out:
                out.write("\n".join(input[i]))

def mates(wc, input):
    if wc.seq_type == "pe":
        return "-mt1 @{}, -mt2 @{} -Pi".format(input[1], input[2])
    return ""

rule dcc:
    input:
        sheet = get_sheets,
        align = get_aligns
    output: "circrna/dcc_{seq_type}/CircRNACount"
    threads: 16
    params:
        ref = config["ref"]["bwa"],
        ann = config["ref"]["annotation"],
        wkdir = "circrna/dcc_{seq_type}",
        mates = mates,
    singularity: "docker://andremrsantos/dcc"
    shell: """
    DCC {input.sheet} \
    -D -F -fg -G -Nr 5 1 {params.mates} \
    -T {threads} \
    -A {params.ref} -an {params.ann} \
    -O {params.wkdir}/ \
    -t {params.wkdir}/_tmp
    """
