# ============================================================
# Step 1: SRA -> paired FASTQ.gz -> FastQC -> fastp trimming
# Raw and trimmed FASTQ are temp() and deleted once downstream rules finish.
# The .sra files are never touched.
# ============================================================

rule sra_to_fastq:
    input:
        sra=f"{SRA_DIR}/{{sample}}/{{sample}}.sra",
    output:
        r1=temp(f"{FASTQ_DIR}/{{sample}}_1.fastq.gz"),
        r2=temp(f"{FASTQ_DIR}/{{sample}}_2.fastq.gz"),
    params:
        tmp=f"{TMP_DIR}/{{sample}}",
    threads: 4
    resources:
        mem_mb=3000,
    log:
        f"{LOG_DIR}/sra_to_fastq/{{sample}}.log",
    shell:
        r"""
        set -euo pipefail
        exec > {log} 2>&1

        rm -rf "{params.tmp}"
        mkdir -p "{params.tmp}/tmp" "$(dirname {output.r1})"

        fasterq-dump "{input.sra}" \
            --split-files \
            --threads {threads} \
            --temp "{params.tmp}/tmp" \
            --outdir "{params.tmp}"

        pigz -p {threads} "{params.tmp}/{wildcards.sample}_1.fastq" \
                          "{params.tmp}/{wildcards.sample}_2.fastq"

        mv "{params.tmp}/{wildcards.sample}_1.fastq.gz" "{output.r1}"
        mv "{params.tmp}/{wildcards.sample}_2.fastq.gz" "{output.r2}"

        rm -rf "{params.tmp}"
        """


rule fastqc_raw:
    input:
        f"{FASTQ_DIR}/{{sample}}_{{read}}.fastq.gz",
    output:
        html=f"{QC_DIR}/{{sample}}_{{read}}_fastqc.html",
        zip=f"{QC_DIR}/{{sample}}_{{read}}_fastqc.zip",
    threads: 1
    resources:
        mem_mb=1000,
    log:
        f"{LOG_DIR}/fastqc/{{sample}}_{{read}}.log",
    shell:
        r"""
        set -euo pipefail
        mkdir -p "{QC_DIR}"
        fastqc --threads {threads} --outdir "{QC_DIR}" {input} > {log} 2>&1
        """


rule fastp:
    input:
        r1=f"{FASTQ_DIR}/{{sample}}_1.fastq.gz",
        r2=f"{FASTQ_DIR}/{{sample}}_2.fastq.gz",
    output:
        r1=temp(f"{TRIM_DIR}/{{sample}}_1.trimmed.fastq.gz"),
        r2=temp(f"{TRIM_DIR}/{{sample}}_2.trimmed.fastq.gz"),
        html=f"{FASTP_DIR}/{{sample}}.fastp.html",
        json=f"{FASTP_DIR}/{{sample}}.fastp.json",
    threads: 4
    resources:
        mem_mb=2000,
    log:
        f"{LOG_DIR}/fastp/{{sample}}.log",
    shell:
        r"""
        set -euo pipefail
        mkdir -p "{TRIM_DIR}" "{FASTP_DIR}"

        fastp \
            --in1 {input.r1} --in2 {input.r2} \
            --out1 {output.r1} --out2 {output.r2} \
            --detect_adapter_for_pe \
            {FASTP_QC} \
            --html {output.html} --json {output.json} \
            --thread {threads} \
            > {log} 2>&1
        """


# Only requested when trimmed_fastqc: true in the config
rule fastqc_trimmed:
    input:
        f"{TRIM_DIR}/{{sample}}_{{read}}.trimmed.fastq.gz",
    output:
        html=f"{QC_TRIM_DIR}/{{sample}}_{{read}}.trimmed_fastqc.html",
        zip=f"{QC_TRIM_DIR}/{{sample}}_{{read}}.trimmed_fastqc.zip",
    threads: 1
    resources:
        mem_mb=1000,
    log:
        f"{LOG_DIR}/fastqc_trimmed/{{sample}}_{{read}}.log",
    shell:
        r"""
        set -euo pipefail
        mkdir -p "{QC_TRIM_DIR}"
        fastqc --threads {threads} --outdir "{QC_TRIM_DIR}" {input} > {log} 2>&1
        """
