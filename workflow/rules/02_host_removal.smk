# ============================================================
# Step 2: host (human, GRCh38) read removal
# The alignment is streamed to samtools and only pairs with BOTH mates
# unmapped are kept (SAM flag 12 = read unmapped + mate unmapped). This
# differs from bowtie2 --un-conc, which also keeps pairs where only one
# mate aligns. --mm memory-maps the index so parallel jobs share it.
# ============================================================

rule bowtie2_index:
    input:
        HOST_FASTA,
    output:
        multiext(HOST_INDEX, ".1.bt2", ".2.bt2", ".3.bt2", ".4.bt2",
                 ".rev.1.bt2", ".rev.2.bt2"),
    threads: 8
    log:
        f"{LOG_DIR}/bowtie2_index.log",
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname {HOST_INDEX})"
        bowtie2-build --threads {threads} {input} {HOST_INDEX} > {log} 2>&1
        """


rule remove_host:
    input:
        r1=f"{TRIM_DIR}/{{sample}}_1.trimmed.fastq.gz",
        r2=f"{TRIM_DIR}/{{sample}}_2.trimmed.fastq.gz",
        idx=multiext(HOST_INDEX, ".1.bt2", ".2.bt2", ".3.bt2", ".4.bt2",
                     ".rev.1.bt2", ".rev.2.bt2"),
    output:
        r1=f"{HOST_DIR}/{{sample}}_1.nohost.fastq.gz",
        r2=f"{HOST_DIR}/{{sample}}_2.nohost.fastq.gz",
        stats=f"{BT2_DIR}/{{sample}}.bowtie2.log",
    threads: 7
    resources:
        mem_mb=8000,
    shell:
        r"""
        set -euo pipefail
        mkdir -p "{HOST_DIR}" "{BT2_DIR}"

        bowtie2 \
            -x {HOST_INDEX} \
            -1 {input.r1} -2 {input.r2} \
            --very-sensitive-local --mm \
            --threads {threads} \
            2> {output.stats} \
          | samtools view -@ 1 -u -f 12 - \
          | samtools fastq -@ 2 -n -0 /dev/null -s /dev/null \
                -1 {output.r1} -2 {output.r2} -

        # Fail loudly if bowtie2 did not finish or nothing was written
        grep -q "overall alignment rate" {output.stats}
        test -s {output.r1}
        test -s {output.r2}
        """
