# ============================================================
# Step 4: VIRGO2 microbiome profiling
#   install -> map (per sample) -> compile -> taxonomy -> gene annotation
#
# Needs: git, git-lfs, samtools, bowtie2, python3.
# VIRGO2 is cloned and indexed by the workflow into virgo2.home. If you already
# installed it by hand, create the flag so the rule does not delete and re-clone it:
#   touch reference/VIRGO2.installed
# ============================================================

rule install_virgo2:
    output:
        touch(VIRGO2_READY),
    threads: 8
    log:
        f"{LOG_DIR}/install_virgo2.log",
    shell:
        r"""
        set -euo pipefail
        exec > {log} 2>&1

        git lfs install
        rm -rf "{VIRGO2_HOME}"
        git clone https://github.com/ravel-lab/VIRGO2.git "{VIRGO2_HOME}"

        cd "{VIRGO2_HOME}"
        git lfs pull

        # 'VIRGO2.py install' unzips the annotation/fasta files and then crashes at
        # its own index step (AttributeError: no attribute 'threads'), so that crash
        # is tolerated here and the index is built explicitly below.
        python3 VIRGO2.py install || echo "VIRGO2.py install exited non-zero (known index-step bug); building the index manually"

        # Make sure the FASTA is unzipped, then build the bowtie2 index
        if [ ! -s FastaFiles/VIRGO2.fa ] && [ -s FastaFiles/VIRGO2.fa.gz ]; then
            gunzip -k FastaFiles/VIRGO2.fa.gz
        fi
        test -s FastaFiles/VIRGO2.fa

        bowtie2-build -f --seed 343 --threads {threads} FastaFiles/VIRGO2.fa Index/VIRGO2

        # Fail loudly if the index was not built
        compgen -G "Index/VIRGO2.1.bt2*" > /dev/null
        """


# VIRGO2 'map' accepts single-end reads only, so both mates are concatenated
# (concatenated gzip is valid gzip) and mapped together.
# NOTE: R1 and R2 reads share the same read names after concatenation; check
# that VIRGO2's map script tolerates this.
rule virgo2_map:
    input:
        r1=f"{HOST_DIR}/{{sample}}_1.nohost.fastq.gz",
        r2=f"{HOST_DIR}/{{sample}}_2.nohost.fastq.gz",
        ready=VIRGO2_READY,
    output:
        flag=touch(f"{VIRGO2_FLAG_DIR}/{{sample}}.done"),
    params:
        tmp=f"{TMP_DIR}/{{sample}}_virgo2",
        prefix=f"{VIRGO2_MAP_DIR}/{{sample}}",
    threads: 4
    resources:
        mem_mb=8000,
    log:
        f"{LOG_DIR}/virgo2_map/{{sample}}.log",
    shell:
        r"""
        set -euo pipefail
        mkdir -p "{VIRGO2_MAP_DIR}" "{VIRGO2_FLAG_DIR}"
        rm -rf "{params.tmp}"; mkdir -p "{params.tmp}"

        cat {input.r1} {input.r2} > "{params.tmp}/{wildcards.sample}.fastq.gz"

        python3 "{VIRGO2_HOME}/VIRGO2.py" map \
            -r "$PWD/{params.tmp}/{wildcards.sample}.fastq.gz" \
            -p {threads} \
            -o "$PWD/{params.prefix}" \
            > {log} 2>&1

        # Fail loudly if nothing was written for this sample
        compgen -G "{params.prefix}*" > /dev/null

        rm -rf "{params.tmp}"
        """


# compile reads EVERY result in VIRGO2_MAP_DIR, so keep other projects' results
# out of that folder.
rule virgo2_compile:
    input:
        expand(f"{VIRGO2_FLAG_DIR}/{{sample}}.done", sample=SAMPLES),
    output:
        nr=f"{VIRGO2_COMP_DIR}/virgo2_compiled.summary.NR.txt",
    log:
        f"{LOG_DIR}/virgo2_compile.log",
    shell:
        r"""
        set -euo pipefail
        mkdir -p "{VIRGO2_COMP_DIR}"

        python3 "{VIRGO2_HOME}/VIRGO2.py" compile \
            -i "$PWD/{VIRGO2_MAP_DIR}" \
            -o "$PWD/{VIRGO2_COMP_DIR}/VIRGO2_output" \
            > {log} 2>&1

        # Copy the compiled summary.NR table to a fixed name for downstream rules
        nr=$(find "{VIRGO2_COMP_DIR}" -maxdepth 1 -name '*summary.NR.txt' \
                  ! -name 'virgo2_compiled*' | head -n 1)
        test -n "$nr"
        cp "$nr" {output.nr}
        """


# Relative abundances and read counts per sample
# (defaults: bacteria only, species gene-detection threshold applied)
rule virgo2_taxonomy:
    input:
        nr=f"{VIRGO2_COMP_DIR}/virgo2_compiled.summary.NR.txt",
    output:
        directory(VIRGO2_TAX_DIR),
    log:
        f"{LOG_DIR}/virgo2_taxonomy.log",
    shell:
        r"""
        set -euo pipefail
        mkdir -p "{output}"

        # relative abundances
        python3 "{VIRGO2_HOME}/VIRGO2.py" taxonomy \
            -i "$PWD/{input.nr}" -o "$PWD/{output}/relabund" \
            > {log} 2>&1

        # read counts
        python3 "{VIRGO2_HOME}/VIRGO2.py" taxonomy \
            -i "$PWD/{input.nr}" -o "$PWD/{output}/readcounts" -r 1 \
            >> {log} 2>&1

        test -n "$(ls -A "{output}")"
        """


# ------------------------------------------------------------
# Gene annotation. Merges the compiled VIRGO2 gene table with the gene annotation
# tables shipped with VIRGO2 (AnnotationTables/): gene length, taxon, eggNOG, PFAM,
# EC, gene product, KEGG, CAZy, AMR, phage and biosynthetic gene cluster annotations.
# Only genes with at least annot_min_total_count reads summed over all samples are
# kept, so the output stays a manageable size.
#
# Output: one row per detected gene = read counts per sample + one column (or a few)
# per annotation table; several annotation values for one gene are joined with ";".
# The summary table shows, per annotation table, how many detected genes matched.
#
# The VIRGO2 README lists the tables but not their exact layout, so for each table
# the rule detects whether there is a header row and which column holds the gene ID
# (the one whose values match the compiled gene IDs), records both in the summary
# table, and stops with an error if no gene ID matches at all.
# ------------------------------------------------------------
rule virgo2_annotate:
    input:
        nr=f"{VIRGO2_COMP_DIR}/virgo2_compiled.summary.NR.txt",
        ready=VIRGO2_READY,
    output:
        table=f"{VIRGO2_ANNOT_DIR}/virgo2_genes_annotated.tsv",
        summary=f"{VIRGO2_ANNOT_DIR}/virgo2_annotation_summary.tsv",
    resources:
        mem_mb=20000,
    run:
        import csv, glob, re
        from itertools import chain, islice

        csv.field_size_limit(10**9)

        # 1. compiled gene x sample table
        all_ids, kept = set(), {}
        with open(input.nr) as fh:
            rd = csv.reader(fh, delimiter="\t")
            first = next(rd)
            try:
                float(first[1])
                samples, rows = [f"sample_{i + 1}" for i in range(len(first) - 1)], chain([first], rd)
            except ValueError:
                samples, rows = first[1:], rd
            for row in rows:
                if not row:
                    continue
                gene = row[0]
                all_ids.add(gene)
                total = 0.0
                for x in row[1:]:
                    try:
                        total += float(x)
                    except ValueError:
                        pass
                if total >= ANNOT_MIN_TOTAL_COUNT:
                    kept[gene] = row[1:]
        if not kept:
            raise ValueError("No gene has reads in the compiled VIRGO2 table")

        # 2. annotation tables shipped with VIRGO2
        found = []
        for p in glob.glob(f"{VIRGO2_HOME}/AnnotationTables/*"):
            m = re.match(r"^(\d+)\.VIRGO2\.([A-Za-z0-9_]+)\.txt(\.gz)?$", os.path.basename(p))
            if m and int(m.group(1)) in ANNOT_TABLE_IDS:
                found.append((int(m.group(1)), m.group(2), p))
        if not found:
            raise ValueError(f"No annotation tables found in {VIRGO2_HOME}/AnnotationTables/ "
                             "(were they unzipped by 'VIRGO2.py install'?)")
        found.sort()

        out_cols, per_table, summary = [], [], []
        for num, kind, path in found:
            opener = gzip.open if path.endswith(".gz") else open
            with opener(path, "rt") as fh:
                rd = csv.reader(fh, delimiter="\t")
                first = next(rd)
                header = not any(f in all_ids for f in first)
                buf = list(islice(rd, 2000))
                sample = ([] if header else [first]) + buf
                ncol = max(len(r) for r in [first] + buf)
                # The gene ID column is the one whose values best match the compiled
                # gene IDs (some tables have an extra row-number column before it).
                hits = [sum(1 for r in sample if i < len(r) and r[i] in all_ids)
                        for i in range(ncol)]
                key = max(range(ncol), key=lambda i: hits[i]) if max(hits) > 0 else 0
                # value columns: everything except the gene ID and any all-integer
                # row-number column to its left
                value_cols = [i for i in range(ncol) if i != key and not (
                    i < key and all(r[i].isdigit() for r in sample if i < len(r)))]
                if header:
                    names = [first[i] if i < len(first) and first[i] else f"col{i + 1}"
                             for i in value_cols]
                    cols = [f"{kind}:{n}" for n in names]
                    key_label = first[key] if key < len(first) and first[key] else f"column {key + 1}"
                else:
                    cols = [kind] if len(value_cols) == 1 else [f"{kind}_{j + 1}" for j in range(len(value_cols))]
                    key_label = f"column {key + 1}"
                ann = {}
                for row in chain([] if header else [first], buf, rd):
                    if key >= len(row) or row[key] not in kept:
                        continue
                    vals = ann.setdefault(row[key], [[] for _ in value_cols])
                    for j, i in enumerate(value_cols):
                        v = row[i] if i < len(row) else ""
                        if v and v not in vals[j]:
                            vals[j].append(v)
            out_cols += cols
            per_table.append((cols, ann))
            summary.append((os.path.basename(path), header, key_label, ",".join(cols), len(ann),
                            100.0 * len(ann) / len(kept)))

        if all(n == 0 for _, _, _, _, n, _ in summary):
            raise ValueError("No gene ID in the compiled table matched any annotation table; "
                             "check the file layouts (see the summary table)")

        # 3. write
        with open(output.table, "w") as out:
            out.write("gene\t" + "\t".join(samples) + "\t" + "\t".join(out_cols) + "\n")
            for gene, counts in kept.items():
                fields = [gene] + counts
                for cols, ann in per_table:
                    vals = ann.get(gene)
                    for i in range(len(cols)):
                        fields.append(";".join(vals[i]) if vals and vals[i] else "NA")
                out.write("\t".join(fields) + "\n")

        with open(output.summary, "w") as out:
            out.write("annotation_table\theader_detected\tgene_id_column\tcolumns\tdetected_genes_annotated\tpercent_of_detected_genes\n")
            out.write(f"compiled_table\t-\t-\t-\t{len(kept)}\t(of {len(all_ids)} genes in the table)\n")
            for name, header, key_label, cols, n, pct in summary:
                out.write(f"{name}\t{header}\t{key_label}\t{cols}\t{n}\t{pct:.1f}\n")
