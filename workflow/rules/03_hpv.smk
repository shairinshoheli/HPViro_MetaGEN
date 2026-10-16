# ============================================================
# Step 3: HPV classification, modelled on Ma et al. 2014 (doi 10.1128/jvi.00093-14)
#
#   1. blast_hpv         BLASTN of each read against the HPV database, E <= max_evalue
#                        (a hit at that E-value makes the read an HPV "candidate")
#   2. hpv_stage1        best hit per read; "strict" = identity >= min_pident and
#                        query coverage >= min_qcov (paper: >=95% and full coverage).
#                        Reads whose best score is shared by 2+ types AND pass the
#                        strict rule are ambiguous and excluded (listed separately).
#   3. hpv_tblastx       TBLASTX of the candidate reads against the same HPV database
#   4. hpv_stage2        strict reads whose TBLASTX best hit agrees with the BLASTN
#                        type are CONFIRMED. All other candidates go to the GenBank check.
#   5. hpv_genbank       BLASTN of those reads against GenBank (NCBI remote BLAST or a
#                        local database). Best hit HPV + definite type + strict
#                        identity/coverage -> that type; best hit HPV otherwise ->
#                        "unclassified"; best hit non-HPV (or no hit) -> discarded.
#   6. hpv_count_matrix  counts per sample (unit set by count_unit) + a "read fate"
#                        table showing how many reads took each route.
#   7. normalisation     per million raw read pairs, per million non-human read pairs,
#                        and per-sample genotype proportions.
#
# Choices made by this workflow, not stated in the paper: the TBLASTX E-value, the
# GenBank E-value, the batch size and wait time for remote BLAST, the rule that tied
# reads are excluded, and the title-matching rules that decide whether a GenBank hit
# is an HPV and which type it is.
# ============================================================

import hashlib


def clean(x):
    return x[4:] if x.startswith("lcl|") else x


def seqhash(seq):
    return hashlib.md5(seq.upper().encode()).hexdigest()


def read_fasta(path):
    """Read a FASTA file with one-line sequences into {name: sequence}."""
    out, name = {}, None
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                name = line[1:].split()[0]
                out[name] = ""
            elif name is not None:
                out[name] += line
    return out


def extract_reads(fastq_gz, wanted, mate, out_fh):
    """Write the reads named in `wanted` (names like SRRxxx.1/1) to out_fh as FASTA."""
    if not wanted:
        return
    proc = subprocess.Popen(["pigz", "-dc", fastq_gz], stdout=subprocess.PIPE, text=True)
    n, name, keep = 0, None, False
    for line in proc.stdout:
        n += 1
        k = n % 4
        if k == 1:
            name = line.split()[0][1:] + "/" + str(mate)
            keep = name in wanted
        elif k == 2 and keep:
            out_fh.write(f">{name}\n{line}")
    if proc.wait() != 0:
        raise RuntimeError(f"pigz failed on {fastq_gz}")


# Title matching for GenBank hits (heuristic, see note above)
_HPV_TITLE = re.compile(
    r"human\s+papillomavirus|\bHPV\b|\b(?:alpha|beta|gamma|mu|nu)papillomavirus", re.I)
_HPV_TYPE = re.compile(r"(?:human\s+papillomavirus|\bHPV)[\s-]*(?:type[\s-]*)?(\d+)", re.I)


def classify_genbank_hit(title, pident, qcov):
    """Return (class, assigned_type) for the best GenBank hit of a read."""
    if not _HPV_TITLE.search(title):
        return "non_hpv", ""
    m = _HPV_TYPE.search(title)
    if m and pident >= MIN_PIDENT and qcov >= MIN_QCOV - 1e-9:
        return "typed", "HPV" + m.group(1)
    return "unclassified", ""


CAND_COLS = "qid\tmate\tref\tpident\tqcov\tbits\tstrict\tambiguous\n"
GB_COLS = ["seqhash", "class", "assigned_type", "pident", "qcov", "bits", "title"]


# ------------------------------------------------------------
# Step 1: BLASTN of host-depleted reads against the HPV database
# (FASTQ is converted to FASTA on the fly; each mate is searched separately)
# ------------------------------------------------------------
rule blast_hpv:
    input:
        reads=f"{HOST_DIR}/{{sample}}_{{read}}.nohost.fastq.gz",
        db=multiext(HPV_DB, ".nhr", ".nin", ".nsq"),
    output:
        f"{BLAST_DIR}/{{sample}}_{{read}}.hpv.blast.tsv.gz",
    threads: 4
    resources:
        mem_mb=2000,
    log:
        f"{LOG_DIR}/blast_hpv/{{sample}}_{{read}}.log",
    shell:
        r"""
        set -euo pipefail
        mkdir -p "{BLAST_DIR}"

        zcat {input.reads} \
            | awk 'NR%4==1{{print ">"substr($1,2)"/{wildcards.read}"}} NR%4==2{{print}}' \
            | blastn \
                -task megablast \
                -query - \
                -db {HPV_DB} \
                -evalue {MAX_EVALUE} \
                -max_target_seqs 1000 \
                -max_hsps 1 \
                -num_threads {threads} \
                -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen" \
                2> {log} \
            | pigz -p 2 > {output}
        """


# ------------------------------------------------------------
# Step 2: best hit per read; write candidate reads as FASTA
# ------------------------------------------------------------
rule hpv_stage1:
    input:
        b1=f"{BLAST_DIR}/{{sample}}_1.hpv.blast.tsv.gz",
        b2=f"{BLAST_DIR}/{{sample}}_2.hpv.blast.tsv.gz",
        r1=f"{HOST_DIR}/{{sample}}_1.nohost.fastq.gz",
        r2=f"{HOST_DIR}/{{sample}}_2.nohost.fastq.gz",
    output:
        stage1=f"{CAND_DIR}/{{sample}}.stage1.tsv",
        fasta=f"{CAND_DIR}/{{sample}}.candidates.fa",
    resources:
        mem_mb=2000,
    run:
        from itertools import groupby

        wanted = {1: set(), 2: set()}
        with open(output.stage1, "w") as out:
            out.write(CAND_COLS)
            for mate, path in ((1, input.b1), (2, input.b2)):
                with gzip.open(path, "rt") as fh:
                    rows = (l.rstrip("\n").split("\t") for l in fh if l.strip())
                    for qid, grp in groupby(rows, key=lambda x: x[0]):
                        hits = {}   # ref -> (pident, qcov, bits)
                        for f in grp:
                            pident, evalue, bits = float(f[2]), float(f[10]), float(f[11])
                            if evalue > MAX_EVALUE or bits <= MIN_BITSCORE:
                                continue
                            qcov = (abs(int(f[7]) - int(f[6])) + 1) / int(f[12])
                            ref = clean(f[1])
                            if ref not in hits or bits > hits[ref][2]:
                                hits[ref] = (pident, qcov, bits)
                        if not hits:
                            continue
                        best = max(h[2] for h in hits.values())
                        tied = sorted(r for r, h in hits.items() if h[2] == best)
                        pident, qcov, bits = hits[tied[0]]
                        strict = any(hits[r][0] >= MIN_PIDENT and hits[r][1] >= MIN_QCOV - 1e-9
                                     for r in tied)
                        ambiguous = strict and len(tied) > 1
                        out.write(f"{qid}\t{mate}\t{','.join(tied)}\t{pident:.2f}\t"
                                  f"{qcov:.3f}\t{bits:.1f}\t{int(strict)}\t{int(ambiguous)}\n")
                        if not ambiguous:
                            wanted[mate].add(qid)

        with open(output.fasta, "w") as fa:
            extract_reads(input.r1, wanted[1], 1, fa)
            extract_reads(input.r2, wanted[2], 2, fa)


# ------------------------------------------------------------
# Step 3: TBLASTX of the candidate reads against the HPV database
# (an empty candidate file gives an empty result instead of an error)
# ------------------------------------------------------------
rule hpv_tblastx:
    input:
        fasta=f"{CAND_DIR}/{{sample}}.candidates.fa",
        db=multiext(HPV_DB, ".nhr", ".nin", ".nsq"),
    output:
        f"{CAND_DIR}/{{sample}}.tblastx.tsv",
    threads: 4
    resources:
        mem_mb=2000,
    log:
        f"{LOG_DIR}/hpv_tblastx/{{sample}}.log",
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname {log})"
        if [ -s {input.fasta} ]; then
            tblastx -query {input.fasta} -db {HPV_DB} \
                -evalue {TBLASTX_EVALUE} -max_target_seqs 1000 -max_hsps 1 \
                -num_threads {threads} \
                -outfmt "6 qseqid sseqid bitscore evalue" \
                > {output} 2> {log}
        else
            : > {output}
        fi
        """


# ------------------------------------------------------------
# Step 4: confirm strict BLASTN assignments with TBLASTX; the rest goes to GenBank
# ------------------------------------------------------------
rule hpv_stage2:
    input:
        stage1=f"{CAND_DIR}/{{sample}}.stage1.tsv",
        tbx=f"{CAND_DIR}/{{sample}}.tblastx.tsv",
        fasta=f"{CAND_DIR}/{{sample}}.candidates.fa",
    output:
        confirmed=f"{CAND_DIR}/{{sample}}.confirmed.tsv",
        pending=f"{CAND_DIR}/{{sample}}.pending.fa",
    run:
        best = {}   # qid -> (best bit score, set of refs)
        with open(input.tbx) as fh:
            for line in fh:
                f = line.rstrip("\n").split("\t")
                if len(f) < 3:
                    continue
                qid, ref, bits = f[0], clean(f[1]), float(f[2])
                b = best.get(qid)
                if b is None or bits > b[0]:
                    best[qid] = (bits, {ref})
                elif bits == b[0]:
                    b[1].add(ref)

        seqs = read_fasta(input.fasta)
        with open(input.stage1) as fh, \
                open(output.confirmed, "w") as conf, open(output.pending, "w") as pend:
            conf.write("qid\tmate\ttype\n")
            fh.readline()
            for line in fh:
                qid, mate, ref, pident, qcov, bits, strict, amb = line.rstrip("\n").split("\t")
                if amb == "1":
                    continue
                if strict == "1" and qid in best and best[qid][1] & set(ref.split(",")):
                    conf.write(f"{qid}\t{mate}\t{ref.split(',')[0]}\n")
                else:
                    pend.write(f">{qid}\n{seqs[qid]}\n")


# ------------------------------------------------------------
# Step 5: GenBank check (all samples together; identical reads are sent once)
#
# genbank.mode = "remote": NCBI remote BLAST. Reads are sent in batches with a pause
#   between them to respect NCBI usage limits. Answers are cached in
#   results/hpv_counts/genbank_cache.tsv, so an interrupted run resumes and a rerun
#   never asks about the same read sequence twice. No API key is used or required.
# genbank.mode = "local": a local BLAST database (genbank.db = its path).
# genbank.mode = "off": skip; those reads are reported as "unverified".
# ------------------------------------------------------------
rule hpv_genbank:
    input:
        pending=expand(f"{CAND_DIR}/{{sample}}.pending.fa", sample=SAMPLES),
    output:
        f"{COUNT_DIR}/hpv_genbank_results.tsv",
    threads: 4
    log:
        f"{LOG_DIR}/hpv_genbank.log",
    run:
        import os, tempfile, time

        cache_path = f"{COUNT_DIR}/genbank_cache.tsv"
        cache = {}
        if os.path.exists(cache_path):
            with open(cache_path) as fh:
                for line in fh:
                    f = line.rstrip("\n").split("\t")
                    if len(f) == len(GB_COLS) and f[0] != "seqhash":
                        cache[f[0]] = f

        seqs = {}   # seqhash -> sequence
        for p in input.pending:
            for name, seq in read_fasta(p).items():
                seqs[seqhash(seq)] = seq.upper()

        def log_msg(msg):
            with open(log[0], "a") as lf:
                lf.write(msg + "\n")
            print(msg, flush=True)

        def query(batch):
            with tempfile.NamedTemporaryFile("w", suffix=".fa", delete=False) as t:
                for h in batch:
                    t.write(f">{h}\n{seqs[h]}\n")
            cmd = ["blastn", "-query", t.name, "-db", GENBANK_DB,
                   "-evalue", str(GENBANK_EVALUE),
                   "-max_target_seqs", "10", "-max_hsps", "1",
                   "-outfmt", "6 qseqid sseqid pident length qstart qend evalue bitscore qlen stitle"]
            cmd += ["-remote"] if GENBANK_MODE == "remote" else ["-num_threads", str(threads)]
            for attempt in range(1, 4):
                res = subprocess.run(cmd, capture_output=True, text=True)
                if res.returncode == 0:
                    os.unlink(t.name)
                    return res.stdout
                log_msg(f"blastn attempt {attempt} failed: {res.stderr.strip()[:300]}")
                time.sleep(30 * attempt)
            os.unlink(t.name)
            raise RuntimeError("GenBank BLAST failed 3 times; see " + log[0])

        todo = [h for h in seqs if h not in cache]
        if GENBANK_MODE != "off" and todo:
            log_msg(f"{len(todo)} unique reads to check against GenBank ({GENBANK_MODE})")
            with open(cache_path, "a") as cf:
                if os.path.getsize(cache_path) == 0:
                    cf.write("\t".join(GB_COLS) + "\n")
                for i in range(0, len(todo), GENBANK_BATCH):
                    batch = todo[i:i + GENBANK_BATCH]
                    best = {}   # seqhash -> (bits, pident, qcov, title)
                    for line in query(batch).splitlines():
                        f = line.split("\t", 9)
                        if len(f) < 10:
                            continue
                        h, pident, qs, qe, bits, qlen = f[0], float(f[2]), int(f[4]), int(f[5]), float(f[7]), int(f[8])
                        qcov = (abs(qe - qs) + 1) / qlen
                        if h not in best or bits > best[h][0]:
                            best[h] = (bits, pident, qcov, f[9])
                    for h in batch:
                        if h in best:
                            bits, pident, qcov, title = best[h]
                            cls, atype = classify_genbank_hit(title, pident, qcov)
                            row = [h, cls, atype, f"{pident:.2f}", f"{qcov:.3f}", f"{bits:.1f}",
                                   title.replace("\t", " ").replace("\n", " ")]
                        else:
                            row = [h, "no_hit", "", "", "", "", ""]
                        cache[h] = row
                        cf.write("\t".join(row) + "\n")
                    cf.flush()
                    log_msg(f"checked {min(i + GENBANK_BATCH, len(todo))}/{len(todo)}")
                    if GENBANK_MODE == "remote" and i + GENBANK_BATCH < len(todo):
                        time.sleep(GENBANK_SLEEP)

        with open(output[0], "w") as out:
            out.write("\t".join(GB_COLS) + "\n")
            for h in seqs:
                row = cache.get(h, [h, "unverified", "", "", "", "", ""])
                out.write("\t".join(row) + "\n")


# ------------------------------------------------------------
# Step 6: raw count matrix (samples x HPV types, plus an "unclassified" column)
#
# Final label of each candidate read:
#   confirmed by BLASTN + TBLASTX            -> its type
#   sent to GenBank: typed                   -> that type
#   sent to GenBank: HPV but not typed       -> "unclassified"
#   sent to GenBank: non-HPV / no hit        -> discarded
#   GenBank check off or not answered        -> "unverified" (not counted)
#   tied between types with a strict hit     -> listed in hpv_ambiguous_fragments.tsv
#
# count_unit:
#   "fragments": a pair counts once; if its mates carry two different types the
#                pair is listed as a conflict instead of being counted.
#   "reads":     every labelled read counts (a pair can count twice).
#   "reads_r1":  only forward reads count (as in the paper).
# Genotypes below min_fragments in a sample are set to 0 in the main matrix; the
# unfiltered matrix is written alongside.
# ------------------------------------------------------------
rule hpv_count_matrix:
    input:
        stage1=expand(f"{CAND_DIR}/{{sample}}.stage1.tsv", sample=SAMPLES),
        confirmed=expand(f"{CAND_DIR}/{{sample}}.confirmed.tsv", sample=SAMPLES),
        pending=expand(f"{CAND_DIR}/{{sample}}.pending.fa", sample=SAMPLES),
        gb=f"{COUNT_DIR}/hpv_genbank_results.tsv",
        fasta=HPV_FASTA,
    output:
        frags=f"{COUNT_DIR}/hpv_fragment_counts.tsv",
        frags_raw=f"{COUNT_DIR}/hpv_fragment_counts_unfiltered.tsv",
        frags_amb=f"{COUNT_DIR}/hpv_ambiguous_fragments.tsv",
        fate=f"{COUNT_DIR}/hpv_read_fate.tsv",
    resources:
        mem_mb=4000,
    run:
        import os
        from collections import Counter

        if COUNT_UNIT not in ("fragments", "reads", "reads_r1"):
            raise ValueError("count_unit must be 'fragments', 'reads' or 'reads_r1'")

        # Columns: every sequence in the HPV FASTA, then anything else found
        refs = []
        with open(input.fasta) as fh:
            for line in fh:
                if line.startswith(">"):
                    refs.append(clean(line[1:].split()[0]))
        refs = list(dict.fromkeys(refs))

        gb = {}
        with open(input.gb) as fh:
            fh.readline()
            for line in fh:
                f = line.rstrip("\n").split("\t")
                gb[f[0]] = (f[1], f[2])

        counts = {s: Counter() for s in SAMPLES}
        amb = Counter()                    # (sample, kind, "A,B") -> n
        fate_cols = ["candidates", "ambiguous_ties", "confirmed_blast_tblastx",
                     "sent_to_genbank", "genbank_typed", "genbank_unclassified",
                     "genbank_non_hpv", "genbank_no_hit", "unverified"]
        fate = {s: Counter() for s in SAMPLES}

        for s in SAMPLES:
            confirmed = {}
            with open(f"{CAND_DIR}/{s}.confirmed.tsv") as fh:
                fh.readline()
                for line in fh:
                    qid, mate, typ = line.rstrip("\n").split("\t")
                    confirmed[qid] = typ
            pending_hash = {n: seqhash(q) for n, q in read_fasta(f"{CAND_DIR}/{s}.pending.fa").items()}

            labels = {}                    # qid -> type or "unclassified"
            with open(f"{CAND_DIR}/{s}.stage1.tsv") as fh:
                fh.readline()
                for line in fh:
                    qid, mate, ref, pident, qcov, bits, strict, ambig = line.rstrip("\n").split("\t")
                    fate[s]["candidates"] += 1
                    if ambig == "1":
                        fate[s]["ambiguous_ties"] += 1
                        amb[(s, "read_tie", ref)] += 1
                        continue
                    if qid in confirmed:
                        fate[s]["confirmed_blast_tblastx"] += 1
                        labels[qid] = confirmed[qid]
                        continue
                    fate[s]["sent_to_genbank"] += 1
                    cls, atype = gb.get(pending_hash[qid], ("unverified", ""))
                    if cls == "typed":
                        fate[s]["genbank_typed"] += 1
                        labels[qid] = atype
                    elif cls == "unclassified":
                        fate[s]["genbank_unclassified"] += 1
                        labels[qid] = "unclassified"
                    elif cls == "non_hpv":
                        fate[s]["genbank_non_hpv"] += 1
                    elif cls == "no_hit":
                        fate[s]["genbank_no_hit"] += 1
                    else:
                        fate[s]["unverified"] += 1

            if COUNT_UNIT == "reads":
                for lab in labels.values():
                    counts[s][lab] += 1
            elif COUNT_UNIT == "reads_r1":
                for qid, lab in labels.items():
                    if qid.endswith("/1"):
                        counts[s][lab] += 1
            else:
                pairs = {}
                for qid, lab in labels.items():
                    pairs.setdefault(qid.rsplit("/", 1)[0], set()).add(lab)
                for pair, labs in pairs.items():
                    typed = {l for l in labs if l != "unclassified"}
                    if len(typed) == 1:
                        counts[s][next(iter(typed))] += 1
                    elif len(typed) > 1:
                        amb[(s, "pair_conflict", ",".join(sorted(typed)))] += 1
                    else:
                        counts[s]["unclassified"] += 1

        found = {c for s in SAMPLES for c in counts[s]} - {"unclassified"}
        extra = sorted(found - set(refs))
        if extra:
            print("WARNING: HPV types not in the HPV FASTA (added as extra columns):", extra)
        cols = refs + extra + ["unclassified"]

        def kept(s, c):
            v = counts[s][c]
            return v if v >= MIN_FRAGMENTS else 0

        os.makedirs(os.path.dirname(output.frags), exist_ok=True)
        with open(output.frags_raw, "w") as out:
            out.write("sample\t" + "\t".join(cols) + "\n")
            for s in SAMPLES:
                out.write(s + "\t" + "\t".join(str(counts[s][c]) for c in cols) + "\n")

        with open(output.frags, "w") as out:
            out.write("sample\t" + "\t".join(cols) + "\n")
            for s in SAMPLES:
                out.write(s + "\t" + "\t".join(str(kept(s, c)) for c in cols) + "\n")

        with open(output.frags_amb, "w") as out:
            out.write("sample\tkind\ttied_or_conflicting_types\tn\n")
            for (s, kind, types), n in sorted(amb.items(), key=lambda x: (x[0][0], -x[1])):
                out.write(f"{s}\t{kind}\t{types}\t{n}\n")

        with open(output.fate, "w") as out:
            out.write("sample\t" + "\t".join(fate_cols) + "\n")
            for s in SAMPLES:
                out.write(s + "\t" + "\t".join(str(fate[s][c]) for c in fate_cols) + "\n")


# ------------------------------------------------------------
# Step 7a: normalisation per million RAW read pairs
# Denominator = raw read pairs from fastp (before trimming and host removal):
# fastp "before_filtering total_reads" counts both mates, so pairs = total_reads / 2.
# ------------------------------------------------------------
rule hpv_normalize:
    input:
        frags=f"{COUNT_DIR}/hpv_fragment_counts.tsv",
        fastp=expand(f"{FASTP_DIR}/{{sample}}.fastp.json", sample=SAMPLES),
    output:
        sizes=f"{COUNT_DIR}/library_sizes.tsv",
        frags=f"{COUNT_DIR}/hpv_fragment_counts_per_million_pairs.tsv",
    run:
        import json

        total_reads = {}
        for s in SAMPLES:
            with open(f"{FASTP_DIR}/{s}.fastp.json") as fh:
                total_reads[s] = json.load(fh)["summary"]["before_filtering"]["total_reads"]
        total_pairs = {s: total_reads[s] / 2 for s in SAMPLES}

        with open(output.sizes, "w") as out:
            out.write("sample\traw_reads\traw_read_pairs\n")
            for s in SAMPLES:
                out.write(f"{s}\t{total_reads[s]}\t{int(total_pairs[s])}\n")

        with open(input.frags) as fh, open(output.frags, "w") as out:
            out.write(fh.readline())          # header
            for line in fh:
                f = line.rstrip("\n").split("\t")
                vals = [float(x) / total_pairs[f[0]] * 1e6 for x in f[1:]]
                out.write(f[0] + "\t" + "\t".join(f"{v:.4f}" for v in vals) + "\n")


# ------------------------------------------------------------
# Step 7b: number of non-human read pairs per sample (pairs left after host
# removal; one line of R1 = one pair), then normalisation per million
# NON-HUMAN read pairs. Only built when hpv.nonhost_norm is true.
# The two normalisations answer different questions: raw pairs reflect sequencing
# depth; non-human pairs reflect depth after removing host reads.
# ------------------------------------------------------------
rule count_nonhost_pairs:
    input:
        r1=f"{HOST_DIR}/{{sample}}_1.nohost.fastq.gz",
    output:
        f"{COUNT_DIR}/nonhost_pairs/{{sample}}.txt",
    threads: 2
    resources:
        mem_mb=500,
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname {output})"
        n=$(pigz -dc -p {threads} {input.r1} | wc -l)
        echo $((n / 4)) > {output}
        """


rule hpv_normalize_nonhost:
    input:
        frags=f"{COUNT_DIR}/hpv_fragment_counts.tsv",
        sizes=expand(f"{COUNT_DIR}/nonhost_pairs/{{sample}}.txt", sample=SAMPLES),
    output:
        sizes=f"{COUNT_DIR}/nonhost_library_sizes.tsv",
        frags=f"{COUNT_DIR}/hpv_fragment_counts_per_million_nonhost_pairs.tsv",
    run:
        nonhost = {}
        for s in SAMPLES:
            with open(f"{COUNT_DIR}/nonhost_pairs/{s}.txt") as fh:
                nonhost[s] = int(fh.read().strip())

        with open(output.sizes, "w") as out:
            out.write("sample\tnonhost_read_pairs\n")
            for s in SAMPLES:
                out.write(f"{s}\t{nonhost[s]}\n")

        with open(input.frags) as fh, open(output.frags, "w") as out:
            out.write(fh.readline())          # header
            for line in fh:
                f = line.rstrip("\n").split("\t")
                n = nonhost[f[0]]
                if n == 0:
                    out.write(f[0] + "\t" + "\t".join("NA" for _ in f[1:]) + "\n")
                else:
                    vals = [float(x) / n * 1e6 for x in f[1:]]
                    out.write(f[0] + "\t" + "\t".join(f"{v:.4f}" for v in vals) + "\n")


# ------------------------------------------------------------
# Step 7c: proportional genotype table (% of HPV fragments per sample)
# proportion = fragments assigned to genotype / all fragments assigned to any HPV
# genotype in that sample x 100. Fragments tied between references are not in the
# denominator. Samples with no HPV fragments get NA. Uses the main matrix, so
# min_fragments applies here too.
# ------------------------------------------------------------
rule hpv_proportions:
    input:
        f"{COUNT_DIR}/hpv_fragment_counts.tsv",
    output:
        f"{COUNT_DIR}/hpv_genotype_proportions_percent.tsv",
    run:
        with open(input[0]) as fh, open(output[0], "w") as out:
            out.write(fh.readline())          # header
            for line in fh:
                f = line.rstrip("\n").split("\t")
                vals = [int(x) for x in f[1:]]
                total = sum(vals)
                if total == 0:
                    out.write(f[0] + "\t" + "\t".join("NA" for _ in vals) + "\n")
                else:
                    out.write(f[0] + "\t" + "\t".join(f"{v / total * 100:.2f}" for v in vals) + "\n")
