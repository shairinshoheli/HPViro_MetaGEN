# ============================================================
# Shared setup: config values, sample detection, directories, helpers.
# Included first by the Snakefile; every later rules file uses these names.
# ============================================================

import os
import glob
import gzip
import hashlib
import re
import subprocess

READS = [1, 2]

# ---------- Inputs ----------
SRA_DIR = config["sra_dir"]
TMP_DIR = config["tmp_dir"]
HOST_FASTA = config["host"]["fasta"]
HOST_INDEX = config["host"]["index"]

# All SRA files found as sra/SRR_ACCESSION/SRR_ACCESSION.sra are used.
SAMPLES = sorted(
    os.path.basename(os.path.dirname(path))
    for path in glob.glob(f"{SRA_DIR}/SRR*/*.sra")
)

if not SAMPLES:
    raise ValueError(
        f"No SRA files found under {SRA_DIR}/. Expected files like "
        f"{SRA_DIR}/SRR1234567/SRR1234567.sra"
    )

print(f"Automatically detected {len(SAMPLES)} SRA samples: {SAMPLES}")

# ---------- Output directories ----------
FASTQ_DIR   = "results/fastq"
QC_DIR      = "results/fastqc"
QC_TRIM_DIR = "results/fastqc_trimmed"
TRIM_DIR    = "results/trimmed"
FASTP_DIR   = "results/fastp"
HOST_DIR    = "results/host_removed"
BT2_DIR     = "results/bowtie2"
BLAST_DIR   = "results/blast_hpv"
COUNT_DIR   = "results/hpv_counts"
CAND_DIR    = "results/hpv_candidates"
LOG_DIR     = "logs"

# ---------- fastp / FastQC ----------
FASTP_QC       = config["fastp_qc"]
TRIMMED_FASTQC = bool(config["trimmed_fastqc"])

# ---------- HPV ----------
_hpv = config["hpv"]
HPV_DB         = _hpv["db"]
HPV_FASTA      = _hpv["fasta"]
MIN_PIDENT     = float(_hpv["min_pident"])
MIN_QCOV       = float(_hpv["min_qcov"])
MIN_BITSCORE   = float(_hpv["min_bitscore"])
MAX_EVALUE     = float(_hpv["max_evalue"])
MIN_FRAGMENTS  = int(_hpv["min_fragments"])
NONHOST_NORM   = bool(_hpv["nonhost_norm"])
COUNT_UNIT     = _hpv["count_unit"]
TBLASTX_EVALUE = float(_hpv["tblastx_evalue"])

_gb = _hpv["genbank"]
GENBANK_MODE   = str(_gb["mode"])
GENBANK_DB     = str(_gb["db"])
GENBANK_EVALUE = float(_gb["evalue"])
GENBANK_BATCH  = int(_gb["batch"])
GENBANK_SLEEP  = int(_gb["sleep"])

# ---------- VIRGO2 / VISTA ----------
VIRGO2_HOME  = config["virgo2"]["home"]
VIRGO2_READY = "reference/VIRGO2.installed"   # flag written when install succeeds
ANNOTATE_GENES        = bool(config["virgo2"]["annotate"])
ANNOT_MIN_TOTAL_COUNT = int(config["virgo2"]["annot_min_total_count"])
ANNOT_TABLE_IDS       = list(config["virgo2"]["annot_table_ids"])

VISTA_HOME     = config["vista"]["home"]
VISTA_READY    = "reference/VISTA.installed"
VISTA_ARCHIVE  = config["vista"]["archive"]
VISTA_DATA_URL = config["vista"]["data_url"]

VIRGO2_DIR      = "results/virgo2"
VIRGO2_MAP_DIR  = f"{VIRGO2_DIR}/map"        # per-sample mapping output (compile reads this whole dir)
VIRGO2_FLAG_DIR = f"{VIRGO2_DIR}/flags"      # kept out of the map dir on purpose
VIRGO2_COMP_DIR = f"{VIRGO2_DIR}/compiled"
VIRGO2_TAX_DIR  = f"{VIRGO2_DIR}/taxonomy"
VIRGO2_ANNOT_DIR = f"{VIRGO2_DIR}/annotation"
VISTA_DIR       = "results/vista"

# ---------- Create input/result/log directories at parse time ----------
# NOT created on purpose: VIRGO2_TAX_DIR and VISTA_DIR (declared as directory()
# outputs; if they already existed Snakemake would treat those rules as done)
# and VIRGO2_HOME / VISTA_HOME (created by git clone).
for _d in [
    SRA_DIR,
    os.path.dirname(HOST_FASTA),
    os.path.dirname(HOST_INDEX),
    os.path.dirname(HPV_DB),
    os.path.dirname(HPV_FASTA),
    FASTQ_DIR, QC_DIR, QC_TRIM_DIR, TRIM_DIR, FASTP_DIR, HOST_DIR, BT2_DIR,
    BLAST_DIR, COUNT_DIR,
    VIRGO2_MAP_DIR, VIRGO2_FLAG_DIR, VIRGO2_COMP_DIR,
    TMP_DIR,
    LOG_DIR,
    *[f"{LOG_DIR}/{_r}" for _r in (
        "sra_to_fastq", "fastqc", "fastp", "fastqc_trimmed", "blast_hpv",
        "virgo2_map",
    )],
]:
    os.makedirs(_d, exist_ok=True)


wildcard_constraints:
    sample=r"SRR\d+",
    read=r"[12]",


# ---------- Optional targets ----------
TRIMMED_FASTQC_TARGETS = (
    expand(f"{QC_TRIM_DIR}/{{sample}}_{{read}}.trimmed_fastqc.html",
           sample=SAMPLES, read=READS)
    if TRIMMED_FASTQC else []
)
ANNOTATE_TARGETS = (
    [f"{VIRGO2_ANNOT_DIR}/virgo2_genes_annotated.tsv"] if ANNOTATE_GENES else []
)
NONHOST_NORM_TARGETS = (
    [f"{COUNT_DIR}/hpv_fragment_counts_per_million_nonhost_pairs.tsv"]
    if NONHOST_NORM else []
)
