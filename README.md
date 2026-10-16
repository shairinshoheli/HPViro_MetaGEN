# HPViro_MetaGEN
A complete pipeline for processing shotgun metagenomic data to look at viruses such as HPV and the composition of the vaginal bacteriome
Snakemake workflow that takes SRA files of paired-end shotgun metagenomes and produces:

- QC reports (FastQC, fastp) and human-read-depleted FASTQ
- HPV genotype counts, using a BLASTN → TBLASTX → GenBank check scheme modelled on Ma et al. 2014 (doi:10.1128/jvi.00093-14)
- VIRGO2 taxonomic abundances, read counts and gene annotation
- VISTA metagenomic subspecies (mgSs) and community state types (mgCST)

## Workflow steps

| Step | Rules file | What it does |
|---|---|---|
| 1 | `workflow/rules/01_fastq_qc.smk` | `fasterq-dump` → FastQC (raw) → fastp trimming |
| 2 | `workflow/rules/02_host_removal.smk` | bowtie2 vs GRCh38; keeps only pairs with both mates unmapped |
| 3 | `workflow/rules/03_hpv.smk` | HPV BLASTN, TBLASTX confirmation, optional GenBank check, count matrices, normalisation |
| 4 | `workflow/rules/04_virgo2.smk` | VIRGO2 install, map, compile, taxonomy, gene annotation |
| 5 | `workflow/rules/05_vista.smk` | VISTA install and mgSs/mgCST classification |

Shared settings and helpers are in `workflow/rules/common.smk`; all user settings are in `config/config.yaml`.

## Environment

```bash
conda env create -f environment.yml    # or: mamba env create -f environment.yml
conda activate hpv-virgo2-vista
```

`environment.yml` installs Snakemake and the command-line tools (sra-tools, FastQC, fastp, bowtie2, samtools, BLAST+, pigz, git-lfs, R). Versions are not pinned; after a working run, record yours with `conda env export --no-builds > environment.lock.yml`.

Not checked here: VIRGO2's Python dependencies and VISTA's R packages. Install whatever their READMEs list (https://github.com/ravel-lab/VIRGO2, https://github.com/JHolm-Lab/VISTA) into this environment.

## Inputs you must supply

Place these under the project root (paths can be changed in `config/config.yaml`; none are included in the repo):

1. **SRA files**: `sra/SRRxxxxxxx/SRRxxxxxxx.sra`. Every sample found there is processed.
2. **Human reference**: `reference/GRCh38/GRCh38.primary_assembly.genome.fa` (the bowtie2 index is built by the workflow).
3. **HPV BLAST database**: The fasta files were take from the PaVE database by NCBI (https://pave.niaid.nih.gov/) and the database was built: `reference/hpvrefs_db/hpv_db.{nhr,nin,nsq}` plus the FASTA it was built from, `hpv_all.fasta` (its sequence names define the count-matrix columns). Build with `makeblastdb -in hpv_all.fasta -dbtype nucl -out reference/hpvrefs_db/hpv_db`; add `-parse_seqids` only if your sequence IDs are compatible.
4. **VISTA data** (only if the automatic download fails): download `VISTA_data.tar.gz` in a browser from https://doi.org/10.6084/m9.figshare.28684934.v10 and save it as `reference/VISTA_data.tar.gz`.

VIRGO2 and VISTA are cloned into `reference/` by the workflow. If you installed either by hand, create its flag so it is not re-cloned: `touch reference/VIRGO2.installed reference/VISTA.installed`.

## Running

Run from the project root, not from `workflow/`:

```bash
snakemake -n --cores 14                                    # dry run
snakemake --cores 14 --resources mem_mb=48000 --rerun-incomplete --keep-going
```

Set `--cores` and `mem_mb` to your machine. Per-rule thread and memory requests are in the rule files. Disk: raw and trimmed FASTQ are deleted automatically (`temp()`), but `tmp_dir` and `results/` need substantial space.

## Key settings (`config/config.yaml`)

- `hpv.min_pident`, `hpv.min_qcov`, `hpv.max_evalue`: thresholds for a "strict" BLASTN assignment (defaults 95 %, full read coverage, 1e-10).
- `hpv.count_unit`: `fragments` (a pair counts once), `reads`, or `reads_r1` (forward reads only, as in the paper).
- `hpv.min_fragments`: minimum fragments per genotype per sample; the unfiltered matrix is always written too.
- `hpv.genbank.mode`: `off` (default; unconfirmed reads are reported as "unverified"), `remote` (NCBI remote BLAST, slow, needs internet, no API key used) or `local` (set `hpv.genbank.db` to your local database).
- `fastp_qc`: fastp quality options (follow the VIRGO2 methods).

The TBLASTX E-value, GenBank E-value, tie-exclusion rule and GenBank title-matching are choices of this workflow, not specified in Ma et al. 2014.

## Main outputs

| Path | Content |
|---|---|
| `results/fastqc/`, `results/fastp/` | QC reports |
| `results/host_removed/` | host-depleted FASTQ |
| `results/hpv_counts/hpv_fragment_counts.tsv` | samples × HPV genotypes (+ `unclassified`) |
| `results/hpv_counts/hpv_fragment_counts_unfiltered.tsv` | same, without `min_fragments` filtering |
| `results/hpv_counts/hpv_read_fate.tsv` | how many reads took each classification route |
| `results/hpv_counts/hpv_ambiguous_fragments.tsv` | tied reads and conflicting mate pairs |
| `results/hpv_counts/*per_million*.tsv`, `*proportions_percent.tsv` | normalised tables |
| `results/virgo2/taxonomy/` | VIRGO2 relative abundances and read counts |
| `results/virgo2/annotation/` | per-gene counts with annotation tables merged |
| `results/vista/` | mgSs / mgCST classifications |

## Known caveats

- VIRGO2 `map` takes single-end input, so R1 and R2 are concatenated; the workflow notes that identical read names across mates should be checked against VIRGO2's map script.
- `VIRGO2.py install` is reported (in the rule comments) to fail at its own index step; the rule builds the bowtie2 index itself.
- `VIRGO2 compile` reads every result in `results/virgo2/map/`; do not mix projects in that folder.

## Data and privacy

No data, credentials or API keys are in this repo (`.gitignore` excludes `sra/`, `reference/`, `results/`, FASTQ, and archives). If any samples come from controlled-access data, do not commit derived per-sample tables without checking the data-use terms.

## Citation

Please cite the tools used (Snakemake, fastp, FastQC, bowtie2, samtools, BLAST+, VIRGO2, VISTA) and Ma et al. 2014 for the HPV scheme. Add a `LICENSE` file before publishing.
