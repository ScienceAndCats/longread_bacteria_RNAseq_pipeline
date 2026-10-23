# Whiteley Bacteria Mapping Pipeline

This repository contains a long-read FASTQ processing pipeline for Oxford Nanopore bacterial sequencing runs. Each FASTQ file is an independent sample, including files whose names contain `_R1`, `_R2`, `_1`, or `_2`. The pipeline trims user-configured adapter sequences, separates reads over 6000 bp after trimming, uses BBSplit to remove reads that map to a decoy/pangenome reference, then maps remaining reads to bacterial and optional host RNA references before genomic DNA references. It discovers references and annotations from folders, counts gene-level alignments, calculates bacterial coverage metrics, and reports the read counts and percentages in CDS, rRNA, other RNA, and intergenic/unannotated regions.

## What the program does

`map_bacteria_with_decoys.sh` runs the analysis in these stages:

1. Finds input FASTQ files matching the configured glob and, when sampling is enabled, randomly selects up to the configured number of reads from each file.
2. Uses Cutadapt to remove user-configured adapter sequences from both read ends, repeating the process for up to `ADAPTER_MAX_TRIMMING_ROUNDS` rounds (default: `3`). Each round removes at most one adapter from the 5′ end and one from the 3′ end, allowing repeated terminal copies to be exposed and removed in later rounds. The minimum read length is applied after all rounds. Reads from `MIN_READ_LENGTH` through **6000 bp inclusive** are kept for alignment; longer reads are saved in a separate FASTQ without further processing and their counts are included in both summary CSVs.
3. Uses `bbsplit.sh` to map eligible trimmed reads to a decoy/pangenome reference and keeps reads that do **not** map to the decoys.
4. Maps decoy-unmapped reads in this order: **bacterial RNA → host RNA (if enabled) → bacterial DNA → host DNA (if enabled)**. Only primary-unmapped reads advance to the next stage, so an RNA match takes priority over a later genomic DNA match. RNA stages without reference files are skipped. A final leftover FASTQ contains reads that failed every configured alignment stage. Every stage uses configurable `k` (default `11`) and fixed `mapmode=pb msa=MultiStateAligner9PacBio fastareadlen=6000`, with `qtrim=f` to preserve adapter-trimmed sequences and `interleaved=f` for single-read FASTQ input.
5. Combines the RNA and DNA alignments for each organism, assigns each mapped read to one exclusive category for the breakdown CSV, and uses `featureCounts` from Subread to count CDS features against the discovered, normalized bacterial and host annotations.
6. Uses `samtools` to create, sort, and index combined bacterial BAM files. Genomic coverage is calculated separately from the DNA-stage alignments only, so transcript coordinates do not inflate the genomic reference length.
7. Uses `gene_position_profile.py` on each sorted combined bacterial BAM to calculate strand-aware, normalized feature-position and aggregate metagene profiles against the combined RNA/DNA annotation.
8. Runs `bacteria_with_decoys_csvConversion.py` to combine cutadapt, length-partition, BBSplit, coverage, and featureCounts outputs into `BACTERIA_<project-directory>.csv`.
9. Writes `BACTERIA_<project-directory>_read_breakdown.csv`, with one sample per row and read counts, percentages of total reads, and the source of every value.

The summary conversion scripts use only the Python standard library. Reference preparation uses Biopython for GenBank annotations, read classification uses pysam, and positional profiling uses pysam and matplotlib.

## Dependencies

The pipeline environment includes these tools:

- Python 3.13.15
- cutadapt 5.2
- BBMap 39.79 (provides BBSplit / `bbsplit.sh`)
- samtools 1.24
- Subread / featureCounts 2.1.1
- pysam
- matplotlib
- biopython (reads `.gbff`/GenBank annotations)

A conda environment file is provided in `environment.yml`.
The shell pipeline also requires `sha256sum` (GNU coreutils) to identify reference files for its index cache.

## Create the conda environment

```bash
conda env create -f environment.yml
conda activate bacteria-RNAseq-pipeline
```

If you are running on an HPC system that requires module loading, load your site-specific conda or mamba module before creating or activating the environment.

## Configure inputs and references

Edit `config.env` before running the pipeline. The key settings are:

| Setting | Purpose |
| --- | --- |
| `FASTQ_DIR` | Directory containing input FASTQ files; each file is an independent sample. |
| `FASTQ_GLOB` | Shell glob for FASTQ inputs, such as `*.fastq`. |
| `SLURM_MAX_SAMPLES` | Maximum samples accepted per submission (default: `40`); larger input folders produce an error until this limit is increased. |
| `SLURM_MAX_CONCURRENT` | Maximum sample array tasks running at once (default: `8`). |
| `SLURM_ACCOUNT` | Required PACE allocation/account for every submitted job. |
| `SLURM_MAIL_USER` | Your notification email address; leave empty to disable job emails. |
| `SLURM_MAIL_TYPE` | Slurm notification events (default: `END,FAIL`); used when an email address is set. |
| `CONDA_SH` | Path to your conda installation's `etc/profile.d/conda.sh`. |
| `CONDA_ENV` | Conda environment name or path to activate on compute nodes (default: `bacteria-RNAseq-pipeline`). |
| `PIPELINE_SKIP_CONDA` | Set `true` only when your site's job setup already activates the required environment (default: `false`). |
| `OUTPUT_DIR` | Directory where output files should be written. |
| `SAMPLE_READS` | Set to `true` to analyze a random subset of each input, or `false` to analyze every read. |
| `SAMPLE_SIZE` | Maximum reads sampled from each FASTQ (default: `5000`). |
| `SAMPLE_SEED` | Seed used to make random sampling reproducible. |
| `DECOY_REFERENCE_DIR` | Folder containing the decoy/pangenome FASTA files. |
| `BACTERIA_REFERENCE_DIR` | Folder containing bacterial genomic DNA, CDS/RNA FASTA files, and their annotations, such as `PAO1`. |
| `HOST_REFERENCE_DIR` | Optional host folder containing DNA and RNA references, including nested folders and RNA shards; leave empty to disable host mapping and counting. |
| `BBSPLIT_K` | BBSplit k-mer length: integer `1`–`15` (default: `11`). All stages use fixed `mapmode=pb msa=MultiStateAligner9PacBio fastareadlen=6000 qtrim=f interleaved=f`. |
| `ADAPTER_SEQUENCES` | Adapter DNA sequences separated by whitespace, on one or multiple lines. Default: `TTTCTGTTGGTGCTGATATTGC`. Set to `""` to skip adapter matching. |
| `ADAPTER_REVERSE_COMPLEMENTS` | Also match the reverse complement of each configured adapter (default: `true`); read orientation remains unchanged. |
| `ADAPTER_MAX_TRIMMING_ROUNDS` | Maximum rounds of trimming at both ends (default: `3`); each round removes at most one adapter per end. |
| `CUTADAPT_ERROR_RATE` | Maximum adapter alignment error rate (default: `0.15`); mismatches and indels are allowed. |
| `CUTADAPT_MIN_OVERLAP` | Minimum adapter overlap, in nucleotides (default: `12`). |
| `MIN_READ_LENGTH` | Minimum retained length **after all trimming**, from `1` through `6000` (default: `20`). Existing project configs with explicit values such as `22` keep that value. The alignment maximum is fixed at `6000` bp. |
| `THREADS` | Number of threads used by every multithreaded step (default: `16`). |
| `GENE_POSITION_BINS` | Number of equal normalized 5'-to-3' bins (default: `100`). |
| `METAGENE_MIN_FEATURE_READS` | Assigned reads required for a feature to contribute to aggregate profiles (default: `10`). |
| `CSV_CONVERSION_SCRIPT` | Path to `bacteria_with_decoys_csvConversion.py`. |

### Quick sampling mode

Set `SAMPLE_READS="true"` in `config.env` for a quick exploratory run. Before trimming or mapping, the pipeline uses reservoir sampling to select up to `SAMPLE_SIZE` complete long-read records from each input. Files with 5,000 reads or fewer are used in full with the default setting. The temporary sampled inputs are written under `OUTPUT_DIR/.bacteria_sampled_fastq`; original FASTQ files are never modified. Sampling is reproducible for the same input paths and `SAMPLE_SEED`. Set `SAMPLE_READS="false"` for a full analysis.

### Reference folders

Point each setting to a folder, without a filename or accession root. For the example reference tree:

```bash
DECOY_REFERENCE_DIR="/storage/home/hcoda1/2/hwilms3/r-mwhiteley3-0/referencegenomes/decoy"
BACTERIA_REFERENCE_DIR="/storage/home/hcoda1/2/hwilms3/r-mwhiteley3-0/referencegenomes/PAO1"
HOST_REFERENCE_DIR="/storage/home/hcoda1/2/hwilms3/r-mwhiteley3-0/referencegenomes/human"
# Use the mouse folder instead for mouse samples, or "" to disable host mapping.
```

`prepare_references.py` searches these folders recursively. FASTA files may end in `.fa`, `.fasta`, or `.fna`, with or without `.gz`; annotations may be GFF/GFF3 or GenBank/GBFF, also plain or gzip-compressed. Existing Bowtie2 and minimap2 indexes are ignored because BBSplit builds its own indexes.

- A flat bacterial folder can contain `*_genomic.fna.gz`, `*_genomic.gbff.gz` and/or `*_genomic.gff.gz`, `*_cds_from_genomic.fna.gz`, and `*_rna_from_genomic.fna.gz`.
- A host folder can contain `human_DNA/` and `human_RNA/`, or the corresponding mouse folders. All RNA shards, such as `human.1.rna.fna.gz` through `human.15.rna.fna.gz` with their `.rna.gbff.gz` annotations, are collected; no single shard needs to be selected in the config.
- Whole genomic FASTA files form the DNA reference. `cds_from_genomic` sequences form coding RNA references alongside `rna_from_genomic` sequences; they are excluded from the genomic DNA reference.

Reference preparation writes combined FASTA files, normalized GFF annotations with a `gene_id` attribute, and `bbsplit_alignments/references/manifest.json`. RNA sequence IDs receive an `RNA__` prefix so transcript coordinates remain distinct from genomic coordinates in combined alignment files. Reference FASTA and annotation records must describe the same sequences and genome version.

BBSplit indexes are cached under `OUTPUT_DIR/bbsplit_alignments/indexes/<stage>/k<BBSPLIT_K>/<reference-file-SHA256>`. Separate RNA and DNA stages, reference checksums, and the k-mer length prevent reuse of an index built from different sequences or parameters. Recreate or update the conda environment from `environment.yml` when upgrading an existing installation, including Biopython for GBFF support.

Existing configurations using `DECOY_BBSPLIT_REFERENCE`, `BACTERIA_BBSPLIT_REFERENCE`, and optional `HOST_BBSPLIT_REFERENCE` basenames remain supported as a fallback. New configurations should use the folder settings above. A legacy basename with only a genomic FASTA and GFF runs without an RNA stage.

The 6000-bp partition is applied after all adapter trimming and the minimum-length filter, so a read originally over 6000 bp can still be aligned if trimming brings it within range. The final Cutadapt log and JSON report the minimum-length filter before this partition. The per-sample `*.length_filter.tsv` records the eligible and oversized counts separately; oversized reads are excluded from every alignment stage and from final leftover reads. Samples with no eligible reads still appear in both summary CSVs with zero alignment and coverage counts.

### Adapter trimming

List the adapter sequences present in your Nanopore library. Each sequence is searched at both exposed read ends; an internal adapter-like motif is retained. For example:

```bash
ADAPTER_SEQUENCES="TTTCTGTTGGTGCTGATATTGC"
ADAPTER_REVERSE_COMPLEMENTS="true"
ADAPTER_MAX_TRIMMING_ROUNDS="3"
CUTADAPT_ERROR_RATE="0.15"
CUTADAPT_MIN_OVERLAP="12"
MIN_READ_LENGTH="20"
```

For multiple adapters, put the DNA sequences in the same whitespace-separated string; line breaks are allowed:

```bash
ADAPTER_SEQUENCES="TTTCTGTTGGTGCTGATATTGC
ACGTTGCAACGTTGCA"
```

The second sequence above is illustrative; replace it with an adapter used in your library. With `ADAPTER_REVERSE_COMPLEMENTS="true"`, the pipeline adds each sequence's reverse complement to the matching list automatically. It does not reverse-complement or reorient reads. Headers and base qualities remain associated with the original read, with trimmed bases and their qualities removed together.

In each round, Cutadapt first trims at most one 5′ adapter, then at most one 3′ adapter. A second or third round can remove another adapter copy exposed by the preceding round. `ADAPTER_MAX_TRIMMING_ROUNDS="3"` therefore permits up to three adapter removals per end. Increase the limit if the library contains more repeated terminal copies. The minimum-length filter runs only after all rounds, followed by the fixed 6000-bp alignment limit. Reads with no adapter match are retained if they pass the length filters. Setting `ADAPTER_SEQUENCES=""` skips matching and still applies those filters.

Each sample uses the `<sample>_<length>bp` output prefix:

| File | Contents |
| --- | --- |
| `*.adapter_roundN.5p.log` | Native Cutadapt report for the 5′ pass of round `N`. |
| `*.adapter_roundN.3p.log` | Native Cutadapt report for the 3′ pass of round `N`. |
| `*.cutadapt_log.txt` and `*.cutadapt.json` | Native Cutadapt reports for the final minimum-length filter. Adapter matches are recorded in the round logs, not accumulated in these final filter reports. |

To trim a FASTQ without any reference mapping, run from the repository:

```bash
bash trim_fastq.sh config.env reads.fastq.gz trimmed.fastq sample_20bp
```

When the optional local `tests/` fixtures are present, run them in the pipeline environment:

```bash
python3 -m unittest discover -s tests -v
```

The local checks cover trimming at both ends, multiple configured adapters, reverse-complement adapter matching without reorienting reads, repeated terminal copies, internal motifs, empty adapter lists, and final length filtering. Length-partition and CSV checks cover the 6000/6001-bp boundary, gzip output, eligible/oversized counts, and read-count conservation. When dependencies are installed, integration checks run BBSplit through feature counting, coverage, and positional profiling with host mapping enabled and disabled, empty samples, reference folders and RNA shards, GenBank annotations, RNA-before-DNA routing, and exclusive region counts and percentages.

To enable host mapping and feature counting, set `HOST_REFERENCE_DIR` to the `human` or `mouse` folder. Setting it to `""` skips host reference preparation, alignment, and counting. Host RNA is searched after bacterial RNA but before bacterial DNA; host DNA is searched last.

## Run the pipeline

To process the configured FASTQ files sequentially, run:

```bash
bash map_bacteria_with_decoys.sh config.env
```

You can keep multiple config files and pass the one for the current run:

```bash
bash map_bacteria_with_decoys.sh configs/project_a.env
```

### Submit sample jobs on Slurm

The included `run_bacteria_pipeline.sbatch` follows the example Slurm template:
one node, one task, 16 CPUs, 64 GB RAM, and a four-hour time limit. Set your
conda path, environment, PACE account, and notification email in `config.env`.
The submission helper passes account/email settings to Slurm for every job;
compute jobs source the saved configuration before activating conda. No edits
to the sbatch file or environment-variable exports are required for these settings.
The helper uses Python 3 on the login node before the compute jobs activate conda.

Keep the input/reference/output folders and cluster settings in `config.env`:

```bash
FASTQ_DIR="/path/to/reads"
FASTQ_GLOB="*.fastq*"
OUTPUT_DIR="/path/to/results"
SLURM_MAX_SAMPLES="40"
SLURM_MAX_CONCURRENT="8"
SLURM_ACCOUNT="YOUR_PACE_ACCOUNT"
SLURM_MAIL_USER="YOUR_USERNAME@example.edu"
SLURM_MAIL_TYPE="END,FAIL"
CONDA_SH="/path/to/YOUR_USERNAME/miniforge3/etc/profile.d/conda.sh"
CONDA_ENV="bacteria-RNAseq-pipeline"
PIPELINE_SKIP_CONDA="false"
```

Leave `SLURM_MAIL_USER` empty to disable emails. If your site's module/job setup
already activates the environment, set `PIPELINE_SKIP_CONDA="true"` in the config.
Otherwise, the helper checks the conda script path and environment setting before
submitting any jobs. The saved configuration also fixes these settings for queued
jobs if you later edit your project config.

Submit from your Slurm login node:

```bash
bash submit_bacteria_pipeline.sh config.env
```

The helper submits all three stages automatically and prints their job IDs:

1. A preparation job combines the reference shards and builds shared BBSplit
   indexes once.
2. A sample array depends on successful preparation. Each `SLURM_ARRAY_TASK_ID`
   selects one immutable manifest entry containing exactly one FASTQ file.
   With 40 samples, the array is `0-39%8`; smaller folders submit exactly their
   sample count. Every discovered sample is included. More than 40 samples is
   an error by default; increase `SLURM_MAX_SAMPLES` to accept a larger run.
3. A separate final job depends on successful completion of the entire array
   (`afterok`). It verifies every sample's completion marker and combines all
   sample rows and featureCounts columns into the final CSVs. It does not rerun
   alignments. A failed sample prevents the final job from publishing a complete
   result. Invalid dependencies are cancelled by Slurm.

The array's CPU allocation follows `THREADS`; the preparation job uses the same
allocation. Each sample writes to its own directory and reads the prepared
references/indexes without rebuilding or modifying them. The final job requests
one CPU and 8 GB RAM. Adjust memory/time in the template for your reference sizes
and cluster. Arrays may run fewer than eight tasks when cluster resources are
busy; eight is the concurrency limit.

Each FASTQ file creates its own sample job. Filenames such as `sample_R1.fastq.gz`, `sample_R2.fastq.gz`, `sample_1.fastq`, and `sample_2.fastq` remain independent samples with their original name stems; the helper does not combine or rename them. Separate lanes or chunks also remain separate sample jobs when their names differ. Sampling applies independently to each file, and `SAMPLE_SIZE` is its maximum number of reads. For `.fq`/`.fq.gz` files, set `FASTQ_GLOB="*.fq*"`.

Each submission creates a fresh directory, so results from different runs do
not mix:

```text
OUTPUT_DIR/slurm_runs/<run-id>/
├── run.json                          # frozen sample list, settings, job IDs
├── config.env                        # configuration snapshot
├── logs/                             # preparation, array, final job logs
├── shared/bbsplit_alignments/        # references and indexes shared by tasks
├── samples/<sample>/                 # per-sample reports, alignments, profiles
├── featurecounts_BACTERIA_*.txt/.csv  # all samples' feature columns
├── host/featurecounts_HOST_*.txt/.csv # optional combined host feature counts
├── BACTERIA_<run-id>.csv
├── BACTERIA_<run-id>_read_breakdown.csv
└── aggregation.completed.json
```

The two final summary CSVs contain one row per manifest sample, including empty
samples. Percentages retain the total-input-read denominator. References are
prepared from the submitted configuration; changing config files later does not
change queued jobs. Input paths, sizes, and modification times are recorded and
checked to detect changed input files. To start another run after a failure,
correct its cause and submit again; the new submission has its own run directory.

`--prepare-only`, `--prepared-dir`, and `--skip-summary` on
`map_bacteria_with_decoys.sh` are used internally by the Slurm worker. The normal
sequential command remains available without those flags.

### Local files and version control

`.gitignore` excludes local `work/` and `tests/` files, Python caches, personal
config copies, reads, references, alignments, scheduler logs, and pipeline output.
Existing work and test files have been removed from Git tracking and remain on
disk locally. The distributed `config.env`, pipeline source, and Slurm template
remain tracked. Tests can still be run locally with the command above when the
ignored `tests/` directory is present.

## Important outputs

- `*_<MIN_READ_LENGTH>bp.trim.fastq` — adapter-trimmed FASTQ files with final lengths from `MIN_READ_LENGTH` through 6000 bp inclusive, used as the alignment input.
- `*_<MIN_READ_LENGTH>bp.length_filter.tsv` — post-trimming partition counts: `input`, `eligible`, `oversized`, and `max_read_length`; input is the number passing the minimum-length filter and the maximum is fixed at `6000`.
- `bbsplit_alignments/oversized_reads/*_<MIN_READ_LENGTH>bp_over_6000bp.fastq.gz` — reads longer than 6000 bp after all adapter trimming, saved without alignment or other downstream processing.
- `bbsplit_alignments/decoy/` — decoy SAM files, BBSplit reports and diagnostic logs, and decoy-unmapped reads.
- `bbsplit_alignments/references/` — prepared RNA/DNA FASTA files, normalized annotations, and the reference manifest.
- `bbsplit_alignments/bacteria/` — combined bacterial SAM/BAM files, aggregate BBSplit reports, category reports, and coverage reports.
- `bbsplit_alignments/{bacteria,host}/{rna,dna}/` — separate molecular-stage alignments, reports, diagnostic logs, and unmapped FASTQ files.
- `bbsplit_alignments/host/` — optional host SAM files, BBSplit reports and diagnostic logs, and host featureCounts results.
- `*.bbsplit.txt` and `*.bbsplit.stderr.txt` in each alignment stage — stable TSV `input`/`aligned`/`unmapped` counts and raw BBSplit diagnostic logs, respectively. RNA/DNA stage reports describe the actual sequential inputs. The combined organism reports summarize both molecular stages: bacterial logical input is the decoy-unmapped total, and host logical input subtracts all bacterial alignments.
- `bbsplit_alignments/{bacteria,host}/*.categories.tsv` — exclusive `CDS`, `rRNA`, `other_RNA`, and `intergenic_unannotated` counts; their sum equals `aligned`.
- `bbsplit_alignments/leftover_reads/*_leftover.fastq.gz` — final eligible reads that mapped to none of the configured decoy, bacterial, or host references. These files are always created; when host mapping is disabled, they contain the bacterial-unmapped reads. Oversized reads are excluded.
- `featurecounts_BACTERIA_summary.txt` and `featurecounts_BACTERIA_summary.csv` — featureCounts results.
- `bbsplit_alignments/host/featurecounts_HOST_summary.txt` and `.csv` — optional, separate host featureCounts results.
- `bbsplit_alignments/bacteria/BACTERIA_*_coverage.txt` — samtools genomic coverage reports calculated from DNA-stage alignments only. Because RNA matches are removed first, these report coverage from fallback DNA reads; they do not project RNA alignments back onto the genome.
- `bbsplit_alignments/bacteria/BACTERIA_*.feature_read_positions.csv` — unique per-read positions and explicitly flagged ambiguous overlaps.
- `bbsplit_alignments/bacteria/BACTERIA_*.feature_position_bins.csv` — read counts and within-feature fractions for every feature and normalized bin.
- `bbsplit_alignments/bacteria/BACTERIA_*.metagene_{CDS,non_CDS}_profile.csv` and `.png` — separate protein-coding and noncoding aggregate profiles.
- `bbsplit_alignments/bacteria/BACTERIA_*.metagene_non_CDS_by_type.csv` and `.png` — subtype profiles (for example rRNA, tRNA, ncRNA, tmRNA, and other).
- `bbsplit_alignments/bacteria/BACTERIA_*.metagene_CDS_vs_non_CDS.png` — comparison of mean normalized profiles.
- `bbsplit_alignments/bacteria/BACTERIA_*.feature_position_summary.csv` — examined, unique, ambiguous, outside-feature, contributing-feature, and non-CDS subtype counts.
- `BACTERIA_<project-directory>.csv` — combined summary table, including eligible reads and reads excluded for exceeding 6000 bp after trimming.
- `BACTERIA_<project-directory>_read_breakdown.csv` — per-sample read disposition, including eligible and oversized counts, and separate bacterial/host CDS, rRNA, other RNA, and intergenic/unannotated counts. Each attribute has a raw count, a percentage of total input reads, and a source column. Existing non-rRNA columns remain available as aligned minus rRNA. Leftover reads passed the minimum- and maximum-length filters but aligned to none of the configured references.

### Read-region categories

Each primary mapped read is counted once for its organism across the RNA and DNA stages. Classification uses strand-compatible overlap with aligned reference blocks; skipped CIGAR regions do not create a feature overlap. If multiple annotation classes overlap a read, precedence is **rRNA → CDS → other RNA → intergenic/unannotated**. The four exclusive category counts sum to the organism's aligned count. The existing non-rRNA count includes CDS, other RNA, and intergenic/unannotated reads and therefore should not be added to those counts.

Other RNA includes annotated noncoding RNAs other than rRNA, such as tRNA, tmRNA, ncRNA, and small RNAs. An RNA transcript alignment outside its annotated CDS is also other RNA. A generic genomic `gene` or `mRNA` record alone does not make a genomic interval CDS or noncoding RNA: genomic reads with no matching CDS or RNA feature, including intronic or UTR sequence without an explicit qualifying feature, fall in intergenic/unannotated.

For genomic noncoding RNA transcripts with annotated exon children, only their exons contribute RNA overlaps; their intervening introns remain intergenic/unannotated. Feature counting runs separately for each sample and combines the resulting columns, including zero columns for empty samples and absent feature types.

Percentages use the original reads reported by Cutadapt as the denominator, or the sampled input when sampling is enabled. They are not percentages of organism-aligned reads. Samples with no reads report `0.00%`; disabled host categories report zero. Category reports are required to regenerate the expanded breakdown from an existing output directory; rerun the updated pipeline to produce them.

## Notes

- Input files must contain one Nanopore molecule per FASTQ record and end in `.fastq`, `.fastq.gz`, `.fq`, or `.fq.gz`; use `FASTQ_GLOB` to narrow which files are selected. Every selected file remains an independent sample regardless of numeric or `_R1`/`_R2` filename suffixes.
- Ensure the configured bacterial reference and annotation use compatible genome versions.

## Normalized feature-position analysis

Feature-position analysis includes the combined bacterial RNA and DNA alignments. RNA transcript records retain their own coordinate system, identified by the `RNA__` sequence prefix, and genomic DNA records retain genomic coordinates. These are separate feature observations rather than transcript-to-genome projections.

The post-alignment analysis collapses differently sized annotated features onto a common biological coordinate: **0% is the 5' end and 100% is the 3' end**. It follows genomic start-to-end on `+` features and reverses genomic coordinates on `-` features. Each long-read observation uses the midpoint of its aligned portion. Unmapped, secondary, and supplementary alignments are ignored; duplicates are not filtered.

Protein-coding loci (`CDS`) use the complete parent gene interval when the annotation hierarchy provides it (falling back to the combined CDS extent). Ordinary `gene` records are never treated as noncoding merely because they lack a CDS type. Non-CDS loci come from biologically informative RNA annotations such as rRNA, tRNA, tmRNA, ncRNA, sRNA, misc_RNA, and snRNA, with RNA children preferred over their generic gene parents. Subtype aggregates prevent abundant rRNA from silently dominating every noncoding comparison.

Each feature is normalized independently: its bin fraction is `reads in bin / reads assigned to that feature`. Aggregate CSVs and plots report the mean and median of those per-feature fractions, rather than pooling raw reads. Features below `METAGENE_MIN_FEATURE_READS` remain visible in detailed and binned files but are excluded from aggregates. A midpoint overlapping multiple features in the same class is flagged in the detailed CSV, counted as ambiguous in the summary, and excluded from positional profiles rather than assigned arbitrarily. CDS and non-CDS assignment are evaluated separately, so biologically overlapping classes remain independently measurable.
