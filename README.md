# Whiteley Bacteria Mapping Pipeline

This repository contains a long-read FASTQ processing pipeline for bacterial sequencing runs, with Oxford Nanopore as the default platform. The pipeline trims a configurable long-read adapter, removes reads that map to a decoy/pangenome reference, maps the remaining reads to a configurable bacterial reference genome, counts gene-level alignments, calculates bacterial coverage metrics, and writes a project-level CSV summary.

## What the program does

`map_bacteria_with_decoys.sh` runs the analysis in these stages:

1. Finds input FASTQ files matching the configured glob and, when sampling is enabled, randomly selects up to the configured number of reads from each file.
2. Uses Cutadapt to trim adapters in the configured `legacy` or `neb_e7330` mode, then applies the minimum read length. NEB mode extracts both flanks in either orientation and cleans repeated terminal adapters. Reads lacking adapters or having only one confirmed flank are retained if they pass the final length filter.
3. Uses `minimap2` to map trimmed reads to a decoy/pangenome index and keeps reads that do **not** map to the decoys.
4. Uses `minimap2` with the `map-ont` preset again to map decoy-unmapped reads to the bacterial reference index, retaining the reads that also fail this second alignment. When `HOST_MINIMAP2_REFERENCE` is set, only those reads that mapped to neither the decoy nor the bacterium are mapped to the host index. It saves a final leftover FASTQ containing reads that failed every configured alignment stage, whether or not host mapping is enabled.
5. Uses `featureCounts` from Subread to assign aligned reads to CDS features in the bacterial GFF annotation sharing the reference basename and, when host mapping is enabled, independently counts host alignments against the matching host annotation.
6. Uses `samtools` to create, sort, index, and calculate coverage from BAM files.
7. Uses `gene_position_profile.py` on each sorted bacterial BAM to calculate strand-aware, normalized feature-position and aggregate metagene profiles.
8. Runs `bacteria_with_decoys_csvConversion.py` to combine cutadapt, minimap2, coverage, and featureCounts outputs into `BACTERIA_<project-directory>.csv`.
9. Writes `BACTERIA_<project-directory>_read_breakdown.csv`, with one sample per row and read counts, percentages of total reads, and the source of every value.

The summary conversion scripts use only the Python standard library. Positional profiling uses pysam and matplotlib.

## Dependencies

The pipeline environment pins these command-line tools:

- Python 3.13.15
- cutadapt 5.2 (NEB mode uses its rightmost 3' adapter matching)
- minimap2 2.30
- samtools 1.24
- Subread / featureCounts 2.1.1
- pysam
- matplotlib

A conda environment file is provided in `environment.yml`.

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
| `FASTQ_DIR` | Directory containing input FASTQ files. |
| `FASTQ_GLOB` | Shell glob for FASTQ inputs, such as `*.fastq`. |
| `OUTPUT_DIR` | Directory where output files should be written. |
| `SAMPLE_READS` | Set to `true` to analyze a random subset of each input, or `false` to analyze every read. |
| `SAMPLE_SIZE` | Maximum reads sampled from each FASTQ (default: `5000`). |
| `SAMPLE_SEED` | Seed used to make random sampling reproducible. |
| `DECOY_MINIMAP2_REFERENCE` | Minimap2 reference basename for the decoy/pangenome reference. |
| `BACTERIA_MINIMAP2_REFERENCE` | Shared basename for the bacterial minimap2 index (or FASTA) and GFF annotation. |
| `HOST_MINIMAP2_REFERENCE` | Optional shared basename for the host index (or FASTA) and GFF annotation; leave empty to disable host mapping and counting. |
| `MINIMAP2_PRESET` | Long-read minimap2 preset: `map-ont` (default), `map-hifi`, or `map-pb`. Other presets are rejected. |
| `ADAPTER_TRIMMING_MODE` | `legacy` (default) or `neb_e7330` for NEBNext E7330S Illumina libraries sequenced with Nanopore. |
| `CUTADAPT_ERROR_RATE` | Maximum adapter alignment error rate (default: `0.15`); mismatches and indels are allowed. |
| `CUTADAPT_MIN_OVERLAP` | Minimum adapter overlap, in nucleotides (default: `12`). |
| `ADAPTER_MAX_TRIMMING_ROUNDS` | NEB mode round limit (default: `10`), including initial flank extraction. |
| `ADAPTER_NEB_5P`, `ADAPTER_NEB_3P` | NEB E7330S 5' and 3' flanks in the insert's forward orientation. |
| `TRIM_NANOPORE_ADAPTERS` | Enable Nanopore adapter removal (default: `true`); NEB mode limits it to terminal matches on incompletely bounded reads. |
| `ADAPTER_NANOPORE` | Long-read adapter passed to cutadapt; defaults to the Oxford Nanopore ligation adapter and may be changed for another library preparation. |
| `ADAPTER_NANOPORE_RC` | Its reverse complement; leave empty to compute it in NEB mode when customizing the adapter. |
| `TRIM_TRUSEQ_SMALL_RNA` | NEB mode opt-in for terminal TruSeq small-RNA adapter removal (default: `false`). |
| `ADAPTER_TRUSEQ_SMALL_RNA` | Optional other-kit adapter `TGGAATTCTCGGGTGCCAAGG`. |
| `TRIM_ILLUMINA_ADAPTERS` | Legacy mode only: optional second Cutadapt pass for standard Illumina and small-RNA kit adapters; unmatched reads are retained. |
| `ADAPTER_ILLUMINA` | Standard Illumina 3' adapter used by the optional second trimming pass. |
| `ADAPTER_ILLUMINA_SMALL_RNA` | Illumina small-RNA kit 3' adapter used by the optional second trimming pass. |
| `MIN_READ_LENGTH` | Minimum retained length **after all trimming** (default: `20`). Existing project configs with explicit values such as `22` keep that value. |
| `THREADS` | Number of threads used by every multithreaded step (default: `16`). |
| `GENE_POSITION_BINS` | Number of equal normalized 5'-to-3' bins (default: `100`). |
| `METAGENE_MIN_FEATURE_READS` | Assigned reads required for a feature to contribute to aggregate profiles (default: `10`). |
| `CSV_CONVERSION_SCRIPT` | Path to `bacteria_with_decoys_csvConversion.py`. |

### Quick sampling mode

Set `SAMPLE_READS="true"` in `config.env` for a quick exploratory run. Before trimming or mapping, the pipeline uses reservoir sampling to select up to `SAMPLE_SIZE` complete long-read records from each input. Files with 5,000 reads or fewer are used in full with the default setting. The temporary sampled inputs are written under `OUTPUT_DIR/.bacteria_sampled_fastq`; original FASTQ files are never modified. Sampling is reproducible for the same input paths and `SAMPLE_SEED`. Set `SAMPLE_READS="false"` for a full analysis.

The minimap2 reference settings are basenames, not individual `.mmi` files. For example, configure `/refs/bacteria_reference`; the pipeline reuses `/refs/bacteria_reference.mmi`, or builds it from `/refs/bacteria_reference.fa`, `.fasta`, or `.fna` (optionally gzip-compressed). Bacterial and host annotations are discovered from the same basename using `.gff*`. Exactly one matching annotation must exist. For feature counting, the pipeline uses the annotation's `locus` attribute, falling back to `locus_tag` and then `gene`.

In `legacy` mode, set `TRIM_ILLUMINA_ADAPTERS="true"` for libraries that may also contain Illumina-derived adapters. Cutadapt first removes the configured Nanopore 3' adapter, then searches for the standard Illumina and Illumina small-RNA adapters in up to two rounds. The minimum length is applied after both passes. Adapter choices and unmatched-read retention remain compatible with existing configurations; the new error/overlap defaults are `0.15`/`12` (set `0.1`/`3` to recover the old matching tolerance). The Nanopore and optional Illumina reports are `*_nanopore_adapter_log.txt` and `*_illumina_adapter_log.txt`; `*_cutadapt_log.txt` and `*.cutadapt.json` report the final length-filter pass.

### NEBNext E7330S libraries sequenced using Nanopore

For NEBNext Small RNA Library Prep Set for Illumina (E7330S) libraries, including Plasmidsaurus Premium PCR sequencing, configure:

```bash
ADAPTER_TRIMMING_MODE="neb_e7330"
ADAPTER_NEB_5P="GTTCAGAGTTCTACAGTCCGACGATC"
ADAPTER_NEB_3P="AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"
CUTADAPT_ERROR_RATE="0.15"
CUTADAPT_MIN_OVERLAP="12"
ADAPTER_MAX_TRIMMING_ROUNDS="10"
MIN_READ_LENGTH="20"  # use 22 if required by your analysis
TRIM_NANOPORE_ADAPTERS="true"
TRIM_TRUSEQ_SMALL_RNA="false"
```

`trim_neb_adapters.py` streams FASTQ through Cutadapt's existing adapter parser,
`AdapterCutter`, and `ReverseComplementer`, which implement `-g`/`-a`, `--times`,
and `--revcomp`. It needs no dependencies beyond those installed with Cutadapt.
Matching is performed in one process; `THREADS` still controls the legacy
Cutadapt commands and downstream multithreaded tools.

1. Find a linked `-g '5P...3P;rightmost'` match requiring both flanks, testing
   the read and its reverse complement with Cutadapt's `--revcomp` engine.
   Keep the span between the leftmost best 5' match and rightmost best 3' match.
   This single linked extraction removes both outer flanks, including any
   preceding/following PCR primer, index, P5/P7 or Nanopore sequence. Read
   orientation is normalized to the supplied forward NEB adapters.
2. Clean exposed ends using non-internal `-g X5P` and `-a 3PX` matching with
   `--times 2` per cleanup round. The literal `X` disallows internal matches
   while allowing partial adapters at the read ends. Consecutive 5' and/or 3'
   copies are removed until no terminal matches remain or the total round
   budget is reached. This avoids a broad repeated search through the insert.
3. Without a linked match, retain reads and try terminal NEB matches, using
   `--revcomp` until a NEB flank establishes orientation. Optional terminal
   Nanopore cleanup recognizes both `ADAPTER_NANOPORE` and its reverse
   complement at either end; these matches never determine insert orientation.
   Complete NEB constructs receive no additional Nanopore pass, protecting
   Nanopore-like motifs within the bounded insert.
4. Apply `MIN_READ_LENGTH` only after these operations. IDs and header comments
   are unchanged; qualities are sliced with their bases and reversed when a
   read is reverse-complemented. Plain and gzipped FASTQ input are supported.

A single confirmed flank is labeled `5p_only` or `3p_only`, **not** a fully
bounded biological insert. If the only NEB motif is internal and its outer
sequence cannot be bounded, the read is left intact and its remaining motif
is reported. Terminal cleanup does not guess which unknown bases are adapters.
TruSeq small-RNA sequence is excluded by default; enable
`TRIM_TRUSEQ_SMALL_RNA` explicitly for terminal matching from that kit.
`TRIM_ILLUMINA_ADAPTERS` applies only to legacy mode.

An internal 3'-adapter followed by a 5'-adapter is flagged
`ambiguous_concatemer=1`. The bounded span is retained unsplit, with its internal
junction, for review; it still undergoes the final length filter. These reads
and partial/unbounded reads remain in the FASTQ supplied to the unchanged
mapping pipeline. Neither flag certifies a clean insert. A biological sequence
identical to a terminal adapter (or to a complete library layout) cannot be
distinguished from that adapter using sequence matching alone.

NEB mode writes these files per sample, using the same `<sample>_<length>bp` prefix:

| File | Contents |
| --- | --- |
| `*.cutadapt_log.txt` | Cutadapt version/settings and whole-operation counts, including final length filtering; compatible with existing summary readers. |
| `*.adapter_stats.csv` | One row: total/detected/trimmed reads, both/single flanks, reverse complements, repeated adapters, multiple rounds, length discards, remaining motifs, round-limit hits, ambiguous concatemers, retained count and retained length minimum/maximum/mean. |
| `*.adapter_reads.csv` | Per-read evidence and retention flags, including original header, lengths, flank match counts and ambiguity flags; read numbers disambiguate duplicate IDs. |
| `*.read_lengths.csv` | Exact retained-read length distribution (`length,retained_reads`). |

The NEB reports are generated from observed per-read Cutadapt match objects
and final lengths, rather than estimates from marginal per-adapter totals.
Both-flank counts require accepted matches to each NEB flank. Repeated-adapter
counts require multiple accepted matches to the same adapter type; the initial
linked extraction counts as one round removing two flanks, and each cleanup
round removes up to two adapters. Counts include reads subsequently discarded
by length, except fields explicitly labeled retained. `reads_with_adapters_detected`
also includes diagnostic interior motif matches; `trimmed_reads` counts actual
length changes. Remaining-motif counts include internal biological lookalikes
and do not by themselves prove adapter contamination. A separate terminal
remaining count identifies incomplete cleanup at the round limit.

To trim a FASTQ without any reference mapping, run from the repository:

```bash
bash trim_fastq.sh config.env reads.fastq.gz trimmed.fastq sample_20bp
```

Run the synthetic regression suite in the pipeline environment:

```bash
python3 -m unittest discover -s tests -v
```

It covers 22-nt inserts in both orientations; two/three repeated copies on
either/both sides; P5/P7, index and Nanopore outer flanks; substitutions,
insertions and deletions; terminal partial adapters; dimers; reads without
adapters; internal motifs; ambiguous concatemers; round limits; short/long
inserts; exact qualities/headers; gzip/plain FASTQ; and legacy final filtering.

To enable host mapping and feature counting, set `HOST_MINIMAP2_REFERENCE` to the shared host reference basename and provide its FASTA/index and matching `.gff*` annotation. Setting it to `""` skips host index preparation, alignment, and counting. Host input consists exclusively of reads that did not align to either the decoy or bacterial reference.

## Run the pipeline

From this repository, run:

```bash
bash map_bacteria_with_decoys.sh config.env
```

You can keep multiple config files and pass the one for the current run:

```bash
bash map_bacteria_with_decoys.sh configs/project_a.env
```

## Important outputs

- `*_<MIN_READ_LENGTH>bp.trim.fastq` — adapter-trimmed FASTQ files after the final length filter.
- `minimap2_alignments/decoy/` — decoy SAM files, minimap2 reports and diagnostic logs, and decoy-unmapped reads.
- `minimap2_alignments/bacteria/` — bacterial SAM/BAM files, minimap2 logs, and coverage reports.
- `minimap2_alignments/bacteria/*_unmapped_to_bacteria.fastq.gz` — reads that mapped to neither the decoy nor bacterial reference and are used as the optional host-alignment input.
- `minimap2_alignments/host/` — optional host SAM files, minimap2 logs, and host featureCounts results.
- `minimap2_alignments/leftover_reads/*_leftover.fastq.gz` — final reads that mapped to none of the configured decoy, bacterial, or host references. These files are always created; when host mapping is disabled, they contain the bacterial-unmapped reads.
- `featurecounts_BACTERIA_summary.txt` and `featurecounts_BACTERIA_summary.csv` — featureCounts results.
- `minimap2_alignments/host/featurecounts_HOST_summary.txt` and `.csv` — optional, separate host featureCounts results.
- `minimap2_alignments/bacteria/BACTERIA_*_coverage.txt` — samtools coverage reports.
- `minimap2_alignments/bacteria/BACTERIA_*.feature_read_positions.csv` — unique per-read positions and explicitly flagged ambiguous overlaps.
- `minimap2_alignments/bacteria/BACTERIA_*.feature_position_bins.csv` — read counts and within-feature fractions for every feature and normalized bin.
- `minimap2_alignments/bacteria/BACTERIA_*.metagene_{CDS,non_CDS}_profile.csv` and `.png` — separate protein-coding and noncoding aggregate profiles.
- `minimap2_alignments/bacteria/BACTERIA_*.metagene_non_CDS_by_type.csv` and `.png` — subtype profiles (for example rRNA, tRNA, ncRNA, tmRNA, and other).
- `minimap2_alignments/bacteria/BACTERIA_*.metagene_CDS_vs_non_CDS.png` — comparison of mean normalized profiles.
- `minimap2_alignments/bacteria/BACTERIA_*.feature_position_summary.csv` — examined, unique, ambiguous, outside-feature, contributing-feature, and non-CDS subtype counts.
- `BACTERIA_<project-directory>.csv` — combined summary table.
- `BACTERIA_<project-directory>_read_breakdown.csv` — per-sample read disposition. Each attribute has a raw count, a percentage of total input reads, and a source column. Bacterial and host non-rRNA values are the corresponding aligned counts minus reads assigned to an `rRNA` feature. Leftover reads passed cutadapt but aligned to none of the decoy, bacterial, or host references.

## Notes

- Input files must contain one long molecule per FASTQ record and end in `.fastq` or `.fastq.gz`; use `FASTQ_GLOB` to narrow which files are selected.
- Ensure the configured bacterial reference and annotation use compatible genome versions.

## Normalized feature-position analysis

The post-alignment analysis collapses differently sized annotated features onto a common biological coordinate: **0% is the 5' end and 100% is the 3' end**. It follows genomic start-to-end on `+` features and reverses genomic coordinates on `-` features. Each long-read observation uses the midpoint of its aligned portion. Unmapped, secondary, and supplementary alignments are ignored; duplicates are not filtered.

Protein-coding loci (`CDS`) use the complete parent gene interval when the annotation hierarchy provides it (falling back to the combined CDS extent). Ordinary `gene` records are never treated as noncoding merely because they lack a CDS type. Non-CDS loci come from biologically informative RNA annotations such as rRNA, tRNA, tmRNA, ncRNA, sRNA, misc_RNA, and snRNA, with RNA children preferred over their generic gene parents. Subtype aggregates prevent abundant rRNA from silently dominating every noncoding comparison.

Each feature is normalized independently: its bin fraction is `reads in bin / reads assigned to that feature`. Aggregate CSVs and plots report the mean and median of those per-feature fractions, rather than pooling raw reads. Features below `METAGENE_MIN_FEATURE_READS` remain visible in detailed and binned files but are excluded from aggregates. A midpoint overlapping multiple features in the same class is flagged in the detailed CSV, counted as ambiguous in the summary, and excluded from positional profiles rather than assigned arbitrarily. CDS and non-CDS assignment are evaluated separately, so biologically overlapping classes remain independently measurable.
