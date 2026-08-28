# Whiteley Bacteria Mapping Pipeline

This repository contains an Oxford Nanopore long-read FASTQ processing pipeline (with paired-read compatibility) for bacterial sequencing runs. The pipeline trims adapters, removes reads that map to a decoy/pangenome reference, maps the remaining reads to a configurable bacterial reference genome, counts gene-level alignments, calculates bacterial coverage metrics, and writes a project-level CSV summary.

## What the program does

`map_bacteria_with_decoys.sh` runs the analysis in these stages:

1. Finds input FASTQ files matching the configured glob and, when sampling is enabled, randomly selects up to the configured number of reads from each file.
2. Uses `cutadapt` to remove configured Illumina and Oxford Nanopore adapter sequences and discard reads shorter than the configured minimum length.
3. Uses `minimap2` to map trimmed reads to a decoy/pangenome index and keeps reads that do **not** map to the decoys.
4. Uses `minimap2` with the `map-ont` preset again to map decoy-unmapped reads to the bacterial reference index, retaining the reads that also fail this second alignment. When `HOST_MINIMAP2_REFERENCE` is set, only those reads that mapped to neither the decoy nor the bacterium are mapped to the host index.
5. Uses `featureCounts` from Subread to assign aligned reads to CDS features in the bacterial GFF annotation sharing the reference basename and, when host mapping is enabled, independently counts host alignments against the matching host annotation.
6. Uses `samtools` to create, sort, index, and calculate coverage from BAM files.
7. Uses `gene_position_profile.py` on each sorted bacterial BAM to calculate strand-aware, normalized feature-position and aggregate metagene profiles.
8. Runs `bacteria_with_decoys_csvConversion.py` to combine cutadapt, minimap2, coverage, and featureCounts outputs into `BACTERIA_<project-directory>.csv`.
9. Writes `BACTERIA_<project-directory>_read_breakdown.csv`, with one sample per row and read counts, percentages of total reads, and the source of every value.

The summary conversion scripts use only the Python standard library. Positional profiling uses pysam and matplotlib.

## Dependencies

The pipeline environment pins these command-line tools:

- Python 3.13.5
- cutadapt 5.1
- minimap2 2.30
- samtools 1.22.1
- Subread / featureCounts 2.1.1
- pysam
- matplotlib

A conda environment file is provided in `environment.yml`.

## Create the conda environment

```bash
conda env create -f environment.yml
conda activate whiteley-bacteria-pipeline
```

If you are running on an HPC system that requires module loading, load your site-specific conda or mamba module before creating or activating the environment.

## Configure inputs and references

Edit `config.env` before running the pipeline. The key settings are:

| Setting | Purpose |
| --- | --- |
| `FASTQ_DIR` | Directory containing input FASTQ files. |
| `FASTQ_GLOB` | Shell glob for FASTQ inputs, such as `*.fastq`. |
| `READ_LAYOUT` | `single` (default) or `paired`; paired files may be named `<sample>_R1.fastq[.gz]` and `<sample>_R2.fastq[.gz]`, or include the same suffix after the read marker (for example, `_R1_001` and `_R2_001`). |
| `OUTPUT_DIR` | Directory where output files should be written. |
| `SAMPLE_READS` | Set to `true` to analyze a random subset of each input, or `false` to analyze every read. |
| `SAMPLE_SIZE` | Maximum reads sampled from each FASTQ (default: `5000`). |
| `SAMPLE_SEED` | Seed used to make random sampling reproducible. |
| `DECOY_MINIMAP2_REFERENCE` | Minimap2 reference basename for the decoy/pangenome reference. |
| `BACTERIA_MINIMAP2_REFERENCE` | Shared basename for the bacterial minimap2 index (or FASTA) and GFF annotation. |
| `HOST_MINIMAP2_REFERENCE` | Optional shared basename for the host index (or FASTA) and GFF annotation; leave empty to disable host mapping and counting. |
| `MINIMAP2_PRESET` | Minimap2 preset, defaulting to `map-ont` for Oxford Nanopore reads. |
| `ADAPTER_SINGLE`, `ADAPTER_NANOPORE` | Illumina and Nanopore adapters both passed to cutadapt for single long reads. |
| `ADAPTER_R1`, `ADAPTER_R2` | Mate-specific paired-end adapters (both default to `CTGTCTCTTATACACATCT`). |
| `MIN_READ_LENGTH` | Minimum read length retained by cutadapt. |
| `THREADS` | Number of threads used by every multithreaded step (default: `16`). |
| `GENE_POSITION_BINS` | Number of equal normalized 5'-to-3' bins (default: `100`). |
| `METAGENE_MIN_FEATURE_READS` | Assigned reads/fragments required for a feature to contribute to aggregate profiles (default: `10`). |
| `CSV_CONVERSION_SCRIPT` | Path to `bacteria_with_decoys_csvConversion.py`. |

### Quick sampling mode

Set `SAMPLE_READS="true"` in `config.env` for a quick exploratory run. Before trimming or mapping, the pipeline uses reservoir sampling to select up to `SAMPLE_SIZE` complete read records. For paired input it selects the same record positions from both mates and verifies that their record counts agree. Files with 5,000 reads or fewer are used in full with the default setting. The temporary sampled inputs are written under `OUTPUT_DIR/.bacteria_sampled_fastq`; original FASTQ files are never modified. Sampling is reproducible for the same input paths and `SAMPLE_SEED`. Set `SAMPLE_READS="false"` for a full analysis.

### Paired-end input

Nanopore data should normally use the default `READ_LAYOUT="single"`. Legacy paired input remains supported: set `READ_LAYOUT="paired"`, set `MINIMAP2_PRESET="sr"`, and provide matching `_R1`/`_R2` files. Cutadapt receives both mates, each minimap2 stage receives both FASTQs, and featureCounts counts fragments rather than individual mates.

The minimap2 reference settings are basenames, not individual `.mmi` files. For example, configure `/refs/bacteria_reference`; the pipeline reuses `/refs/bacteria_reference.mmi`, or builds it from `/refs/bacteria_reference.fa`, `.fasta`, or `.fna` (optionally gzip-compressed). Bacterial and host annotations are discovered from the same basename using `.gff*`. Exactly one matching annotation must exist. For feature counting, the pipeline uses the annotation's `locus` attribute, falling back to `locus_tag` and then `gene`.

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

- `*_22bp.trim.fastq` (single) or `*_R{1,2}_22bp.trim.fastq` (paired) — adapter-trimmed FASTQ files.
- `minimap2_alignments/decoy/` — decoy SAM files, minimap2 reports and diagnostic logs, and decoy-unmapped reads.
- `minimap2_alignments/bacteria/` — bacterial SAM/BAM files, minimap2 logs, and coverage reports.
- `minimap2_alignments/bacteria/*_unmapped_to_bacteria.fastq.gz` — reads that mapped to neither the decoy nor bacterial reference and are used as the optional host-alignment input.
- `minimap2_alignments/host/` — optional host SAM files, minimap2 logs, and host featureCounts results.
- `featurecounts_BACTERIA_summary.txt` and `featurecounts_BACTERIA_summary.csv` — featureCounts results.
- `minimap2_alignments/host/featurecounts_HOST_summary.txt` and `.csv` — optional, separate host featureCounts results.
- `minimap2_alignments/bacteria/BACTERIA_*_coverage.txt` — samtools coverage reports.
- `minimap2_alignments/bacteria/BACTERIA_*.feature_read_positions.csv` — unique per-read/per-fragment positions and explicitly flagged ambiguous overlaps.
- `minimap2_alignments/bacteria/BACTERIA_*.feature_position_bins.csv` — read counts and within-feature fractions for every feature and normalized bin.
- `minimap2_alignments/bacteria/BACTERIA_*.metagene_{CDS,non_CDS}_profile.csv` and `.png` — separate protein-coding and noncoding aggregate profiles.
- `minimap2_alignments/bacteria/BACTERIA_*.metagene_non_CDS_by_type.csv` and `.png` — subtype profiles (for example rRNA, tRNA, ncRNA, tmRNA, and other).
- `minimap2_alignments/bacteria/BACTERIA_*.metagene_CDS_vs_non_CDS.png` — comparison of mean normalized profiles.
- `minimap2_alignments/bacteria/BACTERIA_*.feature_position_summary.csv` — examined, unique, ambiguous, outside-feature, contributing-feature, and non-CDS subtype counts.
- `BACTERIA_<project-directory>.csv` — combined summary table.
- `BACTERIA_<project-directory>_read_breakdown.csv` — per-sample read disposition. Each attribute has a raw count, a percentage of total input reads, and a source column. Bacterial and host non-rRNA values are the corresponding aligned counts minus reads assigned to an `rRNA` feature. Leftover reads passed cutadapt but aligned to none of the decoy, bacterial, or host references.

## Notes

- Input files must end in `.fastq` or `.fastq.gz`; use `FASTQ_GLOB` to narrow which files are selected.
- Ensure the configured bacterial reference and annotation use compatible genome versions.

## Normalized feature-position analysis

The post-alignment analysis collapses differently sized annotated features onto a common biological coordinate: **0% is the 5' end and 100% is the 3' end**. It follows genomic start-to-end on `+` features and reverses genomic coordinates on `-` features. Single-end observations use the midpoint of the aligned portion. In paired mode, the complete aligned fragment midpoint is used and the pair is counted once. Unmapped, secondary, and supplementary alignments are ignored; duplicates are not filtered.

Protein-coding loci (`CDS`) use the complete parent gene interval when the annotation hierarchy provides it (falling back to the combined CDS extent). Ordinary `gene` records are never treated as noncoding merely because they lack a CDS type. Non-CDS loci come from biologically informative RNA annotations such as rRNA, tRNA, tmRNA, ncRNA, sRNA, misc_RNA, and snRNA, with RNA children preferred over their generic gene parents. Subtype aggregates prevent abundant rRNA from silently dominating every noncoding comparison.

Each feature is normalized independently: its bin fraction is `reads in bin / reads assigned to that feature`. Aggregate CSVs and plots report the mean and median of those per-feature fractions, rather than pooling raw reads. Features below `METAGENE_MIN_FEATURE_READS` remain visible in detailed and binned files but are excluded from aggregates. A midpoint overlapping multiple features in the same class is flagged in the detailed CSV, counted as ambiguous in the summary, and excluded from positional profiles rather than assigned arbitrarily. CDS and non-CDS assignment are evaluated separately, so biologically overlapping classes remain independently measurable.
