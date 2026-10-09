"""Optional end-to-end tests using the installed pipeline tools, without mocks.

Run with the pipeline environment active and bbsplit.sh on PATH. The tests
skip when dependencies are unavailable and never download or install tools.
"""

import csv
import gzip
import importlib.util
import os
import random
import shlex
import shutil
import subprocess
import tempfile
import unittest
from pathlib import Path


REPO = Path(__file__).resolve().parents[1]
ADAPTER = "TTTCTGTTGGTGCTGATATTGC"
ELIGIBLE = "Reads eligible for alignment (<=6000 bp)"
OVERSIZED = "Reads over 6000 bp (not processed)"
MISSING = [command for command in ("bbsplit.sh", "cutadapt", "samtools", "featureCounts")
           if shutil.which(command) is None]
MISSING.extend(module for module in ("pysam", "matplotlib")
               if importlib.util.find_spec(module) is None)


def random_sequence(rng, length):
    return "".join(rng.choice("ACGT") for _ in range(length))


def write_fastq(path, records):
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "wt") as handle:
        for name, sequence in records:
            handle.write(f"@{name}\n{sequence}\n+\n{'I' * len(sequence)}\n")


def read_fastq(path):
    opener = gzip.open if path.suffix == ".gz" else open
    records = []
    with opener(path, "rt") as handle:
        while True:
            header = handle.readline()
            if not header:
                break
            sequence, plus, quality = [handle.readline().rstrip("\r\n") for _ in range(3)]
            if not header.startswith("@") or not plus.startswith("+") or len(sequence) != len(quality):
                raise AssertionError(f"Malformed pipeline FASTQ: {path}")
            records.append((header[1:].strip(), sequence, quality))
    return records


def read_csv(path):
    with path.open(newline="") as handle:
        return [{key.strip(): value for key, value in row.items()} for row in csv.DictReader(handle)]


def read_metrics(path):
    with path.open(newline="") as handle:
        return {row["metric"]: int(row["count"]) for row in csv.DictReader(handle, delimiter="\t")}


@unittest.skipIf(bool(MISSING), "Requires real pipeline dependencies: " + ", ".join(MISSING))
class BBSplitPipelineTests(unittest.TestCase):
    def make_fixture(self, directory):
        rng = random.Random(284731)
        bacteria = random_sequence(rng, 15000)
        decoy = random_sequence(rng, 1200)
        host = random_sequence(rng, 1200)
        unrelated = random_sequence(rng, 100)
        reference_dir = directory / "reference genomes"
        reference_dir.mkdir()
        for name, sequence in (("bacteria", bacteria), ("decoy", decoy), ("host", host)):
            # A description after the scaffold ID catches SAM/GFF name mismatches.
            fasta = reference_dir / f"{name}.fasta"
            fasta.write_text(f">{name}_ref synthetic reference with description\n" +
                             "\n".join(sequence[start:start + 80] for start in range(0, len(sequence), 80)) + "\n")
        (reference_dir / "bacteria.gff3").write_text(
            "##gff-version 3\n"
            "bacteria_ref\ttest\tCDS\t501\t10500\t.\t+\t0\tID=bacteria_cds;locus=bacteria_cds\n"
            "bacteria_ref\ttest\trRNA\t11501\t11700\t.\t+\t.\tID=bacteria_rrna;locus=bacteria_rrna\n"
        )
        (reference_dir / "host.gff3").write_text(
            "##gff-version 3\n"
            "host_ref\ttest\tCDS\t301\t700\t.\t+\t0\tID=host_cds;locus=host_cds\n"
            "host_ref\ttest\trRNA\t901\t1100\t.\t+\t.\tID=host_rrna;locus=host_rrna\n"
        )
        expected = {
            "bacteria_5999": bacteria[1000:6999],
            "bacteria_6000_after_adapter": bacteria[4000:10000],
            "bacteria_rrna": bacteria[11520:11620],
            "decoy": decoy[300:400],
            "host": host[400:500],
            "unrelated": unrelated,
        }
        excluded = bacteria[2000:8001]
        input_dir = directory / "input reads"
        input_dir.mkdir()
        records = [(name, sequence + (ADAPTER if name == "bacteria_6000_after_adapter" else ""))
                   for name, sequence in expected.items()]
        records.extend((("over_6000_after_adapter", excluded + ADAPTER), ("short_19", unrelated[:19])))
        write_fastq(input_dir / "mixed.fastq.gz", records)
        write_fastq(input_dir / "all_excluded.fastq", [("excluded_only", excluded + ADAPTER)])
        write_fastq(input_dir / "empty.fastq.gz", [])
        return reference_dir, input_dir, expected, excluded

    def pipeline(self, directory, *, host_enabled, k=None):
        reference_dir, input_dir, expected, excluded = self.make_fixture(directory)
        output_dir = directory / "pipeline output"
        config = directory / "config.env"
        settings = {
            "FASTQ_DIR": str(input_dir), "FASTQ_GLOB": "*.fastq*", "OUTPUT_DIR": str(output_dir),
            "DECOY_BBSPLIT_REFERENCE": str(reference_dir / "decoy"),
            "BACTERIA_BBSPLIT_REFERENCE": str(reference_dir / "bacteria"),
            "HOST_BBSPLIT_REFERENCE": str(reference_dir / "host") if host_enabled else "",
            "THREADS": "1", "MIN_READ_LENGTH": "20", "ADAPTER_TRIMMING_MODE": "legacy",
            "GENE_POSITION_BINS": "10", "METAGENE_MIN_FEATURE_READS": "1",
        }
        if k is not None:
            settings["BBSPLIT_K"] = str(k)
        config.write_text("".join(f"{key}={shlex.quote(value)}\n" for key, value in settings.items()))
        # Keep Java's memory budget bounded on CI hosts with large reported RAM.
        # This wrapper runs the real mapper and changes only its heap allocation.
        wrapper_dir = directory / "tool wrappers"
        wrapper_dir.mkdir()
        mapper = shutil.which("bbsplit.sh")
        wrapper = wrapper_dir / "bbsplit.sh"
        wrapper.write_text(f"#!/usr/bin/env bash\nexec {shlex.quote(mapper)} -Xmx1g \"$@\"\n")
        wrapper.chmod(0o755)
        environment = os.environ.copy()
        environment["PATH"] = str(wrapper_dir) + os.pathsep + environment["PATH"]
        result = subprocess.run(["bash", str(REPO / "map_bacteria_with_decoys.sh"), str(config)],
                                cwd=directory, env=environment, text=True, capture_output=True, timeout=240)
        self.assertEqual(result.returncode, 0,
                         f"Pipeline stdout:\n{result.stdout[-5000:]}\nPipeline stderr:\n{result.stderr[-14000:]}")
        return output_dir, expected, excluded

    def verify_pipeline(self, output, expected, excluded, *, host_enabled, k):
        alignments = output / "bbsplit_alignments"
        trimmed = read_fastq(output / "mixed_20bp.trim.fastq")
        self.assertEqual([(name, sequence) for name, sequence, _ in trimmed], list(expected.items()))
        self.assertEqual([len(sequence) for _, sequence, _ in trimmed], [5999, 6000, 100, 100, 100, 100])
        self.assertTrue(all(quality == "I" * len(sequence) for _, sequence, quality in trimmed))
        self.assertEqual(read_fastq(alignments / "oversized_reads/mixed_20bp_over_6000bp.fastq.gz"),
                         [("over_6000_after_adapter", excluded, "I" * 6001)])
        self.assertEqual(read_fastq(alignments / "oversized_reads/all_excluded_20bp_over_6000bp.fastq.gz"),
                         [("excluded_only", excluded, "I" * 6001)])
        self.assertEqual(read_metrics(output / "mixed_20bp.length_filter.tsv"),
                         {"input": 7, "eligible": 6, "oversized": 1, "max_read_length": 6000})
        for sample in ("all_excluded", "empty"):
            self.assertEqual(read_fastq(output / f"{sample}_20bp.trim.fastq"), [])
            self.assertEqual(read_fastq(alignments / f"leftover_reads/{sample}_20bp_leftover.fastq.gz"), [])
        leftover = read_fastq(alignments / "leftover_reads/mixed_20bp_leftover.fastq.gz")
        leftovers = {name: sequence for name, sequence, _ in leftover}
        self.assertEqual(leftovers, {name: expected[name] for name in
                                    (("unrelated",) if host_enabled else ("host", "unrelated"))})
        for stage, report_name, counts in (
                ("decoy", "mixed_20bp.mapped_to_other_bugs.bbsplit.txt", (6, 1, 5)),
                ("bacteria", "BACTERIA_mixed_20bp.bbsplit.txt", (5, 3, 2))):
            self.assertEqual(read_metrics(alignments / stage / report_name),
                             dict(zip(("input", "aligned", "unmapped"), counts)))
        if host_enabled:
            self.assertEqual(read_metrics(alignments / "host/HOST_mixed_20bp.bbsplit.txt"),
                             {"input": 2, "aligned": 1, "unmapped": 1})
        else:
            self.assertFalse((alignments / "host").exists())
        for stage in (("decoy", "bacteria", "host") if host_enabled else ("decoy", "bacteria")):
            stage_index = alignments / f"indexes/{stage}/k{k}"
            self.assertTrue(stage_index.is_dir())
            self.assertTrue(list(stage_index.rglob(f"*_index_k{k}_*.block")), f"Missing actual {stage} index")
        self.assertEqual({path.name for path in alignments.glob("indexes/*/k*")}, {f"k{k}"})

        main_rows = {row["Sample"]: row for row in read_csv(output / f"BACTERIA_{output.name}.csv")}
        breakdown_rows = {row["Sample"]: row for row in read_csv(output / f"BACTERIA_{output.name}_read_breakdown.csv")}
        self.assertEqual(set(main_rows), {"mixed", "all_excluded", "empty"})
        self.assertEqual(set(breakdown_rows), set(main_rows))
        for sample, raw, passing, eligible, oversized in (("mixed", 8, 7, 6, 1),
                                                          ("all_excluded", 1, 1, 0, 1),
                                                          ("empty", 0, 0, 0, 0)):
            main = main_rows[sample]
            breakdown = breakdown_rows[sample]
            for column, value in (("Raw Reads", raw), ("Reads Written", passing),
                                  (ELIGIBLE, eligible), (OVERSIZED, oversized)):
                self.assertEqual(int(main[column]), value, (sample, column))
            for column, value in (("Total reads", raw), ("Reads passing cutadapt length filter", passing),
                                  (ELIGIBLE, eligible), (OVERSIZED, oversized)):
                self.assertEqual(int(breakdown[f"{column} (reads)"]), value, (sample, column))
        mixed = breakdown_rows["mixed"]
        for column, value in (("Decoy aligned", 1), ("Bacteria aligned", 3), ("Bacteria non-rRNA", 2),
                              ("Bacteria rRNA", 1), ("Host aligned", int(host_enabled)),
                              ("Leftover reads", 1 if host_enabled else 2)):
            self.assertEqual(int(mixed[f"{column} (reads)"]), value, column)
        self.assertEqual(mixed[f"{OVERSIZED} (% of total reads)"], "12.50%")
        self.assertEqual(main_rows["mixed"]["Overall Allignment"], "3")
        self.assertEqual(main_rows["mixed"]["Reads Assigned to features"], "2")
        self.assertEqual(main_rows["mixed"]["Features with non-zero reads"], "1")
        self.assertAlmostEqual(float(main_rows["mixed"]["Percent Covered"].rstrip("%")),
                               9100 / 15000 * 100, delta=.2)
        self.assertGreater(float(main_rows["mixed"]["Mean Coverage"].rstrip("x")), 0)
        self.assertEqual(main_rows["all_excluded"]["Overall Allignment"], "0")
        self.assertEqual(main_rows["empty"]["Overall Allignment"], "0")

        prefix = alignments / "bacteria/BACTERIA_mixed_20bp"
        positions = read_csv(Path(str(prefix) + ".feature_read_positions.csv"))
        self.assertEqual({row["read_id"] for row in positions},
                         {"bacteria_5999", "bacteria_6000_after_adapter", "bacteria_rrna"})
        self.assertEqual({row["reference"] for row in positions}, {"bacteria_ref"})
        self.assertEqual({row["feature_id"] for row in positions}, {"bacteria_cds", "bacteria_rrna"})
        self.assertTrue(all(row["assignment_status"] == "unique" for row in positions))
        summary = {row["metric"]: int(row["value"]) for row in
                   read_csv(Path(str(prefix) + ".feature_position_summary.csv"))}
        self.assertEqual(summary["total_mapped_reads_examined"], 3)
        self.assertEqual(summary["CDS.uniquely_assigned"], 2)
        self.assertEqual(summary["non_CDS.rRNA.uniquely_assigned"], 1)
        for suffix in (".sorted.bam", ".sorted.bam.bai", ".feature_position_bins.csv",
                       ".metagene_CDS_profile.csv", ".metagene_non_CDS_profile.csv",
                       ".metagene_CDS_profile.png", ".metagene_non_CDS_profile.png"):
            path = Path(str(prefix) + suffix)
            self.assertGreater(path.stat().st_size, 0, path)

    def test_default_k_with_host(self):
        with tempfile.TemporaryDirectory() as temp:
            output, expected, excluded = self.pipeline(Path(temp), host_enabled=True)
            self.verify_pipeline(output, expected, excluded, host_enabled=True, k=11)

    def test_configured_k_without_host(self):
        with tempfile.TemporaryDirectory() as temp:
            output, expected, excluded = self.pipeline(Path(temp), host_enabled=False, k=9)
            self.verify_pipeline(output, expected, excluded, host_enabled=False, k=9)


if __name__ == "__main__":
    unittest.main()
