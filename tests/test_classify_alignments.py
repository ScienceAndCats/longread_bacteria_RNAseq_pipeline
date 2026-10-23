"""Exclusive feature assignments across sequential RNA and DNA references."""

import csv
import importlib.util
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest


REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))
HAS_PYSAM = importlib.util.find_spec("pysam") is not None
if HAS_PYSAM:
    import classify_alignments as classifier


def annotation(path, features):
    path.write_text("##gff-version 3\n" + "".join(
        f"ref\ttest\t{kind}\t{start}\t{end}\t.\t{strand}\t.\tgene_id=g{i}\n"
        for i, (kind, start, end, strand) in enumerate(features)))


def sam_file(path, records=(), length=1000):
    lines = [f"@HD\tVN:1.6\tSO:unsorted\n@SQ\tSN:ref\tLN:{length}\n"]
    for name, flag, position, cigar in records:
        if flag & 4:
            lines.append(f"{name}\t{flag}\t*\t0\t0\t*\t*\t0\t0\t*\t*\n")
        else:
            lines.append(f"{name}\t{flag}\tref\t{position}\t60\t{cigar}\t*\t0\t0\t*\t*\n")
    path.write_text("".join(lines))


@unittest.skipUnless(HAS_PYSAM, "pysam is required")
class ClassifierTests(unittest.TestCase):
    def test_priority_and_single_count_despite_overlapping_features_and_duplicate_names(self):
        with tempfile.TemporaryDirectory() as temp:
            path = Path(temp)
            sam, gff = path / "reads.sam", path / "ref.gff"
            annotation(gff, [("CDS", 10, 50, "+"), ("CDS", 20, 55, "+"),
                             ("rRNA", 30, 45, "+"), ("tRNA", 40, 80, "+"),
                             ("gene", 1, 1000, "+")])
            sam_file(sam, [("same_name", 0, 12, "5M"), ("same_name", 0, 22, "5M"),
                           ("rrna_priority", 0, 40, "5M"), ("trna", 0, 65, "5M"),
                           ("intergenic", 0, 90, "5M"), ("secondary", 256, 40, "5M"),
                           ("supplementary", 2048, 40, "5M"), ("unmapped", 4, 0, "*")])
            self.assertEqual(classifier.classify(sam, gff), {
                "aligned": 5, "CDS": 2, "rRNA": 1, "other_RNA": 1, "intergenic_unannotated": 1})

    def test_strands_and_feature_boundaries(self):
        with tempfile.TemporaryDirectory() as temp:
            path = Path(temp)
            sam, gff = path / "reads.sam", path / "ref.gff"
            annotation(gff, [("CDS", 10, 20, "+"), ("rRNA", 10, 20, "-"),
                             ("tRNA", 40, 50, "-"), ("ncRNA", 60, 70, ".")])
            sam_file(sam, [("plus", 0, 10, "11M"), ("minus", 16, 10, "11M"),
                           ("wrong_strand", 0, 40, "5M"), ("right_strand", 16, 40, "5M"),
                           ("before", 0, 5, "5M"), ("after", 0, 21, "5M"),
                           ("strandless", 0, 60, "5M")])
            self.assertEqual(classifier.classify(sam, gff), {
                "aligned": 7, "CDS": 1, "rRNA": 1, "other_RNA": 1, "intergenic_unannotated": 4})

    def test_skipped_and_deleted_reference_bases_do_not_overlap_features(self):
        with tempfile.TemporaryDirectory() as temp:
            path = Path(temp)
            sam, gff = path / "reads.sam", path / "ref.gff"
            annotation(gff, [("CDS", 20, 30, "+"), ("rRNA", 60, 70, "+")])
            sam_file(sam, [("splice", 0, 10, "5M30N5M"),
                           ("deletion", 0, 55, "5M11D5M"),
                           ("cds_hit", 0, 10, "5M5N5M")])
            self.assertEqual(classifier.classify(sam, gff), {
                "aligned": 3, "CDS": 1, "rRNA": 0, "other_RNA": 0, "intergenic_unannotated": 2})

    def test_rna_transcript_fallback_and_dna_mrna_exclusion(self):
        with tempfile.TemporaryDirectory() as temp:
            path = Path(temp)
            dna_sam, dna_gff = path / "dna.sam", path / "dna.gff"
            rna_sam, rna_gff = path / "rna.sam", path / "rna.gff"
            features = [("mRNA", 1, 100, "+"), ("CDS", 20, 80, "+")]
            annotation(dna_gff, features)
            annotation(rna_gff, features)
            sam_file(dna_sam, [("dna_utr", 0, 5, "5M"), ("dna_cds", 0, 30, "5M")])
            sam_file(rna_sam, [("rna_utr", 0, 5, "5M"), ("rna_cds", 0, 30, "5M")])
            self.assertEqual(classifier.classify(dna_sam, dna_gff, rna_sam, rna_gff), {
                "aligned": 4, "CDS": 2, "rRNA": 0, "other_RNA": 1, "intergenic_unannotated": 1})

    def test_genomic_noncoding_exons_exclude_introns_and_preserve_strands(self):
        with tempfile.TemporaryDirectory() as temp:
            path = Path(temp)
            sam, gff = path / "reads.sam", path / "ref.gff"
            # A percent-encoded comma belongs to the ID; it is not a Parent
            # list delimiter. An intermediary exercises ancestor traversal.
            gff.write_text(
                "##gff-version 3\n"
                "ref\ttest\tgene\t10\t100\t.\t+\t.\tID=gene_nc\n"
                "ref\ttest\tlnc_RNA\t10\t100\t.\t+\t.\tID=rna%3Anc%2C1;Parent=gene_nc\n"
                "ref\ttest\ttranscript\t10\t100\t.\t+\t.\tID=bridge;Parent=rna%3Anc%2C1\n"
                "ref\ttest\texon\t10\t20\t.\t+\t.\tParent=bridge\n"
                "ref\ttest\texon\t80\t90\t.\t+\t.\tParent=rna%3Anc%2C1\n"
                "ref\ttest\tncRNA\t200\t220\t.\t+\t.\tID=bacterial_nc\n"
                "ref\ttest\tncRNA\t300\t380\t.\t-\t.\tID=minus_nc\n"
                "ref\ttest\texon\t310\t320\t.\t-\t.\tParent=minus_nc\n"
                "ref\ttest\texon\t360\t370\t.\t-\t.\tParent=minus_nc\n")
            sam_file(sam, [("exon", 0, 10, "5M"), ("exon_antisense", 16, 10, "5M"),
                           ("intron", 0, 30, "5M"), ("junction", 0, 15, "5M60N5M"),
                           ("single_span", 0, 205, "5M"), ("minus_exon", 16, 315, "5M"),
                           ("minus_antisense", 0, 315, "5M"), ("minus_intron", 16, 340, "5M")])
            self.assertEqual(classifier.classify(sam, gff), {
                "aligned": 8, "CDS": 0, "rRNA": 0, "other_RNA": 4, "intergenic_unannotated": 4})

    def test_rna_reference_coordinates_keep_full_transcript_extent(self):
        with tempfile.TemporaryDirectory() as temp:
            path = Path(temp)
            sam, gff = path / "reads.sam", path / "ref.gff"
            gff.write_text(
                "ref\ttest\tncRNA\t1\t100\t.\t+\t.\tID=nc\n"
                "ref\ttest\texon\t1\t20\t.\t+\t.\tParent=nc\n"
                "ref\ttest\texon\t80\t100\t.\t+\t.\tParent=nc\n")
            sam_file(sam, [("rna_middle", 0, 40, "5M")])
            self.assertEqual(classifier.classify_file(sam, gff, rna_reference=True), {
                "CDS": 0, "rRNA": 0, "other_RNA": 1, "intergenic_unannotated": 0})

    def test_empty_samples_and_unannotated_contigs(self):
        with tempfile.TemporaryDirectory() as temp:
            path = Path(temp)
            sam, gff = path / "reads.sam", path / "ref.gff"
            annotation(gff, [])
            sam_file(sam)
            self.assertEqual(sum(classifier.classify(sam, gff).values()), 0)
            sam_file(sam, [("no_annotations", 0, 500, "5M")])
            self.assertEqual(classifier.classify(sam, gff)["intergenic_unannotated"], 1)

    def test_invalid_annotation_references_and_coordinates_are_rejected(self):
        with tempfile.TemporaryDirectory() as temp:
            path = Path(temp)
            sam, gff = path / "reads.sam", path / "ref.gff"
            sam_file(sam)
            for record in ("other\ttest\tCDS\t1\t10\t.\t+\t.\tgene_id=g\n",
                           "ref\ttest\tCDS\t0\t10\t.\t+\t.\tgene_id=g\n",
                           "ref\ttest\tCDS\t1\t1001\t.\t+\t.\tgene_id=g\n"):
                with self.subTest(record=record):
                    gff.write_text(record)
                    with self.assertRaises(ValueError):
                        classifier.classify(sam, gff)

    def test_cli_writes_conserved_metric_report(self):
        with tempfile.TemporaryDirectory() as temp:
            path = Path(temp)
            sam, gff, report = path / "reads.sam", path / "ref.gff", path / "counts.tsv"
            sam_file(sam, [("cds", 0, 10, "5M")])
            annotation(gff, [("CDS", 1, 20, "+")])
            command = [sys.executable, str(REPO / "classify_alignments.py"),
                       "--dna-sam", str(sam), "--dna-gff", str(gff), "--output", str(report)]
            subprocess.run(command, check=True, capture_output=True, text=True)
            with report.open(newline="") as source:
                self.assertEqual(list(csv.reader(source, delimiter="\t")), [
                    ["metric", "count"], ["aligned", "1"], ["CDS", "1"], ["rRNA", "0"],
                    ["other_RNA", "0"], ["intergenic_unannotated", "0"]])
            invalid = subprocess.run(command + ["--rna-sam", str(sam)], capture_output=True, text=True)
            self.assertNotEqual(invalid.returncode, 0)
            self.assertIn("must be supplied together", invalid.stderr)


if __name__ == "__main__":
    unittest.main()
