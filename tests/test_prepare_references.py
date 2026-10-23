"""Folder discovery and transcript/genomic coordinate regression fixtures."""

import gzip
import importlib.util
import json
import tempfile
import unittest
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from prepare_references import prepare_references


def write(path, text):
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.name.endswith(".gz"):
        with gzip.open(path, "wt") as handle:
            handle.write(text)
    else:
        path.write_text(text)


def gff(ref, feature_type="CDS", start=1, end=10, attrs="ID=gene1;locus_tag=gene1"):
    return f"##gff-version 3\n{ref}\ttest\t{feature_type}\t{start}\t{end}\t.\t+\t.\t{attrs}\n"


class ReferenceDiscoveryTests(unittest.TestCase):
    def fixture(self, directory):
        decoy, bacteria = directory / "decoy", directory / "PAO1"
        write(decoy / "pangenome.fasta", ">decoy description\nACGTACGTACGT\n")
        write(decoy / "pangenome.1.bt2l", "index, not FASTA")
        write(decoy / "pangenome.mmi", "another index")
        write(bacteria / "GCF_genomic.fna.gz", ">NC_001.1 genomic sequence\n" + "ACGT" * 25 + "\n")
        write(bacteria / "GCF_genomic.gff.gz", gff("NC_001.1", end=60))
        return decoy, bacteria

    def test_flat_bacterial_extracted_transcripts_and_multiple_decoy_files(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            decoy, bacteria = self.fixture(root)
            write(decoy / "extra.fna.gz", ">decoy2\nGGGTTT\n")
            write(bacteria / "GCF_cds_from_genomic.fna.gz",
                  ">lcl|NC_001.1_cds_1 [locus_tag=gene1] [protein_id=NP_001] [location=complement(1..60)]\n" + "ACGT" * 15 + "\n")
            write(bacteria / "GCF_rna_from_genomic.fna.gz",
                  ">lcl|NC_001.1_rna_1 [locus_tag=rrs] [gbkey=rRNA] [product=16S ribosomal RNA]\nACGTACGT\n"
                  ">lcl|NC_001.1_rna_2 [locus_tag=trna] [gbkey=tRNA]\nCCGGTT\n")
            manifest = prepare_references(str(decoy), str(bacteria), "", root / "prepared")
            self.assertIsNone(manifest["host"])
            dna = Path(manifest["bacteria"]["dna"]["fasta"]).read_text()
            rna = Path(manifest["bacteria"]["rna"]["fasta"]).read_text()
            annotations = Path(manifest["bacteria"]["rna"]["gff"]).read_text()
            self.assertEqual(dna.count(">"), 1)
            self.assertEqual(rna.count(">"), 3)
            self.assertIn(">RNA__lcl|NC_001.1_cds_1", rna)
            self.assertIn("\tCDS\t1\t60\t.\t+\t0\tgene_id=gene1", annotations)
            self.assertIn("\trRNA\t1\t8\t", annotations)
            self.assertIn("\ttRNA\t1\t6\t", annotations)
            self.assertEqual(Path(manifest["decoy"]["fasta"]).read_text().count(">"), 2)
            self.assertEqual(json.loads((root / "prepared/manifest.json").read_text()), manifest)

    def test_recursive_host_gff_shards_preserve_transcript_cds_coordinates(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            decoy, bacteria = self.fixture(root)
            host = root / "human"
            write(host / "human_DNA/GCF_genomic.fna.gz", ">NC_host.1\n" + "G" * 100 + "\n")
            write(host / "human_DNA/GCF_genomic.gff.gz", gff("NC_host.1", end=50))
            write(host / "human_RNA/human.1.rna.fna.gz", ">NM_001.1\n" + "G" * 40 + "\n")
            write(host / "human_RNA/human.1.rna.gff3.gz", gff("NM_001.1", start=10, end=30))
            write(host / "human_RNA/human.2.rna.fna.gz", ">NR_001.1\n" + "A" * 20 + "\n")
            write(host / "human_RNA/human.2.rna.gff3.gz", gff("NR_001.1", "ncRNA", end=20))
            manifest = prepare_references(str(decoy), str(bacteria), str(host), root / "prepared")
            rna = manifest["host"]["rna"]
            self.assertEqual(Path(rna["fasta"]).read_text().count(">"), 2)
            text = Path(rna["gff"]).read_text()
            self.assertIn("RNA__NM_001.1\treference\tCDS\t10\t30", text)
            self.assertIn("RNA__NM_001.1\treference\ttranscript\t1\t40", text)
            self.assertIn("RNA__NR_001.1\treference\tncRNA\t1\t20", text)

    @unittest.skipIf(importlib.util.find_spec("Bio") is None, "Biopython not installed")
    def test_genbank_shards_and_gff_preference(self):
        from Bio import SeqIO
        from Bio.Seq import Seq
        from Bio.SeqFeature import SeqFeature, SimpleLocation
        from Bio.SeqRecord import SeqRecord

        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            decoy, bacteria = self.fixture(root)
            # A malformed redundant GBFF must not be read when matching GFF exists.
            write(bacteria / "GCF_genomic.gbff.gz", "redundant and invalid")
            host = root / "mouse"
            write(host / "mouse_DNA/assembly.fna.gz", ">NC_mouse.1\n" + "G" * 100 + "\n")
            record = SeqRecord(Seq("G" * 100), id="NC_mouse.1", annotations={"molecule_type": "DNA"})
            record.features = [SeqFeature(SimpleLocation(9, 50, strand=1), type="CDS", qualifiers={"locus_tag": ["m_gene"]})]
            target = host / "mouse_DNA/assembly.gbff.gz"
            with gzip.open(target, "wt") as handle:
                SeqIO.write(record, handle, "genbank")
            write(host / "mouse_RNA/mouse.1.rna.fna.gz", ">NM_mouse.1\n" + "G" * 40 + "\n")
            record = SeqRecord(Seq("G" * 40), id="NM_mouse.1", annotations={"molecule_type": "mRNA"})
            record.features = [SeqFeature(SimpleLocation(5, 30, strand=1), type="CDS", qualifiers={"gene": ["m_gene"]})]
            target = host / "mouse_RNA/mouse.1.rna.gbff.gz"
            with gzip.open(target, "wt") as handle:
                SeqIO.write(record, handle, "genbank")
            manifest = prepare_references(str(decoy), str(bacteria), str(host), root / "prepared")
            text = Path(manifest["host"]["rna"]["gff"]).read_text()
            self.assertIn("RNA__NM_mouse.1\treference\tCDS\t6\t30", text)
            self.assertIn("gene_id=m_gene", text)

    def test_legacy_basename_inputs_without_rna(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            decoy, bacteria = self.fixture(root)
            manifest = prepare_references(str(decoy / "pangenome"), str(bacteria / "GCF_genomic"), "", root / "prepared")
            self.assertIsNone(manifest["bacteria"]["rna"])

    def test_rejects_annotation_reference_mismatch_and_conflicting_ids(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            decoy, bacteria = self.fixture(root)
            write(bacteria / "GCF_genomic.gff.gz", gff("wrong_reference"))
            with self.assertRaisesRegex(ValueError, "no matching FASTA"):
                prepare_references(str(decoy), str(bacteria), "", root / "prepared")
            write(bacteria / "GCF_genomic.gff.gz", gff("NC_001.1"))
            write(bacteria / "other_genomic.fa", ">NC_001.1\nTTTT\n")
            with self.assertRaisesRegex(ValueError, "Conflicting reference sequences"):
                prepare_references(str(decoy), str(bacteria), "", root / "prepared")

    def test_rejects_transcript_coordinates_from_genome(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            decoy, bacteria = self.fixture(root)
            write(bacteria / "bacteria.rna.fna", ">transcript\nACGT\n")
            write(bacteria / "bacteria.rna.gff", gff("transcript", end=100))
            with self.assertRaisesRegex(ValueError, "outside FASTA"):
                prepare_references(str(decoy), str(bacteria), "", root / "prepared")

    def test_genomic_ncrna_exons_preserve_parent_hierarchy(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            decoy, bacteria = self.fixture(root)
            write(bacteria / "GCF_genomic.gff.gz",
                  "##gff-version 3\n"
                  "NC_001.1\ttest\tgene\t1\t90\t.\t+\t.\tID=gene-nc;gene=lnc1\n"
                  "NC_001.1\ttest\tlnc_RNA\t1\t90\t.\t+\t.\tID=rna-nc;Parent=gene-nc;gene=lnc1\n"
                  "NC_001.1\ttest\texon\t1\t20\t.\t+\t.\tID=exon-1;Parent=rna-nc\n"
                  "NC_001.1\ttest\texon\t70\t90\t.\t+\t.\tID=exon-2;Parent=rna-nc\n")
            manifest = prepare_references(str(decoy), str(bacteria), "", root / "prepared")
            annotations = Path(manifest["bacteria"]["dna"]["gff"]).read_text()
            self.assertEqual(annotations.count("\texon\t"), 2)
            self.assertIn("ID=exon-1;Parent=rna-nc", annotations)
            self.assertIn("ID=rna-nc;Parent=gene-nc", annotations)

    def test_streamed_fasta_wrapping_deduplicates_identical_sequences(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            decoy, bacteria = self.fixture(root)
            write(bacteria / "identical_genomic.fasta", ">NC_001.1 duplicate\n" + "\n".join(["ACGT"] * 25) + "\n")
            manifest = prepare_references(str(decoy), str(bacteria), "", root / "prepared")
            lines = Path(manifest["bacteria"]["dna"]["fasta"]).read_text().splitlines()
            self.assertEqual(lines, [">NC_001.1", "ACGT" * 20, "ACGT" * 5])


if __name__ == "__main__":
    unittest.main()
