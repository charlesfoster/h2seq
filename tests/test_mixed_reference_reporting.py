import csv
import gzip
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


PROJECT_DIR = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_DIR / "bin"))


class MixedReferenceTests(unittest.TestCase):
    @unittest.skipUnless(shutil.which("samtools"), "samtools is required for the synthetic BAM regression test")
    def test_synthetic_bam_preserves_known_mixture_proportions(self):
        samtools = shutil.which("samtools")
        sequence = "A" * 20
        qualities = "I" * 20
        records = []
        for index in range(90):
            records.append(f"three_a_{index}\t0\t3a_ref\t{index + 1}\t60\t20M\t*\t0\t0\t{sequence}\t{qualities}")
        for index in range(10):
            records.append(f"one_a_{index}\t0\t1a_ref\t{index + 1}\t60\t20M\t*\t0\t0\t{sequence}\t{qualities}")
        for index in range(5):
            records.append(f"ambiguous_{index}\t0\t3a_ref\t{index + 201}\t0\t20M\t*\t0\t0\t{sequence}\t{qualities}")
        for index in range(5):
            records.append(f"unmapped_{index}\t4\t*\t0\t0\t*\t*\t0\t0\t{sequence}\t{qualities}")

        sam = "\n".join(
            [
                "@HD\tVN:1.6\tSO:unsorted",
                "@SQ\tSN:1a_ref\tLN:1000",
                "@SQ\tSN:3a_ref\tLN:1000",
                *records,
                "",
            ]
        )

        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            source_sam = root / "synthetic.sam"
            source_bam = root / "synthetic.bam"
            filtered_sam = root / "filtered.sam"
            filtered_bam = root / "synthetic.mixed_filtered.bam"
            assignment = root / "synthetic/short_reads/coverage/synthetic.mixed_assignment.tsv"
            assignment.parent.mkdir(parents=True)
            source_sam.write_text(sam)
            subprocess.run([samtools, "view", "-b", "-o", source_bam, source_sam], check=True)

            collate = subprocess.Popen([samtools, "collate", "-O", "-u", source_bam], stdout=subprocess.PIPE)
            view = subprocess.Popen([samtools, "view", "-h", "-"], stdin=collate.stdout, stdout=subprocess.PIPE, text=True)
            collate.stdout.close()
            awk = subprocess.run(
                [
                    "awk",
                    "-v",
                    "sample_id=synthetic",
                    "-v",
                    "read_type=short",
                    "-v",
                    "min_mapq=10",
                    "-v",
                    f"summary_output={assignment}",
                    "-f",
                    str(PROJECT_DIR / "bin/filter_mixed_reference_alignments.awk"),
                ],
                stdin=view.stdout,
                text=True,
                capture_output=True,
                check=True,
            )
            view.stdout.close()
            self.assertEqual(view.wait(), 0)
            self.assertEqual(collate.wait(), 0)
            filtered_sam.write_text(awk.stdout)
            subprocess.run([samtools, "view", "-b", "-o", filtered_bam, filtered_sam], check=True)
            subprocess.run([samtools, "sort", "-o", root / "sorted.bam", filtered_bam], check=True)
            filtered_bam = root / "sorted.bam"
            subprocess.run([samtools, "index", "-c", filtered_bam], check=True)

            self.assertEqual(
                subprocess.check_output([samtools, "view", "-c", filtered_bam], text=True).strip(), "100"
            )
            self.assertEqual(
                subprocess.check_output([samtools, "view", "-c", filtered_bam, "3a_ref"], text=True).strip(), "90"
            )
            self.assertEqual(
                subprocess.check_output([samtools, "view", "-c", filtered_bam, "1a_ref"], text=True).strip(), "10"
            )

            with assignment.open(newline="") as handle:
                assignment_rows = list(csv.DictReader(handle, delimiter="\t"))
            by_key = {(row["assignment"], row["reference_name"]): row for row in assignment_rows}
            self.assertEqual(by_key[("assigned", "3a_ref")]["fragments"], "90")
            self.assertEqual(by_key[("assigned", "1a_ref")]["fragments"], "10")
            self.assertEqual(by_key[("ambiguous", "")]["fragments"], "5")
            self.assertEqual(by_key[("unassigned", "")]["fragments"], "5")

            selection_dir = root / "synthetic/short_reads/reference_selection"
            selection_dir.mkdir(parents=True)
            (selection_dir / "synthetic.best_reference.tsv").write_text(
                "sample_id\tgenotype\tsubtype\tbest_ref\tmixed_infection\tsecondary_genotypes\tselection_status\tqc_fail_reason\n"
                "synthetic\t3\t3a\t3a_ref\ttrue\t1a:0.1\tpass\t\n"
            )
            (assignment.parent / "synthetic.coverage_summary.tsv").write_text(
                "sample_id\tread_type\treference_name\treference_length\tpositions_covered\tgenome_coverage_pct\tmean_depth\n"
                "synthetic\tshort\t1a_ref\t1000\t800\t80.0\t10.0\n"
                "synthetic\tshort\t3a_ref\t1000\t900\t90.0\t90.0\n"
            )
            combined = root / "combined.csv"
            components = root / "components.csv"
            subprocess.run(
                [
                    sys.executable,
                    str(PROJECT_DIR / "bin/build_run_summary.py"),
                    "--outdir",
                    str(root),
                    "--pipeline-version",
                    "test",
                    "--output",
                    str(combined),
                    "--component-output",
                    str(components),
                ],
                check=True,
            )
            with combined.open(newline="") as handle:
                summary = next(csv.DictReader(handle))
            self.assertEqual(summary["component_fractions"], "3a:0.9000;1a:0.1000")
            self.assertEqual(summary["ambiguous_fragments"], "5")

    def test_fragment_filter_separates_assigned_and_ambiguous(self):
        sam = """@HD\tVN:1.6\tSO:queryname
@SQ\tSN:1a_ref\tLN:100
@SQ\tSN:3a_ref\tLN:100
read1\t99\t3a_ref\t1\t60\t50M\t=\t1\t50\tAAAA\tIIII
read1\t147\t3a_ref\t1\t60\t50M\t=\t1\t-50\tTTTT\tIIII
read1\t355\t1a_ref\t1\t0\t50M\t=\t1\t50\tAAAA\tIIII
read2\t99\t1a_ref\t1\t50\t50M\t=\t1\t50\tAAAA\tIIII
read2\t147\t1a_ref\t1\t50\t50M\t=\t1\t-50\tTTTT\tIIII
read3\t99\t3a_ref\t1\t0\t50M\t=\t1\t50\tAAAA\tIIII
read3\t147\t3a_ref\t1\t60\t50M\t=\t1\t-50\tTTTT\tIIII
read4\t99\t3a_ref\t1\t60\t50M\t1a_ref\t1\t0\tAAAA\tIIII
read4\t147\t1a_ref\t1\t60\t50M\t3a_ref\t1\t0\tTTTT\tIIII
read5\t77\t*\t0\t0\t*\t*\t0\t0\tAAAA\tIIII
read5\t141\t*\t0\t0\t*\t*\t0\t0\tTTTT\tIIII
"""
        with tempfile.TemporaryDirectory() as tmpdir:
            summary = Path(tmpdir) / "assignment.tsv"
            command = [
                "awk",
                "-v",
                "sample_id=sample",
                "-v",
                "read_type=short",
                "-v",
                "min_mapq=10",
                "-v",
                f"summary_output={summary}",
                "-f",
                str(PROJECT_DIR / "bin/filter_mixed_reference_alignments.awk"),
            ]
            result = subprocess.run(command, input=sam, text=True, capture_output=True, check=True)
            self.assertIn("read1", result.stdout)
            self.assertIn("read2", result.stdout)
            self.assertNotIn("\t355\t", result.stdout)
            self.assertNotIn("read3", result.stdout)
            self.assertNotIn("read4", result.stdout)
            self.assertNotIn("read5", result.stdout)

            with summary.open(newline="") as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            by_key = {(row["assignment"], row["reference_name"]): row for row in rows}
            self.assertEqual(by_key[("assigned", "1a_ref")]["fragments"], "1")
            self.assertEqual(by_key[("assigned", "3a_ref")]["fragments"], "1")
            self.assertEqual(by_key[("ambiguous", "")]["fragments"], "2")
            self.assertEqual(by_key[("unassigned", "")]["fragments"], "1")

    def test_single_reference_bam_is_passed_through(self):
        sam = """@HD\tVN:1.6\tSO:queryname
@SQ\tSN:3a_ref\tLN:100
mapped\t0\t3a_ref\t1\t0\t4M\t*\t0\t0\tAAAA\tIIII
mapped\t256\t3a_ref\t2\t0\t4M\t*\t0\t0\tAAAA\tIIII
unmapped\t4\t*\t0\t0\t*\t*\t0\t0\tTTTT\tIIII
"""
        with tempfile.TemporaryDirectory() as tmpdir:
            summary = Path(tmpdir) / "assignment.tsv"
            command = [
                "awk",
                "-v",
                "sample_id=sample",
                "-v",
                "read_type=long",
                "-v",
                "min_mapq=10",
                "-v",
                f"summary_output={summary}",
                "-f",
                str(PROJECT_DIR / "bin/filter_mixed_reference_alignments.awk"),
            ]
            result = subprocess.run(command, input=sam, text=True, capture_output=True, check=True)
            self.assertEqual(result.stdout, sam)

    def test_multi_reference_coverage_and_sample_summary_use_main_reference(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            coverage_dir = root / "sample/short_reads/coverage"
            selection_dir = root / "sample/short_reads/reference_selection"
            coverage_dir.mkdir(parents=True)
            selection_dir.mkdir(parents=True)

            regions = coverage_dir / "sample.regions.bed.gz"
            per_base = coverage_dir / "sample.per-base.bed.gz"
            with gzip.open(regions, "wt") as handle:
                handle.write("1a_ref\t0\t10\twhole_genome\t5\n")
                handle.write("3a_ref\t0\t10\twhole_genome\t20\n")
            with gzip.open(per_base, "wt") as handle:
                handle.write("1a_ref\t0\t10\t5\n")
                handle.write("3a_ref\t0\t10\t20\n")

            coverage = coverage_dir / "sample.coverage_summary.tsv"
            subprocess.run(
                [
                    sys.executable,
                    str(PROJECT_DIR / "bin/summarise_mosdepth_coverage.py"),
                    "--sample-id",
                    "sample",
                    "--read-type",
                    "short",
                    "--reference-name",
                    "3a_ref",
                    "--genome-bed-gz",
                    str(regions),
                    "--per-base-bed-gz",
                    str(per_base),
                    "--min-depth",
                    "15",
                    "--summary-output",
                    str(coverage),
                ],
                check=True,
            )

            best_ref = selection_dir / "sample.best_reference.tsv"
            best_ref.write_text(
                "sample_id\tgenotype\tsubtype\tbest_ref\tmixed_infection\tsecondary_genotypes\tselection_status\tqc_fail_reason\n"
                "sample\t3\t3a\t3a_ref\ttrue\t1a:0.1\tpass\t\n"
            )
            assignment = coverage_dir / "sample.mixed_assignment.tsv"
            assignment.write_text(
                "sample_id\tread_type\tassignment\treference_name\tfragments\tprimary_alignment_records\tmin_mapq\treference_count\n"
                "sample\tshort\tassigned\t1a_ref\t10\t20\t10\t2\n"
                "sample\tshort\tassigned\t3a_ref\t90\t180\t10\t2\n"
                "sample\tshort\tambiguous\t\t2\t4\t10\t2\n"
            )
            for reference_name, polyprotein_coverage in [("3a_ref", "80.0"), ("1a_ref", "95.0")]:
                hcv_coverage = coverage_dir / f"sample.{reference_name}.hcv_glue_coverage.tsv"
                hcv_coverage.write_text(
                    "sample_id\tread_type\treference_name\tfeature\tcoverage_pct\n"
                    f"sample\tshort\t{reference_name}\tPolyprotein\t{polyprotein_coverage}\n"
                    f"sample\tshort\t{reference_name}\tCore\t100.0\n"
                )

            combined = root / "combined.csv"
            components = root / "components.csv"
            subprocess.run(
                [
                    sys.executable,
                    str(PROJECT_DIR / "bin/build_run_summary.py"),
                    "--outdir",
                    str(root),
                    "--pipeline-version",
                    "test",
                    "--output",
                    str(combined),
                    "--component-output",
                    str(components),
                ],
                check=True,
            )

            with combined.open(newline="") as handle:
                sample_rows = list(csv.DictReader(handle))
            self.assertEqual(len(sample_rows), 1)
            self.assertEqual(sample_rows[0]["selected_reference"], "3a_ref")
            self.assertEqual(float(sample_rows[0]["genome_coverage"]), 100.0)
            self.assertEqual(float(sample_rows[0]["mean_depth"]), 20.0)
            self.assertEqual(sample_rows[0]["mixed_infection"], "true")
            self.assertEqual(sample_rows[0]["secondary_subtypes"], "1a")
            self.assertEqual(sample_rows[0]["component_fractions"], "3a:0.9000;1a:0.1000")
            self.assertEqual(sample_rows[0]["ambiguous_fragments"], "2")
            self.assertEqual(sample_rows[0]["polyprotein_coverage"], "80.0")
            self.assertNotIn("secondary_genotypes", sample_rows[0])

            with components.open(newline="") as handle:
                component_rows = list(csv.DictReader(handle))
            self.assertEqual(len(component_rows), 2)
            by_role = {row["component_role"]: row for row in component_rows}
            self.assertEqual(by_role["main"]["reference_name"], "3a_ref")
            self.assertEqual(by_role["secondary"]["reference_name"], "1a_ref")
            self.assertEqual(by_role["main"]["ambiguous_fragments"], "2")
            self.assertEqual(float(by_role["main"]["assigned_fraction_of_assigned"]), 0.9)
            self.assertEqual(by_role["main"]["polyprotein_coverage"], "80.0")
            self.assertEqual(by_role["secondary"]["polyprotein_coverage"], "95.0")

            from render_hcv_report import read_assignment_stats

            assignment_stats = read_assignment_stats(assignment, "1a_ref")
            self.assertEqual(assignment_stats["assigned_fragments"], 10)
            self.assertEqual(assignment_stats["assigned_fraction"], 0.1)
            self.assertEqual(assignment_stats["ambiguous_fragments"], 2)

            from build_multiqc_sections import build_coverage_section, build_region_coverage_section

            coverage_payload = build_coverage_section(root)
            self.assertEqual(
                set(coverage_payload["data"]),
                {"sample (short) [main: 3a_ref]", "sample (short) [secondary: 1a_ref]"},
            )
            self.assertEqual(coverage_payload["data"]["sample (short) [main: 3a_ref]"]["coverage_10x_pct"], 100.0)

            region_payload = build_region_coverage_section(root)
            self.assertEqual(
                set(region_payload["ycats"]),
                {"sample (short) [main: 3a_ref]", "sample (short) [secondary: 1a_ref]"},
            )


if __name__ == "__main__":
    unittest.main()
