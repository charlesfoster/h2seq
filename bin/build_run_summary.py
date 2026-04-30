#!/usr/bin/env python3
import argparse
import csv
import json
from datetime import date
from pathlib import Path


def parse_args():
    parser = argparse.ArgumentParser(description="Build a run-level sample summary CSV.")
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--pipeline-version", default="")
    parser.add_argument("--min-reference-coverage-pct", type=float, default=50.0)
    parser.add_argument("--output", required=True)
    parser.add_argument("--multiqc-output")
    return parser.parse_args()


def row_key(sample_id, read_type):
    return (sample_id, read_type)


def infer_sample_and_read_type(path):
    parts = Path(path).parts
    sample_id = None
    read_type = None
    for i, part in enumerate(parts):
        if part in {"long_reads", "short_reads"} and i > 0:
            sample_id = parts[i - 1]
            read_type = "long" if part == "long_reads" else "short"
            break
    if sample_id is None:
        name = Path(path).name
        sample_id = name.split(".")[0]
        read_type = "unknown"
    return sample_id, read_type


def ensure_row(rows, sample_id, read_type):
    key = row_key(sample_id, read_type)
    if key not in rows:
        rows[key] = {
            "sample_id": sample_id,
            "read_type": read_type,
            "designated_genotype": "",
            "designated_subtype": "",
            "selected_reference": "",
            "reference_length": "",
            "genome_coverage": "",
            "polyprotein_coverage": "",
            "mean_depth": "",
            "total_reads": "",
            "reads_passing_qc": "",
            "qc_status": "",
            "qc_fail_reason": "",
            "genotype_subtype_status": "",
            "pipeline_version": "",
            "analysis_date": "",
        }
    return rows[key]


def load_tsv_rows(path):
    with open(path, newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def extract_seqkit_count(records, path):
    counts = []
    for record in records:
        value = (record.get("num_seqs", "") or "").strip()
        if not value:
            raise ValueError(f"Missing num_seqs value in {path}")
        counts.append(value)

    if not counts:
        raise ValueError(f"No seqkit rows found in {path}")
    if len(set(counts)) != 1:
        raise ValueError(f"Inconsistent num_seqs values in {path}: {counts}")

    return counts[0]


def to_float(value):
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def add_qc_fail(row, reason):
    reasons = [item for item in row.get("qc_fail_reason", "").split(";") if item]
    if reason not in reasons:
        reasons.append(reason)
    row["qc_status"] = "qc_fail"
    row["qc_fail_reason"] = ";".join(reasons)


def finalize_qc(rows, min_reference_coverage_pct):
    for row in rows.values():
        coverage_pct = to_float(row.get("genome_coverage"))
        if coverage_pct is not None and coverage_pct < min_reference_coverage_pct:
            add_qc_fail(row, "low_ref_coverage")

        if row.get("qc_status") == "qc_fail":
            row["genotype_subtype_status"] = "unreliable_qc_fail"
            row["designated_genotype"] = ""
            row["designated_subtype"] = ""
        elif row.get("designated_genotype") or row.get("designated_subtype"):
            row["qc_status"] = row.get("qc_status") or "pass"
            row["genotype_subtype_status"] = "assigned"
        else:
            row["qc_status"] = row.get("qc_status") or ""
            row["genotype_subtype_status"] = ""


def main():
    args = parse_args()
    outdir = Path(args.outdir)
    rows = {}

    for path in sorted(outdir.rglob("*.coverage_summary.tsv")):
        for record in load_tsv_rows(path):
            row = ensure_row(rows, record["sample_id"], record["read_type"])
            row["selected_reference"] = row["selected_reference"] or record.get("reference_name", "")
            row["reference_length"] = record.get("reference_length", "")
            row["genome_coverage"] = record.get("genome_coverage_pct", "")
            row["mean_depth"] = record.get("mean_depth", "")

    for path in sorted(outdir.rglob("*.hcv_glue_coverage.tsv")):
        for record in load_tsv_rows(path):
            row = ensure_row(rows, record["sample_id"], record["read_type"])
            if record["feature"] == "Polyprotein":
                row["polyprotein_coverage"] = record.get("coverage_pct", "")

    for path in sorted(outdir.rglob("*.best_reference.tsv")):
        for record in load_tsv_rows(path):
            sample_id, read_type = infer_sample_and_read_type(path)
            row = ensure_row(rows, sample_id, read_type)
            row["designated_genotype"] = record.get("genotype", "")
            row["designated_subtype"] = record.get("subtype", "")
            row["selected_reference"] = record.get("best_ref", row["selected_reference"])
            row["qc_status"] = record.get("selection_status", row["qc_status"])
            row["qc_fail_reason"] = record.get("qc_fail_reason", row["qc_fail_reason"])

    for path in sorted(outdir.rglob("*.empty.fasta")):
        if "consensus_main.empty.fasta" not in path.name:
            continue
        sample_id, read_type = infer_sample_and_read_type(path)
        row = ensure_row(rows, sample_id, read_type)
        add_qc_fail(row, "genome_insufficient_for_hcv-glue")

    for path in sorted(outdir.rglob("*.fastp.json")):
        sample_id, read_type = infer_sample_and_read_type(path)
        row = ensure_row(rows, sample_id, read_type)
        with open(path) as handle:
            payload = json.load(handle)
        before = payload.get("summary", {}).get("before_filtering", {})
        after = payload.get("summary", {}).get("after_filtering", {})
        filtering = payload.get("filtering_result", {})
        row["total_reads"] = before.get("total_reads", "")
        row["reads_passing_qc"] = filtering.get("passed_filter_reads", after.get("total_reads", ""))

    seqkit_paths = (
        sorted(outdir.rglob("*.raw_long.tsv"))
        + sorted(outdir.rglob("*.clean_long.tsv"))
        + sorted(outdir.rglob("*.raw_short.tsv"))
        + sorted(outdir.rglob("*.clean_short.tsv"))
    )

    for path in seqkit_paths:
        records = load_tsv_rows(path)
        if not records:
            continue
        sample_id, read_type = infer_sample_and_read_type(path)
        row = ensure_row(rows, sample_id, read_type)
        count = extract_seqkit_count(records, path)
        name = Path(path).name
        if ".raw_long." in name or ".raw_short." in name:
            row["total_reads"] = count
        elif ".clean_long." in name or ".clean_short." in name:
            row["reads_passing_qc"] = count

    finalize_qc(rows, args.min_reference_coverage_pct)

    for row in rows.values():
        row["pipeline_version"] = args.pipeline_version
        row["analysis_date"] = date.today().isoformat()

    fieldnames = [
        "sample_id",
        "read_type",
        "designated_genotype",
        "designated_subtype",
        "selected_reference",
        "reference_length",
        "genome_coverage",
        "polyprotein_coverage",
        "mean_depth",
        "total_reads",
        "reads_passing_qc",
        "qc_status",
        "qc_fail_reason",
        "genotype_subtype_status",
        "pipeline_version",
        "analysis_date",
    ]

    with open(args.output, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for key in sorted(rows):
            writer.writerow(rows[key])

    if args.multiqc_output:
        table_headers = {
            "read_type": {"title": "Read Type"},
            "designated_genotype": {"title": "Genotype"},
            "designated_subtype": {"title": "Subtype"},
            "selected_reference": {"title": "Reference"},
            "reference_length": {"title": "Reference Length"},
            "genome_coverage": {"title": "Genome Cov %", "format": "{:,.1f}"},
            "polyprotein_coverage": {"title": "Polyprotein Cov %", "format": "{:,.1f}"},
            "mean_depth": {"title": "Mean Depth", "format": "{:,.2f}"},
            "total_reads": {"title": "Reads Analysed"},
            "reads_passing_qc": {"title": "Reads Passing QC"},
            "qc_status": {"title": "QC Status"},
            "qc_fail_reason": {"title": "QC Fail Reason"},
            "genotype_subtype_status": {"title": "Genotype/Subtype Status"},
            "pipeline_version": {"title": "Pipeline Version"},
            "analysis_date": {"title": "Analysis Date"},
        }

        table_data = {}
        for key in sorted(rows):
            row = rows[key]
            sample_label = row["sample_id"]
            if row["read_type"] not in {"", "unknown"}:
                sample_label = f"{sample_label} ({row['read_type']})"
            table_data[sample_label] = {column: row.get(column, "") for column in table_headers}

        multiqc_payload = {
            "id": "h2seq_run_summary",
            "section_name": "h2seq Run Summary",
            "description": "Per-sample summary of QC, mapping, coverage and consensus outputs.",
            "plot_type": "table",
            "pconfig": {
                "id": "h2seq_run_summary_table",
                "title": "h2seq Run Summary",
            },
            "headers": table_headers,
            "data": table_data,
        }

        with open(args.multiqc_output, "w") as handle:
            json.dump(multiqc_payload, handle, indent=2)


if __name__ == "__main__":
    main()
