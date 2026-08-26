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
    parser.add_argument("--component-output", required=True)
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
            "mixed_infection": "",
            "secondary_subtypes": "",
            "component_fractions": "",
            "ambiguous_fragments": "",
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


def reference_genotype_subtype(reference_name):
    subtype = reference_name.split("_", 1)[0] if reference_name else ""
    genotype = "".join(char for char in subtype if char.isdigit())
    return genotype, subtype


def apply_final_component_assignments(rows, assignment_by_key):
    for key, assignment_rows in assignment_by_key.items():
        row = ensure_row(rows, *key)
        assigned = [record for record in assignment_rows if record.get("assignment") == "assigned"]
        assigned_total = sum(int(float(record.get("fragments") or 0)) for record in assigned)
        ambiguous = sum(
            int(float(record.get("fragments") or 0))
            for record in assignment_rows
            if record.get("assignment") == "ambiguous"
        )
        row["ambiguous_fragments"] = ambiguous
        if not assigned_total:
            continue

        selected_reference = row.get("selected_reference", "")
        components = []
        for record in assigned:
            reference_name = record.get("reference_name", "")
            _genotype, subtype = reference_genotype_subtype(reference_name)
            fragments = int(float(record.get("fragments") or 0))
            components.append(
                {
                    "reference_name": reference_name,
                    "subtype": subtype,
                    "fraction": fragments / assigned_total,
                }
            )
        components.sort(key=lambda component: (component["reference_name"] != selected_reference, -component["fraction"]))
        row["mixed_infection"] = "true" if len(components) > 1 else "false"
        row["secondary_subtypes"] = ";".join(
            dict.fromkeys(
                component["subtype"]
                for component in components
                if component["reference_name"] != selected_reference
            )
        )
        row["component_fractions"] = ";".join(
            f"{component['subtype']}:{component['fraction']:.4f}" for component in components
        )


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
    coverage_by_key = {}
    assignment_by_key = {}
    hcv_coverage_by_key = {}

    for path in sorted(outdir.rglob("*.coverage_summary.tsv")):
        for record in load_tsv_rows(path):
            key = row_key(record["sample_id"], record["read_type"])
            coverage_by_key.setdefault(key, {})[record.get("reference_name", "")] = record
            ensure_row(rows, *key)

    for path in sorted(outdir.rglob("*.mixed_assignment.tsv")):
        for record in load_tsv_rows(path):
            key = row_key(record["sample_id"], record["read_type"])
            assignment_by_key.setdefault(key, []).append(record)

    for path in sorted(outdir.rglob("*.hcv_glue_coverage.tsv")):
        for record in load_tsv_rows(path):
            if record["feature"] != "Polyprotein":
                continue
            key = row_key(record["sample_id"], record["read_type"])
            hcv_coverage_by_key.setdefault(key, {})[record.get("reference_name", "")] = record.get(
                "coverage_pct", ""
            )
            ensure_row(rows, *key)

    for path in sorted(outdir.rglob("*.best_reference.tsv")):
        for record in load_tsv_rows(path):
            sample_id, read_type = infer_sample_and_read_type(path)
            row = ensure_row(rows, sample_id, read_type)
            row["designated_genotype"] = record.get("genotype", "")
            row["designated_subtype"] = record.get("subtype", "")
            row["mixed_infection"] = record.get("mixed_infection", "")
            row["selected_reference"] = record.get("best_ref", row["selected_reference"])
            row["qc_status"] = record.get("selection_status", row["qc_status"])
            row["qc_fail_reason"] = record.get("qc_fail_reason", row["qc_fail_reason"])

    for key, reference_records in coverage_by_key.items():
        row = ensure_row(rows, *key)
        selected_reference = row.get("selected_reference", "")
        if not selected_reference and len(reference_records) == 1:
            selected_reference = next(iter(reference_records))
            row["selected_reference"] = selected_reference
        if selected_reference not in reference_records:
            raise ValueError(
                f"Selected reference {selected_reference!r} has no coverage row for {key}; "
                f"available references: {sorted(reference_records)}"
            )
        record = reference_records[selected_reference]
        row["reference_length"] = record.get("reference_length", "")
        row["genome_coverage"] = record.get("genome_coverage_pct", "")
        row["mean_depth"] = record.get("mean_depth", "")
        row["polyprotein_coverage"] = hcv_coverage_by_key.get(key, {}).get(selected_reference, "")

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

    apply_final_component_assignments(rows, assignment_by_key)
    finalize_qc(rows, args.min_reference_coverage_pct)

    for row in rows.values():
        row["pipeline_version"] = args.pipeline_version
        row["analysis_date"] = date.today().isoformat()

    fieldnames = [
        "sample_id",
        "read_type",
        "designated_genotype",
        "designated_subtype",
        "mixed_infection",
        "secondary_subtypes",
        "component_fractions",
        "ambiguous_fragments",
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

    component_fieldnames = [
        "sample_id",
        "read_type",
        "component_role",
        "genotype",
        "subtype",
        "reference_name",
        "reference_length",
        "positions_covered",
        "genome_coverage",
        "mean_depth",
        "polyprotein_coverage",
        "assigned_fragments",
        "assigned_alignment_records",
        "assigned_fraction_of_assigned",
        "ambiguous_fragments",
        "unassigned_fragments",
        "assignment_min_mapq",
        "pipeline_version",
        "analysis_date",
    ]
    with open(args.component_output, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=component_fieldnames)
        writer.writeheader()
        for key in sorted(coverage_by_key):
            sample_row = rows[key]
            assignment_rows = assignment_by_key.get(key, [])
            assigned_by_reference = {
                record.get("reference_name", ""): record
                for record in assignment_rows
                if record.get("assignment") == "assigned"
            }
            assigned_total = sum(int(float(record.get("fragments") or 0)) for record in assigned_by_reference.values())
            ambiguous_fragments = sum(
                int(float(record.get("fragments") or 0))
                for record in assignment_rows
                if record.get("assignment") == "ambiguous"
            )
            unassigned_fragments = sum(
                int(float(record.get("fragments") or 0))
                for record in assignment_rows
                if record.get("assignment") == "unassigned"
            )
            assignment_min_mapq = next(
                (record.get("min_mapq", "") for record in assignment_rows if record.get("min_mapq", "") != ""),
                "",
            )
            component_records = sorted(
                coverage_by_key[key].items(),
                key=lambda item: (item[0] != sample_row["selected_reference"], item[0]),
            )
            for reference_name, coverage in component_records:
                assignment = assigned_by_reference.get(reference_name, {})
                assigned_fragments = int(float(assignment.get("fragments") or 0))
                genotype, subtype = reference_genotype_subtype(reference_name)
                writer.writerow(
                    {
                        "sample_id": key[0],
                        "read_type": key[1],
                        "component_role": "main" if reference_name == sample_row["selected_reference"] else "secondary",
                        "genotype": genotype,
                        "subtype": subtype,
                        "reference_name": reference_name,
                        "reference_length": coverage.get("reference_length", ""),
                        "positions_covered": coverage.get("positions_covered", ""),
                        "genome_coverage": coverage.get("genome_coverage_pct", ""),
                        "mean_depth": coverage.get("mean_depth", ""),
                        "polyprotein_coverage": hcv_coverage_by_key.get(key, {}).get(reference_name, ""),
                        "assigned_fragments": assigned_fragments,
                        "assigned_alignment_records": assignment.get("primary_alignment_records", ""),
                        "assigned_fraction_of_assigned": (assigned_fragments / assigned_total) if assigned_total else "",
                        "ambiguous_fragments": ambiguous_fragments,
                        "unassigned_fragments": unassigned_fragments,
                        "assignment_min_mapq": assignment_min_mapq,
                        "pipeline_version": args.pipeline_version,
                        "analysis_date": date.today().isoformat(),
                    }
                )

    if args.multiqc_output:
        table_headers = {
            "read_type": {"title": "Read Type"},
            "designated_genotype": {"title": "Genotype"},
            "designated_subtype": {"title": "Subtype"},
            "mixed_infection": {"title": "Mixed Infection"},
            "secondary_subtypes": {"title": "Secondary Subtypes"},
            "component_fractions": {"title": "Final Component Fractions"},
            "ambiguous_fragments": {"title": "Ambiguous Fragments"},
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
