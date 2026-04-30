#!/usr/bin/env python3
import argparse
import csv
import gzip
import re
import sys
from collections import OrderedDict
from pathlib import Path


CIGAR_RE = re.compile(r"(\d+)([MIDNSHP=X])")


def parse_args():
    parser = argparse.ArgumentParser(
        description="Select the best reference from competitive minimap2 alignments."
    )
    parser.add_argument("--primary-sam", required=True, help="Headered SAM or SAM.GZ of primary mapped alignments")
    parser.add_argument("--bam", required=True, help="Competitive BAM used to derive the SAM")
    parser.add_argument("--panel-fasta", required=True, help="Reference-panel FASTA")
    parser.add_argument("--sample-name", required=True, help="Sample name")
    parser.add_argument("--read-type", choices=["long", "short"], required=True, help="Read type being ranked")
    parser.add_argument("--output", required=True, help="Compatible best-reference TSV")
    parser.add_argument("--best-ref-txt", required=True, help="TXT containing the selected reference ID")
    parser.add_argument(
        "--alternate-subtype-txt",
        required=True,
        help="Pattern file used by downstream reference extraction",
    )
    parser.add_argument("--ranking-output", required=True, help="Full ranked reference-selection TSV")
    parser.add_argument(
        "--use-coverage-breadth",
        action="store_true",
        help="Use reference coverage breadth as an early tie-breaker",
    )
    return parser.parse_args()


def open_text(path):
    path = str(path)
    if path.endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8")
    return open(path, "r", encoding="utf-8")


def parse_fasta_lengths(path):
    lengths = OrderedDict()
    current_name = None
    current_length = 0

    with open_text(path) as handle:
        for line in handle:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith(">"):
                if current_name is not None:
                    lengths[current_name] = current_length
                current_name = line[1:].split()[0]
                current_length = 0
            elif current_name is not None:
                current_length += len(line.strip())

    if current_name is not None:
        lengths[current_name] = current_length

    return lengths


def empty_stats(reference_length):
    return {
        "reference_length": reference_length,
        "mapped_reads": 0,
        "mapped_bases": 0,
        "aligned_bases": 0,
        "mapq_sum": 0,
        "identity_weighted_sum": 0.0,
        "identity_weight_total": 0,
        "intervals": [],
    }


def parse_sq_header(line):
    name = None
    length = None
    for field in line.rstrip("\n").split("\t")[1:]:
        if field.startswith("SN:"):
            name = field[3:]
        elif field.startswith("LN:"):
            length = int(field[3:])
    return name, length


def consumes_query(op):
    return op in {"M", "I", "S", "=", "X"}


def aligned_query_bases(cigar):
    total = 0
    for length, op in CIGAR_RE.findall(cigar):
        if op in {"M", "I", "=", "X"}:
            total += int(length)
    return total


def inferred_query_length(cigar):
    total = 0
    for length, op in CIGAR_RE.findall(cigar):
        if consumes_query(op):
            total += int(length)
    return total


def reference_intervals(pos_1based, cigar):
    intervals = []
    ref_pos = pos_1based - 1
    for length_text, op in CIGAR_RE.findall(cigar):
        length = int(length_text)
        if op in {"M", "=", "X"}:
            intervals.append((ref_pos, ref_pos + length))
            ref_pos += length
        elif op in {"D", "N"}:
            ref_pos += length
    return intervals


def get_nm(optional_fields):
    for field in optional_fields:
        if field.startswith("NM:i:"):
            return int(field[5:])
    return None


def covered_bases(intervals):
    if not intervals:
        return 0

    merged_total = 0
    current_start = None
    current_end = None
    for start, end in sorted(intervals):
        if current_start is None:
            current_start = start
            current_end = end
        elif start <= current_end:
            current_end = max(current_end, end)
        else:
            merged_total += current_end - current_start
            current_start = start
            current_end = end

    if current_start is not None:
        merged_total += current_end - current_start

    return merged_total


def fmt_float(value):
    if value is None:
        return ""
    return f"{value:.6f}"


def extract_genotype_subtype(reference_id):
    parts = reference_id.split("_")
    if not parts or not parts[0]:
        raise ValueError(f"Unable to parse genotype/subtype from reference ID: {reference_id}")

    subtype = parts[0]
    genotype = "".join(char for char in subtype if char.isdigit())
    return genotype, subtype


def load_alignment_stats(primary_sam, panel_fasta):
    fasta_lengths = parse_fasta_lengths(panel_fasta)
    stats = OrderedDict((name, empty_stats(length)) for name, length in fasta_lengths.items())

    with open_text(primary_sam) as handle:
        for line in handle:
            if not line.strip():
                continue
            if line.startswith("@SQ"):
                name, length = parse_sq_header(line)
                if name is not None and name not in stats:
                    stats[name] = empty_stats(length or 0)
                elif name is not None and length is not None:
                    stats[name]["reference_length"] = length
                continue
            if line.startswith("@"):
                continue

            fields = line.rstrip("\n").split("\t")
            if len(fields) < 11:
                continue

            flag = int(fields[1])
            if flag & 0x904:
                continue

            reference_name = fields[2]
            if reference_name == "*":
                continue

            pos = int(fields[3])
            mapq = 0 if fields[4] == "*" else int(fields[4])
            cigar = fields[5]
            sequence = fields[9]
            optional_fields = fields[11:]

            if reference_name not in stats:
                stats[reference_name] = empty_stats(0)

            reference_stats = stats[reference_name]
            query_length = len(sequence) if sequence != "*" else inferred_query_length(cigar)
            aligned_bases = aligned_query_bases(cigar)

            reference_stats["mapped_reads"] += 1
            reference_stats["mapped_bases"] += query_length
            reference_stats["aligned_bases"] += aligned_bases
            reference_stats["mapq_sum"] += mapq
            reference_stats["intervals"].extend(reference_intervals(pos, cigar))

            nm = get_nm(optional_fields)
            if nm is not None and aligned_bases > 0:
                identity = max(0.0, 1.0 - (nm / aligned_bases))
                reference_stats["identity_weighted_sum"] += identity * aligned_bases
                reference_stats["identity_weight_total"] += aligned_bases

    return stats


def summarize_rows(args, stats):
    rows = []
    for reference_name, reference_stats in stats.items():
        mapped_reads = reference_stats["mapped_reads"]
        mapped_bases = reference_stats["mapped_bases"]
        aligned_bases = reference_stats["aligned_bases"]
        reference_length = reference_stats["reference_length"]
        covered_reference_bases = covered_bases(reference_stats["intervals"])

        aligned_fraction = (aligned_bases / mapped_bases) if mapped_bases else None
        mean_mapq = (reference_stats["mapq_sum"] / mapped_reads) if mapped_reads else None
        mean_identity = None
        if reference_stats["identity_weight_total"]:
            mean_identity = reference_stats["identity_weighted_sum"] / reference_stats["identity_weight_total"]

        reference_coverage_breadth = None
        if reference_length:
            reference_coverage_breadth = covered_reference_bases / reference_length

        rows.append(
            {
                "sample_id": args.sample_name,
                "read_type": args.read_type,
                "source_panel_fasta": str(Path(args.panel_fasta).resolve()),
                "competitive_bam": str(Path(args.bam).resolve()),
                "reference_id": reference_name,
                "reference_length": reference_length,
                "mapped_reads": mapped_reads,
                "mapped_bases": mapped_bases,
                "aligned_bases": aligned_bases,
                "aligned_fraction": aligned_fraction,
                "mean_identity": mean_identity,
                "mean_mapq": mean_mapq,
                "covered_reference_bases": covered_reference_bases,
                "reference_coverage_breadth": reference_coverage_breadth,
            }
        )

    return rows


def sort_rows(rows, use_coverage_breadth):
    def descending_float(value):
        return value if value is not None else -1.0

    if use_coverage_breadth:
        return sorted(
            rows,
            key=lambda row: (
                -row["aligned_bases"],
                -descending_float(row["reference_coverage_breadth"]),
                -row["mapped_bases"],
                -row["mapped_reads"],
                -descending_float(row["mean_identity"]),
                -descending_float(row["mean_mapq"]),
                row["reference_id"],
            ),
        )

    return sorted(
        rows,
        key=lambda row: (
            -row["aligned_bases"],
            -row["mapped_bases"],
            -row["mapped_reads"],
            -descending_float(row["mean_identity"]),
            -descending_float(row["mean_mapq"]),
            row["reference_id"],
        ),
    )


def write_ranking(path, sorted_rows, selection_rule):
    columns = [
        "sample_id",
        "read_type",
        "source_panel_fasta",
        "competitive_bam",
        "reference_id",
        "reference_length",
        "mapped_reads",
        "mapped_bases",
        "aligned_bases",
        "aligned_fraction",
        "mean_identity",
        "mean_mapq",
        "covered_reference_bases",
        "reference_coverage_breadth",
        "selection_rank",
        "selected",
        "selection_rule",
    ]

    with open(path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t")
        writer.writeheader()
        for rank, row in enumerate(sorted_rows, start=1):
            output_row = dict(row)
            output_row["aligned_fraction"] = fmt_float(row["aligned_fraction"])
            output_row["mean_identity"] = fmt_float(row["mean_identity"])
            output_row["mean_mapq"] = fmt_float(row["mean_mapq"])
            output_row["reference_coverage_breadth"] = fmt_float(row["reference_coverage_breadth"])
            output_row["selection_rank"] = rank
            output_row["selected"] = "true" if rank == 1 else "false"
            output_row["selection_rule"] = selection_rule
            writer.writerow(output_row)


def write_compatible_outputs(args, sorted_rows, selection_rule):
    winner = sorted_rows[0]
    best_genotype, best_subtype = extract_genotype_subtype(winner["reference_id"])
    close_hits = [
        row["reference_id"]
        for row in sorted_rows[1:]
        if row["aligned_bases"] == winner["aligned_bases"] and row["aligned_bases"] > 0
    ]
    other_subtypes = []
    for reference_id in close_hits:
        _genotype, subtype = extract_genotype_subtype(reference_id)
        if subtype != best_subtype:
            other_subtypes.append(reference_id)

    with open(args.best_ref_txt, "w", encoding="utf-8") as handle:
        handle.write(f"{winner['reference_id']}\n")

    with open(args.alternate_subtype_txt, "w", encoding="utf-8") as handle:
        handle.write(f"{winner['reference_id']}\n")

    columns = [
        "sample_id",
        "genotype",
        "subtype",
        "best_ref",
        "close_hits",
        "other_potential_subtypes",
        "selection_method",
        "selection_rule",
        "reference_length",
        "mapped_reads",
        "mapped_bases",
        "aligned_bases",
        "aligned_fraction",
        "mean_identity",
        "mean_mapq",
        "covered_reference_bases",
        "reference_coverage_breadth",
    ]

    with open(args.output, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t")
        writer.writeheader()
        writer.writerow(
            {
                "sample_id": args.sample_name,
                "genotype": best_genotype,
                "subtype": best_subtype,
                "best_ref": winner["reference_id"],
                "close_hits": ";".join(close_hits),
                "other_potential_subtypes": ";".join(dict.fromkeys(other_subtypes)),
                "selection_method": "minimap2_competitive_bam",
                "selection_rule": selection_rule,
                "reference_length": winner["reference_length"],
                "mapped_reads": winner["mapped_reads"],
                "mapped_bases": winner["mapped_bases"],
                "aligned_bases": winner["aligned_bases"],
                "aligned_fraction": fmt_float(winner["aligned_fraction"]),
                "mean_identity": fmt_float(winner["mean_identity"]),
                "mean_mapq": fmt_float(winner["mean_mapq"]),
                "covered_reference_bases": winner["covered_reference_bases"],
                "reference_coverage_breadth": fmt_float(winner["reference_coverage_breadth"]),
            }
        )


def main():
    args = parse_args()
    stats = load_alignment_stats(args.primary_sam, args.panel_fasta)
    rows = summarize_rows(args, stats)

    if not rows:
        raise SystemExit(f"No references were found in {args.primary_sam} or {args.panel_fasta}.")
    if not any(row["mapped_reads"] > 0 for row in rows):
        raise SystemExit(f"No primary mapped reads were found in {args.bam} for sample {args.sample_name}.")

    sorted_rows = sort_rows(rows, args.use_coverage_breadth)
    if args.use_coverage_breadth:
        selection_rule = (
            "aligned_bases desc; reference_coverage_breadth desc; mapped_bases desc; "
            "mapped_reads desc; mean_identity desc; mean_mapq desc; reference_id asc"
        )
    else:
        selection_rule = (
            "aligned_bases desc; mapped_bases desc; mapped_reads desc; "
            "mean_identity desc; mean_mapq desc; reference_id asc"
        )

    write_ranking(args.ranking_output, sorted_rows, selection_rule)
    write_compatible_outputs(args, sorted_rows, selection_rule)


if __name__ == "__main__":
    try:
        main()
    except ValueError as exc:
        print(str(exc), file=sys.stderr)
        sys.exit(1)
