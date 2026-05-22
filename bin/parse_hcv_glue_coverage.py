#!/usr/bin/env python3
import argparse
import csv
import re
from pathlib import Path


FEATURES = [
    "Polyprotein",
    "Core",
    "E1",
    "E2",
    "p7",
    "NS2",
    "NS3",
    "NS4A",
    "NS4B",
    "NS5A",
    "NS5B",
]
FEATURE_REMAP = {
    "Precursor polyprotein": "Polyprotein",
}


def parse_args():
    parser = argparse.ArgumentParser(description="Parse HCV-GLUE HTML feature coverage into a TSV.")
    parser.add_argument("--sample-id", required=True)
    parser.add_argument("--read-type", required=True, choices=["long", "short"])
    parser.add_argument("--input", required=True)
    parser.add_argument("--output", required=True)
    return parser.parse_args()


def extract_feature_table(html_text):
    match = re.search(r'<table id="featureCoverageTable">(.*?)</table>', html_text, re.IGNORECASE | re.DOTALL)
    if not match:
        return {}
    table_html = match.group(1)
    rows = re.findall(r"<tr>(.*?)</tr>", table_html, re.IGNORECASE | re.DOTALL)
    coverage = {}
    for row_html in rows:
        cells = re.findall(r"<t[dh][^>]*>(.*?)</t[dh]>", row_html, re.IGNORECASE | re.DOTALL)
        if len(cells) < 2:
            continue
        feature = re.sub(r"<.*?>", "", cells[0]).strip()
        value = re.sub(r"<.*?>", "", cells[1]).strip().rstrip("%")
        if not feature or not value:
            continue
        feature = FEATURE_REMAP.get(feature, feature)
        coverage[feature] = value
    return coverage


def main():
    args = parse_args()
    html_text = Path(args.input).read_text(errors="ignore")
    coverage = extract_feature_table(html_text)

    with open(args.output, "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["sample_id", "read_type", "feature", "coverage_pct"])
        for feature in FEATURES:
            writer.writerow([args.sample_id, args.read_type, feature, coverage.get(feature, "")])


if __name__ == "__main__":
    main()
