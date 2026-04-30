#!/usr/bin/env python3
import argparse
import gzip


def parse_args():
    parser = argparse.ArgumentParser(description="Create a BED mask for low-depth positions from mosdepth per-base output.")
    parser.add_argument("--input", required=True)
    parser.add_argument("--min-depth", required=True, type=float)
    parser.add_argument("--output", required=True)
    return parser.parse_args()


def main():
    args = parse_args()
    current = None

    with gzip.open(args.input, "rt") as handle, open(args.output, "w") as out_handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            chrom, start, end, depth = line.split("\t")[:4]
            start = int(start)
            end = int(end)
            depth = float(depth)

            if depth < args.min_depth:
                if current and current[0] == chrom and current[2] == start:
                    current[2] = end
                else:
                    if current:
                        out_handle.write(f"{current[0]}\t{current[1]}\t{current[2]}\n")
                    current = [chrom, start, end]
            elif current:
                out_handle.write(f"{current[0]}\t{current[1]}\t{current[2]}\n")
                current = None

        if current:
            out_handle.write(f"{current[0]}\t{current[1]}\t{current[2]}\n")


if __name__ == "__main__":
    main()
