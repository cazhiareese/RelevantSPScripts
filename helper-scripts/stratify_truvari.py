#!/usr/bin/env python3
"""
stratify_truvari.py — Break Truvari bench results into per-size-bin, per-SV-type
precision / recall / F1.

Truvari writes, in each <graph>.<method>/ output directory:
    tp-base.vcf.gz   true positives, baseline (truth) representation
    fn.vcf.gz        false negatives  (truth not matched)        -> recall
    tp-comp.vcf.gz   true positives, comparison (callset) representation
    fp.vcf.gz        false positives  (callset not matched)      -> precision

Recall    per bin/type = TP_base / (TP_base + FN)
Precision per bin/type = TP_comp / (TP_comp + FP)
F1                      = 2*P*R / (P+R)

This scans BENCH_DIR for subdirectories containing those files and emits one
combined TSV with both per-(type x bin) rows and an ALL-types row per run.

Usage:
    python stratify_truvari.py --bench-dir graphs/truvari \\
        --out graphs/truvari/stratified_summary.tsv
"""

import argparse
import gzip
import sys
from collections import defaultdict
from pathlib import Path


# Size bins — six-bin A-F scheme, boundaries shared with the SV size-bin
# reporting (sv_accession_analysis.py). The top bin F is left OPEN-ended here
# (>=250 kb) instead of capping at 1 Mb, so no truth/callset variant is dropped
# from the precision/recall accounting.
# NOTE: keep in sync with generate_sv_truth.py / parse_syri_output.py.
BIN_ORDER = ["A_50-150bp", "B_151-500bp", "C_501bp-5kb",
             "D_5-50kb", "E_50-250kb", "F_0.25-1Mb"]
SV_TYPES = ["DEL", "INS", "INV", "DUP", "BND"]   # BND = intra-chr translocation


def size_bin_label(size: int) -> str:
    size = abs(size)
    if size < 151:
        return "A_50-150bp"
    elif size < 501:
        return "B_151-500bp"
    elif size < 5_001:
        return "C_501bp-5kb"
    elif size < 50_001:
        return "D_5-50kb"
    elif size < 250_001:
        return "E_50-250kb"
    return "F_0.25-1Mb"


def _info_field(info: str, key: str):
    """Return the value of an INFO key, or None."""
    for field in info.split(";"):
        if field.startswith(key + "="):
            return field[len(key) + 1:]
    return None


def variant_type_and_size(ref: str, alt: str, info: str):
    """
    Determine (sv_type, size) for a VCF record from INFO, falling back to
    REF/ALT lengths for sequence-resolved callset records.
    """
    svtype = _info_field(info, "SVTYPE")
    svlen = _info_field(info, "SVLEN")

    size = None
    if svlen is not None:
        try:
            size = abs(int(svlen.split(",")[0]))
        except ValueError:
            size = None

    if size is None:
        # sequence-resolved: size from REF/ALT length difference, or symbolic END
        if alt.startswith("<"):
            end = _info_field(info, "END")
            # symbolic without SVLEN — can't size reliably; skip
            size = None
        else:
            size = abs(len(alt) - len(ref))
            if size == 0:
                size = max(len(ref), len(alt))

    if svtype is None:
        # infer from sequence-resolved alleles
        if alt.startswith("<") and alt.endswith(">"):
            svtype = alt.strip("<>").split(":")[0]
        elif len(alt) > len(ref):
            svtype = "INS"
        elif len(alt) < len(ref):
            svtype = "DEL"
        else:
            svtype = "OTHER"

    # NOTE: DUP is kept as DUP (no fold to INS). Recall is read from the truth
    # side (tp-base/fn keep SVTYPE=DUP), so this gives a real per-type DUP recall.
    # --dup-to-ins still governs *matching* on the truvari side; this only affects
    # how rows are labelled in the report.
    return svtype, size


def count_by_bin_type(vcf_path: Path):
    """Return {(sv_type, size_bin): count} for a (possibly gzipped) VCF."""
    counts = defaultdict(int)
    if not vcf_path.exists():
        return counts
    opener = gzip.open if str(vcf_path).endswith(".gz") else open
    with opener(vcf_path, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 8:
                continue
            ref, alt, info = cols[3], cols[4], cols[7]
            svtype, size = variant_type_and_size(ref, alt, info)
            if size is None or size < 50:
                continue
            counts[(svtype, size_bin_label(size))] += 1
    return counts


def f1(p, r):
    return (2 * p * r / (p + r)) if (p + r) > 0 else 0.0


def process_run(run_dir: Path):
    """Yield metric rows for one <graph>.<method> truvari directory."""
    tp_base = count_by_bin_type(run_dir / "tp-base.vcf.gz")
    fn      = count_by_bin_type(run_dir / "fn.vcf.gz")
    tp_comp = count_by_bin_type(run_dir / "tp-comp.vcf.gz")
    fp      = count_by_bin_type(run_dir / "fp.vcf.gz")

    name = run_dir.name
    if "." in name:
        graph, method = name.split(".", 1)
    else:
        graph, method = name, ""

    keys = set(tp_base) | set(fn) | set(tp_comp) | set(fp)
    # roll-up accumulators for the ALL row
    agg = defaultdict(lambda: [0, 0, 0, 0])  # bin -> [tpb, fn, tpc, fp]

    # Iterate EVERY SV type present (not just DEL/INS/INV/DUP) so the FP/TP
    # totals reconcile with truvari's summary.json. Graph callers emit
    # inversions as large substitutions and other complex events that land in
    # an OTHER bucket — dropping them would understate false positives.
    present_types = {t for (t, _b) in keys}
    type_order = [t for t in SV_TYPES if t in present_types] + \
                 sorted(present_types - set(SV_TYPES))

    rows = []
    for sv_type in type_order:
        # per-type roll-up (this type, summed across all size bins)
        type_agg = [0, 0, 0, 0]  # [tpb, fn, tpc, fp]
        for size_bin in BIN_ORDER:
            k = (sv_type, size_bin)
            if k not in keys:
                continue
            tpb = tp_base.get(k, 0)
            fnn = fn.get(k, 0)
            tpc = tp_comp.get(k, 0)
            fpp = fp.get(k, 0)
            recall    = tpb / (tpb + fnn) if (tpb + fnn) else 0.0
            precision = tpc / (tpc + fpp) if (tpc + fpp) else 0.0
            rows.append([graph, method, sv_type, size_bin,
                         tpb, fnn, tpc, fpp,
                         f"{recall:.4f}", f"{precision:.4f}", f"{f1(precision, recall):.4f}"])
            a = agg[size_bin]
            a[0] += tpb; a[1] += fnn; a[2] += tpc; a[3] += fpp
            type_agg[0] += tpb; type_agg[1] += fnn; type_agg[2] += tpc; type_agg[3] += fpp

        # this-type ALL-bins summary row
        tpb, fnn, tpc, fpp = type_agg
        recall    = tpb / (tpb + fnn) if (tpb + fnn) else 0.0
        precision = tpc / (tpc + fpp) if (tpc + fpp) else 0.0
        rows.append([graph, method, sv_type, "ALL",
                     tpb, fnn, tpc, fpp,
                     f"{recall:.4f}", f"{precision:.4f}", f"{f1(precision, recall):.4f}"])

    # ALL-types per bin
    for size_bin in BIN_ORDER:
        if size_bin not in agg:
            continue
        tpb, fnn, tpc, fpp = agg[size_bin]
        recall    = tpb / (tpb + fnn) if (tpb + fnn) else 0.0
        precision = tpc / (tpc + fpp) if (tpc + fpp) else 0.0
        rows.append([graph, method, "ALL", size_bin,
                     tpb, fnn, tpc, fpp,
                     f"{recall:.4f}", f"{precision:.4f}", f"{f1(precision, recall):.4f}"])

    # ALL-types ALL-bins
    tpb = sum(a[0] for a in agg.values()); fnn = sum(a[1] for a in agg.values())
    tpc = sum(a[2] for a in agg.values()); fpp = sum(a[3] for a in agg.values())
    recall    = tpb / (tpb + fnn) if (tpb + fnn) else 0.0
    precision = tpc / (tpc + fpp) if (tpc + fpp) else 0.0
    rows.append([graph, method, "ALL", "ALL",
                 tpb, fnn, tpc, fpp,
                 f"{recall:.4f}", f"{precision:.4f}", f"{f1(precision, recall):.4f}"])
    return rows


COLUMNS = ["graph", "method", "sv_type", "size_bin",
           "TP_base", "FN", "TP_comp", "FP",
           "recall", "precision", "f1"]


def main():
    ap = argparse.ArgumentParser(description="Stratify Truvari results by size bin and SV type.")
    ap.add_argument("--bench-dir", required=True,
                    help="Directory containing <graph>.<method>/ truvari output subdirs.")
    ap.add_argument("--out", required=True, help="Output TSV path.")
    args = ap.parse_args()

    bench_dir = Path(args.bench_dir)
    # Accept either a single truvari run dir (fn.vcf.gz directly inside) or a
    # parent dir holding several <graph>.<method>/ run subdirectories.
    if (bench_dir / "fn.vcf.gz").exists():
        run_dirs = [bench_dir]
    else:
        run_dirs = sorted(d for d in bench_dir.iterdir()
                          if d.is_dir() and (d / "fn.vcf.gz").exists())

    if not run_dirs:
        print(f"[ERROR] No truvari run dirs (with fn.vcf.gz) found in {bench_dir}. "
              f"Run 'make -f Makefile.benchmark bench-<graph>' first.", file=sys.stderr)
        sys.exit(1)

    all_rows = []
    for d in run_dirs:
        print(f"  processing {d.name}", file=sys.stderr)
        all_rows.extend(process_run(d))

    with open(args.out, "w") as f:
        f.write("\t".join(COLUMNS) + "\n")
        for row in all_rows:
            f.write("\t".join(str(x) for x in row) + "\n")

    # Echo the ALL/ALL lines so the headline numbers are visible immediately.
    print("\nHeadline (sv_type=ALL, size_bin=ALL):", file=sys.stderr)
    for row in all_rows:
        if row[2] == "ALL" and row[3] == "ALL":
            print(f"  {row[0]:<8} {row[1]:<8} "
                  f"recall={row[8]} precision={row[9]} f1={row[10]}", file=sys.stderr)


if __name__ == "__main__":
    main()
