#!/usr/bin/env python3
"""
make_ground_truth_vcf.py — Build a ground-truth (silver-standard) VCF from the
manually dotplot-validated SVs.

Input:
  - a validation CSV (Azucena.csv / IR64.csv) whose cells contain validated SV
    IDs of the form TYPE_id (e.g. DEL_110, DUP_212, TRA_504, INS_58, INV_155);
    all columns are scanned (bins A–D, Complex, Round 2, Round 3).
  - the SyRI output directory + prefix, so each TYPE_id can be mapped back to its
    coordinates (the id is the global index assigned in parse_syri_output).

Output:
  - a VCF (Nipponbare reference coordinates) containing ONLY the validated SVs:
      DEL/INS/INV/DUP -> symbolic <DEL>/<INS>/<INV>/<DUP> with SVLEN/END/SIZEBIN
      TRA/invTRA      -> SVTYPE=BND at the source breakpoint (MATEPOS = query dest)

Usage:
    python make_ground_truth_vcf.py \\
        --syri-dir bio_align/syri_azucena --prefix azucena \\
        --csv bio_align/Azucena.csv --sample Azucena \\
        --out bio_align/ground_truth/azucena_ground_truth.vcf
"""

import argparse
import csv
import re
import sys
from datetime import datetime
from pathlib import Path

# reuse the SyRI parser (same directory)
from parse_syri_output import load_syri_svs

_ID_RE = re.compile(r"^[A-Za-z]+_\d+$")


def read_validated_ids(csv_path: Path) -> set:
    """Scan every cell; return the set of (sv_type, id) validated in the CSV."""
    ids = set()
    with open(csv_path, newline="") as f:
        for row in csv.reader(f):
            for cell in row:
                c = (cell or "").strip()
                if _ID_RE.match(c):
                    t, n = c.rsplit("_", 1)
                    ids.add((t, int(n)))
    return ids


VCF_HEADER = """\
##fileformat=VCFv4.2
##fileDate={date}
##source=make_ground_truth_vcf.py (dotplot-validated SyRI SVs)
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">
##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Length of the SV">
##INFO=<ID=END,Number=1,Type=Integer,Description="End position">
##INFO=<ID=SIZEBIN,Number=1,Type=String,Description="Size bin category">
##INFO=<ID=MATEPOS,Number=1,Type=Integer,Description="Translocation destination in the query assembly">
##INFO=<ID=SYRISV,Number=1,Type=String,Description="Original SyRI SV type label">
##ALT=<ID=DEL,Description="Deletion">
##ALT=<ID=INS,Description="Insertion">
##ALT=<ID=INV,Description="Inversion">
##ALT=<ID=DUP,Description="Duplication">
##ALT=<ID=BND,Description="Breakend (translocation)">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample}
"""

_SYMBOLIC = {"DEL": "<DEL>", "INS": "<INS>", "INV": "<INV>", "DUP": "<DUP>"}
_SVLEN_SIGN = {"DEL": -1, "INS": 1, "INV": 1, "DUP": 1}


def write_ground_truth(svs, sample, path: Path):
    path.parent.mkdir(parents=True, exist_ok=True)
    n = 0
    with open(path, "w") as f:
        f.write(VCF_HEADER.format(date=datetime.now().strftime("%Y%m%d"),
                                  sample=sample.upper()))
        for sv in sorted(svs, key=lambda x: x["ref_start"]):
            svtype = sv["sv_type"]
            chrom  = sv["ref_chr"]
            pos    = sv["ref_start"]
            sv_id  = f"{svtype}_{sv['id']}"

            if svtype in ("TRA", "invTRA"):
                # translocation -> BND at the source breakpoint (Nipponbare frame);
                # the destination locus in the query is recorded as MATEPOS.
                mate = sv.get("qry_start", pos)
                alt  = f"N[{chrom}:{sv['ref_end']}["
                info = (f"SVTYPE=BND;SVLEN={sv['size']};END={sv['ref_end']};"
                        f"MATEPOS={mate};SIZEBIN={sv['size_bin']};SYRISV={svtype}")
            else:
                sign  = _SVLEN_SIGN.get(svtype, 1)
                svlen = sign * sv["size"]
                end   = sv["ref_end"] if sv["ref_end"] != pos else pos + sv["size"]
                alt   = _SYMBOLIC.get(svtype, "<BND>")
                info  = (f"SVTYPE={svtype};SVLEN={svlen};END={end};"
                         f"SIZEBIN={sv['size_bin']};SYRISV={svtype}")

            f.write(f"{chrom}\t{pos}\t{sv_id}\tN\t{alt}\t.\tPASS\t"
                    f"{info}\tGT\t1/1\n")
            n += 1
    return n


def main():
    ap = argparse.ArgumentParser(description="Build ground-truth VCF from validated SVs.")
    ap.add_argument("--syri-dir", required=True)
    ap.add_argument("--prefix",   required=True, help="SyRI prefix (e.g. azucena).")
    ap.add_argument("--csv",      required=True, help="Validation CSV with TYPE_id cells.")
    ap.add_argument("--sample",   required=True, help="Accession name (Azucena / IR64).")
    ap.add_argument("--min-size", type=int, default=50)
    ap.add_argument("--out",      required=True)
    args = ap.parse_args()

    validated = read_validated_ids(Path(args.csv))
    print(f"Validated IDs in {args.csv}: {len(validated)}", file=sys.stderr)

    all_svs = load_syri_svs(Path(args.syri_dir), args.prefix, args.min_size)
    by_key = {(s["sv_type"], s["id"]): s for s in all_svs}

    keep, missing = [], []
    for key in validated:
        if key in by_key:
            keep.append(by_key[key])
        else:
            missing.append(key)

    from collections import Counter
    kept_counts = Counter(s["sv_type"] for s in keep)
    print(f"Matched {len(keep)} / {len(validated)} validated SVs to SyRI coords.",
          file=sys.stderr)
    print("  by type: " + "  ".join(f"{t}:{kept_counts[t]}"
          for t in sorted(kept_counts)), file=sys.stderr)
    if missing:
        print(f"[WARN] {len(missing)} validated IDs had no SyRI match "
              f"(check prefix / SyRI files): {sorted(missing)[:10]}...", file=sys.stderr)

    n = write_ground_truth(keep, args.sample, Path(args.out))
    print(f"✔ Wrote {n} SVs -> {args.out}", file=sys.stderr)


if __name__ == "__main__":
    main()
