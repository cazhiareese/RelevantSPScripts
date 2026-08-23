#!/usr/bin/env python3
"""
sample_svs_for_validation.py — Randomly sample SVs from an Assemblytics BED
for manual dotplot validation.

For each sampled SV the script:
  1. Writes a sampled_svs.tsv with coordinates and size-bin labels
  2. Writes a run_dotplots.sh that extracts flanking regions from both
     assemblies, re-runs nucmer on each region, and generates a mummerplot PNG

Requires mummer_env (nucmer, mummerplot, samtools) for run_dotplots.sh.

Usage:
    python simulation/sample_svs_for_validation.py \\
        --bed    bio_align/ir64_asm.Assemblytics_structural_variants.bed \\
        --ref-fasta   bio_align/Nipponbare_chr1.fasta \\
        --query-fasta bio_align/IR64_chr1.fasta \\
        --n 5 \\
        --outdir bio_align/validation_samples/ir64/
"""

import argparse
import random
import sys
from collections import defaultdict
from pathlib import Path


# =============================================================================
# SV TYPE REMAPPING
# =============================================================================

# Assemblytics uses 7 internal type names. Map them to the 4 canonical types
# used in the pangenome benchmark (DEL, INS, INV, DUP) so sampling is
# stratified by the types that actually matter for Truvari benchmarking.
ASSEMBLYTICS_TYPE_MAP = {
    "Deletion":           "DEL",
    "Tandem_contraction": "DEL",
    "Repeat_contraction": "DEL",
    "Insertion":          "INS",
    "Inversion":          "INV",
    "Tandem_expansion":   "DUP",
    "Repeat_expansion":   "DUP",
}


# =============================================================================
# ASSEMBLYTICS BED COLUMN INDICES (0-based)
# =============================================================================

# Standard Assemblytics_structural_variants.bed columns:
# ref_chr  ref_start  ref_end  name  size  strand  sv_type
# ref_gap_size  query_gap_size  query_coordinates  method

COL_REF_CHR    = 0
COL_REF_START  = 1
COL_REF_END    = 2
COL_NAME       = 3
COL_SIZE       = 4
COL_STRAND     = 5
COL_TYPE       = 6
COL_REF_GAP    = 7
COL_QRY_GAP    = 8
COL_QRY_COORDS = 9
COL_METHOD     = 10


# =============================================================================
# SIZE BIN LABELLING (matches generate_sv_truth.py and parse_sv_responses.py)
# =============================================================================

def size_bin_label(size: int) -> str:
    size = abs(size)
    if size < 500:
        return "A_50-500bp"
    elif size < 5_000:
        return "B_500bp-5kb"
    elif size < 50_000:
        return "C_5-50kb"
    elif size <= 100_000:
        return "D_50-100kb"
    else:
        return "E_>100kb"


# =============================================================================
# PARSE ASSEMBLYTICS BED
# =============================================================================

def parse_assemblytics_bed(path: str, min_size: int = 50) -> list[dict]:
    """
    Parse an Assemblytics_structural_variants.bed file.
    Skips header lines (starting with #) and SVs smaller than min_size.
    """
    svs = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 7:
                continue
            try:
                size = abs(int(parts[COL_SIZE]))
            except ValueError:
                continue
            if size < min_size:
                continue
            raw_type = parts[COL_TYPE]
            sv_type  = ASSEMBLYTICS_TYPE_MAP.get(raw_type, raw_type)
            svs.append({
                "ref_chr":      parts[COL_REF_CHR],
                "ref_start":    int(parts[COL_REF_START]),
                "ref_end":      int(parts[COL_REF_END]),
                "name":         parts[COL_NAME],
                "size":         size,
                "strand":       parts[COL_STRAND] if len(parts) > COL_STRAND else ".",
                "sv_type":      sv_type,
                "assemblytics_type": raw_type,
                "query_coords": parts[COL_QRY_COORDS] if len(parts) > COL_QRY_COORDS else ".",
                "method":       parts[COL_METHOD]     if len(parts) > COL_METHOD     else ".",
                "size_bin":     size_bin_label(size),
            })
    return svs


# =============================================================================
# RANDOM SAMPLING
# =============================================================================

def sample_svs(svs: list[dict], n: int, seed: int) -> list[dict]:
    """
    Sample up to n SVs per SV type. Types with fewer than n entries
    contribute all their SVs. Sampling is stratified so every type is
    represented, regardless of how common it is.
    """
    rng = random.Random(seed)
    by_type: dict[str, list] = defaultdict(list)
    for sv in svs:
        by_type[sv["sv_type"]].append(sv)

    sampled = []
    for sv_type, group in sorted(by_type.items()):
        k = min(n, len(group))
        chosen = rng.sample(group, k)
        sampled.extend(chosen)
        print(f"  {sv_type:<30} sampled {k:>3} / {len(group)}", file=sys.stderr)

    return sampled


# =============================================================================
# QUERY COORDINATE PARSER
# =============================================================================

def parse_query_coords(qry_str: str) -> tuple[str, int, int] | None:
    """
    Parse the query_coordinates field from Assemblytics BED.
    Handles formats:  chr:start-end   and   chr:start+size
    Returns (chrom, start, end) or None on failure.
    """
    if not qry_str or qry_str in (".", ""):
        return None
    try:
        chrom, coords = qry_str.rsplit(":", 1)
        if "-" in coords:
            start, end = coords.split("-", 1)
            return chrom, int(start), int(end)
        elif "+" in coords:
            start, size = coords.split("+", 1)
            start = int(start)
            return chrom, start, start + int(size)
    except (ValueError, AttributeError):
        pass
    return None


# =============================================================================
# OUTPUT: SAMPLED SVs TSV
# =============================================================================

TSV_COLUMNS = [
    "sv_type", "assemblytics_type", "size_bin", "name",
    "ref_chr", "ref_start", "ref_end", "size",
    "query_coords", "method",
]


def write_sampled_tsv(svs: list[dict], path: Path):
    with open(path, "w") as f:
        f.write("\t".join(TSV_COLUMNS) + "\n")
        for sv in svs:
            f.write("\t".join(str(sv.get(c, "")) for c in TSV_COLUMNS) + "\n")


# =============================================================================
# OUTPUT: DOTPLOT SHELL SCRIPT
# =============================================================================

def write_dotplot_script(
    svs: list[dict],
    ref_fasta: Path,
    query_fasta: Path,
    flank: int,
    outdir: Path,
    path: Path,
):
    """
    Write a self-contained bash script that generates one mummerplot PNG
    per sampled SV.

    Strategy per SV:
      1. Extract ref_region (SV ± flank) from the reference FASTA
      2. Extract qry_region (query_coords ± flank) from the query FASTA
      3. Run nucmer on those two small sequences
      4. Run mummerplot --png on the resulting delta
    The output PNG shows a local dotplot centred on the SV, making it easy
    to visually confirm or reject the call.
    """
    with open(path, "w") as f:
        f.write("#!/usr/bin/env bash\n")
        f.write("# Auto-generated by sample_svs_for_validation.py\n")
        f.write("# Requires: conda environment mummer_env\n")
        f.write("#   (nucmer, mummerplot, samtools must be in that env)\n\n")
        f.write("set -euo pipefail\n\n")
        f.write(f"REF=\"{ref_fasta}\"\n")
        f.write(f"QUERY=\"{query_fasta}\"\n")
        f.write(f"OUTDIR=\"{outdir}\"\n")
        f.write(f'CRUN="conda run -n mummer_env"\n\n')
        f.write('mkdir -p "$OUTDIR"\n\n')

        for sv in svs:
            # Sanitise name for use in filenames
            safe_name = sv["name"].replace(":", "_").replace("|", "_").replace("/", "_")
            sv_type   = sv["sv_type"]
            ref_chr   = sv["ref_chr"]

            # Reference region (SV ± flank); clamp start to 0
            ref_s = max(0, sv["ref_start"] - flank)
            ref_e = sv["ref_end"] + flank
            ref_region = f"{ref_chr}:{ref_s}-{ref_e}"

            # Query region (query_coords ± flank)
            qry = parse_query_coords(sv["query_coords"])
            if qry is None:
                # Assemblytics did not record query coordinates for this SV
                # (common for Between_alignments events at contig boundaries).
                # Cannot build a valid dotplot without knowing the query locus —
                # skip rather than querying the wrong FASTA with ref coordinates.
                f.write(
                    f"# SKIPPED  {sv_type}  {ref_chr}:{sv['ref_start']}-{sv['ref_end']}"
                    f"  — no query_coords in BED (Between_alignments boundary?)\n\n"
                )
                continue

            qry_chr, qry_s, qry_e = qry
            qry_s = max(0, qry_s - flank)
            qry_e = qry_e + flank
            qry_region = f"{qry_chr}:{qry_s}-{qry_e}"

            tmp_ref   = f"/tmp/{safe_name}_ref.fa"
            tmp_qry   = f"/tmp/{safe_name}_qry.fa"
            tmp_delta = f"/tmp/{safe_name}"
            plot_out  = f'"$OUTDIR"/{safe_name}_{sv_type}'

            f.write(
                f"# {sv_type}  {ref_chr}:{sv['ref_start']}-{sv['ref_end']}"
                f"  size={sv['size']} bp  ({sv['size_bin']})\n"
            )
            f.write(f'$CRUN samtools faidx "$REF"   "{ref_region}" > {tmp_ref}\n')
            f.write(f'$CRUN samtools faidx "$QUERY" "{qry_region}"  > {tmp_qry}\n')
            f.write(f"$CRUN nucmer --prefix {tmp_delta} {tmp_ref} {tmp_qry}\n")
            f.write(
                f'$CRUN mummerplot --png --large {tmp_delta}.delta'
                f' -p {plot_out}\n'
            )
            f.write(f'echo "✔  {safe_name}_{sv_type}.png"\n\n')

    path.chmod(0o755)


# =============================================================================
# CLI AND MAIN
# =============================================================================

def parse_args():
    p = argparse.ArgumentParser(
        description="Randomly sample Assemblytics SVs and generate mummerplot validation scripts."
    )
    p.add_argument("--bed", required=True,
                   help="Assemblytics_structural_variants.bed file.")
    p.add_argument("--ref-fasta", required=True,
                   help="Reference FASTA (Nipponbare_chr1.fasta).")
    p.add_argument("--query-fasta", required=True,
                   help="Query FASTA (IR64_chr1.fasta or Azucena_chr1.fasta).")
    p.add_argument("--n", type=int, default=5,
                   help="SVs to sample per SV type (default: 5).")
    p.add_argument("--min-size", type=int, default=50,
                   help="Minimum SV size in bp (default: 50).")
    p.add_argument("--flank", type=int, default=20_000,
                   help="Flanking bp around each SV for dotplot context "
                        "(default: 20000). Larger flank = more alignment "
                        "context but slower nucmer.")
    p.add_argument("--seed", type=int, default=42,
                   help="Random seed for reproducibility (default: 42).")
    p.add_argument("--outdir", required=True,
                   help="Output directory for sampled_svs.tsv and run_dotplots.sh.")
    return p.parse_args()


def main():
    args = parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    print(f"Parsing {args.bed} ...", file=sys.stderr)
    svs = parse_assemblytics_bed(args.bed, min_size=args.min_size)
    print(f"  {len(svs)} SVs ≥ {args.min_size} bp loaded", file=sys.stderr)

    if not svs:
        print("[ERROR] No SVs passed the size filter. "
              "Check --bed and --min-size.", file=sys.stderr)
        sys.exit(1)

    print(f"\nSampling {args.n} per SV type (seed={args.seed}):", file=sys.stderr)
    sampled = sample_svs(svs, args.n, args.seed)
    print(f"\n  Total sampled: {len(sampled)}", file=sys.stderr)

    tsv_path = outdir / "sampled_svs.tsv"
    write_sampled_tsv(sampled, tsv_path)
    print(f"\n✔ Sampled SVs   → {tsv_path}")

    sh_path = outdir / "run_dotplots.sh"
    write_dotplot_script(
        sampled,
        ref_fasta=Path(args.ref_fasta).resolve(),
        query_fasta=Path(args.query_fasta).resolve(),
        flank=args.flank,
        outdir=outdir.resolve(),
        path=sh_path,
    )
    print(f"✔ Dotplot script → {sh_path}")
    print(f"\nNext: bash {sh_path}")


if __name__ == "__main__":
    main()
