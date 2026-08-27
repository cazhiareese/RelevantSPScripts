#!/usr/bin/env python3
"""
swave_gene_impact.py — chr1 pggb-swave typed calls → gene-impact annotation.

- Decomposes compound swave SVs (SVTYPE like DEL+DEL+invDUP) into typed
  sub-events via the BKPS INFO field, each with reference coordinates.
- Overlaps each sub-event with IRGSP-1.0 gene/exon/CDS/promoter tracks.
- Assigns a gene-impact category (deleted / truncated / duplicated /
  inverted / inversion_breakpoint / insertion_in_gene / promoter).
- Flags each parent call `validated` if its footprint overlaps a truvari
  tp-comp call (within 500 bp) for that sample.

chr1 only. Per-sample (a call is kept for a cultivar only if that sample
carries the ALT). Sub-events < 50 bp are dropped (benchmark sizemin).
"""
import re, sys, gzip
from pathlib import Path

ANN = Path("/Users/cazhia/Desktop/SP-Scripts/annotation")
VCF = Path("/Users/cazhia/Desktop/SP-Scripts/benchmarking/real/pggb-swave/swave.sample_level.split.vcf")
TPCOMP = {
    "Azucena": "/Users/cazhia/Desktop/SP-Scripts/benchmarking/real/pggb-swave/bench_out/Azucena/truvari/tp-comp.vcf.gz",
    "IR64":    "/Users/cazhia/Desktop/SP-Scripts/benchmarking/real/pggb-swave/bench_out/IR64/truvari/tp-comp.vcf.gz",
}
VALID_TOL = 500
MIN_SIZE = 50

# ---- load feature tracks (chr '1'; BED 0-based start) -------------------
def load_named(bed):
    out = []
    for l in open(ANN / bed):
        c = l.rstrip("\n").split("\t")
        if c[0] != "1":
            continue
        out.append((int(c[1]), int(c[2]), c[3]))  # start0, end, id|name
    return out

def load_plain(bed):
    out = []
    for l in open(ANN / bed):
        c = l.rstrip("\n").split("\t")
        if c[0] != "1":
            continue
        out.append((int(c[1]), int(c[2])))
    return out

genes    = load_named("genes.bed")
promoters= load_named("promoter.bed")
exons    = load_plain("exon.bed")
cds      = load_plain("cds.bed")

def overlaps_any(s0, e, intervals):
    return any(s0 < ie and e > is0 for is0, ie in intervals)

def genes_hit(s0, e, named):
    return [(is0, ie, nm) for is0, ie, nm in named if s0 < ie and e > is0]

# ---- tp-comp footprints per sample (for validated flag) ----------------
def load_tpcomp(path):
    ivs = []
    with gzip.open(path, "rt") as f:
        for l in f:
            if l.startswith("#"):
                continue
            c = l.split("\t")
            chrom = c[0].replace("chr", "")
            if chrom != "1":
                continue
            pos = int(c[1]); ref = c[3]
            s0 = pos - 1 - VALID_TOL
            e = pos - 1 + max(1, len(ref)) + VALID_TOL
            ivs.append((s0, e))
    return ivs

tpcomp = {s: load_tpcomp(p) for s, p in TPCOMP.items()}

# ---- BKPS sub-event parser ---------------------------------------------
# e.g. DEL_2162_Nipponbare#1#chr1_697025_699186
#      invDUP_1411_Nipponbare#1#chr1_708725_710135_Nipponbare#1#chr1_702996_702997
#      hyperCPX_DUP_50191_Nipponbare#1#chr1_752472_802662
BK = re.compile(
    r"^(?P<type>[A-Za-z]+(?:_[A-Za-z]+)*)_(?P<size>\d+)_"
    r"(?P<chrom>[^_]+)_(?P<start>\d+)_(?P<end>\d+)"
    r"(?:_(?P<c2>[^_]+)_(?P<s2>\d+)_(?P<e2>\d+))?$"
)

def coarse(t):
    tl = t.lower()
    if "dup" in tl:  return "DUP"
    if "inv" in tl:  return "INV"
    if "del" in tl:  return "DEL"
    if "ins" in tl:  return "INS"
    return "COMPLEX"

def impact(coarse_t, fclass, s0, e, ghits):
    """Return list of (gene, gene_impact), feature-aware, over the gene hits."""
    res = []
    for gs0, ge, nm in ghits:
        full = gs0 >= s0 and ge <= e
        coding = fclass in ("CDS(coding)", "exon/UTR")
        if coarse_t == "DEL":
            gi = "deleted_gene" if full else ("truncated_gene" if fclass == "CDS(coding)"
                 else "utr_deletion" if fclass == "exon/UTR"
                 else "promoter_deletion" if fclass == "promoter" else "intronic_deletion")
        elif coarse_t == "DUP":
            gi = ("whole_gene_duplication" if full else "duplicated_gene") if fclass != "promoter" else "promoter_duplication"
        elif coarse_t == "INV":
            gi = "inverted_gene" if full else ("inversion_breakpoint" if coding
                 else "promoter_inversion" if fclass == "promoter" else "intronic_inversion")
        elif coarse_t == "INS":
            gi = "coding_insertion" if coding else ("promoter_insertion" if fclass == "promoter" else "intronic_insertion")
        else:
            gi = "complex_" + (fclass.split("(")[0] if fclass else "genic")
        res.append((nm, gi))
    return res

def feature_class(s0, e, ghits):
    if overlaps_any(s0, e, cds):   return "CDS(coding)"
    if overlaps_any(s0, e, exons): return "exon/UTR"
    if ghits:                      return "intron"
    return None  # caller checks promoter next

# ---- parse VCF ----------------------------------------------------------
def info_get(info, key):
    m = re.search(rf"(?:^|;){key}=([^;]+)", info)
    return m.group(1) if m else None

samples = None
rows = []
for l in open(VCF):
    if l.startswith("##"):
        continue
    if l.startswith("#CHROM"):
        cols = l.rstrip("\n").split("\t")
        samples = cols[9:]  # Azucena, IR64, Nipponbare
        idx = {s: 9 + i for i, s in enumerate(samples)}
        continue
    c = l.rstrip("\n").split("\t")
    chrom = c[0].replace("Nipponbare#1#chr", "").replace("chr", "")
    if chrom != "1":
        continue
    pos = int(c[1]); info = c[7]
    end = info_get(info, "END"); bkps = info_get(info, "BKPS")
    parent_type = info_get(info, "SVTYPE") or "."
    endi = int(end) if end else pos
    for sample in ("Azucena", "IR64"):
        gt = c[idx[sample]].split(":")[0]
        if "1" not in gt:          # sample must carry the ALT
            continue
        validated = overlaps_any(pos - 1, endi, tpcomp[sample])
        # decompose BKPS into sub-events (fallback: whole record)
        subs = []
        if bkps:
            for tok in bkps.split(","):
                m = BK.match(tok)
                if not m:
                    continue
                st, en = int(m["start"]), int(m["end"])
                subs.append((m["type"], st, en))
        if not subs:
            subs = [(parent_type, pos, endi)]
        for stype, st, en in subs:
            if en - st < MIN_SIZE:
                continue
            s0 = st - 1
            ct = coarse(stype)
            ghits = genes_hit(s0, en, genes)
            fclass = feature_class(s0, en, ghits)
            if fclass is None:
                phits = genes_hit(s0, en, promoters)
                if phits:
                    fclass = "promoter"; ghits_used = phits
                else:
                    continue  # intergenic — skip
            else:
                ghits_used = ghits
            imp = impact(ct, fclass, s0, en, ghits_used)
            for (idnm, gi) in (imp if imp else [(".", ".")]):
                rows.append([sample, "1", st, en, stype, ct, parent_type, pos,
                             fclass, gi, idnm, "yes" if validated else "no"])

# ---- write ----------------------------------------------------------------
hdr = ["sample","chrom","sub_start","sub_end","sub_type","coarse_type",
       "parent_svtype","parent_pos","feature_class","gene_impact","gene","validated"]
rows.sort(key=lambda x: (x[0], x[2]))
out_all = ANN / "swave_chr1_gene_impact.tsv"          # all chr1 candidates
out_val = ANN / "swave_chr1_gene_impact.validated.tsv"  # validated subset
with open(out_all, "w") as w:
    w.write("\t".join(hdr) + "\n")
    for r in rows:
        w.write("\t".join(map(str, r)) + "\n")
with open(out_val, "w") as w:
    w.write("\t".join(hdr) + "\n")
    for r in rows:
        if r[11] == "yes":
            w.write("\t".join(map(str, r)) + "\n")

from collections import Counter
val = [r for r in rows if r[11] == "yes"]
print(f"wrote {out_all} ({len(rows)} rows)  and  {out_val} ({len(val)} validated rows)\n")
for label, data in (("VALIDATED (tp-comp confirmed)", val), ("ALL chr1 candidates", rows)):
    print(f"########## {label} ##########")
    for sample in ("Azucena", "IR64"):
        sub = [r for r in data if r[0] == sample]
        genes_n = len({r[10] for r in sub})
        print(f"  {sample}: {len(sub)} rows, {genes_n} genes | impact={dict(Counter(r[9] for r in sub))}")
    print()
