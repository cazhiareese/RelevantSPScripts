#!/usr/bin/env python3
"""
03_annotate_pggb-swave_wholegenome_gene_impact.py

Whole-genome (chr1-12) pggb-swave typed calls -> gene-impact annotation.

Same logic as the chr1 script (02_...), extended genome-wide:
  - decompose compound swave SVs (SVTYPE = DEL+DEL+invDUP etc.) into typed
    sub-events via BKPS reference coordinates,
  - overlap each sub-event with IRGSP-1.0 gene/exon/CDS/promoter tracks
    (bedtools, for genome-wide speed),
  - assign a feature class and a gene-impact category.

VALIDATION STATUS (important):
  - chr1 calls were benchmarked (truvari vs curated truth); they get
    validation_status = 'validated' (overlaps a tp-comp call) or
    'chr1_unmatched'.
  - chr2-12 were NOT benchmarked -> validation_status = 'candidate'
    (computationally predicted, unvalidated). Do not report chr2-12 counts
    with the same confidence as chr1.
"""
import re, subprocess, gzip, sys
from pathlib import Path
from collections import Counter, defaultdict

BASE   = Path("/Users/cazhia/Desktop/SP-Scripts/0_annotation")
BT     = "/Users/cazhia/miniforge3/bin/bedtools"
VCF    = BASE / "whole_genome" / "swave.sample_level.split.vcf"
TRACKS = BASE / "reference_tracks"
GENES  = TRACKS / "IRGSP_genes.bed"
PROMO  = TRACKS / "IRGSP_promoters_2kb_upstream.bed"
EXON   = TRACKS / "IRGSP_exons_merged.bed"
CDS    = TRACKS / "IRGSP_CDS_merged.bed"
TPCOMP = {
    "Azucena": BASE / "inputs" / "pggb-swave_Azucena_chr1_tp-comp.vcf.gz",
    "IR64":    BASE / "inputs" / "pggb-swave_IR64_chr1_tp-comp.vcf.gz",
}
OUTDIR = BASE / "whole_genome"
CHROMS = {str(i) for i in range(1, 13)}
MIN_SIZE, VALID_TOL = 50, 500

# ---- chr1 tp-comp footprints (validated flag; chr1 only) ------------------
def load_tpcomp(path):
    ivs = []
    with gzip.open(path, "rt") as f:
        for l in f:
            if l.startswith("#"):
                continue
            c = l.split("\t")
            if c[0].replace("chr", "") != "1":
                continue
            pos, ref = int(c[1]), c[3]
            ivs.append((pos - 1 - VALID_TOL, pos - 1 + max(1, len(ref)) + VALID_TOL))
    return ivs
tpcomp = {s: load_tpcomp(p) for s, p in TPCOMP.items()}
def validated_chr1(s0, e, sample):
    return any(s0 < ie and e > is0 for is0, ie in tpcomp[sample])

# ---- BKPS sub-event parser (same as chr1 script) --------------------------
BK = re.compile(r"^(?P<type>[A-Za-z]+(?:_[A-Za-z]+)*)_(?P<size>\d+)_"
                r"(?P<chrom>[^_]+)_(?P<start>\d+)_(?P<end>\d+)"
                r"(?:_(?P<c2>[^_]+)_(?P<s2>\d+)_(?P<e2>\d+))?$")
def coarse(t):
    tl = t.lower()
    for k, v in (("dup","DUP"),("inv","INV"),("del","DEL"),("ins","INS")):
        if k in tl: return v
    return "COMPLEX"
def info_get(info, key):
    m = re.search(rf"(?:^|;){key}=([^;]+)", info); return m.group(1) if m else None

# ---- parse VCF -> sub-event list -----------------------------------------
subs = []   # each: dict with chrom,start(1b),end,subtype,coarse,parent_pos,parent_svtype,sample,vstatus
with open(VCF) as fh:
    idx = None
    for l in fh:
        if l.startswith("##"): continue
        if l.startswith("#CHROM"):
            cols = l.rstrip("\n").split("\t"); idx = {s: 9+i for i,s in enumerate(cols[9:])}; continue
        c = l.rstrip("\n").split("\t")
        chrom = c[0].replace("Nipponbare#1#chr", "").replace("chr", "")
        if chrom not in CHROMS: continue
        pos, info = int(c[1]), c[7]
        endi = int(info_get(info, "END") or pos)
        ptype = info_get(info, "SVTYPE") or "."
        bkps = info_get(info, "BKPS")
        for sample in ("Azucena", "IR64"):
            if "1" not in c[idx[sample]].split(":")[0]:  # sample must carry ALT
                continue
            if chrom == "1":
                vstatus = "validated" if validated_chr1(pos-1, endi, sample) else "chr1_unmatched"
            else:
                vstatus = "candidate"
            events = []
            if bkps:
                for tok in bkps.split(","):
                    m = BK.match(tok)
                    if m: events.append((m["type"], int(m["start"]), int(m["end"])))
            if not events:
                events = [(ptype, pos, endi)]
            for st, s, e in events:
                if e - s < MIN_SIZE: continue
                subs.append(dict(chrom=chrom, start=s, end=e, subtype=st, coarse=coarse(st),
                                 ppos=pos, ptype=ptype, sample=sample, vstatus=vstatus))

# ---- write sub-events BED, run bedtools overlaps --------------------------
sub_bed = OUTDIR / "_subevents.bed"
with open(sub_bed, "w") as w:
    for i, s in enumerate(subs):
        w.write(f"{s['chrom']}\t{s['start']-1}\t{s['end']}\t{i}\n")

def run(args):
    return subprocess.run(args, capture_output=True, text=True, check=True).stdout

def hits_named(track):  # id -> set(gene "id|name")
    out = defaultdict(set)
    res = run([BT, "intersect", "-a", str(sub_bed), "-b", str(track), "-wa", "-wb"])
    for line in res.splitlines():
        a = line.split("\t"); out[int(a[3])].add(a[7])
    return out
def hits_flag(track):   # set(id) overlapping
    res = run([BT, "intersect", "-a", str(sub_bed), "-b", str(track), "-u"])
    return {int(l.split("\t")[3]) for l in res.splitlines()}

gene_hits = hits_named(GENES)
prom_hits = hits_named(PROMO)
in_cds    = hits_flag(CDS)
in_exon   = hits_flag(EXON)

# ---- classify + write -----------------------------------------------------
def feature_class(i, has_gene):
    if i in in_cds:  return "CDS(coding)"
    if i in in_exon: return "exon/UTR"
    if has_gene:     return "intron"
    return None
def impact(ct, fclass, full):
    if ct == "DEL":
        return ("deleted_gene" if full else "truncated_gene" if fclass=="CDS(coding)"
                else "utr_deletion" if fclass=="exon/UTR"
                else "promoter_deletion" if fclass=="promoter" else "intronic_deletion")
    if ct == "DUP":
        return "promoter_duplication" if fclass=="promoter" else ("whole_gene_duplication" if full else "duplicated_gene")
    if ct == "INV":
        return ("inverted_gene" if full else "inversion_breakpoint" if fclass in ("CDS(coding)","exon/UTR")
                else "promoter_inversion" if fclass=="promoter" else "intronic_inversion")
    if ct == "INS":
        return ("coding_insertion" if fclass in ("CDS(coding)","exon/UTR")
                else "promoter_insertion" if fclass=="promoter" else "intronic_insertion")
    return "complex_" + (fclass.split("(")[0] if fclass else "genic")

hdr = ["sample","chrom","sub_start","sub_end","sub_type","coarse_type","parent_svtype",
       "parent_pos","feature_class","gene_impact","gene","validation_status"]
rows = []
for i, s in enumerate(subs):
    gh = gene_hits.get(i, set())
    fclass = feature_class(i, bool(gh))
    if fclass is None:
        ph = prom_hits.get(i, set())
        if not ph: continue           # intergenic
        fclass, used = "promoter", ph
    else:
        used = gh
    s0, e = s["start"]-1, s["end"]
    for g in sorted(used):
        gs, ge = None, None
        # gene interval for full-containment test (parse from genes track name lookup omitted; approx via coords)
        rows.append([s["sample"], s["chrom"], s["start"], s["end"], s["subtype"], s["coarse"],
                     s["ptype"], s["ppos"], fclass,
                     impact(s["coarse"], fclass, False), g, s["vstatus"]])

# full-containment refinement for deleted/inverted: recompute using gene coords
gene_iv = {}
for l in open(GENES):
    c = l.rstrip("\n").split("\t")
    gene_iv[c[3]] = (int(c[1]), int(c[2]))
for r in rows:
    ct = r[5]; g = r[10]; s0 = r[2]-1; e = r[3]
    if g in gene_iv and ct in ("DEL","DUP","INV"):
        gs, ge = gene_iv[g]
        full = gs >= s0 and ge <= e
        r[9] = impact(ct, r[8], full)

rows.sort(key=lambda x: (x[0], int(x[1]) if x[1].isdigit() else 99, x[2]))
out_all = OUTDIR / "pggb-swave_wholegenome_gene_impact.tsv"
with open(out_all, "w") as w:
    w.write("\t".join(hdr) + "\n")
    for r in rows: w.write("\t".join(map(str, r)) + "\n")
sub_bed.unlink(missing_ok=True)

print(f"wrote {out_all}  ({len(rows)} genic/promoter sub-event x gene rows)\n")
for sample in ("Azucena","IR64"):
    ss = [r for r in rows if r[0]==sample]
    print(f"== {sample}: {len(ss)} rows, {len({r[10] for r in ss})} genes ==")
    print("   by validation_status:", dict(Counter(r[11] for r in ss)))
    print("   chr1 validated genes :", len({r[10] for r in ss if r[11]=='validated'}))
    print("   gene_impact          :", dict(Counter(r[9] for r in ss)))
