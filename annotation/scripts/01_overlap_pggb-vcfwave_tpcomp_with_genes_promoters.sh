#!/usr/bin/env bash
# Overlap tp-comp SVs with IRGSP-1.0 gene annotation.
set -euo pipefail
cd /Users/cazhia/Desktop/SP-Scripts/annotation
BT=/Users/cazhia/miniforge3/bin/bedtools
GFF=Oryza_sativa.IRGSP-1.0.63.gff3

# --- One-time: derive feature BEDs from GFF (seqid is bare "1".."12") ---
if [ ! -s genes.bed ]; then
  awk -F'\t' '$3=="gene"||$3=="ncRNA_gene"{
      id=$9; sub(/.*ID=gene:/,"",id); sub(/;.*/,"",id);
      nm=$9; if(nm ~ /Name=/){sub(/.*Name=/,"",nm); sub(/;.*/,"",nm)} else nm=".";
      print $1"\t"$4-1"\t"$5"\t"id"|"nm"\t.\t"$7
  }' "$GFF" | sort -k1,1 -k2,2n > genes.bed
  awk -F'\t' '$3=="exon"{print $1"\t"$4-1"\t"$5}' "$GFF" | sort -k1,1 -k2,2n | $BT merge -i - > exon.bed
  awk -F'\t' '$3=="CDS"{print $1"\t"$4-1"\t"$5}'  "$GFF" | sort -k1,1 -k2,2n | $BT merge -i - > cds.bed
  # Promoter: 2 kb upstream of TSS, strand-aware, keeps gene id|name
  awk -F'\t' '$3=="gene"||$3=="ncRNA_gene"{
      id=$9; sub(/.*ID=gene:/,"",id); sub(/;.*/,"",id);
      nm=$9; if(nm ~ /Name=/){sub(/.*Name=/,"",nm); sub(/;.*/,"",nm)} else nm=".";
      if($7=="+"){ s=$4-1-2000; e=$4-1 } else { s=$5; e=$5+2000 }
      if(s<0)s=0; if(e>s) print $1"\t"s"\t"e"\t"id"|"nm"\t.\t"$7
  }' "$GFF" | sort -k1,1 -k2,2n > promoter.bed
fi

sample="$1"          # Azucena | IR-64
vcf="${sample}-tp-comp.vcf.gz"
out="${sample}_sv_gene_overlap"

# SV reference footprint -> BED, remap chrN -> N
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%ID\n' "$vcf" | awk -F'\t' '{
    chr=$1; sub(/^chr/,"",chr);
    reflen=length($3);
    start=$2-1; end=$2-1+reflen;
    svlen=length($4)-reflen; kind=(svlen>0?"INS":(svlen<0?"DEL":"OTHER"));
    print chr"\t"start"\t"end"\t"$5"\t"(svlen<0?-svlen:svlen)"\t"kind
}' | sort -k1,1 -k2,2n > "${out}.svs.bed"

# Overlap against genes / exon / cds
$BT intersect -a "${out}.svs.bed" -b genes.bed    -wa -wb > "${out}.genes.txt"    || true
$BT intersect -a "${out}.svs.bed" -b promoter.bed -wa -wb > "${out}.promoter.txt" || true
$BT intersect -a "${out}.svs.bed" -b cds.bed  -u > "${out}.in_cds.bed"  || true
$BT intersect -a "${out}.svs.bed" -b exon.bed -u > "${out}.in_exon.bed" || true

# Per-SV classification table
awk -F'\t' '
  FNR==NR{if(FILENAME ~ /in_cds/)  cds[$4]=1;  next}
' /dev/null >/dev/null 2>&1 || true

python3 - "$out" <<'PY'
import sys
out=sys.argv[1]
def load_ids(f):
    s=set()
    try:
        for l in open(f):
            c=l.split('\t'); s.add(c[3])
    except FileNotFoundError: pass
    return s
cds=load_ids(f"{out}.in_cds.bed")
exon=load_ids(f"{out}.in_exon.bed")
def load_gene_map(f):
    d={}
    try:
        for l in open(f):
            c=l.rstrip('\n').split('\t')
            d.setdefault(c[3],set()).add(c[9])  # svid -> id|name
    except FileNotFoundError: pass
    return d
genes=load_gene_map(f"{out}.genes.txt")
prom=load_gene_map(f"{out}.promoter.txt")
rows=[]
for l in open(f"{out}.svs.bed"):
    c=l.rstrip('\n').split('\t')
    chrom,st,en,svid,size,kind=c[0],c[1],c[2],c[3],c[4],c[5]
    gl=genes.get(svid,set()); pl=prom.get(svid,set())
    if svid in cds: cls,gg="CDS(coding)",gl
    elif svid in exon: cls,gg="exon/UTR",gl
    elif gl: cls,gg="intron",gl
    elif pl: cls,gg="promoter",pl
    else: cls,gg="intergenic",set()
    rows.append((chrom,int(st)+1,int(en),svid,kind,size,cls,";".join(sorted(gg)) if gg else "."))
rows.sort(key=lambda r:(r[0],r[1]))
hdr=["chrom","pos","end","sv_id","type","size","overlap_class","genes"]
with open(f"{out}.tsv","w") as w:
    w.write("\t".join(hdr)+"\n")
    for r in rows: w.write("\t".join(map(str,r))+"\n")
from collections import Counter
cc=Counter(r[6] for r in rows)
print(f"== {out}  (n={len(rows)} SVs) ==")
for k in ["CDS(coding)","exon/UTR","intron","promoter","intergenic"]:
    print(f"  {k:12s} {cc.get(k,0)}")
PY
