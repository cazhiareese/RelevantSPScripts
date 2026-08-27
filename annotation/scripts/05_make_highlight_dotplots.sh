#!/usr/bin/env bash
# Generate slide dotplots for the two highlighted loci:
#   1. SLB1/SLB2 58 kb deletion in IR64 (chr1:28.99-29.04 Mb)  -> generated fresh
#   2. Gn1a intronic deletion in Azucena (validated, = DEL_34)  -> reuse existing
# Run: conda run -n mummer_env bash 05_make_highlight_dotplots.sh   (or just bash, it wraps with CRUN)
set -euo pipefail
cd /Users/cazhia/Desktop/SP-Scripts
CRUN="conda run -n mummer_env"
BA=bio_align
OUT=0_annotation/figures
mkdir -p "$OUT" /tmp/dp

REF="$BA/Nipponbare_chr1.fasta"; REFCHR="Nipponbare#1#chr1"

make_del_dotplot () {         # $1 name  $2 queryfasta  $3 qchr  $4 sv_start  $5 sv_end  $6 flank
  local name=$1 qfa=$2 qchr=$3 s=$4 e=$5 fl=$6
  local rs=$((s-fl)) re=$((e+fl))
  echo ">> $name : ref $REFCHR:$rs-$re  (SV $s-$e)"
  $CRUN samtools faidx "$REF" "$REFCHR:${rs}-${re}" > /tmp/dp/${name}_ref.fa
  # pass1: locate homologous region in the query assembly
  $CRUN nucmer -l 100 -c 1000 --prefix /tmp/dp/${name}_p1 /tmp/dp/${name}_ref.fa "$qfa"
  $CRUN show-coords -T -r -l -c /tmp/dp/${name}_p1.delta > /tmp/dp/${name}_p1.coords
  read qlo qhi < <(awk 'NR>4 && $5>2000{lo=($3<$4?$3:$4); hi=($3>$4?$3:$4);
        if(min==""||lo<min)min=lo; if(hi>max)max=hi} END{print min-5000, max+5000}' /tmp/dp/${name}_p1.coords)
  echo "   query $qchr:$qlo-$qhi"
  $CRUN samtools faidx "$qfa" "${qchr}:${qlo}-${qhi}" > /tmp/dp/${name}_qry.fa
  # pass2: window-vs-window alignment + plot
  $CRUN nucmer -l 100 -c 200 --prefix /tmp/dp/${name} /tmp/dp/${name}_ref.fa /tmp/dp/${name}_qry.fa
  $CRUN mummerplot --png --large -p "$OUT/${name}" /tmp/dp/${name}.delta
  echo "   wrote $OUT/${name}.png"
}

# 1. SLB1/SLB2 deletion in IR64 (58 kb region), 20 kb flanks
make_del_dotplot "SLB1_SLB2_IR64_deletion" "$BA/IR64_chr1.fasta" "IR64#1#chr1" 28986552 29044703 20000

# 2. Gn1a — reuse the validated Azucena DEL_34 dotplot
cp "$BA/validation_samples/azucena/Azucena_DEL_34.png" "$OUT/Gn1a_Azucena_validated_deletion.png"
echo ">> Gn1a : copied validated Azucena_DEL_34.png -> $OUT/Gn1a_Azucena_validated_deletion.png"

ls -la "$OUT"
