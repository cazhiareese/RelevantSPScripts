# =============================================================================
# Makefile — generic SV benchmarking of every method against its ground truth.
#   truvari bench per (method x sample), logs land in <method>/bench_out/<sample>/
# =============================================================================
#
# One makefile, two branches (choose with BRANCH=):
#
#   BRANCH=real       pangenome callsets are multi-sample (IR64, Azucena[, Nipponbare]).
#                     Each pangenome method is filtered to a sample and scored
#                     against that accession's SILVER truth:
#                        real/ground_truth/ir64_ground_truth.vcf
#                        real/ground_truth/azucena_ground_truth.vcf
#                     svim-asm is single-sample (Nipponbare-vs-<acc>), scored directly.
#
#   BRANCH=simulated  every method is filtered to the injected "Simulated" hap and
#                     scored against  simulated/pggb-swave/sv_truth.vcf.
#
# Methods discovered per branch:
#   <graph>-swave/swave.sample_level.split.vcf   (graph = cactus|minigraph|pggb)
#   <graph>-vcfwave/<graph>.vcfwave.vcf.gz
#   svim-asm variants.vcf (single sample)
#
# Output: <method-folder>/bench_out/<sample>/{truvari/, truth.vcf.gz, calls.vcf.gz,
#         truvari.stratified.tsv, bench.log}
#
# SILVER truth (real branch): only validated SVs are in the truth, so a callset SV
# absent from it may be a REAL SV you didn't validate, NOT a false positive.
#   => RECALL is the trustworthy metric; precision is a LOWER BOUND.
#
# Quick start:
#   make BRANCH=real                       # all real-branch benchmarks
#   make BRANCH=real bench-IR64-pggb-vcfwave   # one run
#   make BRANCH=real summary
#   make BRANCH=simulated                  # all simulated-branch benchmarks
#   make BRANCH=simulated summary
# =============================================================================

SHELL          := /bin/bash
.SHELLFLAGS    := -o pipefail -c
.DELETE_ON_ERROR:

# ---- Config (override on the command line) ----------------------------------
BRANCH         ?=
GRAPHS         ?= cactus minigraph pggb
WAVE_METHODS   ?= swave vcfwave

# SV types to DROP from the metric, per method family. The graph callers (swave,
# vcfwave) decompose translocations into local INS/DEL and never emit SVTYPE=BND,
# so scoring BND only sinks their recall with unrepresentable truth. We drop BND
# for svim-asm too, so every method is scored against the SAME non-BND truth set
# (a fair comparison); under these lenient settings svim-asm matches no BND truth
# anyway, so dropping them removes only guaranteed-FN events from its denominator.
#   GRAPH_DROP_TYPES=      -> also score BND for the graph methods
#   GRAPH_DROP_TYPES="BND INV" -> drop more types
GRAPH_DROP_TYPES   ?= BND
SVIMASM_DROP_TYPES ?= BND

# Common contig name that truth + calls are renamed to before matching, and its
# length (Nipponbare chr1) — truvari requires a length in the ##contig header, and
# the real ground-truth VCFs ship without a contig line, so we inject one.
BENCH_CHROM    ?= chr1
CHROM_LEN      ?= 43270923

# Tools. truvari comes from the local venv; bcftools/tabix must be on PATH
# (override e.g. BCFTOOLS="conda run -n vgsv bcftools").
VENV           ?= benchmark_venv
TRUVARI        ?= $(VENV)/bin/truvari
PYTHON         ?= $(VENV)/bin/python
BCFTOOLS       ?= bcftools
TABIX          ?= tabix
STRATIFY       ?= $(PYTHON) ../helper-scripts/stratify_truvari.py

# Truvari matching. Symbolic truth -> --pctseq 0 (no reference FASTA needed).
# --bnddist matches translocation (BND) breakends by proximity; --dup-to-ins folds
# <DUP>; --typeignore tolerates caller type labels (needed so the graph callers'
# INS/DEL decompositions match DUP/other truth). Lenient (see recall caveat).
TRUVARI_OPTS   ?= --pctseq 0 --pctsize 0.5 --refdist 1000 --sizemin 50 \
                  --sizemax 125000 --dup-to-ins  --bnddist 500 --typeignore

# ---- Per-branch wiring ------------------------------------------------------
# TRUTH_<sample> / SVIMASM_<sample> are looked up by computed name in the rules.
ifeq ($(BRANCH),real)
ROOT             := real
SAMPLES          := IR64 Azucena
TRUTH_IR64       := real/ground_truth/ir64_ground_truth.vcf
TRUTH_Azucena    := real/ground_truth/azucena_ground_truth.vcf
SVIMASM_IR64     := real/svimasm/IR64/variants.vcf
SVIMASM_Azucena  := real/svimasm/Azucena/Azucena/variants.vcf
SVIMASM_OUTROOT  := real/svimasm
endif
ifeq ($(BRANCH),simulated)
ROOT               := simulated
SAMPLES            := Simulated
TRUTH_Simulated    := simulated/pggb-swave/sv_truth.vcf
SVIMASM_Simulated  := simulated/svimasm_out/variants.vcf
SVIMASM_OUTROOT    := simulated/svimasm_out
endif

# Callset path for a pangenome (graph,method).  $(call CALLSET,<graph>,<method>)
CALLSET = $(if $(filter swave,$(2)),$(ROOT)/$(1)-swave/swave.sample_level.converted.vcf,$(ROOT)/$(1)-vcfwave/$(1).vcfwave.vcf.gz)

# Build a bcftools -e expression that drops records of the given SVTYPEs.
#   $(call EXCL_EXPR,BND INV) -> INFO/SVTYPE="BND" || INFO/SVTYPE="INV"
empty :=
space := $(empty) $(empty)
EXCL_EXPR = $(subst $(space), || ,$(strip $(foreach t,$(1),INFO/SVTYPE="$(t)")))

# `make BRANCH=<b>` (no explicit target) runs the benchmarks; bare `make` shows help.
ifeq ($(strip $(BRANCH)),)
.DEFAULT_GOAL := help
else
.DEFAULT_GOAL := bench
endif

# ----------------------------------------------------------------------------
# One benchmark job.
#   $(call BENCHJOB,<slug>,<callset>,<truth>,<filter-sample>,<outdir>,<drop-types>)
#   <filter-sample> empty => no sample filter (single-sample svim-asm callset).
#   <drop-types>    empty => keep all types; else drop those SVTYPEs from BOTH
#                            truth and calls before matching (e.g. BND).
# Emitted via $(eval) so each becomes a real `bench-<slug>` target.
# awk reheads CHROM (and any ##contig, injecting one if absent) to $(BENCH_CHROM).
# ----------------------------------------------------------------------------
define REHEAD_AWK
awk 'BEGIN{FS=OFS="\t";c=0;st=0;sl=0;en=0;mp=0} /^##contig=<ID=/{sub(/ID=[^,>]+/,"ID=$(BENCH_CHROM)");print;c=1;next} /^##INFO=<ID=SVTYPE[,>]/{st=1} /^##INFO=<ID=SVLEN[,>]/{sl=1} /^##INFO=<ID=END[,>]/{en=1} /^##INFO=<ID=MATEPOS[,>]/{mp=1} /^#CHROM/{if(!c)print "##contig=<ID=$(BENCH_CHROM),length=$(CHROM_LEN)>"; if(!st)print "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"SV type\">"; if(!sl)print "##INFO=<ID=SVLEN,Number=.,Type=Integer,Description=\"SV length\">"; if(!en)print "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position\">"; if(!mp)print "##INFO=<ID=MATEPOS,Number=1,Type=Integer,Description=\"Mate position\">"; print;next} /^#/{print;next} {$$$$1="$(BENCH_CHROM)";print}'
endef

define BENCHJOB
BENCH_JOBS     += $(1)
SUMMARY_ITEMS  += $(1)=$(5)
.PHONY: bench-$(1)
bench-$(1):
	@mkdir -p $(5)
	@if [ ! -s "$(2)" ]; then echo "SKIP [$(1)]: callset not found: $(2)"; exit 0; fi; \
	 if [ ! -s "$(3)" ]; then echo "SKIP [$(1)]: truth not found: $(3)"; exit 0; fi; \
	 echo ">>> [$(1)] normalize truth -> $(BENCH_CHROM)$(if $(6), (drop types: $(6)),)"; \
	 zcat -f "$(3)" | $(REHEAD_AWK) \
	   $(if $(6),| $(BCFTOOLS) view -e '$(call EXCL_EXPR,$(6))' -Ov -,) \
	   | $(BCFTOOLS) sort -Oz -o "$(5)/truth.vcf.gz" -; \
	 $(TABIX) -f -p vcf "$(5)/truth.vcf.gz"; \
	 echo ">>> [$(1)] normalize calls$(if $(4), (filter sample=$(4)),)"; \
	 zcat -f "$(2)" | $(REHEAD_AWK) \
	   | $(BCFTOOLS) norm -m -any -Ov - 2>/dev/null \
	   $(if $(4),| $(BCFTOOLS) view -s $(4) --min-ac 1 -Ov -,) \
	   $(if $(6),| $(BCFTOOLS) view -e '$(call EXCL_EXPR,$(6))' -Ov -,) \
	   | $(BCFTOOLS) sort -Oz -o "$(5)/calls.vcf.gz" -; \
	 $(TABIX) -f -p vcf "$(5)/calls.vcf.gz"; \
	 echo ">>> [$(1)] truvari bench"; rm -rf "$(5)/truvari"; \
	 $(TRUVARI) bench -b "$(5)/truth.vcf.gz" -c "$(5)/calls.vcf.gz" -o "$(5)/truvari" $(TRUVARI_OPTS); \
	 $(STRATIFY) --bench-dir "$(5)/truvari" --out "$(5)/truvari.stratified.tsv" >/dev/null 2>&1 || true; \
	 { \
	   echo "==================================================================="; \
	   echo " benchmark [$(BRANCH)] : $(1)"; \
	   echo " date    : $$$$(date '+%F %T')"; \
	   echo " truth   : $(3)"; \
	   echo " callset : $(2)$(if $(4),   (filtered to sample $(4)),)"; \
	   echo " opts    : $(TRUVARI_OPTS)"; \
	   $(if $(6),echo " dropped : SVTYPE in [$(6)] excluded from truth+calls";,) \
	   $(if $(filter real,$(BRANCH)),echo " NOTE    : SILVER truth -> RECALL is the metric; precision is a lower bound.";,) \
	   echo "==================================================================="; \
	   echo; echo "-- summary.json --"; cat "$(5)/truvari/summary.json" 2>/dev/null; echo; \
	   echo "-- Stratified (size-bin x type; BND=translocation) --"; \
	   column -t "$(5)/truvari.stratified.tsv" 2>/dev/null || echo "(unavailable)"; \
	 } | tee "$(5)/bench.log"; \
	 echo ">>> log: $(5)/bench.log"
endef

# ---- Generate jobs: pangenome (sample-filtered) + svim-asm (single-sample) --
$(foreach s,$(SAMPLES),$(foreach g,$(GRAPHS),$(foreach m,$(WAVE_METHODS),\
  $(eval $(call BENCHJOB,$(s)-$(g)-$(m),$(call CALLSET,$(g),$(m)),$(TRUTH_$(s)),$(s),$(ROOT)/$(g)-$(m)/bench_out/$(s),$(GRAPH_DROP_TYPES))))))

$(foreach s,$(SAMPLES),\
  $(eval $(call BENCHJOB,$(s)-svimasm,$(SVIMASM_$(s)),$(TRUTH_$(s)),,$(SVIMASM_OUTROOT)/bench_out/$(s),$(SVIMASM_DROP_TYPES))))

# ----------------------------------------------------------------------------
# Aggregates (defined AFTER the evals so $(BENCH_JOBS) is populated).
# ----------------------------------------------------------------------------
.PHONY: bench summary help clean
bench: $(addprefix bench-,$(BENCH_JOBS))
	@test -n "$(strip $(BENCH_JOBS))" \
	  || { echo "ERROR: set BRANCH=real or BRANCH=simulated (got '$(BRANCH)')" >&2; exit 1; }
	@echo ">>> all [$(BRANCH)] benchmarks done — logs in <method>/bench_out/<sample>/bench.log"

summary:
	@test -n "$(strip $(SUMMARY_ITEMS))" \
	  || { echo "ERROR: set BRANCH=real or BRANCH=simulated (got '$(BRANCH)')" >&2; exit 1; }
	@printf '%-30s %8s %10s %6s\n' "job" "recall" "precision" "f1"
	@for item in $(SUMMARY_ITEMS); do \
	  slug="$${item%%=*}"; dir="$${item##*=}"; s="$$dir/truvari/summary.json"; \
	  if [ -f "$$s" ]; then \
	    printf '%-30s ' "$$slug"; \
	    $(PYTHON) -c "import json;d=json.load(open('$$s'));print('%8.4f %10.4f %6.4f'%(d['recall'],d['precision'],d['f1']))"; \
	  fi; \
	done
	@echo "(real = SILVER truth: recall is the metric, precision is a lower bound)"

help:
	@echo "Generic SV benchmark — pick a branch with BRANCH=real | simulated."
	@echo ""
	@echo "  make BRANCH=<b>                      run every method x sample"
	@echo "  make BRANCH=<b> bench-<slug>         one run (see slugs below)"
	@echo "  make BRANCH=<b> summary              recall/precision/f1 table"
	@echo "  make BRANCH=<b> clean                remove this branch's bench_out/"
	@echo ""
	@echo "Branch = '$(BRANCH)'"
ifneq ($(strip $(BENCH_JOBS)),)
	@echo "  ROOT    = $(ROOT)"
	@echo "  SAMPLES = $(SAMPLES)"
	@echo "  jobs    :"
	@for j in $(BENCH_JOBS); do echo "     bench-$$j"; done
else
	@echo "  (set BRANCH=real or BRANCH=simulated to list jobs)"
endif
	@echo ""
	@echo "  TRUVARI      = $(TRUVARI)"
	@echo "  TRUVARI_OPTS = $(TRUVARI_OPTS)"
	@echo "  GRAPH_DROP_TYPES   = $(GRAPH_DROP_TYPES)   (SVTYPEs dropped for swave/vcfwave; graph callers emit no BND)"
	@echo "  SVIMASM_DROP_TYPES = $(SVIMASM_DROP_TYPES)   (BND dropped so svim-asm is scored on the same non-BND truth)"

clean:
	@test -n "$(ROOT)" || { echo "ERROR: set BRANCH=real or BRANCH=simulated" >&2; exit 1; }
	rm -rf $(ROOT)/*/bench_out $(SVIMASM_OUTROOT)/bench_out
	@echo ">>> removed $(BRANCH) bench_out/ dirs"
