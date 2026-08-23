# =============================================================================
# Makefile — Build Chr1 pangenome graphs for the REAL-SV arm (NO simulated hap).
#            minigraph, minigraph-cactus, pggb on Nipponbare + IR64 + Azucena.
# =============================================================================
#
# This is the real-variation arm: the graph contains only the three real rice
# cultivars (no Simulated haplotype), so the SVs called from it are the genuine
# IR64/Azucena-vs-Nipponbare structural variants — to be validated against the
# dotplot-confirmed SyRI silver truth.
#
# Self-contained (config inlined — does NOT include config.mk), so it can be
# dropped into real_sv_arm/ and run directly.
#
# Inputs (lowercase, in $(INPUT_DIR)):
#     nipponbare_chr1.fna
#     ir64_chr1.fna
#     azucena_chr1.fna
# The first cultivar in $(CULTIVARS) is the reference backbone / --reference.
#
# Usage:
#     make all                 # build all three graphs
#     make minigraph | cactus | pggb
#     make inputs              # just stage/index the PanSN FASTAs
#     make THREADS=64 all
#     make CACTUS_PANGENOME="sudo docker run ... cactus-pangenome" cactus
#     make clean | help
# =============================================================================

SHELL          := /bin/bash
.SHELLFLAGS    := -o pipefail -c
.DELETE_ON_ERROR:

# ---- Configuration (override on the command line) ---------------------------
# Real cultivars only — Simulated is intentionally absent. FIRST = reference.
CULTIVARS      ?= Nipponbare IR64 Azucena 
REFERENCE      ?= $(firstword $(CULTIVARS))

THREADS        ?= 32

INPUT_DIR      ?= .
OUT_DIR        ?= graphs
PANSN_DIR      ?= $(OUT_DIR)/pansn
PREFIX         ?= chr1_pangenome
EXT            ?= fna

# Tool commands (override to use containers / modules / conda).
MINIGRAPH        ?= minigraph
CACTUS_PANGENOME ?= cactus-pangenome
PGGB             ?= pggb
SAMTOOLS         ?= samtools
BGZIP            ?= bgzip

# ---- Graph-build parameters ----
MINIGRAPH_OPTS ?= -cxggs
PGGB_P         ?= 90
PGGB_S         ?= 5000
PGGB_K         ?= 23
# pggb -V spec: REF:DELIM[:LEN]. REF = reference sample, DELIM = PanSN separator.
# Makes pggb run `vg deconstruct` and emit a VCF against the reference.
PGGB_VCF_SPEC  ?= $(REFERENCE):#
CACTUS_OPTS    ?= --gfa clip full --vcf

# ---- Derived values ---------------------------------------------------------
# Lowercase a cultivar name to get its input filename stem: $(call lc,IR64)=ir64
lc = $(shell printf '%s' '$(1)' | tr '[:upper:]' '[:lower:]')

N_HAPS         := $(words $(CULTIVARS))
OTHERS         := $(filter-out $(REFERENCE),$(CULTIVARS))
MINIGRAPH_ORDER:= $(REFERENCE) $(OTHERS)

RAW_FASTAS     := $(foreach c,$(CULTIVARS),$(INPUT_DIR)/$(call lc,$(c))_chr1.$(EXT))
PANSN_FASTAS   := $(foreach c,$(CULTIVARS),$(PANSN_DIR)/$(c)_chr1.$(EXT))

COMBINED       := $(OUT_DIR)/$(PREFIX).combined.fa.gz
SEQFILE        := $(OUT_DIR)/$(PREFIX).seqfile.txt

MINIGRAPH_GFA  := $(OUT_DIR)/minigraph/$(PREFIX).minigraph.gfa
CACTUS_DIR     := $(OUT_DIR)/minigraph-cactus
CACTUS_JOBSTORE:= $(OUT_DIR)/mc_jobstore
CACTUS_GFA     := $(CACTUS_DIR)/$(PREFIX).gfa.gz
PGGB_DIR       := $(OUT_DIR)/pggb
PGGB_STAMP     := $(PGGB_DIR)/.pggb.done

# ----------------------------------------------------------------------------
# Top-level targets
# ----------------------------------------------------------------------------
.PHONY: all minigraph cactus pggb inputs help clean

all: minigraph cactus pggb

inputs: $(PANSN_FASTAS) $(COMBINED) $(SEQFILE)

# ----------------------------------------------------------------------------
# Stage inputs: normalise each header to "<Cultivar>#1#chr1" (PanSN).
# Generated per-cultivar so a lowercase input file (ir64_chr1.fna) maps to a
# cased PanSN sample (IR64 -> >IR64#1#chr1). Safe if inputs are already PanSN.
# $(call STAGE_INPUT,<CasedSample>)
# ----------------------------------------------------------------------------
define STAGE_INPUT
$(PANSN_DIR)/$(1)_chr1.$(EXT): $(INPUT_DIR)/$(call lc,$(1))_chr1.$(EXT)
	@mkdir -p $(PANSN_DIR)
	@echo ">>> PanSN-normalising $$< -> $$@ (sample $(1))"
	@awk -v s="$(1)" 'NR==1{print ">"s"#1#chr1"; next} /^>/{next} {print}' $$< > $$@
	@$(SAMTOOLS) faidx $$@
endef
$(foreach c,$(CULTIVARS),$(eval $(call STAGE_INPUT,$(c))))

# ----------------------------------------------------------------------------
# Combined bgzipped + faidx'd FASTA (pggb input)
# ----------------------------------------------------------------------------
$(COMBINED): $(PANSN_FASTAS)
	@mkdir -p $(OUT_DIR)
	@echo ">>> Concatenating $(N_HAPS) haplotypes -> $(COMBINED)"
	cat $(PANSN_FASTAS) | $(BGZIP) -@ $(THREADS) > $(COMBINED)
	$(SAMTOOLS) faidx $(COMBINED)

# ----------------------------------------------------------------------------
# Seqfile for minigraph-cactus: "<sample>\t<path>", reference first
# ----------------------------------------------------------------------------
$(SEQFILE): $(PANSN_FASTAS)
	@mkdir -p $(OUT_DIR)
	@echo ">>> Writing cactus seqfile -> $(SEQFILE)"
	@rm -f $(SEQFILE)
	@$(foreach c,$(MINIGRAPH_ORDER),\
	   printf '%s\t%s\n' "$(c)" "$(abspath $(PANSN_DIR)/$(c)_chr1.$(EXT))" >> $(SEQFILE);)
	@cat $(SEQFILE)

# ----------------------------------------------------------------------------
# 1) minigraph  (reference-based, incremental)
# ----------------------------------------------------------------------------
minigraph: $(MINIGRAPH_GFA)

$(MINIGRAPH_GFA): $(PANSN_FASTAS)
	@mkdir -p $(dir $(MINIGRAPH_GFA))
	@echo ">>> minigraph: ref=$(REFERENCE)  order=[$(MINIGRAPH_ORDER)]"
	$(MINIGRAPH) $(MINIGRAPH_OPTS) -t $(THREADS) \
	    $(foreach c,$(MINIGRAPH_ORDER),$(PANSN_DIR)/$(c)_chr1.$(EXT)) \
	    > $(MINIGRAPH_GFA)
	@echo ">>> minigraph GFA: $(MINIGRAPH_GFA)"

# ----------------------------------------------------------------------------
# 2) minigraph-cactus  (cactus-pangenome)
# ----------------------------------------------------------------------------
cactus: $(CACTUS_GFA)

$(CACTUS_GFA): $(SEQFILE)
	@mkdir -p $(CACTUS_DIR)
	@echo ">>> minigraph-cactus: ref=$(REFERENCE)"
	$(CACTUS_PANGENOME) $(CACTUS_JOBSTORE) $(SEQFILE) \
	    --outDir $(CACTUS_DIR) \
	    --outName $(PREFIX) \
	    --reference $(REFERENCE) \
	    --mgCores $(THREADS) \
	    --consCores $(THREADS) \
	    $(CACTUS_OPTS)
	@echo ">>> minigraph-cactus output: $(CACTUS_DIR)"

# ----------------------------------------------------------------------------
# 3) pggb  (all-vs-all, reference-free)
# ----------------------------------------------------------------------------
pggb: $(PGGB_STAMP)

$(PGGB_STAMP): $(COMBINED)
	@mkdir -p $(PGGB_DIR)
	@echo ">>> pggb: n=$(N_HAPS) p=$(PGGB_P) s=$(PGGB_S) k=$(PGGB_K)"
	-$(PGGB) -i $(COMBINED) -o $(PGGB_DIR) \
	    -n $(N_HAPS) -t $(THREADS) \
	    -p $(PGGB_P) -s $(PGGB_S) -k $(PGGB_K) -V $(PGGB_VCF_SPEC)
	@gfa=$$(ls $(PGGB_DIR)/*.smooth.final.gfa 2>/dev/null | head -1); \
	  if [ -s "$$gfa" ]; then \
	    touch $(PGGB_STAMP); \
	    echo ">>> pggb graph OK: $$gfa  (a failed odgi layout/viz step is cosmetic-only)"; \
	  else \
	    echo "ERROR: pggb produced no *.smooth.final.gfa in $(PGGB_DIR)" >&2; exit 1; \
	  fi

# ----------------------------------------------------------------------------
# Housekeeping
# ----------------------------------------------------------------------------
help:
	@echo "REAL-SV arm graph build (no Simulated hap): all | minigraph | cactus | pggb | inputs"
	@echo "Other: clean | help"
	@echo ""
	@echo "Config (override as 'make VAR=value target'):"
	@echo "  CULTIVARS  = $(CULTIVARS)"
	@echo "  REFERENCE  = $(REFERENCE)   (first cultivar)"
	@echo "  N_HAPS     = $(N_HAPS)"
	@echo "  THREADS    = $(THREADS)"
	@echo "  INPUT_DIR  = $(INPUT_DIR)"
	@echo "  OUT_DIR    = $(OUT_DIR)"
	@echo "  Inputs expected: $(RAW_FASTAS)"
	@echo ""
	@echo "  MINIGRAPH        = $(MINIGRAPH)"
	@echo "  CACTUS_PANGENOME = $(CACTUS_PANGENOME)"
	@echo "  PGGB             = $(PGGB)"

clean:
	rm -rf $(OUT_DIR)
