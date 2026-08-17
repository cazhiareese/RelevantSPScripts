# =============================================================================
# Makefile — Build Chr1 pangenome graphs with minigraph, minigraph-cactus, pggb
# =============================================================================
#
# Inputs (place in $(INPUT_DIR), one sequence each, PanSN headers
#         "<Cultivar>#1#chr1"):
#     Nipponbare_chr1.fa
#     IR64_chr1.fa
#     Azucena_chr1.fa
#     Simulated_chr1.fa
#
# The first cultivar in $(CULTIVARS) is the reference backbone for minigraph
# and the --reference sample for minigraph-cactus.
#
# Usage on the server:
#     make all                 # build all three graphs
#     make minigraph           # minigraph only
#     make cactus              # minigraph-cactus only
#     make pggb                # pggb only
#     make inputs              # just stage/index the PanSN FASTAs
#     make THREADS=64 all      # override thread count
#     make clean               # remove all generated graphs + staged inputs
#     make help                # list targets and current config
#
# Tool invocations are variables so you can point them at a module load,
# conda env, container, or absolute path without editing the recipes, e.g.:
#     make MINIGRAPH="singularity exec minigraph.sif minigraph" minigraph
#     make PGGB="singularity exec pggb.sif pggb" pggb
#
# The vg deconstruct → SV-VCF step (Swave vs vg-vcfbub-vcfwave) lives in a
# separate makefile that consumes the graphs built here:
#     make -f Makefile.deconstruct callsets
#
# Shared configuration (cultivars, paths, tool names) is in config.mk.
# =============================================================================

include config.mk

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
# ----------------------------------------------------------------------------
# $(call STAGE_INPUT,<CasedSample>)
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
	$(PGGB) -i $(COMBINED) -o $(PGGB_DIR) \
	    -n $(N_HAPS) -t $(THREADS) \
	    -p $(PGGB_P) -s $(PGGB_S) -k $(PGGB_K) -V $(PGGB_VCF_SPEC)
	@touch $(PGGB_STAMP)
	@echo ">>> pggb output: $(PGGB_DIR)"

# ----------------------------------------------------------------------------
# Housekeeping
# ----------------------------------------------------------------------------
help:
	@echo "Graph build : all | minigraph | cactus | pggb | inputs"
	@echo "Other       : clean | help"
	@echo ""
	@echo "  Then derive SV callsets with the separate deconstruct makefile:"
	@echo "    make -f Makefile.deconstruct callsets"
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
