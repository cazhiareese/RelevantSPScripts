# =============================================================================
# Makefile – Tool installation for Rice Pangenome & Structural Variant Study
#
# Study: Graph-Based Rice Pangenome for Structural Variant Identification
#        (Lleva & Jacildo, UPLB CMSC 190)
#
# Tools installed:
#   Pangenome construction : minigraph, minigraph-cactus (cactus), pggb
#   SV calling             : syri, vg, swave
#   Alignment              : minimap2
#   Read simulation        : visor, art_modern
#   SV evaluation          : truvari, manta
#   Sequence utilities     : samtools, bcftools
#   R packages             : regioneR
#
# Prerequisites:
#   - conda / mamba  (for most tools)
#   - pip            (Python packages)
#   - git
#   - wget / curl
#
# Usage:
#   make all             – install everything
#   make pangenome       – install minigraph, cactus, pggb
#   make sv-callers      – install syri, vg, swave
#   make alignment       – install minimap2
#   make utilities       – install samtools, bcftools, truvari, visor
#   make simulation      – install visor, art_modern
#   make evaluation      – install truvari, manta, regioneR
#   make check           – verify installed tools
#   make clean           – remove downloaded build artefacts
# =============================================================================

# ── Versions ──────────────────────────────────────────────────────────────────
MINIGRAPH_VER        := 0.21
MINIMAP2_VER         := 2.30
SAMTOOLS_VER         := 1.23.1
BCFTOOLS_VER         := 1.23.1
VISOR_VER            := 1.1.2.1
ART_MODERN_VER       := 1.4.0
TRUVARI_VER          := 5.4.0
MANTA_VER            := 1.6.0
VG_VER               := 1.73.0
SYRI_VER             := 1.6
SWAVE_VER            := 1.4
PGGB_VER             := 0.7.4
CACTUS_VER           := 3.1.4

# ── Paths ─────────────────────────────────────────────────────────────────────
PREFIX      ?= $(HOME)/clleva/tools
BIN         := $(PREFIX)/bin
SRC         := $(PREFIX)/src

# All tools are installed under ~/clleva/tools/{bin,src}.
# These directories are created automatically by the 'dirs' target.
# Override the prefix at the command line if needed:
#   make all PREFIX=/some/other/path

# ── conda / mamba ─────────────────────────────────────────────────────────────
# make runs non-interactive shells that do NOT source .bashrc, so conda is
# never added to PATH automatically.  Point directly to the miniforge3 binaries.
# miniforge3 was installed at ~/clleva/miniforge3 on this server.
# Override at the command line if your install path differs:
#   make all MINIFORGE_ROOT=/path/to/miniforge3
MINIFORGE_ROOT ?= $(HOME)/clleva/miniforge3
CONDA          := $(MINIFORGE_ROOT)/bin/conda
MAMBA          := $(MINIFORGE_ROOT)/bin/mamba
PIP            := $(MINIFORGE_ROOT)/bin/pip

SHELL := /bin/bash

.PHONY: all pangenome sv-callers alignment utilities simulation evaluation check clean \
        minigraph cactus pggb minimap2 samtools bcftools visor art_modern \
        truvari manta vg syri swave regioneR dirs

# =============================================================================
# Top-level targets
# =============================================================================
all: dirs pangenome sv-callers alignment utilities simulation evaluation
	@echo ""
	@echo "✔  All tools installed. Run 'make check' to verify."

pangenome:  minigraph cactus pggb
sv-callers: syri vg swave
alignment:  minimap2
utilities:  samtools bcftools
simulation: visor art_modern
evaluation: truvari manta regioneR

dirs:
	@mkdir -p $(BIN) $(SRC)

# =============================================================================
# Pangenome construction tools
# =============================================================================

# ── minigraph ─────────────────────────────────────────────────────────────────
minigraph: | dirs
minigraph: $(BIN)/minigraph
$(BIN)/minigraph: | dirs
	@echo ">>> Installing minigraph v$(MINIGRAPH_VER)..."
	cd $(SRC) && \
	  wget -q https://github.com/lh3/minigraph/releases/download/v$(MINIGRAPH_VER)/minigraph-$(MINIGRAPH_VER).tar.bz2 && \
	  tar -xjf minigraph-$(MINIGRAPH_VER).tar.bz2 && \
	  cd minigraph-$(MINIGRAPH_VER) && \
	  make -j$$(nproc) && \
	  cp minigraph $(BIN)/
	@echo "✔  minigraph installed."

# ── minigraph-cactus (cactus-pangenome) ───────────────────────────────────────
cactus: | dirs
cactus: $(BIN)/cactus-pangenome
$(BIN)/cactus-pangenome: | dirs
	@echo ">>> Installing minigraph-cactus v$(CACTUS_VER) via conda..."
	$(MAMBA) create -y -n cactus -c conda-forge -c bioconda \
	    cactus=$(CACTUS_VER) python=3.10 || true
	CACTUS_BIN=$$($(CONDA) run -n cactus which cactus-pangenome 2>/dev/null) && \
	  if [ -n "$$CACTUS_BIN" ]; then ln -sf $$CACTUS_BIN $(BIN)/cactus-pangenome; fi
	@echo "✔  cactus-pangenome installed (conda env: cactus)."
	@echo "   Activate with: $(CONDA) activate cactus"

# ── pggb ──────────────────────────────────────────────────────────────────────
pggb: | dirs
pggb: $(BIN)/pggb
$(BIN)/pggb: | dirs
	@echo ">>> Installing pggb v$(PGGB_VER) via conda..."
	$(MAMBA) create -y -n pggb -c conda-forge -c bioconda \
	    pggb=$(PGGB_VER) python=3.10 || true
	PGGB_BIN=$$($(CONDA) run -n pggb which pggb 2>/dev/null) && \
	  if [ -n "$$PGGB_BIN" ]; then ln -sf $$PGGB_BIN $(BIN)/pggb; fi
	@echo "✔  pggb installed (conda env: pggb)."
	@echo "   Activate with: $(CONDA) activate pggb"

# =============================================================================
# Alignment tools
# =============================================================================

# ── minimap2 ──────────────────────────────────────────────────────────────────
minimap2: | dirs
minimap2: $(BIN)/minimap2
$(BIN)/minimap2: | dirs
	@echo ">>> Installing minimap2 v$(MINIMAP2_VER)..."
	cd $(SRC) && \
	  wget -q https://github.com/lh3/minimap2/releases/download/v$(MINIMAP2_VER)/minimap2-$(MINIMAP2_VER)_x64-linux.tar.bz2 && \
	  tar -xjf minimap2-$(MINIMAP2_VER)_x64-linux.tar.bz2 && \
	  cp minimap2-$(MINIMAP2_VER)_x64-linux/minimap2 $(BIN)/
	@echo "✔  minimap2 installed."

# =============================================================================
# Sequence utilities
# =============================================================================

# ── samtools ──────────────────────────────────────────────────────────────────
samtools: | dirs
samtools: $(BIN)/samtools
$(BIN)/samtools: | dirs
	@echo ">>> Installing samtools v$(SAMTOOLS_VER)..."
	cd $(SRC) && \
	  wget -q https://github.com/samtools/samtools/releases/download/$(SAMTOOLS_VER)/samtools-$(SAMTOOLS_VER).tar.bz2 && \
	  tar -xjf samtools-$(SAMTOOLS_VER).tar.bz2 && \
	  cd samtools-$(SAMTOOLS_VER) && \
	  ./configure --prefix=$(PREFIX) && \
	  make -j$$(nproc) && \
	  make install
	@echo "✔  samtools installed."

# ── bcftools ──────────────────────────────────────────────────────────────────
bcftools: | dirs
bcftools: $(BIN)/bcftools
$(BIN)/bcftools: | dirs
	@echo ">>> Installing bcftools v$(BCFTOOLS_VER)..."
	cd $(SRC) && \
	  wget -q https://github.com/samtools/bcftools/releases/download/$(BCFTOOLS_VER)/bcftools-$(BCFTOOLS_VER).tar.bz2 && \
	  tar -xjf bcftools-$(BCFTOOLS_VER).tar.bz2 && \
	  cd bcftools-$(BCFTOOLS_VER) && \
	  ./configure --prefix=$(PREFIX) && \
	  make -j$$(nproc) && \
	  make install
	@echo "✔  bcftools installed."

# =============================================================================
# SV callers and pangenome-based variant tools
# =============================================================================

# ── vg toolkit ────────────────────────────────────────────────────────────────
vg: | dirs
vg: $(BIN)/vg
$(BIN)/vg: | dirs
	@echo ">>> Installing vg v$(VG_VER) (static binary)..."
	wget -q -O $(BIN)/vg \
	  https://github.com/vgteam/vg/releases/download/v$(VG_VER)/vg
	chmod +x $(BIN)/vg
	@echo "✔  vg installed."

# ── SyRi ──────────────────────────────────────────────────────────────────────
syri: | dirs
syni: $(BIN)/syri
$(BIN)/syri: | dirs
	@echo ">>> Installing syri v$(SYRI_VER) via pip..."
	$(PIP) install syri==$(SYRI_VER)
	SYRI_BIN=$$(which syri 2>/dev/null) && \
	  if [ -n "$$SYRI_BIN" ] && [ "$(BIN)/syri" != "$$SYRI_BIN" ]; then \
	    ln -sf $$SYRI_BIN $(BIN)/syri; \
	  fi
	@echo "✔  syri installed."

# ── SWave ─────────────────────────────────────────────────────────────────────
swave: | dirs
swave: $(BIN)/swave
$(BIN)/swave: | dirs
	@echo ">>> Installing SWave v$(SWAVE_VER) via conda..."
	$(MAMBA) create -y -n swave -c conda-forge -c bioconda \
	    swave=$(SWAVE_VER) python=3.10 || true
	SWAVE_BIN=$$($(CONDA) run -n swave which swave 2>/dev/null) && \
	  if [ -n "$$SWAVE_BIN" ]; then ln -sf $$SWAVE_BIN $(BIN)/swave; fi
	@echo "✔  SWave installed (conda env: swave)."
	@echo "   Activate with: $(CONDA) activate swave"

# =============================================================================
# Simulation tools
# =============================================================================

# ── VISOR ─────────────────────────────────────────────────────────────────────
visor: | dirs
visor: $(BIN)/VISOR
$(BIN)/VISOR: | dirs
	@echo ">>> Installing VISOR v$(VISOR_VER) via pip..."
	$(PIP) install VISOR==$(VISOR_VER)
	VISOR_BIN=$$(which VISOR 2>/dev/null) && \
	  if [ -n "$$VISOR_BIN" ] && [ "$(BIN)/VISOR" != "$$VISOR_BIN" ]; then \
	    ln -sf $$VISOR_BIN $(BIN)/VISOR; \
	  fi
	@echo "✔  VISOR installed."

# ── ART Modern (art_modern) ────────────────────────────────────────────────────
art_modern: | dirs
art_modern: $(BIN)/art_modern
$(BIN)/art_modern: | dirs
	@echo ">>> Installing art_modern v$(ART_MODERN_VER) via conda..."
	$(MAMBA) create -y -n art_modern -c conda-forge -c bioconda \
	    art=$(ART_MODERN_VER) || true
	ART_BIN=$$($(CONDA) run -n art_modern which art_illumina 2>/dev/null) && \
	  if [ -n "$$ART_BIN" ]; then ln -sf $$ART_BIN $(BIN)/art_modern; fi
	@echo "✔  art_modern installed (conda env: art_modern)."
	@echo "   Activate with: $(CONDA) activate art_modern"

# =============================================================================
# Evaluation tools
# =============================================================================

# ── Truvari ───────────────────────────────────────────────────────────────────
truvari: | dirs
truvari: $(BIN)/truvari
$(BIN)/truvari: | dirs
	@echo ">>> Installing Truvari v$(TRUVARI_VER) via pip..."
	$(PIP) install truvari==$(TRUVARI_VER)
	TRUVARI_BIN=$$(which truvari 2>/dev/null) && \
	  if [ -n "$$TRUVARI_BIN" ] && [ "$(BIN)/truvari" != "$$TRUVARI_BIN" ]; then \
	    ln -sf $$TRUVARI_BIN $(BIN)/truvari; \
	  fi
	@echo "✔  Truvari installed."

# ── Manta ─────────────────────────────────────────────────────────────────────
manta: | dirs
manta: $(BIN)/manta
$(BIN)/manta: | dirs
	@echo ">>> Installing Manta v$(MANTA_VER) via conda..."
	$(MAMBA) create -y -n manta -c conda-forge -c bioconda \
	    manta=$(MANTA_VER) || true
	MANTA_BIN=$$($(CONDA) run -n manta which configManta.py 2>/dev/null) && \
	  if [ -n "$$MANTA_BIN" ]; then ln -sf $$MANTA_BIN $(BIN)/manta; fi
	@echo "✔  Manta installed (conda env: manta)."
	@echo "   Activate with: $(CONDA) activate manta"

# ── regioneR (R package) ──────────────────────────────────────────────────────
regioneR:
	@echo ">>> Installing regioneR v$(strip $(subst v,,$(REGIONR_VER))) via R/Bioconductor..."
	Rscript -e "if (!requireNamespace('BiocManager', quietly=TRUE)) install.packages('BiocManager', repos='https://cloud.r-project.org'); BiocManager::install('regioneR', ask=FALSE, update=FALSE)"
	@echo "✔  regioneR installed."

# =============================================================================
# Verification
# =============================================================================
check:
	@echo "=== Tool versions ==="
	@echo -n "minigraph      : "; $(BIN)/minigraph --version 2>&1 | head -1 || echo "NOT FOUND"
	@echo -n "cactus-pg      : "; $(CONDA) run -n cactus cactus-pangenome --version 2>&1 | head -1 || echo "NOT FOUND"
	@echo -n "pggb           : "; $(CONDA) run -n pggb pggb --version 2>&1 | head -1 || echo "NOT FOUND"
	@echo -n "minimap2       : "; $(BIN)/minimap2 --version 2>&1 | head -1 || echo "NOT FOUND"
	@echo -n "samtools       : "; $(BIN)/samtools --version 2>&1 | head -1 || echo "NOT FOUND"
	@echo -n "bcftools       : "; $(BIN)/bcftools --version 2>&1 | head -1 || echo "NOT FOUND"
	@echo -n "vg             : "; $(BIN)/vg version 2>&1 | head -1 || echo "NOT FOUND"
	@echo -n "syri           : "; syri --version 2>&1 | head -1 || echo "NOT FOUND"
	@echo -n "swave          : "; $(CONDA) run -n swave swave --version 2>&1 | head -1 || echo "NOT FOUND"
	@echo -n "VISOR          : "; VISOR --version 2>&1 | head -1 || echo "NOT FOUND"
	@echo -n "truvari        : "; truvari version 2>&1 | head -1 || echo "NOT FOUND"
	@echo -n "manta          : "; $(CONDA) run -n manta configManta.py --version 2>&1 | head -1 || echo "NOT FOUND"
	@echo -n "regioneR (R)   : "; Rscript -e "packageVersion('regioneR')" 2>&1 | head -1 || echo "NOT FOUND"

# =============================================================================
# Clean
# =============================================================================
clean:
	@echo "Removing downloaded source archives..."
	rm -f $(SRC)/minigraph-$(MINIGRAPH_VER).tar.bz2
	rm -f $(SRC)/minimap2-$(MINIMAP2_VER)_x64-linux.tar.bz2
	rm -f $(SRC)/samtools-$(SAMTOOLS_VER).tar.bz2
	rm -f $(SRC)/bcftools-$(BCFTOOLS_VER).tar.bz2
	@echo "✔  Source archives removed (compiled binaries and conda envs kept)."
