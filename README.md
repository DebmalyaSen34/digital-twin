# Mycoplasma genitalium Digital Twin

A whole-cell simulation platform for **Mycoplasma genitalium** (iPS189 metabolic model) that integrates genomic mutation analysis, drug response modeling, and interactive 3D network visualization. Built with COBRApy, Vivarium, Plotly, and Dash.

---

## Overview

This project constructs a computational "digital twin" of *M. genitalium* — the smallest free-living organism — by combining:

1. **Constraint-Based Metabolic Modeling (FBA)** — Flux Balance Analysis on the [iPS189 genome-scale model](src/iPS189.xml) (189 genes, ~260 reactions)
2. **Hybrid Multi-Scale Simulation (Vivarium)** — ODE-based drug diffusion → enzyme binding kinetics → metabolic flux, composed via the [Vivarium framework](https://vivarium-collective.github.io/)
3. **Genomic Mutation Simulation** — Real DNA→protein translation from the [NCBI genome](data/ncbi_dataset/), with SILENT/MISSENSE/NONSENSE classification and physiological penalty modeling
4. **Interactive 3D Network Dashboard** — Force-directed 3D metabolic network with real-time cascade visualization powered by Plotly + Dash

```
┌──────────────────────────────────────────────────────────┐
│               WHOLE-CELL SIMULATION STACK                │
├──────────────────────────────────────────────────────────┤
│                                                          │
│   Genome (FASTA)  ──→  DNA Mutation  ──→  Protein        │
│        │                    │            Impact          │
│        ▼                    ▼               │            │
│   Proteome Map        SILENT / MISSENSE /   │            │
│   (master_map)        NONSENSE              │            │
│        │                    │               │            │
│        ▼                    ▼               ▼            │
│  ┌──────────────────────────────────────────────┐        │
│  │          VIVARIUM COMPOSITE PROCESS          │        │
│  │                                              │        │
│  │  Drug Diffusion (ODE) ──→ Enzyme Kinetics    │        │
│  │       membrane transport    binding/inhibit  │        │
│  │              │                    │          │        │
│  │              ▼                    ▼          │        │
│  │         Dynamic FBA (COBRApy)                │        │
│  │         constrained flux balance             │        │
│  └──────────────────────────────────────────────┘        │
│                         │                                │
│                         ▼                                │
│              3D Network Dashboard (Dash)                 │
│              pathway-grouped visualization               │
│                                                          │
└──────────────────────────────────────────────────────────┘
```

---

## 📁 Project Structure

```
twins/
├── README.md
├── requirements.txt
├── src/
│   ├── iPS189.xml                          # SBML metabolic model (189 genes)
│   ├── config/
│   │   ├── file_paths.py                   # Centralized file paths & constants
│   │   └── drug_db.py                      # Drug target database (MIC, Ki, targets)
│   ├── life/
│   │   ├── v0.py                           # Genesis: genome → proteome → viability
│   │   ├── v1.py                           # Basic mutation simulator
│   │   ├── v2.py                           # Strict evolutionary bottleneck simulator
│   │   ├── v3.py                           # v2 + static 3D HTML visualization
│   │   └── v4.py                           # ⭐ Full interactive Dash dashboard
│   ├── vivarium_model/
│   │   ├── drug_diffusion.py               # ODE membrane transport process
│   │   ├── binding_kinetics.py             # Enzyme inhibition kinetics process
│   │   ├── fba_process.py                  # Dynamic FBA process
│   │   ├── whole_cell_composite.py         # Vivarium composite (orchestrator)
│   │   ├── dashboard.py                    # Drug-only Vivarium dashboard
│   │   └── testing_imports.py              # Import verification utility
│   └── dash_whole_cell_simulation.py       # Earlier 2D Cytoscape dashboard
├── data/
│   ├── master_map.csv                      # Reaction → pathway → metabolite mapping
│   ├── biomass.csv                         # Biomass composition
│   ├── metabolites.csv                     # Metabolite annotations
│   ├── pathway_gene_mapping.csv            # Gene → pathway associations
│   ├── exchange_rxns.csv                   # Exchange reaction bounds
│   ├── rxn_grp.csv                         # Reaction groupings
│   └── ncbi_dataset/                       # M. genitalium genome (GCA_000027325.1)
│       └── data/GCA_000027325.1/
│           └── cds_from_genomic.fna        # CDS FASTA sequences
├── debug/                                  # Generated HTML visualizations & logs
├── lib/                                    # Frontend libraries (vis.js, tom-select)
└── utils/                                  # Utility scripts
```

---

## Getting Started

### Prerequisites

- Python 3.10+
- The *M. genitalium* genome FASTA in `data/ncbi_dataset/`

### Installation

```bash
git clone <repo-url>
cd twins
pip install -r requirements.txt
```

### Key Dependencies

| Package | Purpose |
|---------|---------|
| `cobra` | Constraint-based metabolic modeling (FBA) |
| `vivarium-core` | Multi-scale process composition (optional) |
| `dash` | Interactive web dashboard framework |
| `plotly` | 3D network visualization |
| `biopython` | Genome parsing & protein translation |
| `pandas` / `numpy` | Data handling & numerical computation |

> **Note:** Vivarium is optional. If not installed, the system falls back to pure FBA simulation automatically.

---

## Running the Dashboards

### Full Interactive Dashboard (recommended)

The primary dashboard ([`src/life/v4.py`](src/life/v4.py)) supports drug response, mutation simulation, and combined perturbations:

```bash
python src/life/v4.py
```

Open [http://localhost:8051](http://localhost:8051) in your browser.

**Features:**
- **4 simulation modes:** Drug Response, Random Mutation, Targeted Mutation, Drug + Mutation
- **3D force-directed network** with 1000+ nodes grouped by metabolic pathway
- **Real-time cascade visualization** showing direct and ripple effects
- **Time-course plots** for drug diffusion, enzyme activity, and growth rate
- **Flux comparison charts** (wild-type vs perturbed)
- **Click-to-inspect** any gene, reaction, or metabolite

### Vivarium Drug Dashboard

The Vivarium-specific dashboard ([`src/vivarium_model/dashboard.py`](src/vivarium_model/dashboard.py)) focuses on multi-scale drug simulation:

```bash
python src/vivarium_model/dashboard.py
```

Open [http://localhost:8050](http://localhost:8050).

### Static Visualization (no server)

Generate standalone HTML reports:

```bash
python src/life/v3.py
```

Outputs interactive HTML files to `debug/`:
- `mutation_dashboard.html` — Combined single-page report
- `mutation_network_3d.html` — 3D network only
- `mutation_flux_comparison.html` — Flux bar charts
- `mutation_pathway_impact.html` — Pathway-level analysis

---

## Simulation Modes

### Drug Response

Simulates antibiotic effects on *M. genitalium* metabolism using a 6-drug panel:

| Drug | Target Gene(s) | Mechanism | MIC (µM) |
|------|----------------|-----------|----------|
| Trimethoprim | MG228 | DHFR competitive inhibition | 2.0 |
| Methotrexate | MG228, MG006 | DHFR + thymidylate synthase | 0.5 |
| Fosmidomycin | MG066 | DXR slow tight-binding | 10.0 |
| Cerulenin | MG212, MG114 | Fatty acid synthase (covalent) | 5.0 |
| Mupirocin | MG345 | Isoleucyl-tRNA synthetase | 0.1 |
| Glycolysis Inhibitor | MG041, MG429 | Glucose import (cooperative) | 20.0 |

**When Vivarium is available**, the simulation runs three coupled processes:

1. **[`DrugDiffusion`](src/vivarium_model/drug_diffusion.py)** — Models passive membrane transport using Fick's law with *M. genitalium*-specific geometry (0.3 µm diameter, A/V ≈ 2×10⁵ cm⁻¹)
2. **[`EnzymeKinetics`](src/vivarium_model/binding_kinetics.py)** — Computes enzyme activity using competitive/tight-binding inhibition models with literature Ki values
3. **[`DynamicFBA`](src/vivarium_model/fba_process.py)** — Runs FBA with enzyme activity-constrained reaction bounds

These are composed in [`whole_cell_composite.py`](src/vivarium_model/whole_cell_composite.py).

### Mutation Simulation

1. Selects a metabolically active gene (carrying flux under baseline conditions)
2. Loads the actual CDS DNA sequence from the NCBI genome FASTA
3. Introduces a point mutation (biased 40% toward premature stop codons)
4. Translates mutant DNA → protein using Mycoplasma codon table (table 4)
5. Classifies the impact:
   - **SILENT** — No amino acid change → no metabolic effect
   - **MISSENSE** — Amino acid substitution → enzyme limited to 10% capacity
   - **NONSENSE** — Premature stop codon → complete gene knockout
6. Runs FBA with the constrained/knocked-out gene to compute fitness loss

### Combined Mode

Applies both drug and mutation perturbations simultaneously. Fitness is computed multiplicatively:

$$\text{Combined Fitness} = \frac{\text{Growth}_{\text{drug}}}{\text{Growth}_{\text{WT}}} \times \frac{\text{Growth}_{\text{mutation}}}{\text{Growth}_{\text{WT}}} \times \text{Growth}_{\text{WT}}$$

---

## 3D Network Visualization

The metabolic network is rendered as an interactive 3D scatter plot with force-directed layout:

- **Nodes:** Genes (◆ diamond), Reactions (■ square), Metabolites (● circle)
- **Layout:** Pathway-grouped on a Fibonacci sphere with spring-based relaxation (30 iterations)
- **Hub filtering:** Cofactors connected to >15 reactions are excluded to reduce clutter
- **Cascade highlighting:** Direct targets → 1st-order neighbors → 2nd-order ripple effects

### Node Color Key

| Color | Meaning |
|-------|---------|
| 🟣 `#d946ef` | Target gene (drug target or mutated gene) |
| 🟠 `#f97316` | Directly affected reaction |
| 🟡 `#eab308` | Directly affected metabolite |
| 🔴 `#c2410c` | Cascade reaction (secondary) |
| 🟤 `#a16207` | Cascade metabolite (secondary) |
| ⚫ `#64748b` | Unaffected (dimmed) |

---

## Data Sources

| File | Description |
|------|-------------|
| [`src/iPS189.xml`](src/iPS189.xml) | SBML model from Suthers et al. — 189 genes, ~260 reactions, ~270 metabolites |
| [`data/master_map.csv`](data/master_map.csv) | Reaction-to-pathway-to-metabolite mapping for annotation |
| [`data/ncbi_dataset/`](data/ncbi_dataset/) | NCBI RefSeq genome assembly GCA_000027325.1 |

---

## Architecture

### Evolutionary History (src/life/)

| Version | File | Description |
|---------|------|-------------|
| v0 | [`v0.py`](src/life/v0.py) | Genesis: Genome → Proteome → FBA viability check |
| v1 | [`v1.py`](src/life/v1.py) | Basic mutation: random gene → translate → FBA |
| v2 | [`v2.py`](src/life/v2.py) | Strict evolution: active-gene targeting, flux-based penalties |
| v3 | [`v3.py`](src/life/v3.py) | Added 3D Plotly network + static HTML output |
| v4 | [`v4.py`](src/life/v4.py) | **Full Dash dashboard** with drug + mutation + combined modes |

### Vivarium Processes ([`src/vivarium_model/`](src/vivarium_model/))

| Module | Class | Role |
|--------|-------|------|
| [`drug_diffusion.py`](src/vivarium_model/drug_diffusion.py) | [`DrugDiffusion`](src/vivarium_model/drug_diffusion.py) | ODE membrane transport with literature permeability coefficients |
| [`binding_kinetics.py`](src/vivarium_model/binding_kinetics.py) | [`EnzymeKinetics`](src/vivarium_model/binding_kinetics.py) | Competitive/tight-binding enzyme inhibition with drug-specific Ki |
| [`fba_process.py`](src/vivarium_model/fba_process.py) | [`DynamicFBA`](src/vivarium_model/fba_process.py) | Activity-constrained FBA with dynamic bound scaling |
| [`whole_cell_composite.py`](src/vivarium_model/whole_cell_composite.py) | Composer | Orchestrates all three processes with shared state |
| [`dashboard.py`](src/vivarium_model/dashboard.py) | Dash app | Drug-focused hybrid simulation dashboard (port 8050) |

---

## Configuration

Centralized settings are in [`src/config/file_paths.py`](src/config/file_paths.py):

```python
GENOME_FASTA = "data/ncbi_dataset/data/GCA_000027325.1/cds_from_genomic.fna"
MASTER_MAP_FILE = "data/master_map.csv"
MODEL_FILE = "src/iPS189.xml"
HUB_THRESHOLD = 15      # Metabolites with degree > 15 excluded from network
CLR_BG = "#0a0a1a"      # Dashboard background color
```

Drug targets and pharmacokinetic parameters are in [`src/config/drug_db.py`](src/config/drug_db.py).

---

## Key Algorithms

### Hub Metabolite Filtering

Cofactors like ATP, NAD⁺, H₂O are connected to many reactions and create visual noise. Metabolites with degree > `HUB_THRESHOLD` (15) are excluded from the network graph:

$$\text{Hub} = \{ m \mid \deg(m) > 15 \}$$

### Force-Directed 3D Layout

1. Pathways placed on a **Fibonacci sphere** (golden angle spacing) at radius 80
2. Nodes jittered around pathway centers (genes: σ=12, reactions: σ=18, metabolites: σ=25)
3. **30 iterations** of spring relaxation:
   - Attraction: edges pull connected nodes toward distance 8
   - Repulsion: nodes within distance 30 repel with force ∝ $1/d^2$
   - Anchor: nodes pulled toward pathway center (k=0.005)
   - Damping: linearly decreasing from 0.8 to 0

### Mutation Bias

The [`mutate_dna`](src/life/v4.py) function introduces biologically interesting mutations:
- **40% chance:** Targets a codon in the first 80% of the gene and attempts a single-base change to create a premature stop codon (NONSENSE)
- **60% chance:** Random single-base substitution (typically SILENT or MISSENSE)