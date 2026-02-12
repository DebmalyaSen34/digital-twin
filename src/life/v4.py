import dash
from dash import html, dcc, Input, Output, State, ctx, no_update
import cobra
import pandas as pd
import numpy as np
import hashlib
import random
import plotly.graph_objects as go
import plotly.io as pio
import sys
import os

from Bio import SeqIO
from Bio.Seq import Seq

# Add src to path for process imports
project_root = os.path.join(os.path.dirname(__file__), "../..")
sys.path.insert(0, project_root)

# --- CONFIGURATION ---
GENOME_FASTA = "data/ncbi_dataset/data/GCA_000027325.1/cds_from_genomic.fna"
MODEL_FILE = "src/iPS189.xml"
HUB_THRESHOLD = 15
CLR_BG = "#0a0a1a"

# --- Try loading Vivarium ---
VIVARIUM_AVAILABLE = False
try:
    from src.vivarium_model.whole_cell_composite import run_hybrid_simulation
    from src.vivarium_model.binding_kinetics import DRUG_TARGET_PARAMS
    from src.vivarium_model.drug_diffusion import DRUG_PERMEABILITY
    VIVARIUM_AVAILABLE = True
    print("✅ Vivarium modules loaded.")
except ImportError as e:
    print(f"⚠️ Vivarium unavailable (FBA fallback). {e}")

# =====================================================================
# DRUG DATABASE
# =====================================================================
DRUG_DB = {
    "Control (No Treatment)": {
        "targets": [], "desc": "Baseline metabolism.", "mic_uM": 0, "typical_dose_uM": 0,
    },
    "Trimethoprim": {
        "targets": ["MG228"],
        "desc": "Inhibits DHFR (Folate synthesis). Competitive inhibitor with Ki ≈ 5 nM.",
        "mic_uM": 2.0, "typical_dose_uM": 10.0,
    },
    "Methotrexate": {
        "targets": ["MG228", "MG006"],
        "desc": "Blocks dTMP production. Tight-binding DHFR inhibitor (Ki ≈ 1 nM) + thymidylate synthase.",
        "mic_uM": 0.5, "typical_dose_uM": 5.0,
    },
    "Fosmidomycin": {
        "targets": ["MG066"],
        "desc": "Inhibits DXR in isoprenoid synthesis (MEP pathway). Slow tight-binding inhibitor.",
        "mic_uM": 10.0, "typical_dose_uM": 50.0,
    },
    "Cerulenin": {
        "targets": ["MG212", "MG114"],
        "desc": "Irreversible (covalent) inhibitor of fatty acid synthase. Time-dependent killing.",
        "mic_uM": 5.0, "typical_dose_uM": 20.0,
    },
    "Mupirocin": {
        "targets": ["MG345"],
        "desc": "Inhibits isoleucyl-tRNA synthetase. Competitive with isoleucine (Ki ≈ 20 nM).",
        "mic_uM": 0.1, "typical_dose_uM": 1.0,
    },
    "Generic Glycolysis Inhibitor": {
        "targets": ["MG041", "MG429"],
        "desc": "Blocks glucose import/phosphorylation. Cooperative inhibition.",
        "mic_uM": 20.0, "typical_dose_uM": 100.0,
    },
}

# =====================================================================
# LOAD MODEL & GENOME
# =====================================================================
print("Loading metabolic model...")
model = cobra.io.read_sbml_model(MODEL_FILE)
print(f"Model: {len(model.genes)} genes, {len(model.reactions)} reactions, {len(model.metabolites)} metabolites")

print("Loading genome records...")
try:
    GENOME_RECORDS = list(SeqIO.parse(GENOME_FASTA, "fasta"))
    print(f"Genome: {len(GENOME_RECORDS)} CDS records loaded.")
except Exception as e:
    GENOME_RECORDS = []
    print(f"⚠️ Genome FASTA not found: {e}")

try:
    rxn_df = pd.read_csv("data/master_map.csv")
except Exception:
    rxn_df = pd.DataFrame()

# --- Identify hub metabolites ---
_met_degree = {}
for rxn in model.reactions:
    for met in rxn.metabolites:
        _met_degree[met.id] = _met_degree.get(met.id, 0) + 1
HUB_METABOLITES = {mid for mid, deg in _met_degree.items() if deg > HUB_THRESHOLD}
print(f"Hub metabolites excluded (degree > {HUB_THRESHOLD}): {len(HUB_METABOLITES)}")

# --- Pathway & metabolite name caches ---
_pathway_cache = {}
_met_name_cache = {}
if not rxn_df.empty:
    for _, row in rxn_df.iterrows():
        rn = row.get("Reaction_Name", "")
        if rn and rn not in _pathway_cache:
            _pathway_cache[rn] = row.get("Pathway", "Unknown")
        mid = row.get("Metabolite_ID", "")
        if mid and mid not in _met_name_cache:
            _met_name_cache[mid] = row.get("name", mid)


def get_pathway(rxn_name):
    clean = rxn_name.replace("R_", "").replace("_", " ")
    return _pathway_cache.get(clean, _pathway_cache.get(rxn_name, "Unknown"))


def get_met_name(met_id):
    return _met_name_cache.get(met_id, met_id)


# --- Pre-compute mappings ---
_gene_rxn_map = {}
for gene in model.genes:
    _gene_rxn_map[gene.id] = [rxn.id for rxn in gene.reactions]

_rxn_met_map = {}
_rxn_gene_map = {}
_rxn_pathway_map = {}
for rxn in model.reactions:
    _rxn_met_map[rxn.id] = [m.id for m in rxn.metabolites if m.id not in HUB_METABOLITES]
    _rxn_gene_map[rxn.id] = [g.id for g in rxn.genes]
    _rxn_pathway_map[rxn.id] = get_pathway(rxn.name)

_met_rxn_map = {}
for rxn_id, mets in _rxn_met_map.items():
    for mid in mets:
        if mid not in _met_rxn_map:
            _met_rxn_map[mid] = []
        _met_rxn_map[mid].append(rxn_id)

# --- Find active genes for mutation targets ---
print("Computing baseline FBA for active gene list...")
_baseline_model = model.copy()
_baseline_model.objective = "Biomass"
for rxn in _baseline_model.exchanges:
    rxn.lower_bound = -10.0
_baseline_sol = _baseline_model.optimize()
_baseline_growth = _baseline_sol.objective_value
_baseline_fluxes = _baseline_sol.fluxes.abs()

ACTIVE_GENES = []
for gene in _baseline_model.genes:
    for rxn in gene.reactions:
        if _baseline_fluxes[rxn.id] > 0.0001:
            ACTIVE_GENES.append(gene.id)
            break

print(f"Active genes (mutation targets): {len(ACTIVE_GENES)}")

# =====================================================================
# 3D LAYOUT (computed once at startup)
# =====================================================================

def _stable_hash(s, mod=2**31):
    return int(hashlib.md5(s.encode()).hexdigest(), 16) % mod


def compute_3d_layout():
    """Pathway-grouped force-directed 3D layout."""
    np.random.seed(42)

    node_ids = []
    node_types = {}
    node_pathway = {}

    for gene in model.genes:
        node_ids.append(gene.id)
        node_types[gene.id] = "gene"
        rxn_ids = _gene_rxn_map.get(gene.id, [])
        node_pathway[gene.id] = _rxn_pathway_map.get(rxn_ids[0], "Unknown") if rxn_ids else "Unknown"

    seen_mets = set()
    for rxn in model.reactions:
        if not rxn.genes:
            continue
        node_ids.append(rxn.id)
        node_types[rxn.id] = "reaction"
        node_pathway[rxn.id] = _rxn_pathway_map.get(rxn.id, "Unknown")

        for met in rxn.metabolites:
            if met.id in HUB_METABOLITES or met.id in seen_mets:
                continue
            node_ids.append(met.id)
            node_types[met.id] = "metabolite"
            first_rxn = _met_rxn_map.get(met.id, [None])[0]
            node_pathway[met.id] = _rxn_pathway_map.get(first_rxn, "Unknown") if first_rxn else "Unknown"
            seen_mets.add(met.id)

    # Fibonacci sphere for pathway centers
    pathways = sorted(set(node_pathway.values()))
    pw_centers = {}
    n_pw = len(pathways)
    golden_angle = np.pi * (3.0 - np.sqrt(5.0))
    sphere_radius = 80.0

    for i, pw in enumerate(pathways):
        y = 1.0 - (2.0 * i / max(n_pw - 1, 1))
        r = np.sqrt(max(0, 1.0 - y * y))
        theta = golden_angle * i
        pw_centers[pw] = np.array([
            sphere_radius * r * np.cos(theta),
            sphere_radius * y,
            sphere_radius * r * np.sin(theta),
        ])

    # Place nodes with jitter
    positions = {}
    for nid in node_ids:
        pw = node_pathway.get(nid, "Unknown")
        center = pw_centers.get(pw, np.zeros(3))
        nt = node_types.get(nid, "metabolite")
        np.random.seed(_stable_hash(nid) % (2**31))
        spread = {"gene": 12.0, "reaction": 18.0}.get(nt, 25.0)
        positions[nid] = center + np.random.randn(3) * spread

    # Spring-based relaxation
    node_list = list(positions.keys())
    pos_array = np.array([positions[n] for n in node_list])
    node_index = {n: i for i, n in enumerate(node_list)}

    edge_pairs = []
    for rxn in model.reactions:
        if rxn.id not in node_index:
            continue
        ri = node_index[rxn.id]
        for gene in rxn.genes:
            if gene.id in node_index:
                edge_pairs.append((node_index[gene.id], ri))
        for met in rxn.metabolites:
            if met.id in node_index:
                edge_pairs.append((ri, node_index[met.id]))

    n_nodes = len(node_list)
    for iteration in range(30):
        forces = np.zeros_like(pos_array)

        for i, j in edge_pairs:
            diff = pos_array[j] - pos_array[i]
            dist = np.linalg.norm(diff) + 1e-6
            force = diff / dist * (dist - 8.0) * 0.01
            forces[i] += force
            forces[j] -= force

        if n_nodes < 500:
            for i in range(n_nodes):
                for j in range(i + 1, n_nodes):
                    diff = pos_array[j] - pos_array[i]
                    dist_sq = np.dot(diff, diff) + 1e-6
                    if dist_sq < 900:
                        force = diff / np.sqrt(dist_sq) * (-50.0 / dist_sq)
                        forces[i] += force
                        forces[j] -= force
        else:
            n_samples = min(n_nodes * 10, 50000)
            idx_i = np.random.randint(0, n_nodes, n_samples)
            idx_j = np.random.randint(0, n_nodes, n_samples)
            for si in range(n_samples):
                i, j = idx_i[si], idx_j[si]
                if i == j:
                    continue
                diff = pos_array[j] - pos_array[i]
                dist_sq = np.dot(diff, diff) + 1e-6
                if dist_sq < 900:
                    force = diff / np.sqrt(dist_sq) * (-50.0 / dist_sq)
                    forces[i] += force
                    forces[j] -= force

        for idx, nid in enumerate(node_list):
            pw = node_pathway.get(nid, "Unknown")
            center = pw_centers.get(pw, np.zeros(3))
            forces[idx] += (center - pos_array[idx]) * 0.005

        pos_array += forces * 0.8 * (1.0 - iteration / 30.0)

    for idx, nid in enumerate(node_list):
        positions[nid] = pos_array[idx]

    return positions, node_types, node_pathway, pw_centers


print("Computing 3D layout (one-time)...")
NODE_POS, NODE_TYPES, NODE_PATHWAYS, PW_CENTERS = compute_3d_layout()
print(f"3D layout: {len(NODE_POS)} nodes")

# =====================================================================
# SIMULATION BACKENDS
# =====================================================================

_sim_cache = {}


def run_fba_legacy(drug_name, efficacy):
    """Pure-FBA drug simulation — mirrors dashboard.py logic exactly."""
    m = model.copy()
    m.objective = "Biomass"
    for rxn in m.exchanges:
        rxn.lower_bound = -10.0

    sol_wt = m.optimize()
    wt_growth = sol_wt.objective_value if sol_wt.status == "optimal" else 0.0
    wt_fluxes = {r.id: float(sol_wt.fluxes[r.id]) for r in m.reactions} if sol_wt.status == "optimal" else {r.id: 0.0 for r in m.reactions}

    drug_data = DRUG_DB[drug_name]

    # Apply drug effect on a fresh copy (same as dashboard.py)
    m2 = model.copy()
    m2.objective = "Biomass"
    for rxn in m2.exchanges:
        rxn.lower_bound = -10.0

    for tid in drug_data["targets"]:
        if tid in _gene_rxn_map:
            for rxn_id in _gene_rxn_map[tid]:
                try:
                    rxn = m2.reactions.get_by_id(rxn_id)
                except KeyError:
                    continue
                lo, hi = rxn.bounds
                rxn.bounds = (lo * (1.0 - efficacy), hi * (1.0 - efficacy))

    sol_drug = m2.optimize()
    if sol_drug.status == "optimal":
        drug_growth = sol_drug.objective_value
        drug_fluxes = {r.id: float(sol_drug.fluxes[r.id]) for r in m2.reactions}
    else:
        drug_growth = 0.0
        drug_fluxes = {r.id: 0.0 for r in m2.reactions}

    loss = (wt_growth - drug_growth) / wt_growth * 100 if wt_growth > 0 else 0.0

    return {
        "wt_growth": float(wt_growth), "drug_growth": float(drug_growth),
        "wt_fluxes": wt_fluxes, "drug_fluxes": drug_fluxes,
        "time": [0.0], "drug_intracellular": [0.0],
        "enzyme_activities": {}, "growth_rate": [float(drug_growth)],
        "simulation_mode": "legacy_fba",
        "loss_pct": float(loss),
    }


def run_drug_simulation(drug_name, concentration_uM, sim_time=300.0):
    """Run drug simulation with caching."""
    key = ("drug", drug_name, round(concentration_uM, 2), round(sim_time, 1))
    if key not in _sim_cache:
        if VIVARIUM_AVAILABLE and drug_name != "Control (No Treatment)":
            try:
                result = run_hybrid_simulation(
                    drug_name=drug_name,
                    extracellular_concentration=concentration_uM,
                    total_time=sim_time,
                    emit_step=max(sim_time / 30.0, 1.0),
                )
                result["simulation_mode"] = "vivarium_hybrid"

                # Vivarium returns growth_rate as a time series — extract scalar values
                # wt_growth comes from baseline FBA, drug_growth from final time point
                if "wt_growth" not in result or result.get("wt_growth", 0) == 0:
                    # Compute WT growth from a clean FBA
                    m_wt = model.copy()
                    m_wt.objective = "Biomass"
                    for rxn in m_wt.exchanges:
                        rxn.lower_bound = -10.0
                    sol_wt = m_wt.optimize()
                    result["wt_growth"] = float(sol_wt.objective_value) if sol_wt.status == "optimal" else 0.0

                if "drug_growth" not in result or result.get("drug_growth", 0) == 0:
                    # Use last growth rate from time series if available
                    gr = result.get("growth_rate", [])
                    if gr and len(gr) > 0:
                        result["drug_growth"] = float(gr[-1])

                wt = result.get("wt_growth", 0)
                dg = result.get("drug_growth", 0)
                result["loss_pct"] = (wt - dg) / wt * 100 if wt > 0 else 0.0

                # Also ensure wt_fluxes / drug_fluxes exist for the flux comparison chart
                if "wt_fluxes" not in result or not result["wt_fluxes"]:
                    m_wt = model.copy()
                    m_wt.objective = "Biomass"
                    for rxn in m_wt.exchanges:
                        rxn.lower_bound = -10.0
                    sol_wt = m_wt.optimize()
                    if sol_wt.status == "optimal":
                        result["wt_fluxes"] = {r.id: float(sol_wt.fluxes[r.id]) for r in m_wt.reactions}
                    else:
                        result["wt_fluxes"] = {}

                if "drug_fluxes" not in result or not result["drug_fluxes"]:
                    # Run FBA with the drug efficacy to get flux snapshot
                    mic = DRUG_DB[drug_name].get("mic_uM", 10.0)
                    eff = concentration_uM / (mic + concentration_uM) if mic > 0 else 0.0
                    fba_result = run_fba_legacy(drug_name, eff)
                    result["drug_fluxes"] = fba_result["drug_fluxes"]

            except Exception as e:
                print(f"⚠️ Vivarium simulation failed, falling back to FBA: {e}")
                mic = DRUG_DB[drug_name].get("mic_uM", 10.0)
                efficacy = concentration_uM / (mic + concentration_uM) if mic > 0 and concentration_uM > 0 else 0.0
                result = run_fba_legacy(drug_name, efficacy)
        else:
            mic = DRUG_DB[drug_name].get("mic_uM", 10.0)
            efficacy = concentration_uM / (mic + concentration_uM) if mic > 0 and concentration_uM > 0 else 0.0
            result = run_fba_legacy(drug_name, efficacy)
        _sim_cache[key] = result
        if len(_sim_cache) > 50:
            _sim_cache.pop(next(iter(_sim_cache)))
    return _sim_cache[key]


def mutate_dna(dna_seq):
    """Introduce a single point mutation."""
    bases = ['A', 'T', 'C', 'G']
    seq_list = list(dna_seq)
    pos = random.randint(0, len(seq_list) - 1)
    original = seq_list[pos]
    new_base = random.choice([b for b in bases if b != original])
    seq_list[pos] = new_base
    return Seq("".join(seq_list)), f"Position {pos}: {original} → {new_base}"


def analyze_mutation_impact(original_prot, mutated_prot):
    if original_prot == mutated_prot:
        return "SILENT"
    elif "*" in str(mutated_prot)[:-1]:
        return "NONSENSE"
    return "MISSENSE"


def run_mutation_simulation(gene_id):
    """Run mutation simulation on a specific gene, return results dict."""
    key = ("mutation", gene_id)
    if key in _sim_cache:
        return _sim_cache[key]

    m = model.copy()
    m.objective = "Biomass"
    for rxn in m.exchanges:
        rxn.lower_bound = -10.0

    sol_wt = m.optimize()
    wt_growth = sol_wt.objective_value if sol_wt.status == "optimal" else 0.0
    wt_fluxes = {r.id: float(sol_wt.fluxes[r.id]) for r in m.reactions} if sol_wt.status == "optimal" else {r.id: 0.0 for r in m.reactions}
    abs_fluxes = sol_wt.fluxes.abs() if sol_wt.status == "optimal" else None

    # Find genome record
    cleaned_id = gene_id[:2] + "_" + gene_id[2:] if not gene_id.startswith("MG_") else gene_id
    candidate_record = None
    for rec in GENOME_RECORDS:
        if cleaned_id in rec.description or gene_id in rec.description:
            candidate_record = rec
            break

    impact = "SILENT"
    mutation_log = "No FASTA record found"
    if candidate_record:
        original_dna = candidate_record.seq
        mutated_dna, mutation_log = mutate_dna(original_dna)
        orig_prot = original_dna.translate(table=4, to_stop=True)
        mut_prot = mutated_dna.translate(table=4, to_stop=True)
        impact = analyze_mutation_impact(orig_prot, mut_prot)

    # Apply penalty on a fresh copy
    m2 = model.copy()
    m2.objective = "Biomass"
    for rxn in m2.exchanges:
        rxn.lower_bound = -10.0

    try:
        target_gene = m2.genes.get_by_id(gene_id)
    except KeyError:
        target_gene = None

    if target_gene and abs_fluxes is not None:
        if impact == "MISSENSE":
            max_flux = 0
            for rxn in target_gene.reactions:
                if abs_fluxes[rxn.id] > max_flux:
                    max_flux = abs_fluxes[rxn.id]
            new_limit = max(max_flux * 0.1, 1e-9)
            for rxn in target_gene.reactions:
                if rxn.upper_bound > 0:
                    rxn.upper_bound = min(rxn.upper_bound, new_limit)
                if rxn.lower_bound < 0:
                    rxn.lower_bound = max(rxn.lower_bound, -new_limit)
        elif impact == "NONSENSE":
            target_gene.knock_out()

    mut_sol = m2.optimize()
    if mut_sol.status == "optimal":
        new_growth = mut_sol.objective_value
        mut_fluxes = {r.id: float(mut_sol.fluxes[r.id]) for r in m2.reactions}
    else:
        new_growth = 0.0
        mut_fluxes = {r.id: 0.0 for r in m2.reactions}

    loss = (wt_growth - new_growth) / wt_growth * 100 if wt_growth > 0 else 0.0

    result = {
        "wt_growth": float(wt_growth),
        "drug_growth": float(new_growth),
        "wt_fluxes": wt_fluxes,
        "drug_fluxes": mut_fluxes,
        "time": [0.0],
        "drug_intracellular": [0.0],
        "enzyme_activities": {},
        "growth_rate": [float(new_growth)],
        "simulation_mode": "mutation_fba",
        "mutation_impact": impact,
        "mutation_log": mutation_log,
        "loss_pct": float(loss),
        "target_gene": gene_id,
    }
    _sim_cache[key] = result
    if len(_sim_cache) > 50:
        _sim_cache.pop(next(iter(_sim_cache)))
    return result


# =====================================================================
# HIGHLIGHT & CASCADE COMPUTATION
# =====================================================================

def compute_cascade(target_gene_ids, cascade_depth=1):
    """Compute highlight sets from a list of gene targets."""
    highlighted_genes = set(target_gene_ids)
    highlighted_rxns = set()
    highlighted_mets = set()
    secondary_rxns = set()
    secondary_mets = set()
    hit_pathways = set()

    frontier_mets = set()
    for tid in target_gene_ids:
        for rxn_id in _gene_rxn_map.get(tid, []):
            highlighted_rxns.add(rxn_id)
            hit_pathways.add(_rxn_pathway_map.get(rxn_id, "Unknown"))
            for mid in _rxn_met_map.get(rxn_id, []):
                highlighted_mets.add(mid)
                frontier_mets.add(mid)

    visited = set(highlighted_rxns)
    for _ in range(cascade_depth):
        nf = set()
        for mid in frontier_mets:
            for rxn_id in _met_rxn_map.get(mid, []):
                if rxn_id in visited:
                    continue
                secondary_rxns.add(rxn_id)
                visited.add(rxn_id)
                hit_pathways.add(_rxn_pathway_map.get(rxn_id, "Unknown"))
                for m2 in _rxn_met_map.get(rxn_id, []):
                    if m2 not in highlighted_mets:
                        secondary_mets.add(m2)
                    nf.add(m2)
        frontier_mets = nf

    return {
        "highlighted_genes": highlighted_genes,
        "highlighted_rxns": highlighted_rxns,
        "highlighted_mets": highlighted_mets,
        "secondary_rxns": secondary_rxns,
        "secondary_mets": secondary_mets,
        "hit_pathways": hit_pathways,
    }


# =====================================================================
# BUILD 3D FIGURE
# =====================================================================

def build_3d_figure(cascade=None, enzyme_activities=None, title_text=None):
    """Build Plotly 3D scatter for the metabolic network."""
    if cascade is None:
        cascade = {k: set() for k in [
            "highlighted_genes", "highlighted_rxns", "highlighted_mets",
            "secondary_rxns", "secondary_mets", "hit_pathways",
        ]}
    if enzyme_activities is None:
        enzyme_activities = {}

    hg = cascade["highlighted_genes"]
    hr = cascade["highlighted_rxns"]
    hm = cascade["highlighted_mets"]
    sr = cascade["secondary_rxns"]
    sm = cascade["secondary_mets"]
    hp = cascade["hit_pathways"]

    # EDGES
    edge_x, edge_y, edge_z = [], [], []
    edge_hit_x, edge_hit_y, edge_hit_z = [], [], []
    edge_rip_x, edge_rip_y, edge_rip_z = [], [], []

    for rxn in model.reactions:
        if rxn.id not in NODE_POS:
            continue
        rx, ry, rz = NODE_POS[rxn.id]

        for gene in rxn.genes:
            if gene.id not in NODE_POS:
                continue
            gx, gy, gz = NODE_POS[gene.id]
            is_hit = gene.id in hg and rxn.id in hr
            is_rip = rxn.id in sr
            t = (edge_hit_x, edge_hit_y, edge_hit_z) if is_hit else \
                (edge_rip_x, edge_rip_y, edge_rip_z) if is_rip else \
                (edge_x, edge_y, edge_z)
            t[0].extend([gx, rx, None])
            t[1].extend([gy, ry, None])
            t[2].extend([gz, rz, None])

        for met in rxn.metabolites:
            if met.id not in NODE_POS:
                continue
            mx, my, mz = NODE_POS[met.id]
            is_hit = rxn.id in hr and met.id in hm
            is_rip = rxn.id in sr or met.id in sm
            t = (edge_hit_x, edge_hit_y, edge_hit_z) if is_hit else \
                (edge_rip_x, edge_rip_y, edge_rip_z) if is_rip else \
                (edge_x, edge_y, edge_z)
            t[0].extend([rx, mx, None])
            t[1].extend([ry, my, None])
            t[2].extend([rz, mz, None])

    traces = []
    traces.append(go.Scatter3d(
        x=edge_x, y=edge_y, z=edge_z, mode="lines",
        line=dict(color="rgba(30,41,59,0.25)", width=1),
        hoverinfo="skip", showlegend=False,
    ))
    if edge_hit_x:
        traces.append(go.Scatter3d(
            x=edge_hit_x, y=edge_hit_y, z=edge_hit_z, mode="lines",
            line=dict(color="rgba(251,146,60,0.8)", width=3),
            hoverinfo="skip", showlegend=False,
        ))
    if edge_rip_x:
        traces.append(go.Scatter3d(
            x=edge_rip_x, y=edge_rip_y, z=edge_rip_z, mode="lines",
            line=dict(color="rgba(202,138,4,0.4)", width=1.5),
            hoverinfo="skip", showlegend=False,
        ))

    # NODE HELPER
    def _add(ids, color, size, symbol, name, opacity=1.0, border_color=None, border_width=0):
        if not ids:
            return
        xs = [NODE_POS[n][0] for n in ids if n in NODE_POS]
        ys = [NODE_POS[n][1] for n in ids if n in NODE_POS]
        zs = [NODE_POS[n][2] for n in ids if n in NODE_POS]
        labs = [n for n in ids if n in NODE_POS]

        hovers = []
        for nid in labs:
            nt = NODE_TYPES.get(nid, "")
            pw = NODE_PATHWAYS.get(nid, "Unknown")
            if nt == "gene":
                nr = len(_gene_rxn_map.get(nid, []))
                enz = ""
                if nid in enzyme_activities:
                    acts = enzyme_activities[nid]
                    enz = f"<br>Enzyme activity: {acts[-1]:.1%}" if acts else ""
                hovers.append(f"<b>🧬 {nid}</b><br>Gene | Rxns: {nr}<br>Pathway: {pw}{enz}")
            elif nt == "reaction":
                try:
                    rname = model.reactions.get_by_id(nid).name
                except Exception:
                    rname = nid
                hovers.append(f"<b>⚗️ {rname}</b><br>ID: {nid}<br>Pathway: {pw}")
            elif nt == "metabolite":
                mname = get_met_name(nid)
                nc = len(_met_rxn_map.get(nid, []))
                hovers.append(f"<b>🧪 {mname}</b><br>ID: {nid}<br>Connected: {nc}")
            else:
                hovers.append(nid)

        mk = dict(size=size, color=color, opacity=opacity, symbol=symbol)
        if border_color:
            mk["line"] = dict(color=border_color, width=border_width)

        traces.append(go.Scatter3d(
            x=xs, y=ys, z=zs,
            mode="markers+text" if size >= 10 else "markers",
            text=labs if size >= 10 else None,
            textposition="top center" if size >= 10 else None,
            textfont=dict(size=max(7, min(12, size // 2)), color=color) if size >= 10 else None,
            name=name, marker=mk, hovertext=hovers, hoverinfo="text",
            customdata=[[nid, NODE_TYPES.get(nid, "")] for nid in labs],
        ))

    # Categorize
    dg, dr, dm = [], [], []
    tg, hrl, hml = [], [], []
    srl, sml = [], []

    for nid, nt in NODE_TYPES.items():
        if nt == "gene":
            (tg if nid in hg else dg).append(nid)
        elif nt == "reaction":
            (hrl if nid in hr else srl if nid in sr else dr).append(nid)
        elif nt == "metabolite":
            (hml if nid in hm else sml if nid in sm else dm).append(nid)

    _add(dg, "#0ea5e9", 3, "diamond", "Genes (dim)", 0.15)
    _add(dr, "#64748b", 2, "square", "Reactions (dim)", 0.1)
    _add(dm, "#9ca3af", 1.5, "circle", "Metabolites (dim)", 0.08)
    _add(srl, "#c2410c", 5, "square", "Cascade reactions", 0.8)
    _add(sml, "#a16207", 3.5, "circle", "Cascade metabolites", 0.6)
    _add(hml, "#eab308", 7, "circle", "Direct metabolites", 1.0, "#fef08a", 1)
    _add(hrl, "#f97316", 10, "square", "Direct reactions", 1.0, "#ffedd5", 2)
    _add(tg, "#d946ef", 16, "diamond", "🎯 Target genes", 1.0, "#fdf4ff", 3)

    # Pathway annotations
    anns = []
    for pw, c in PW_CENTERS.items():
        if pw == "Unknown":
            continue
        hit = pw in hp
        anns.append(dict(
            x=c[0], y=c[1], z=c[2],
            text=f"<b>{pw}</b>" if hit else pw, showarrow=False,
            font=dict(size=13 if hit else 10, color="#c7d2fe" if hit else "#4a4a7a"),
            bgcolor="rgba(49,46,129,0.6)" if hit else "rgba(20,20,50,0.4)", borderpad=4,
        ))

    ax = dict(showbackground=True, backgroundcolor="#0a0a1a", gridcolor="#111133",
              showticklabels=False, title="", showspikes=False, zeroline=False, showline=False)

    fig = go.Figure(data=traces)
    fig.update_layout(
        title=dict(text=title_text or "", font=dict(size=14, color="#ddd"), x=0.5) if title_text else {},
        scene=dict(xaxis=ax, yaxis=ax, zaxis=ax, bgcolor="#0a0a1a", annotations=anns,
                   camera=dict(eye=dict(x=1.6, y=1.6, z=1.0), up=dict(x=0, y=1, z=0)),
                   aspectmode="data", dragmode="orbit"),
        paper_bgcolor="#0a0a1a", plot_bgcolor="#0a0a1a",
        margin=dict(l=0, r=0, t=40 if title_text else 0, b=0),
        showlegend=True,
        legend=dict(x=0.01, y=0.99, bgcolor="rgba(10,10,26,0.8)",
                    font=dict(color="#888", size=10), bordercolor="#222255", borderwidth=1,
                    itemsizing="constant"),
        hovermode="closest", uirevision="keep",
    )
    return fig


print("Building initial 3D figure...")
INITIAL_FIG = build_3d_figure()
print("Ready.")

# =====================================================================
# DASH APP
# =====================================================================


def _legend_row(color, symbol, label):
    return html.Div(style={"display": "flex", "alignItems": "center", "gap": "6px", "marginBottom": "3px"}, children=[
        html.Span(symbol, style={"color": color, "fontSize": "12px", "width": "14px", "textAlign": "center"}),
        html.Span(label, style={"color": "#666"}),
    ])


def _metric_card(title, value, subtitle=None, accent="#aaa"):
    ch = [
        html.Div(title, style={"color": "#556", "fontSize": "9px", "textTransform": "uppercase", "letterSpacing": "0.5px"}),
        html.Div(value, style={"color": accent, "fontSize": "18px", "fontWeight": "bold", "marginTop": "1px"}),
    ]
    if subtitle:
        ch.append(html.Div(subtitle, style={"color": "#445", "fontSize": "10px", "marginTop": "1px"}))
    return html.Div(style={
        "backgroundColor": "#0d0d28", "borderRadius": "6px",
        "padding": "10px 12px", "border": f"1px solid {accent}15",
    }, children=ch)


app = dash.Dash(__name__, title="iPS189 Whole-Cell Lab — Mutation & Drug Simulator", update_title=None)

app.layout = html.Div(style={
    "backgroundColor": CLR_BG, "minHeight": "100vh", "fontFamily": "'Segoe UI', sans-serif", "overflow": "hidden",
}, children=[
    dcc.Store(id="sim-store", data={}),

    # HEADER
    html.Div(style={
        "background": "linear-gradient(135deg, #0a0a1a 0%, #111144 100%)",
        "padding": "15px 30px", "borderBottom": "1px solid #222255",
        "display": "flex", "alignItems": "center", "justifyContent": "space-between",
    }, children=[
        html.Div([
            html.H1("🧬 iPS189 Whole-Cell Lab", style={"color": "#fff", "margin": 0, "fontSize": "24px"}),
            html.P([
                "Mutation Simulator + Drug Response | 3D Interactive Network",
                html.Span(f"  ({'Vivarium' if VIVARIUM_AVAILABLE else 'FBA fallback'})",
                          style={"color": "#22c55e" if VIVARIUM_AVAILABLE else "#ef4444", "fontWeight": "bold"}),
            ], style={"color": "#556", "margin": "2px 0 0 0", "fontSize": "12px"}),
        ]),
        html.Div(id="header-stats", style={"display": "flex", "gap": "30px", "alignItems": "center"}),
    ]),

    html.Div(style={"display": "flex", "height": "calc(100vh - 70px)"}, children=[

        # LEFT PANEL
        html.Div(style={
            "width": "320px", "flexShrink": 0, "padding": "16px",
            "backgroundColor": "#0d0d24", "overflowY": "auto", "borderRight": "1px solid #1a1a3e",
        }, children=[

            # MODE SELECTOR
            html.Div("⚙️ SIMULATION MODE", style={
                "color": "#818cf8", "fontWeight": "bold", "fontSize": "11px",
                "letterSpacing": "1px", "marginBottom": "8px",
            }),
            dcc.RadioItems(
                id="sim-mode",
                options=[
                    {"label": " 💊 Drug Response", "value": "drug"},
                    {"label": " 🧬 Random Mutation", "value": "mutation_random"},
                    {"label": " 🎯 Targeted Mutation", "value": "mutation_target"},
                    {"label": " 💊+🧬 Drug + Mutation", "value": "combined"},
                ],
                value="drug",
                style={"marginBottom": "16px"},
                inputStyle={"marginRight": "6px", "accentColor": "#818cf8"},
                labelStyle={
                    "display": "block", "marginBottom": "6px", "cursor": "pointer",
                    "color": "#e2e8f0", "fontSize": "12px",
                },
            ),

            html.Hr(style={"borderColor": "#1a1a3e", "margin": "12px 0"}),

            # DRUG CONTROLS
            html.Div(id="drug-controls-section", children=[
                html.H3("💊 Drug Controls", style={"color": "#ddd", "marginTop": 0, "fontSize": "14px"}),
                html.Label("Drug Treatment:", style={"color": "#888", "fontSize": "11px"}),
                dcc.Dropdown(
                    id="drug-select",
                    options=[{"label": k, "value": k} for k in DRUG_DB],
                    value="Control (No Treatment)",
                    style={"marginBottom": "12px"},
                ),
                html.Label("Concentration (µM):", style={"color": "#888", "fontSize": "11px"}),
                html.Div(id="conc-val", style={"color": "#ff6600", "fontSize": "20px", "fontWeight": "bold"}, children="10.0 µM"),
                dcc.Slider(id="conc-slider", min=0, max=200, step=1, value=10,
                           marks={0: {"label": "0", "style": {"color": "#444"}},
                                  50: {"label": "50", "style": {"color": "#444"}},
                                  100: {"label": "100", "style": {"color": "#444"}},
                                  200: {"label": "200 µM", "style": {"color": "#444"}}}),
                html.Label("Sim Time (s):", style={"color": "#888", "fontSize": "11px", "marginTop": "8px"}),
                dcc.Slider(id="time-slider", min=30, max=600, step=30, value=300,
                           marks={30: {"label": "30s", "style": {"color": "#444"}},
                                  300: {"label": "5m", "style": {"color": "#444"}},
                                  600: {"label": "10m", "style": {"color": "#444"}}}),
            ]),

            html.Hr(style={"borderColor": "#1a1a3e", "margin": "12px 0"}),

            # MUTATION CONTROLS
            html.Div(id="mutation-controls-section", children=[
                html.H3("🧬 Mutation Controls", style={"color": "#ddd", "marginTop": 0, "fontSize": "14px"}),
                html.Label("Target Gene:", style={"color": "#888", "fontSize": "11px"}),
                dcc.Dropdown(
                    id="mutation-gene-select",
                    options=[{"label": g, "value": g} for g in sorted(ACTIVE_GENES)],
                    value=ACTIVE_GENES[0] if ACTIVE_GENES else None,
                    placeholder="Select gene to mutate...",
                    style={"marginBottom": "12px"},
                ),
                html.Button("🎲 Roll Random Mutation", id="random-mutation-btn",
                            style={
                                "width": "100%", "padding": "8px", "backgroundColor": "#4c1d95",
                                "color": "#ddd", "border": "1px solid #6d28d9", "borderRadius": "6px",
                                "cursor": "pointer", "fontSize": "12px", "fontWeight": "bold",
                                "marginBottom": "12px",
                            }),
            ]),

            html.Hr(style={"borderColor": "#1a1a3e", "margin": "12px 0"}),

            html.Label("Cascade Depth:", style={"color": "#888", "fontSize": "11px"}),
            dcc.Slider(id="cascade-slider", min=0, max=3, step=1, value=1,
                       marks={i: {"label": str(i), "style": {"color": "#444"}} for i in range(4)}),

            html.Hr(style={"borderColor": "#1a1a3e", "margin": "12px 0"}),

            html.Div(id="info-panel", style={
                "padding": "12px", "backgroundColor": "#111133",
                "borderRadius": "8px", "border": "1px solid #222255", "fontSize": "12px",
            }),

            html.Hr(style={"borderColor": "#1a1a3e", "margin": "12px 0"}),

            html.Div(id="kinetics-panel", style={
                "padding": "12px", "backgroundColor": "#111133",
                "borderRadius": "8px", "border": "1px solid #222255",
            }),

            html.Hr(style={"borderColor": "#1a1a3e", "margin": "12px 0"}),

            html.Div(id="metrics-panel"),

            html.Hr(style={"borderColor": "#1a1a3e", "margin": "12px 0"}),

            # LEGEND
            html.Div(style={"fontSize": "10px", "color": "#555"}, children=[
                html.Div("LEGEND", style={"color": "#444", "fontWeight": "bold", "letterSpacing": "1px", "marginBottom": "6px"}),
                _legend_row("#d946ef", "◆", "Target Gene"),
                _legend_row("#f97316", "■", "Direct Reaction"),
                _legend_row("#eab308", "●", "Direct Metabolite"),
                _legend_row("#c2410c", "■", "Cascade Reaction"),
                _legend_row("#a16207", "●", "Cascade Metabolite"),
                _legend_row("#334155", "●", "Unaffected"),
                html.Div("─" * 20, style={"color": "#333", "margin": "4px 0"}),
                _legend_row("#818cf8", "🖱", "Drag to rotate"),
                _legend_row("#818cf8", "⚙", "Scroll to zoom"),
                _legend_row("#818cf8", "⇧", "Shift+drag to pan"),
            ]),
        ]),

        # CENTER: 3D Network
        html.Div(style={"flex": 1, "position": "relative"}, children=[
            dcc.Graph(
                id="network-3d", figure=INITIAL_FIG,
                style={"width": "100%", "height": "100%"},
                config={"displayModeBar": True, "displaylogo": False, "scrollZoom": True,
                        "modeBarButtonsToRemove": ["toImage", "resetCameraLastSave3d"]},
            ),
            html.Div(id="node-tooltip", style={
                "position": "absolute", "bottom": "12px", "left": "12px",
                "color": "#ddd", "padding": "10px 14px",
                "backgroundColor": "rgba(10,10,30,0.95)", "borderRadius": "8px",
                "border": "1px solid #333366", "fontSize": "12px", "maxWidth": "360px",
                "backdropFilter": "blur(10px)", "pointerEvents": "none",
            }),
            html.Div("🖱 Drag to rotate · Scroll to zoom · Shift+drag to pan", style={
                "position": "absolute", "top": "10px", "right": "15px", "color": "#444",
                "fontSize": "10px", "backgroundColor": "rgba(10,10,26,0.7)",
                "padding": "4px 10px", "borderRadius": "4px",
            }),
        ]),

        # RIGHT: Time-course / Summary
        html.Div(id="timecourse-panel", style={
            "width": "300px", "flexShrink": 0, "padding": "12px",
            "backgroundColor": "#0d0d24", "overflowY": "auto", "borderLeft": "1px solid #1a1a3e",
        }),
    ]),
])


# =====================================================================
# CALLBACKS
# =====================================================================

@app.callback(
    Output("drug-controls-section", "style"),
    Output("mutation-controls-section", "style"),
    Input("sim-mode", "value"),
)
def toggle_controls(mode):
    show = {"display": "block"}
    hide = {"display": "none"}
    if mode == "drug":
        return show, hide
    elif mode in ("mutation_random", "mutation_target"):
        return hide, show
    else:  # combined
        return show, show


@app.callback(
    Output("conc-val", "children"),
    Input("conc-slider", "value"),
)
def update_conc_label(val):
    return f"{val:.1f} µM"


@app.callback(
    Output("mutation-gene-select", "value"),
    Input("random-mutation-btn", "n_clicks"),
    prevent_initial_call=True,
)
def roll_random_gene(_):
    return random.choice(ACTIVE_GENES) if ACTIVE_GENES else no_update


@app.callback(
    Output("sim-store", "data"),
    Input("sim-mode", "value"),
    Input("drug-select", "value"),
    Input("conc-slider", "value"),
    Input("time-slider", "value"),
    Input("mutation-gene-select", "value"),
)
def compute_simulation(mode, drug_name, conc_uM, sim_time, mutation_gene):
    """Run appropriate simulation based on mode."""
    if mode == "drug":
        result = run_drug_simulation(drug_name, float(conc_uM), float(sim_time))
        result["_mode"] = "drug"
        result["_targets"] = DRUG_DB[drug_name]["targets"]
        return result

    elif mode in ("mutation_random", "mutation_target"):
        gene_id = mutation_gene
        if not gene_id:
            gene_id = random.choice(ACTIVE_GENES) if ACTIVE_GENES else None
        if not gene_id:
            return {"_mode": "mutation", "_targets": []}
        result = run_mutation_simulation(gene_id)
        result["_mode"] = "mutation"
        result["_targets"] = [gene_id]
        return result

    elif mode == "combined":
        # Drug simulation
        drug_result = run_drug_simulation(drug_name, float(conc_uM), float(sim_time))
        # Mutation simulation
        gene_id = mutation_gene or (random.choice(ACTIVE_GENES) if ACTIVE_GENES else None)
        if gene_id:
            mut_result = run_mutation_simulation(gene_id)
        else:
            mut_result = {
                "drug_growth": drug_result.get("wt_growth", 0),
                "mutation_impact": "NONE", "mutation_log": "—", "loss_pct": 0.0,
            }

        wt = drug_result.get("wt_growth", 0)
        # Combined fitness: multiply the relative fitness of each perturbation
        drug_fitness = drug_result.get("drug_growth", 0) / wt if wt > 0 else 0
        mut_fitness = mut_result.get("drug_growth", 0) / wt if wt > 0 else 0
        combined_growth = wt * drug_fitness * mut_fitness

        targets = list(DRUG_DB[drug_name]["targets"])
        if gene_id:
            targets.append(gene_id)

        combined_loss = (wt - combined_growth) / wt * 100 if wt > 0 else 0.0

        return {
            "wt_growth": float(wt),
            "drug_growth": float(combined_growth),
            "wt_fluxes": drug_result.get("wt_fluxes", {}),
            "drug_fluxes": drug_result.get("drug_fluxes", {}),
            "time": drug_result.get("time", [0.0]),
            "drug_intracellular": drug_result.get("drug_intracellular", [0.0]),
            "enzyme_activities": drug_result.get("enzyme_activities", {}),
            "growth_rate": drug_result.get("growth_rate", [float(combined_growth)]),
            "simulation_mode": "combined",
            "mutation_impact": mut_result.get("mutation_impact", "NONE"),
            "mutation_log": mut_result.get("mutation_log", "—"),
            "loss_pct": float(combined_loss),
            "_mode": "combined",
            "_targets": targets,
            "target_gene": gene_id,
        }

    return {}


@app.callback(
    Output("network-3d", "figure"),
    Output("info-panel", "children"),
    Output("metrics-panel", "children"),
    Output("header-stats", "children"),
    Output("kinetics-panel", "children"),
    Output("timecourse-panel", "children"),
    Input("cascade-slider", "value"),
    Input("sim-store", "data"),
)
def update_graph(cascade_depth, sim_data):
    mode = sim_data.get("_mode", "drug")
    targets = sim_data.get("_targets", [])
    wt_growth = sim_data.get("wt_growth", 0)
    drug_growth = sim_data.get("drug_growth", 0)
    sim_mode_label = sim_data.get("simulation_mode", "unknown")

    # Cascade
    cascade = compute_cascade(targets, cascade_depth)
    enzyme_acts = sim_data.get("enzyme_activities", {})

    # Title
    if mode == "drug":
        title = None
    elif mode == "mutation":
        impact = sim_data.get("mutation_impact", "?")
        gene = sim_data.get("target_gene", "?")
        loss = sim_data.get("loss_pct", 0)
        title = f"🧬 Mutation: {gene} ({impact}) — {loss:.1f}% fitness loss"
    elif mode == "combined":
        gene = sim_data.get("target_gene", "?")
        impact = sim_data.get("mutation_impact", "?")
        title = f"💊+🧬 Combined: Drug + {gene} ({impact})"
    else:
        title = None

    fig = build_3d_figure(cascade, enzyme_acts, title)

    # INFO PANEL
    info_children = []
    if mode == "drug" or mode == "combined":
        drug_name = ""
        for k, v in DRUG_DB.items():
            if set(v["targets"]) & set(targets):
                drug_name = k
                break
        if not drug_name and targets:
            drug_name = "Custom"
        drug_data = DRUG_DB.get(drug_name, {"desc": "", "targets": [], "mic_uM": 0})
        badge_color = "#22c55e" if "vivarium" in sim_mode_label else "#ef4444"
        info_children.extend([
            html.Div([
                f"💊 {drug_name}",
                html.Span(f" [{sim_mode_label.replace('_', ' ').upper()}]",
                          style={"color": badge_color, "fontSize": "9px", "fontWeight": "bold"}),
            ], style={"fontWeight": "bold", "color": "#ff6600", "fontSize": "14px", "marginBottom": "6px"}),
            html.Div(drug_data["desc"], style={"color": "#999", "fontSize": "11px"}),
        ])
        if drug_data["targets"]:
            info_children.append(html.Div(style={"marginTop": "6px", "display": "flex", "flexWrap": "wrap", "gap": "4px"}, children=[
                html.Span("Drug Targets: ", style={"color": "#777", "fontSize": "11px"}),
                *[html.Span(t, style={"backgroundColor": "#ff0040", "color": "#fff",
                                      "padding": "1px 8px", "borderRadius": "10px", "fontSize": "10px",
                                      "fontWeight": "bold"}) for t in drug_data["targets"]],
            ]))

    if mode in ("mutation", "combined"):
        gene = sim_data.get("target_gene", "?")
        impact = sim_data.get("mutation_impact", "?")
        mutation_log = sim_data.get("mutation_log", "—")
        impact_color = {"SILENT": "#22c55e", "MISSENSE": "#eab308", "NONSENSE": "#ef4444"}.get(impact, "#888")
        info_children.extend([
            html.Div(style={"marginTop": "8px" if mode == "combined" else "0"}, children=[
                html.Div([
                    f"🧬 Mutation: {gene}",
                    html.Span(f" [{impact}]", style={"color": impact_color, "fontWeight": "bold", "fontSize": "11px"}),
                ], style={"fontWeight": "bold", "color": "#a78bfa", "fontSize": "14px", "marginBottom": "4px"}),
                html.Div(mutation_log, style={"color": "#888", "fontSize": "11px", "fontFamily": "monospace"}),
                html.Div(style={"marginTop": "4px"}, children=[
                    html.Span("Effect: ", style={"color": "#777", "fontSize": "11px"}),
                    html.Span(
                        {"SILENT": "No protein change", "MISSENSE": "Enzyme limited to 10% capacity",
                         "NONSENSE": "Gene knocked out (premature stop)"}.get(impact, "Unknown"),
                        style={"color": impact_color, "fontSize": "11px"},
                    ),
                ]),
            ]),
        ])

    hit_pathways = cascade["hit_pathways"] - {"Unknown"}
    if hit_pathways:
        info_children.append(html.Div(style={"marginTop": "8px"}, children=[
            html.Span("Pathways: ", style={"color": "#777", "fontSize": "11px"}),
            html.Span(", ".join(sorted(hit_pathways)), style={"color": "#ff8844", "fontSize": "11px"}),
        ]))

    # KINETICS PANEL
    kinetics_children = []
    if enzyme_acts:
        kinetics_children.append(html.Div("🔬 ENZYME KINETICS", style={
            "color": "#a855f7", "fontWeight": "bold", "fontSize": "11px",
            "letterSpacing": "1px", "marginBottom": "8px",
        }))
        for gid, acts in enzyme_acts.items():
            final = acts[-1] if acts else 1.0
            bw = max(0, min(100, final * 100))
            bc = "#22c55e" if final > 0.7 else ("#eab308" if final > 0.3 else "#ef4444")
            kinetics_children.append(html.Div(style={"marginBottom": "8px"}, children=[
                html.Div(style={"display": "flex", "justifyContent": "space-between"}, children=[
                    html.Span(gid, style={"color": "#ddd", "fontSize": "11px", "fontWeight": "bold"}),
                    html.Span(f"{final:.1%}", style={"color": bc, "fontSize": "11px"}),
                ]),
                html.Div(style={"width": "100%", "height": "6px", "backgroundColor": "#1a1a3e",
                                "borderRadius": "3px", "overflow": "hidden", "marginTop": "2px"}, children=[
                    html.Div(style={"width": f"{bw}%", "height": "100%", "backgroundColor": bc,
                                    "borderRadius": "3px", "transition": "width 0.5s ease"}),
                ]),
            ]))
    elif mode == "mutation":
        impact = sim_data.get("mutation_impact", "?")
        gene = sim_data.get("target_gene", "?")
        if impact == "MISSENSE":
            kinetics_children.append(html.Div(f"🧬 {gene}: 10% capacity (damaged enzyme)",
                                              style={"color": "#eab308", "fontSize": "11px"}))
        elif impact == "NONSENSE":
            kinetics_children.append(html.Div(f"🧬 {gene}: 0% (knocked out)",
                                              style={"color": "#ef4444", "fontSize": "11px"}))
        elif impact == "SILENT":
            kinetics_children.append(html.Div(f"🧬 {gene}: 100% (silent mutation)",
                                              style={"color": "#22c55e", "fontSize": "11px"}))
        else:
            kinetics_children.append(html.Div("No enzyme data", style={"color": "#444", "fontSize": "11px"}))
    else:
        kinetics_children.append(html.Div("No enzyme targets", style={"color": "#444", "fontSize": "11px"}))

    # METRICS
    fitness = (drug_growth / wt_growth * 100) if wt_growth > 0 else 0
    delta_pct = ((drug_growth - wt_growth) / wt_growth * 100) if wt_growth > 0 else 0

    metrics = html.Div(style={"display": "flex", "flexDirection": "column", "gap": "6px"}, children=[
        _metric_card("Wild-Type Growth", f"{wt_growth:.4f}", None, "#44aaff"),
        _metric_card("Perturbed Growth", f"{drug_growth:.4f}",
                     f"{delta_pct:+.1f}%",
                     "#ff4444" if delta_pct < -1 else ("#44ff44" if delta_pct > 1 else "#888")),
        _metric_card("Relative Fitness", f"{fitness:.1f}%", None,
                     "#ff4444" if fitness < 50 else ("#ffaa00" if fitness < 80 else "#44ff44")),
        _metric_card("Impact",
                     f"{len(cascade['highlighted_rxns'])} + {len(cascade['secondary_rxns'])} rxns",
                     f"{len(cascade['highlighted_mets']) + len(cascade['secondary_mets'])} metabolites",
                     "#ff6600"),
    ])

    # HEADER
    fitness_color = "#ff4444" if fitness < 50 else ("#ffaa00" if fitness < 80 else "#44ff44")
    mode_labels = {"drug": "💊 Drug", "mutation": "🧬 Mutation", "combined": "💊+🧬 Combined"}
    header = [
        html.Div([
            html.Div("FITNESS", style={"color": "#556", "fontSize": "9px", "letterSpacing": "1px"}),
            html.Div(f"{fitness:.1f}%", style={"color": fitness_color, "fontSize": "20px", "fontWeight": "bold"}),
        ]),
        html.Div([
            html.Div("IMPACT", style={"color": "#556", "fontSize": "9px", "letterSpacing": "1px"}),
            html.Div(f"{len(cascade['highlighted_rxns']) + len(cascade['secondary_rxns'])} rxns",
                     style={"color": "#ff6600", "fontSize": "14px"}),
        ]),
        html.Div([
            html.Div("MODE", style={"color": "#556", "fontSize": "9px", "letterSpacing": "1px"}),
            html.Div(mode_labels.get(mode, mode), style={"color": "#818cf8", "fontSize": "14px"}),
        ]),
    ]

    # TIMECOURSE PANEL
    timecourse = _build_timecourse(sim_data)

    return fig, info_children, metrics, header, kinetics_children, timecourse


def _build_timecourse(sim_data):
    """Build the right-side time-course / summary panel."""
    mode = sim_data.get("_mode", "drug")
    times = sim_data.get("time", [0.0])
    drug_intra = sim_data.get("drug_intracellular", [0.0])
    enzyme_acts = sim_data.get("enzyme_activities", {})
    growth_rates = sim_data.get("growth_rate", [0.0])
    wt_growth = sim_data.get("wt_growth", 0)
    drug_growth = sim_data.get("drug_growth", 0)

    children = [
        html.Div("📈 ANALYSIS", style={
            "color": "#556", "fontWeight": "bold", "fontSize": "11px",
            "letterSpacing": "1px", "marginBottom": "12px",
        }),
    ]

    has_timeseries = len(times) > 1

    # Drug diffusion plot (only in drug/combined modes with time series)
    if has_timeseries and mode in ("drug", "combined"):
        peak = max(drug_intra) if drug_intra else 0
        children.append(dcc.Graph(
            figure={
                "data": [{"x": times, "y": drug_intra, "type": "scatter", "mode": "lines",
                          "line": {"color": "#22c55e", "width": 2}, "fill": "tozeroy",
                          "fillcolor": "rgba(34,197,94,0.1)", "name": "[Drug]ᵢₙ"}],
                "layout": {
                    "title": {"text": f"Drug Diffusion (peak: {peak:.2f} µM)", "font": {"size": 11, "color": "#888"}},
                    "xaxis": {"title": "Time (s)", "color": "#555", "gridcolor": "#1a1a3e"},
                    "yaxis": {"title": "[Drug] µM", "color": "#555", "gridcolor": "#1a1a3e"},
                    "paper_bgcolor": "rgba(0,0,0,0)", "plot_bgcolor": "rgba(13,13,36,0.8)",
                    "margin": {"l": 50, "r": 10, "t": 30, "b": 35}, "height": 170, "font": {"color": "#888"},
                },
            }, config={"displayModeBar": False}, style={"marginBottom": "8px"},
        ))

    # Enzyme activity plot
    if has_timeseries and enzyme_acts:
        etraces = []
        colors = ["#d946ef", "#f97316", "#eab308", "#22d3ee", "#a855f7"]
        for i, (gid, acts) in enumerate(enzyme_acts.items()):
            fa = acts[-1] if acts else 1.0
            etraces.append({"x": times[:len(acts)], "y": acts, "type": "scatter", "mode": "lines",
                            "line": {"color": colors[i % len(colors)], "width": 2}, "name": f"{gid} ({fa:.0%})"})
        children.append(dcc.Graph(
            figure={
                "data": etraces,
                "layout": {
                    "title": {"text": "Enzyme Activity", "font": {"size": 11, "color": "#888"}},
                    "xaxis": {"title": "Time (s)", "color": "#555", "gridcolor": "#1a1a3e"},
                    "yaxis": {"title": "Activity", "range": [-0.05, 1.05], "color": "#555", "gridcolor": "#1a1a3e"},
                    "paper_bgcolor": "rgba(0,0,0,0)", "plot_bgcolor": "rgba(13,13,36,0.8)",
                    "margin": {"l": 50, "r": 10, "t": 30, "b": 35}, "height": 170, "font": {"color": "#888"},
                    "showlegend": True, "legend": {"font": {"size": 9, "color": "#888"}, "x": 0.55, "y": 0.95},
                    "shapes": [{"type": "line", "x0": times[0], "x1": times[-1], "y0": 0.5, "y1": 0.5,
                                "line": {"color": "#ef4444", "width": 1, "dash": "dot"}}],
                },
            }, config={"displayModeBar": False}, style={"marginBottom": "8px"},
        ))

    # Growth rate plot
    if has_timeseries:
        children.append(dcc.Graph(
            figure={
                "data": [
                    {"x": times[:len(growth_rates)], "y": growth_rates, "type": "scatter", "mode": "lines",
                     "line": {"color": "#3b82f6", "width": 2}, "fill": "tozeroy",
                     "fillcolor": "rgba(59,130,246,0.1)", "name": "Growth"},
                    {"x": [times[0], times[-1]], "y": [wt_growth, wt_growth], "type": "scatter", "mode": "lines",
                     "line": {"color": "#44aaff", "width": 1, "dash": "dash"}, "name": f"WT ({wt_growth:.2f})"},
                ],
                "layout": {
                    "title": {"text": "Growth Rate", "font": {"size": 11, "color": "#888"}},
                    "xaxis": {"title": "Time (s)", "color": "#555", "gridcolor": "#1a1a3e"},
                    "yaxis": {"title": "Growth (h⁻¹)", "color": "#555", "gridcolor": "#1a1a3e"},
                    "paper_bgcolor": "rgba(0,0,0,0)", "plot_bgcolor": "rgba(13,13,36,0.8)",
                    "margin": {"l": 50, "r": 10, "t": 30, "b": 35}, "height": 170, "font": {"color": "#888"},
                    "showlegend": True, "legend": {"font": {"size": 9, "color": "#888"}, "x": 0.55, "y": 0.95},
                },
            }, config={"displayModeBar": False}, style={"marginBottom": "8px"},
        ))

    # Growth gauge (always shown)
    loss = sim_data.get("loss_pct", None)
    if loss is None:
        loss = (wt_growth - drug_growth) / wt_growth * 100 if wt_growth > 0 else 0.0
    loss = max(0.0, min(100.0, float(loss)))

    gauge_max = max(wt_growth * 1.2, drug_growth * 1.2, 0.01)
    wt_ref = float(wt_growth) if wt_growth > 0 else 0.001

    children.append(dcc.Graph(
        figure=go.Figure(go.Indicator(
            mode="gauge+number+delta",
            value=float(drug_growth),
            delta={"reference": wt_ref, "decreasing": {"color": "#ef4444"}},
            gauge={
                "axis": {"range": [0, gauge_max], "tickcolor": "#888"},
                "bar": {"color": "#ef4444" if loss > 50 else ("#eab308" if loss > 10 else "#22c55e")},
                "bgcolor": "#111133", "bordercolor": "#222255",
                "steps": [
                    {"range": [0, gauge_max * 0.4], "color": "rgba(239,68,68,0.2)"},
                    {"range": [gauge_max * 0.4, gauge_max * 0.83], "color": "rgba(234,179,8,0.2)"},
                ],
                "threshold": {"line": {"color": "#3b82f6", "width": 3}, "thickness": 0.8, "value": wt_ref},
            },
            title={"text": f"Fitness (Loss: {loss:.1f}%)", "font": {"size": 12, "color": "#ddd"}},
            number={"font": {"size": 22, "color": "#ddd"}, "suffix": " h⁻¹"},
        )).update_layout(
            paper_bgcolor="rgba(0,0,0,0)", font=dict(color="#888"),
            height=200, margin=dict(l=20, r=20, t=40, b=10),
        ),
        config={"displayModeBar": False}, style={"marginBottom": "8px"},
    ))

    # Flux comparison bar chart (top changed reactions)
    wt_fluxes = sim_data.get("wt_fluxes", {})
    drug_fluxes = sim_data.get("drug_fluxes", {})
    if wt_fluxes and drug_fluxes:
        changes = []
        for rid in set(wt_fluxes.keys()) & set(drug_fluxes.keys()):
            delta = abs(wt_fluxes[rid] - drug_fluxes[rid])
            if delta > 0.001:
                changes.append((rid, wt_fluxes[rid], drug_fluxes[rid], delta))
        changes.sort(key=lambda x: x[3], reverse=True)
        top = changes[:12]

        if top:
            labels = []
            for rid, _, _, _ in top:
                try:
                    labels.append(model.reactions.get_by_id(rid).name[:20])
                except Exception:
                    labels.append(rid[:20])
            children.append(dcc.Graph(
                figure={
                    "data": [
                        {"x": labels, "y": [t[1] for t in top], "type": "bar", "name": "WT",
                         "marker": {"color": "#3b82f6"}, "opacity": 0.8},
                        {"x": labels, "y": [t[2] for t in top], "type": "bar", "name": "Perturbed",
                         "marker": {"color": "#ef4444"}, "opacity": 0.8},
                    ],
                    "layout": {
                        "title": {"text": "Top Changed Fluxes", "font": {"size": 11, "color": "#888"}},
                        "barmode": "group",
                        "xaxis": {"tickangle": -45, "color": "#555", "gridcolor": "#1a1a3e"},
                        "yaxis": {"title": "Flux", "color": "#555", "gridcolor": "#1a1a3e"},
                        "paper_bgcolor": "rgba(0,0,0,0)", "plot_bgcolor": "rgba(13,13,36,0.8)",
                        "margin": {"l": 50, "r": 10, "t": 30, "b": 100}, "height": 280,
                        "font": {"color": "#888"}, "legend": {"font": {"size": 9, "color": "#888"}},
                    },
                }, config={"displayModeBar": False}, style={"marginBottom": "8px"},
            ))

    # Summary box
    sim_mode_label = sim_data.get("simulation_mode", "unknown")
    summary_loss = sim_data.get("loss_pct", None)
    if summary_loss is None:
        summary_loss = (wt_growth - drug_growth) / wt_growth * 100 if wt_growth > 0 else 0.0
    summary_loss = float(summary_loss)

    children.append(html.Div(style={
        "padding": "8px", "backgroundColor": "#111133", "borderRadius": "6px",
        "border": "1px solid #222255", "fontSize": "10px", "color": "#556", "marginTop": "8px",
    }, children=[
        html.Div("DETAILS", style={"fontWeight": "bold", "letterSpacing": "1px", "marginBottom": "4px"}),
        html.Div(f"🔬 Mode: {sim_mode_label}"),
        html.Div(f"📉 WT Growth: {wt_growth:.4f} h⁻¹"),
        html.Div(f"📉 Perturbed: {drug_growth:.4f} h⁻¹"),
        html.Div(f"📊 Fitness Loss: {summary_loss:.1f}%"),
        *([html.Div(f"🧬 Mutation: {sim_data.get('mutation_impact', '?')} — {sim_data.get('mutation_log', '—')}")]
          if mode in ("mutation", "combined") else []),
    ]))

    return html.Div(children)


@app.callback(
    Output("node-tooltip", "children"),
    Input("network-3d", "clickData"),
    State("sim-store", "data"),
)
def show_node_info(click_data, sim_data):
    if not click_data or "points" not in click_data:
        return html.Span("Click any node for details", style={"color": "#444"})

    point = click_data["points"][0]
    custom = point.get("customdata", [])
    if not custom or len(custom) < 2:
        hover = point.get("hovertext", "")
        if hover:
            return html.Div(dcc.Markdown(
                hover.replace("<b>", "**").replace("</b>", "**").replace("<br>", "\n\n")),
                style={"color": "#ddd"})
        return html.Span("Click a node for details", style={"color": "#444"})

    nid, ntype = custom[0], custom[1]
    wt_fluxes = sim_data.get("wt_fluxes", {})
    drug_fluxes = sim_data.get("drug_fluxes", {})
    enzyme_acts = sim_data.get("enzyme_activities", {})

    if ntype == "gene":
        rxn_ids = _gene_rxn_map.get(nid, [])
        pws = {_rxn_pathway_map.get(r, "Unknown") for r in rxn_ids}
        info = [
            html.B(f"🧬 Gene: {nid}", style={"fontSize": "14px", "color": "#ff6688"}),
            html.P(f"Reactions: {len(rxn_ids)}", style={"margin": "4px 0", "color": "#aaa"}),
            html.P(f"Pathways: {', '.join(sorted(pws))}", style={"color": "#888", "fontSize": "11px"}),
        ]
        if nid in enzyme_acts:
            fa = enzyme_acts[nid][-1] if enzyme_acts[nid] else 1.0
            info.extend([
                html.Hr(style={"borderColor": "#333", "margin": "6px 0"}),
                html.Div("🔬 Enzyme Kinetics:", style={"color": "#a855f7", "fontSize": "11px", "fontWeight": "bold"}),
                html.Div(f"Activity: {fa:.1%}", style={
                    "color": "#ef4444" if fa < 0.3 else "#eab308" if fa < 0.7 else "#22c55e",
                    "fontSize": "12px", "fontWeight": "bold"}),
            ])
        if nid == sim_data.get("target_gene"):
            impact = sim_data.get("mutation_impact", "?")
            info.append(html.Div(f"⚠️ MUTATED ({impact})", style={
                "color": "#ef4444", "fontSize": "12px", "fontWeight": "bold", "marginTop": "4px"}))
        return html.Div(info)

    elif ntype == "reaction":
        try:
            rname = model.reactions.get_by_id(nid).name
        except Exception:
            rname = nid
        pw = _rxn_pathway_map.get(nid, "Unknown")
        fw = wt_fluxes.get(nid, 0)
        fd = drug_fluxes.get(nid, 0)
        pct = ((fd - fw) / fw * 100) if abs(fw) > 1e-9 else 0
        fc = "#ff4444" if pct < -10 else ("#ffaa00" if pct < -1 else "#44ff44")
        return html.Div([
            html.B(f"⚗️ {rname}", style={"fontSize": "13px", "color": "#ffaa44"}),
            html.P(f"Pathway: {pw}", style={"color": "#888", "margin": "3px 0", "fontSize": "11px"}),
            html.Div(style={"display": "flex", "gap": "15px", "marginTop": "4px", "fontSize": "12px"}, children=[
                html.Span(f"WT: {fw:.4f}", style={"color": "#44aaff"}),
                html.Span(f"Perturbed: {fd:.4f}", style={"color": fc}),
                html.Span(f"({pct:+.1f}%)", style={"color": fc, "fontWeight": "bold"}),
            ]),
        ])

    elif ntype == "metabolite":
        mname = get_met_name(nid)
        nc = len(_met_rxn_map.get(nid, []))
        return html.Div([
            html.B(f"🧪 {mname}", style={"fontSize": "13px", "color": "#ffdd55"}),
            html.P(f"ID: {nid}", style={"color": "#666", "margin": "3px 0", "fontSize": "11px"}),
            html.P(f"Connected reactions: {nc}", style={"color": "#888", "fontSize": "11px"}),
        ])

    return html.Span(f"Node: {nid}", style={"color": "#888"})


if __name__ == "__main__":
    app.run(debug=True, port=8051)