import dash
from dash import html, dcc, Input, Output, State, ctx, no_update, callback
import dash_cytoscape as cyto
import cobra
import pandas as pd
import json
import numpy as np
import plotly.graph_objects as go
import sys
import os
import hashlib

# Add src to path for process imports
project_root = os.path.join(os.path.dirname(__file__), "../..")
sys.path.insert(0, project_root)

# --- Load extra layouts for Cytoscape ---
cyto.load_extra_layouts()

VIVARIUM_AVAILABLE = False

try:
    from src.vivarium_model.whole_cell_composite import run_hybrid_simulation
    from src.vivarium_model.binding_kinetics import DRUG_TARGET_PARAMS
    from src.vivarium_model.drug_diffusion import DRUG_PERMEABILITY
    VIVARIUM_AVAILABLE = True
    print("✅ Vivarium modules loaded successfully.")
except ImportError as e:
    print("⚠️ Vivarium modules not found. Please ensure Vivarium is installed and src/vivarium_model is in the path.")
    print(f"Error details: {e}")

# --- DRUG DATABASE (extended with pharmacokinetic data) ---
DRUG_DB = {
    "Control (No Treatment)": {
        "targets": [],
        "desc": "Baseline metabolism.",
        "mic_uM": 0,
        "typical_dose_uM": 0,
    },
    "Trimethoprim": {
        "targets": ["MG228"],
        "desc": "Inhibits DHFR (Folate synthesis). Competitive inhibitor with Ki ≈ 5 nM.",
        "mic_uM": 2.0,
        "typical_dose_uM": 10.0,
    },
    "Methotrexate": {
        "targets": ["MG228", "MG006"],
        "desc": "Blocks dTMP production. Tight-binding DHFR inhibitor (Ki ≈ 1 nM) + thymidylate synthase.",
        "mic_uM": 0.5,
        "typical_dose_uM": 5.0,
    },
    "Fosmidomycin": {
        "targets": ["MG066"],
        "desc": "Inhibits DXR in isoprenoid synthesis (MEP pathway). Slow tight-binding inhibitor.",
        "mic_uM": 10.0,
        "typical_dose_uM": 50.0,
    },
    "Cerulenin": {
        "targets": ["MG212", "MG114"],
        "desc": "Irreversible (covalent) inhibitor of fatty acid synthase. Time-dependent killing.",
        "mic_uM": 5.0,
        "typical_dose_uM": 20.0,
    },
    "Mupirocin": {
        "targets": ["MG345"],
        "desc": "Inhibits isoleucyl-tRNA synthetase. Competitive with isoleucine (Ki ≈ 20 nM).",
        "mic_uM": 0.1,
        "typical_dose_uM": 1.0,
    },
    "Generic Glycolysis Inhibitor": {
        "targets": ["MG041", "MG429"],
        "desc": "Blocks glucose import/phosphorylation. Cooperative inhibition.",
        "mic_uM": 20.0,
        "typical_dose_uM": 100.0,
    },
}

# --- Cofactors / hub metabolites to EXCLUDE from the graph ---
HUB_METABOLITES = set()

# Color palette
CLR_BG = "#0a0a1a"

# --- LOAD MODEL & DATA ---
print("Loading metabolic model...")
model = cobra.io.read_sbml_model("src/iPS189.xml")
print(f"Model loaded: {len(model.genes)} genes, {len(model.reactions)} reactions, {len(model.metabolites)} metabolites")

try:
    rxn_df = pd.read_csv("data/master_map.csv")
except Exception:
    rxn_df = pd.DataFrame()

# --- Identify hub metabolites dynamically (top N most-connected) ---
_met_degree = {}
for rxn in model.reactions:
    for met in rxn.metabolites:
        _met_degree[met.id] = _met_degree.get(met.id, 0) + 1

HUB_THRESHOLD = 15
HUB_METABOLITES = {mid for mid, deg in _met_degree.items() if deg > HUB_THRESHOLD}
print(f"Hub metabolites excluded (degree > {HUB_THRESHOLD}): {len(HUB_METABOLITES)}")

# Pre-index the master_map for fast lookups
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

_pathway_rxns = {}
for rxn_id, pw in _rxn_pathway_map.items():
    if pw not in _pathway_rxns:
        _pathway_rxns[pw] = []
    _pathway_rxns[pw].append(rxn_id)


# =====================================================================
# SIMULATION BACKENDS
# =====================================================================

def run_fba_legacy(drug_name, efficacy):
    """Legacy pure-FBA approach (fallback when Vivarium unavailable)."""
    model.objective = "Biomass"
    sol_wt = model.optimize()
    wt_growth = sol_wt.objective_value
    wt_fluxes = sol_wt.fluxes

    drug_data = DRUG_DB[drug_name]
    with model:
        for tid in drug_data["targets"]:
            if tid in _gene_rxn_map:
                for rxn_id in _gene_rxn_map[tid]:
                    rxn = model.reactions.get_by_id(rxn_id)
                    lo, hi = rxn.bounds
                    rxn.bounds = (lo * (1 - efficacy), hi * (1 - efficacy))
        sol_drug = model.optimize()
        if sol_drug.status == "optimal":
            drug_growth = sol_drug.objective_value
            drug_fluxes = sol_drug.fluxes
        else:
            drug_growth = 0.0
            drug_fluxes = wt_fluxes * 0

    return {
        "wt_growth": float(wt_growth),
        "drug_growth": float(drug_growth),
        "wt_fluxes": {rid: float(wt_fluxes.get(rid, 0)) for rid in wt_fluxes.index},
        "drug_fluxes": {rid: float(drug_fluxes.get(rid, 0)) for rid in drug_fluxes.index},
        "time": [0.0],
        "drug_intracellular": [0.0],
        "enzyme_activities": {},
        "growth_rate": [float(drug_growth)],
        "simulation_mode": "legacy_fba",
    }


_sim_cache = {}


def run_hybrid_cached(drug_name, concentration_uM, sim_time=300.0):
    """Run Vivarium hybrid simulation with caching."""
    key = (drug_name, round(concentration_uM, 2), round(sim_time, 1))
    if key not in _sim_cache:
        if VIVARIUM_AVAILABLE and drug_name != "Control (No Treatment)":
            result = run_hybrid_simulation(
                drug_name=drug_name,
                extracellular_concentration=concentration_uM,
                total_time=sim_time,
                emit_step=max(sim_time / 30.0, 1.0),  # ~30 data points
            )
            result["simulation_mode"] = "vivarium_hybrid"
        else:
            # Fallback: convert concentration to an effective "efficacy"
            mic = DRUG_DB[drug_name].get("mic_uM", 10.0)
            if mic > 0 and concentration_uM > 0:
                efficacy = concentration_uM / (mic + concentration_uM)
            else:
                efficacy = 0.0
            result = run_fba_legacy(drug_name, efficacy)

        _sim_cache[key] = result
        # Keep cache bounded
        if len(_sim_cache) > 30:
            oldest = next(iter(_sim_cache))
            del _sim_cache[oldest]

    return _sim_cache[key]


# =====================================================================
# 3D NETWORK LAYOUT (force-directed in 3D)
# =====================================================================

def _stable_hash(s, mod=2**31):
    """Deterministic hash for reproducible layout."""
    return int(hashlib.md5(s.encode()).hexdigest(), 16) % mod


def compute_3d_layout():
    """
    Compute 3D positions for all nodes using a pathway-grouped force-directed
    approach. Pathways are placed in distinct spherical sectors, nodes within
    each pathway are spread locally.
    """
    np.random.seed(42)

    # Collect all node ids & their pathway membership
    node_ids = []
    node_types = {}
    node_pathway = {}

    for gene in model.genes:
        node_ids.append(gene.id)
        node_types[gene.id] = "gene"
        rxn_ids = _gene_rxn_map.get(gene.id, [])
        if rxn_ids:
            node_pathway[gene.id] = _rxn_pathway_map.get(rxn_ids[0], "Unknown")
        else:
            node_pathway[gene.id] = "Unknown"

    seen_mets = set()
    for rxn in model.reactions:
        if not rxn.genes:
            continue
        node_ids.append(rxn.id)
        node_types[rxn.id] = "reaction"
        node_pathway[rxn.id] = _rxn_pathway_map.get(rxn.id, "Unknown")

        for met in rxn.metabolites:
            if met.id in HUB_METABOLITES:
                continue
            if met.id not in seen_mets:
                node_ids.append(met.id)
                node_types[met.id] = "metabolite"
                first_rxn = _met_rxn_map.get(met.id, [None])[0]
                if first_rxn:
                    node_pathway[met.id] = _rxn_pathway_map.get(first_rxn, "Unknown")
                else:
                    node_pathway[met.id] = "Unknown"
                seen_mets.add(met.id)

    # Assign pathway centers on a sphere
    pathways = sorted(set(node_pathway.values()))
    pw_centers = {}
    n_pw = len(pathways)
    golden_angle = np.pi * (3.0 - np.sqrt(5.0))
    sphere_radius = 80.0

    for i, pw in enumerate(pathways):
        # Fibonacci sphere for even distribution
        y = 1.0 - (2.0 * i / max(n_pw - 1, 1))
        r = np.sqrt(max(0, 1.0 - y * y))
        theta = golden_angle * i
        px = sphere_radius * r * np.cos(theta)
        py = sphere_radius * y
        pz = sphere_radius * r * np.sin(theta)
        pw_centers[pw] = np.array([px, py, pz])

    # Place nodes around their pathway center with jitter based on type
    positions = {}
    pw_node_count = {}
    for nid in node_ids:
        pw = node_pathway.get(nid, "Unknown")
        center = pw_centers.get(pw, np.array([0, 0, 0]))
        nt = node_types.get(nid, "metabolite")

        # Use stable hash for deterministic jitter
        h = _stable_hash(nid)
        np.random.seed(h % (2**31))

        if nt == "gene":
            spread = 12.0
        elif nt == "reaction":
            spread = 18.0
        else:
            spread = 25.0

        jitter = np.random.randn(3) * spread
        positions[nid] = center + jitter

    # Simple spring-based relaxation to reduce overlap
    node_list = list(positions.keys())
    pos_array = np.array([positions[n] for n in node_list])

    # Build edge list for attraction
    edge_pairs = []
    node_index = {n: i for i, n in enumerate(node_list)}

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

    # Run a few iterations of force-directed relaxation
    n_nodes = len(node_list)
    for iteration in range(30):
        forces = np.zeros_like(pos_array)

        # Edge attraction (spring force)
        for i, j in edge_pairs:
            diff = pos_array[j] - pos_array[i]
            dist = np.linalg.norm(diff) + 1e-6
            ideal_len = 8.0
            f = (dist - ideal_len) * 0.01
            force = diff / dist * f
            forces[i] += force
            forces[j] -= force

        # Repulsion (only between nearby nodes for performance)
        # Use a simple approach: sample pairs
        if n_nodes < 500:
            for i in range(n_nodes):
                for j in range(i + 1, n_nodes):
                    diff = pos_array[j] - pos_array[i]
                    dist_sq = np.dot(diff, diff) + 1e-6
                    if dist_sq < 900:  # Only repel if close (< 30 units)
                        f = -50.0 / dist_sq
                        dist = np.sqrt(dist_sq)
                        force = diff / dist * f
                        forces[i] += force
                        forces[j] -= force
        else:
            # For large graphs, sample random pairs
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
                    f = -50.0 / dist_sq
                    dist = np.sqrt(dist_sq)
                    force = diff / dist * f
                    forces[i] += force
                    forces[j] -= force

        # Pathway gravity: pull nodes toward their pathway center
        for idx, nid in enumerate(node_list):
            pw = node_pathway.get(nid, "Unknown")
            center = pw_centers.get(pw, np.array([0, 0, 0]))
            diff = center - pos_array[idx]
            forces[idx] += diff * 0.005

        # Apply forces with damping
        damping = 0.8 * (1.0 - iteration / 30.0)
        pos_array += forces * damping

    # Write back
    for idx, nid in enumerate(node_list):
        positions[nid] = pos_array[idx]

    return positions, node_types, node_pathway, pw_centers


print("Computing 3D layout (one-time)...")
NODE_POSITIONS_3D, NODE_TYPES, NODE_PATHWAYS, PATHWAY_CENTERS = compute_3d_layout()
print(f"3D layout computed for {len(NODE_POSITIONS_3D)} nodes")


# =====================================================================
# BUILD 3D PLOTLY FIGURE
# =====================================================================

def build_3d_figure(highlighted_genes=None, highlighted_rxns=None, highlighted_mets=None,
                    secondary_rxns=None, secondary_mets=None, hit_pathways=None,
                    enzyme_activities=None):
    """Build a Plotly 3D scatter figure for the metabolic network."""
    if highlighted_genes is None:
        highlighted_genes = set()
    if highlighted_rxns is None:
        highlighted_rxns = set()
    if highlighted_mets is None:
        highlighted_mets = set()
    if secondary_rxns is None:
        secondary_rxns = set()
    if secondary_mets is None:
        secondary_mets = set()
    if hit_pathways is None:
        hit_pathways = set()
    if enzyme_activities is None:
        enzyme_activities = {}

    all_affected = highlighted_genes | highlighted_rxns | highlighted_mets | secondary_rxns | secondary_mets

    # === EDGES ===
    edge_x, edge_y, edge_z = [], [], []
    edge_hit_x, edge_hit_y, edge_hit_z = [], [], []
    edge_ripple_x, edge_ripple_y, edge_ripple_z = [], [], []

    for rxn in model.reactions:
        if rxn.id not in NODE_POSITIONS_3D:
            continue
        rx, ry, rz = NODE_POSITIONS_3D[rxn.id]

        for gene in rxn.genes:
            if gene.id not in NODE_POSITIONS_3D:
                continue
            gx, gy, gz = NODE_POSITIONS_3D[gene.id]

            is_hit = (gene.id in highlighted_genes and rxn.id in highlighted_rxns)
            is_ripple = (rxn.id in secondary_rxns)

            if is_hit:
                target = (edge_hit_x, edge_hit_y, edge_hit_z)
            elif is_ripple:
                target = (edge_ripple_x, edge_ripple_y, edge_ripple_z)
            else:
                target = (edge_x, edge_y, edge_z)

            target[0].extend([gx, rx, None])
            target[1].extend([gy, ry, None])
            target[2].extend([gz, rz, None])

        for met in rxn.metabolites:
            if met.id not in NODE_POSITIONS_3D:
                continue
            mx, my, mz = NODE_POSITIONS_3D[met.id]

            is_hit = (rxn.id in highlighted_rxns and met.id in highlighted_mets)
            is_ripple = (rxn.id in secondary_rxns or met.id in secondary_mets)

            if is_hit:
                target = (edge_hit_x, edge_hit_y, edge_hit_z)
            elif is_ripple:
                target = (edge_ripple_x, edge_ripple_y, edge_ripple_z)
            else:
                target = (edge_x, edge_y, edge_z)

            target[0].extend([rx, mx, None])
            target[1].extend([ry, my, None])
            target[2].extend([rz, mz, None])

    traces = []

    # Default edges
    traces.append(go.Scatter3d(
        x=edge_x, y=edge_y, z=edge_z,
        mode="lines", name="Connections",
        line=dict(color="rgba(30,41,59,0.25)", width=1),
        hoverinfo="skip", showlegend=False,
    ))

    # Hit edges
    if edge_hit_x:
        traces.append(go.Scatter3d(
            x=edge_hit_x, y=edge_hit_y, z=edge_hit_z,
            mode="lines", name="Direct impact",
            line=dict(color="rgba(251,146,60,0.8)", width=3),
            hoverinfo="skip", showlegend=False,
        ))

    # Ripple edges
    if edge_ripple_x:
        traces.append(go.Scatter3d(
            x=edge_ripple_x, y=edge_ripple_y, z=edge_ripple_z,
            mode="lines", name="Cascade",
            line=dict(color="rgba(202,138,4,0.4)", width=1.5),
            hoverinfo="skip", showlegend=False,
        ))

    # === NODES by category ===
    # Helper to build node traces
    def _add_node_trace(ids, color, size, symbol, name, opacity=1.0, border_color=None, border_width=0):
        if not ids:
            return
        xs = [NODE_POSITIONS_3D[n][0] for n in ids if n in NODE_POSITIONS_3D]
        ys = [NODE_POSITIONS_3D[n][1] for n in ids if n in NODE_POSITIONS_3D]
        zs = [NODE_POSITIONS_3D[n][2] for n in ids if n in NODE_POSITIONS_3D]
        labels = [n for n in ids if n in NODE_POSITIONS_3D]

        hover_texts = []
        for nid in labels:
            nt = NODE_TYPES.get(nid, "")
            pw = NODE_PATHWAYS.get(nid, "Unknown")
            if nt == "gene":
                n_rxns = len(_gene_rxn_map.get(nid, []))
                enz_info = ""
                if nid in enzyme_activities:
                    acts = enzyme_activities[nid]
                    final = acts[-1] if acts else 1.0
                    enz_info = f"<br>Enzyme activity: {final:.1%}"
                hover_texts.append(f"<b>🧬 {nid}</b><br>Type: Gene<br>Reactions: {n_rxns}<br>Pathway: {pw}{enz_info}")
            elif nt == "reaction":
                rxn_obj = model.reactions.get_by_id(nid) if nid in [r.id for r in model.reactions] else None
                rxn_name = rxn_obj.name if rxn_obj else nid
                hover_texts.append(f"<b>⚗️ {rxn_name}</b><br>ID: {nid}<br>Pathway: {pw}")
            elif nt == "metabolite":
                mname = get_met_name(nid)
                n_conn = len(_met_rxn_map.get(nid, []))
                hover_texts.append(f"<b>🧪 {mname}</b><br>ID: {nid}<br>Connected rxns: {n_conn}")
            else:
                hover_texts.append(nid)

        marker = dict(
            size=size, color=color, opacity=opacity, symbol=symbol,
        )
        if border_color and border_width > 0:
            marker["line"] = dict(color=border_color, width=border_width)

        traces.append(go.Scatter3d(
            x=xs, y=ys, z=zs,
            mode="markers+text" if size >= 10 else "markers",
            text=labels if size >= 10 else None,
            textposition="top center" if size >= 10 else None,
            textfont=dict(size=max(7, min(12, size // 2)), color=color) if size >= 10 else None,
            name=name,
            marker=marker,
            hovertext=hover_texts,
            hoverinfo="text",
            customdata=[[nid, NODE_TYPES.get(nid, "")] for nid in labels],
        ))

    # Categorize all nodes
    dimmed_genes = []
    dimmed_rxns = []
    dimmed_mets = []
    target_genes = []
    hit_rxn_list = []
    hit_met_list = []
    sec_rxn_list = []
    sec_met_list = []

    for nid, ntype in NODE_TYPES.items():
        if ntype == "gene":
            if nid in highlighted_genes:
                target_genes.append(nid)
            else:
                dimmed_genes.append(nid)
        elif ntype == "reaction":
            if nid in highlighted_rxns:
                hit_rxn_list.append(nid)
            elif nid in secondary_rxns:
                sec_rxn_list.append(nid)
            else:
                dimmed_rxns.append(nid)
        elif ntype == "metabolite":
            if nid in highlighted_mets:
                hit_met_list.append(nid)
            elif nid in secondary_mets:
                sec_met_list.append(nid)
            else:
                dimmed_mets.append(nid)

    # Dimmed nodes (background)
    _add_node_trace(dimmed_genes, "#0ea5e9", 3, "diamond", "Genes (dim)", opacity=0.15)
    _add_node_trace(dimmed_rxns, "#64748b", 2, "square", "Reactions (dim)", opacity=0.1)
    _add_node_trace(dimmed_mets, "#9ca3af", 1.5, "circle", "Metabolites (dim)", opacity=0.08)

    # Secondary cascade
    _add_node_trace(sec_rxn_list, "#c2410c", 5, "square", "Cascade reactions", opacity=0.8)
    _add_node_trace(sec_met_list, "#a16207", 3.5, "circle", "Cascade metabolites", opacity=0.6)

    # Direct hits
    _add_node_trace(hit_met_list, "#eab308", 7, "circle", "Direct metabolites",
                    opacity=1.0, border_color="#fef08a", border_width=1)
    _add_node_trace(hit_rxn_list, "#f97316", 10, "square", "Direct reactions",
                    opacity=1.0, border_color="#ffedd5", border_width=2)

    # Drug target genes (largest, most prominent)
    _add_node_trace(target_genes, "#d946ef", 16, "diamond", "🎯 Target genes",
                    opacity=1.0, border_color="#fdf4ff", border_width=3)

    # === PATHWAY LABELS as 3D annotations ===
    annotations = []
    for pw, center in PATHWAY_CENTERS.items():
        if pw == "Unknown":
            continue
        is_hit = pw in hit_pathways
        annotations.append(dict(
            x=center[0], y=center[1], z=center[2],
            text=f"<b>{pw}</b>" if is_hit else pw,
            showarrow=False,
            font=dict(
                size=13 if is_hit else 10,
                color="#c7d2fe" if is_hit else "#4a4a7a",
            ),
            bgcolor="rgba(49,46,129,0.6)" if is_hit else "rgba(20,20,50,0.4)",
            borderpad=4,
        ))

    # === LAYOUT ===
    axis_style = dict(
        showbackground=True,
        backgroundcolor="#0a0a1a",
        gridcolor="#111133",
        showticklabels=False,
        title="",
        showspikes=False,
        zeroline=False,
        showline=False,
    )

    fig = go.Figure(data=traces)

    fig.update_layout(
        scene=dict(
            xaxis=axis_style,
            yaxis=axis_style,
            zaxis=axis_style,
            bgcolor="#0a0a1a",
            annotations=annotations,
            camera=dict(
                eye=dict(x=1.6, y=1.6, z=1.0),
                up=dict(x=0, y=1, z=0),
            ),
            aspectmode="data",
            dragmode="orbit",
        ),
        paper_bgcolor="#0a0a1a",
        plot_bgcolor="#0a0a1a",
        margin=dict(l=0, r=0, t=0, b=0),
        showlegend=True,
        legend=dict(
            x=0.01, y=0.99, bgcolor="rgba(10,10,26,0.8)",
            font=dict(color="#888", size=10),
            bordercolor="#222255", borderwidth=1,
            itemsizing="constant",
        ),
        hovermode="closest",
        uirevision="keep",  # Preserve camera angle across updates
    )

    return fig


print("Building initial 3D figure...")
INITIAL_3D_FIG = build_3d_figure()
print("3D figure ready.")


def _legend_row(color, symbol, label):
    return html.Div(style={"display": "flex", "alignItems": "center", "gap": "6px", "marginBottom": "3px"}, children=[
        html.Span(symbol, style={"color": color, "fontSize": "12px", "width": "14px", "textAlign": "center"}),
        html.Span(label, style={"color": "#666"}),
    ])


# --- DASH APP ---
app = dash.Dash(__name__, title="iPS189 Hybrid Whole-Cell Simulator", update_title=None)

app.layout = html.Div(
    style={"backgroundColor": CLR_BG, "minHeight": "100vh", "fontFamily": "'Segoe UI', sans-serif",
           "overflow": "hidden"},
    children=[
        dcc.Store(id="sim-store", data={}),

        # --- HEADER ---
        html.Div(style={
            "background": "linear-gradient(135deg, #0a0a1a 0%, #111144 100%)",
            "padding": "15px 30px", "borderBottom": "1px solid #222255",
            "display": "flex", "alignItems": "center", "justifyContent": "space-between",
        }, children=[
            html.Div([
                html.H1("🦠 iPS189 Hybrid Whole-Cell Simulation", style={
                    "color": "#fff", "margin": 0, "fontSize": "24px",
                }),
                html.P([
                    "Drug diffusion → Enzyme kinetics → Metabolic FBA",
                    html.Span(
                        f"  ({'Vivarium' if VIVARIUM_AVAILABLE else 'FBA fallback'})",
                        style={"color": "#22c55e" if VIVARIUM_AVAILABLE else "#ef4444",
                               "fontWeight": "bold"},
                    ),
                    html.Span("  |  🌐 3D Interactive Network", style={"color": "#818cf8"}),
                ], style={"color": "#556", "margin": "2px 0 0 0", "fontSize": "12px"}),
            ]),
            html.Div(id="header-stats", style={"display": "flex", "gap": "30px", "alignItems": "center"}),
        ]),

        html.Div(style={"display": "flex", "height": "calc(100vh - 70px)"}, children=[

            # --- LEFT PANEL ---
            html.Div(style={
                "width": "300px", "flexShrink": 0, "padding": "16px",
                "backgroundColor": "#0d0d24", "overflowY": "auto",
                "borderRight": "1px solid #1a1a3e",
            }, children=[
                html.H3("💊 Drug Controls", style={"color": "#ddd", "marginTop": 0, "fontSize": "16px"}),

                html.Label("Drug Treatment:", style={"color": "#888", "fontSize": "11px"}),
                dcc.Dropdown(
                    id="drug-select",
                    options=[{"label": k, "value": k} for k in DRUG_DB],
                    value="Control (No Treatment)",
                    style={"marginBottom": "16px"},
                ),

                html.Label("Extracellular Drug Concentration (µM):", style={"color": "#888", "fontSize": "11px"}),
                html.Div(style={"display": "flex", "alignItems": "baseline", "gap": "8px"}, children=[
                    html.Div(id="conc-val", style={
                        "color": "#ff6600", "fontSize": "22px", "fontWeight": "bold",
                    }, children="10.0 µM"),
                ]),
                dcc.Slider(
                    id="conc-slider", min=0, max=200, step=1, value=10,
                    marks={
                        0: {"label": "0", "style": {"color": "#444"}},
                        50: {"label": "50 µM", "style": {"color": "#444"}},
                        100: {"label": "100 µM", "style": {"color": "#444"}},
                        200: {"label": "200 µM", "style": {"color": "#444"}},
                    },
                ),

                html.Label("Simulation Time (seconds):", style={"color": "#888", "fontSize": "11px", "marginTop": "8px"}),
                dcc.Slider(
                    id="time-slider", min=30, max=600, step=30, value=300,
                    marks={
                        30: {"label": "30s", "style": {"color": "#444"}},
                        300: {"label": "5 min", "style": {"color": "#444"}},
                        600: {"label": "10 min", "style": {"color": "#444"}},
                    },
                ),

                html.Label("Cascade Depth:", style={"color": "#888", "fontSize": "11px", "marginTop": "8px"}),
                dcc.Slider(
                    id="cascade-slider", min=0, max=3, step=1, value=1,
                    marks={i: {"label": str(i), "style": {"color": "#444"}} for i in range(4)},
                ),

                html.Hr(style={"borderColor": "#1a1a3e", "margin": "16px 0"}),

                html.Div(id="drug-info", style={
                    "padding": "12px", "backgroundColor": "#111133",
                    "borderRadius": "8px", "border": "1px solid #222255", "fontSize": "12px",
                }),

                html.Hr(style={"borderColor": "#1a1a3e", "margin": "16px 0"}),

                html.Div(id="kinetics-panel", style={
                    "padding": "12px", "backgroundColor": "#111133",
                    "borderRadius": "8px", "border": "1px solid #222255",
                }),

                html.Hr(style={"borderColor": "#1a1a3e", "margin": "16px 0"}),

                html.Div(id="metrics-panel"),

                html.Hr(style={"borderColor": "#1a1a3e", "margin": "16px 0"}),

                html.Div(style={"fontSize": "10px", "color": "#555"}, children=[
                    html.Div("LEGEND", style={
                        "color": "#444", "fontWeight": "bold",
                        "letterSpacing": "1px", "marginBottom": "6px",
                    }),
                    _legend_row("#d946ef", "◆", "Target Gene"),
                    _legend_row("#f97316", "■", "Direct Reaction"),
                    _legend_row("#eab308", "●", "Direct Metabolite"),
                    _legend_row("#c2410c", "■", "Cascade Reaction"),
                    _legend_row("#a16207", "●", "Cascade Metabolite"),
                    _legend_row("#334155", "●", "Unaffected"),
                    html.Div("─" * 20, style={"color": "#333", "margin": "4px 0"}),
                    html.Div("🔬 Hybrid model layers:", style={"color": "#556", "marginBottom": "3px"}),
                    _legend_row("#22c55e", "↘", "Drug diffusion (ODE)"),
                    _legend_row("#a855f7", "⊘", "Enzyme inhibition (kinetics)"),
                    _legend_row("#3b82f6", "⚡", "Metabolic flux (FBA)"),
                    html.Div("─" * 20, style={"color": "#333", "margin": "4px 0"}),
                    html.Div("🌐 3D Controls:", style={"color": "#556", "marginBottom": "3px"}),
                    _legend_row("#818cf8", "🖱", "Drag to rotate"),
                    _legend_row("#818cf8", "⚙", "Scroll to zoom"),
                    _legend_row("#818cf8", "⇧", "Shift+drag to pan"),
                ]),
            ]),

            # --- CENTER: 3D Network ---
            html.Div(style={"flex": 1, "position": "relative"}, children=[
                dcc.Graph(
                    id="network-3d",
                    figure=INITIAL_3D_FIG,
                    style={"width": "100%", "height": "100%"},
                    config={
                        "displayModeBar": True,
                        "modeBarButtonsToRemove": ["toImage", "resetCameraLastSave3d"],
                        "displaylogo": False,
                        "scrollZoom": True,
                    },
                ),

                html.Div(id="node-tooltip", style={
                    "position": "absolute", "bottom": "12px", "left": "12px",
                    "color": "#ddd", "padding": "10px 14px",
                    "backgroundColor": "rgba(10, 10, 30, 0.95)",
                    "borderRadius": "8px", "border": "1px solid #333366",
                    "fontSize": "12px", "maxWidth": "360px",
                    "backdropFilter": "blur(10px)", "pointerEvents": "none",
                    "transition": "opacity 0.3s ease",
                }),

                html.Div("🖱 Drag to rotate · Scroll to zoom · Shift+drag to pan", style={
                    "position": "absolute", "top": "10px", "right": "15px",
                    "color": "#444", "fontSize": "10px",
                    "backgroundColor": "rgba(10,10,26,0.7)", "padding": "4px 10px",
                    "borderRadius": "4px",
                }),
            ]),

            # --- RIGHT: Time-course graphs ---
            html.Div(id="timecourse-panel", style={
                "width": "300px", "flexShrink": 0, "padding": "12px",
                "backgroundColor": "#0d0d24", "overflowY": "auto",
                "borderLeft": "1px solid #1a1a3e",
            }),
        ]),
    ],
)


@app.callback(
    Output("conc-val", "children"),
    Input("conc-slider", "value"),
)
def update_conc_label(val):
    return f"{val:.1f} µM"


@app.callback(
    Output("sim-store", "data"),
    Input("drug-select", "value"),
    Input("conc-slider", "value"),
    Input("time-slider", "value"),
)
def compute_simulation(drug_name, conc_uM, sim_time):
    """Run hybrid simulation (diffusion → kinetics → FBA)."""
    result = run_hybrid_cached(drug_name, float(conc_uM), float(sim_time))

    # Serialize for dcc.Store
    drug_data = DRUG_DB[drug_name]
    relevant_rxn_ids = set()
    for tid in drug_data["targets"]:
        for rxn_id in _gene_rxn_map.get(tid, []):
            relevant_rxn_ids.add(rxn_id)
            for mid in _rxn_met_map.get(rxn_id, []):
                for r2 in _met_rxn_map.get(mid, []):
                    relevant_rxn_ids.add(r2)

    wt_fluxes = result.get("wt_fluxes", {})
    drug_fluxes = result.get("fluxes", result.get("drug_fluxes", {}))
    wt_subset = {rid: float(wt_fluxes.get(rid, 0)) for rid in relevant_rxn_ids}
    drug_subset = {rid: float(drug_fluxes.get(rid, 0)) for rid in relevant_rxn_ids}

    growth_series = result.get("growth_rate", [0.0])
    final_growth = growth_series[-1] if growth_series else 0.0

    return {
        "wt_growth": float(result.get("wt_growth", 0)),
        "drug_growth": float(final_growth),
        "wt_fluxes": wt_subset,
        "drug_fluxes": drug_subset,
        "time": result.get("time", [0.0]),
        "drug_intracellular": result.get("drug_intracellular", [0.0]),
        "enzyme_activities": result.get("enzyme_activities", {}),
        "growth_rate": growth_series,
        "simulation_mode": result.get("simulation_mode", "unknown"),
    }


@app.callback(
    Output("network-3d", "figure"),
    Output("drug-info", "children"),
    Output("metrics-panel", "children"),
    Output("header-stats", "children"),
    Output("kinetics-panel", "children"),
    Output("timecourse-panel", "children"),
    Input("drug-select", "value"),
    Input("cascade-slider", "value"),
    Input("sim-store", "data"),
)
def update_graph(drug_name, cascade_depth, sim_data):
    drug_data = DRUG_DB[drug_name]

    wt_growth = sim_data.get("wt_growth", 0)
    drug_growth = sim_data.get("drug_growth", 0)
    sim_mode = sim_data.get("simulation_mode", "unknown")

    # --- Compute highlight sets ---
    highlighted_genes = set(drug_data["targets"])
    highlighted_rxns = set()
    highlighted_mets = set()
    secondary_rxns = set()
    secondary_mets = set()
    hit_pathways = set()

    frontier_mets = set()
    for tid in drug_data["targets"]:
        for rxn_id in _gene_rxn_map.get(tid, []):
            highlighted_rxns.add(rxn_id)
            pw = _rxn_pathway_map.get(rxn_id, "Unknown")
            hit_pathways.add(pw)
            for mid in _rxn_met_map.get(rxn_id, []):
                highlighted_mets.add(mid)
                frontier_mets.add(mid)

    visited_rxns = set(highlighted_rxns)
    for _ in range(cascade_depth):
        next_frontier = set()
        for mid in frontier_mets:
            for rxn_id in _met_rxn_map.get(mid, []):
                if rxn_id in visited_rxns:
                    continue
                secondary_rxns.add(rxn_id)
                visited_rxns.add(rxn_id)
                pw = _rxn_pathway_map.get(rxn_id, "Unknown")
                hit_pathways.add(pw)
                for m2 in _rxn_met_map.get(rxn_id, []):
                    if m2 not in highlighted_mets:
                        secondary_mets.add(m2)
                    next_frontier.add(m2)
        frontier_mets = next_frontier

    enzyme_activities = sim_data.get("enzyme_activities", {})

    # --- Build 3D figure ---
    fig = build_3d_figure(
        highlighted_genes=highlighted_genes,
        highlighted_rxns=highlighted_rxns,
        highlighted_mets=highlighted_mets,
        secondary_rxns=secondary_rxns,
        secondary_mets=secondary_mets,
        hit_pathways=hit_pathways,
        enzyme_activities=enzyme_activities,
    )

    # --- Drug info ---
    mode_badge = html.Span(
        f" [{sim_mode.replace('_', ' ').upper()}]",
        style={
            "color": "#22c55e" if "vivarium" in sim_mode else "#ef4444",
            "fontSize": "9px", "fontWeight": "bold",
        }
    )

    info_children = [
        html.Div([f"💊 {drug_name}", mode_badge], style={
            "fontWeight": "bold", "color": "#ff6600", "fontSize": "14px", "marginBottom": "6px",
        }),
        html.Div(drug_data['desc'], style={"color": "#999", "fontSize": "11px"}),
    ]
    if drug_data["targets"]:
        info_children.append(
            html.Div(style={"marginTop": "8px", "display": "flex", "flexWrap": "wrap", "gap": "4px"}, children=[
                html.Span("Targets: ", style={"color": "#777", "fontSize": "11px"}),
                *[html.Span(t, style={
                    "backgroundColor": "#ff0040", "color": "#fff",
                    "padding": "1px 8px", "borderRadius": "10px",
                    "fontSize": "10px", "fontWeight": "bold",
                }) for t in drug_data["targets"]]
            ])
        )
        mic = drug_data.get("mic_uM", 0)
        if mic > 0:
            info_children.append(
                html.Div(f"MIC: {mic} µM", style={
                    "color": "#888", "fontSize": "11px", "marginTop": "4px",
                })
            )
    if hit_pathways - {"Unknown"}:
        info_children.append(
            html.Div(style={"marginTop": "8px"}, children=[
                html.Span("Pathways: ", style={"color": "#777", "fontSize": "11px"}),
                html.Span(", ".join(sorted(hit_pathways - {"Unknown"})),
                          style={"color": "#ff8844", "fontSize": "11px"}),
            ])
        )

    # --- Kinetics panel ---
    kinetics_children = []
    if enzyme_activities:
        kinetics_children.append(
            html.Div("🔬 ENZYME KINETICS", style={
                "color": "#a855f7", "fontWeight": "bold", "fontSize": "11px",
                "letterSpacing": "1px", "marginBottom": "8px",
            })
        )
        for gene_id, activities in enzyme_activities.items():
            final_activity = activities[-1] if activities else 1.0
            bar_width = max(0, min(100, final_activity * 100))
            bar_color = "#22c55e" if final_activity > 0.7 else ("#eab308" if final_activity > 0.3 else "#ef4444")
            kinetics_children.append(html.Div(style={"marginBottom": "8px"}, children=[
                html.Div(style={"display": "flex", "justifyContent": "space-between"}, children=[
                    html.Span(f"{gene_id}", style={"color": "#ddd", "fontSize": "11px", "fontWeight": "bold"}),
                    html.Span(f"{final_activity:.1%}", style={"color": bar_color, "fontSize": "11px"}),
                ]),
                html.Div(style={
                    "width": "100%", "height": "6px", "backgroundColor": "#1a1a3e",
                    "borderRadius": "3px", "overflow": "hidden", "marginTop": "2px",
                }, children=[
                    html.Div(style={
                        "width": f"{bar_width}%", "height": "100%", "backgroundColor": bar_color,
                        "borderRadius": "3px", "transition": "width 0.5s ease",
                    })
                ]),
            ]))
    else:
        kinetics_children.append(
            html.Div("No enzyme targets", style={"color": "#444", "fontSize": "11px"})
        )

    # --- Metrics ---
    fitness = (drug_growth / wt_growth * 100) if wt_growth > 0 else 0
    delta_pct = ((drug_growth - wt_growth) / wt_growth * 100) if wt_growth > 0 else 0

    metrics = html.Div(style={"display": "flex", "flexDirection": "column", "gap": "6px"}, children=[
        _metric_card("Wild-Type Growth", f"{wt_growth:.4f}", None, "#44aaff"),
        _metric_card("Drug Growth", f"{drug_growth:.4f}",
                     f"{delta_pct:+.1f}%",
                     "#ff4444" if delta_pct < -1 else ("#44ff44" if delta_pct > 1 else "#888")),
        _metric_card("Relative Fitness", f"{fitness:.1f}%", None,
                     "#ff4444" if fitness < 50 else ("#ffaa00" if fitness < 80 else "#44ff44")),
        _metric_card("Impact",
                     f"{len(highlighted_rxns)} + {len(secondary_rxns)} rxns",
                     f"{len(highlighted_mets) + len(secondary_mets)} metabolites",
                     "#ff6600"),
    ])

    # --- Header stats ---
    header = []
    if drug_name != "Control (No Treatment)":
        fitness_color = "#ff4444" if fitness < 50 else ("#ffaa00" if fitness < 80 else "#44ff44")
        header = [
            html.Div([
                html.Div("FITNESS", style={"color": "#556", "fontSize": "9px", "letterSpacing": "1px"}),
                html.Div(f"{fitness:.1f}%", style={"color": fitness_color, "fontSize": "20px", "fontWeight": "bold"}),
            ]),
            html.Div([
                html.Div("IMPACT", style={"color": "#556", "fontSize": "9px", "letterSpacing": "1px"}),
                html.Div(f"{len(highlighted_rxns) + len(secondary_rxns)} reactions",
                         style={"color": "#ff6600", "fontSize": "14px"}),
            ]),
            html.Div([
                html.Div("MODEL", style={"color": "#556", "fontSize": "9px", "letterSpacing": "1px"}),
                html.Div("Hybrid" if "vivarium" in sim_mode else "FBA",
                         style={"color": "#22c55e" if "vivarium" in sim_mode else "#ef4444", "fontSize": "14px"}),
            ]),
        ]

    timecourse = _build_timecourse_panel(sim_data, drug_name)

    return fig, info_children, metrics, header, kinetics_children, timecourse


def _build_timecourse_panel(sim_data, drug_name):
    """Build the right-side time-course visualization panel."""
    times = sim_data.get("time", [0.0])
    drug_intra = sim_data.get("drug_intracellular", [0.0])
    enzyme_acts = sim_data.get("enzyme_activities", {})
    growth_rates = sim_data.get("growth_rate", [0.0])

    if drug_name == "Control (No Treatment)" or len(times) < 2:
        return html.Div([
            html.Div("📈 TIME COURSE", style={
                "color": "#556", "fontWeight": "bold", "fontSize": "11px",
                "letterSpacing": "1px", "marginBottom": "12px",
            }),
            html.Div("Select a drug to see dynamics", style={"color": "#444", "fontSize": "12px"}),
        ])

    children = [
        html.Div("📈 TIME COURSE", style={
            "color": "#556", "fontWeight": "bold", "fontSize": "11px",
            "letterSpacing": "1px", "marginBottom": "12px",
        }),
    ]

    # 1) Drug diffusion plot
    peak_drug = max(drug_intra) if drug_intra else 0
    children.append(dcc.Graph(
        figure={
            "data": [{
                "x": times, "y": drug_intra,
                "type": "scatter", "mode": "lines",
                "line": {"color": "#22c55e", "width": 2},
                "name": "Intracellular [Drug]",
                "fill": "tozeroy",
                "fillcolor": "rgba(34, 197, 94, 0.1)",
            }],
            "layout": {
                "title": {"text": f"Drug Diffusion (peak: {peak_drug:.2f} µM)", "font": {"size": 11, "color": "#888"}},
                "xaxis": {"title": "Time (s)", "color": "#555", "gridcolor": "#1a1a3e"},
                "yaxis": {"title": "[Drug] µM", "color": "#555", "gridcolor": "#1a1a3e"},
                "paper_bgcolor": "rgba(0,0,0,0)",
                "plot_bgcolor": "rgba(13,13,36,0.8)",
                "margin": {"l": 50, "r": 10, "t": 30, "b": 35},
                "height": 180,
                "font": {"color": "#888"},
            },
        },
        config={"displayModeBar": False},
        style={"marginBottom": "8px"},
    ))

    # 2) Enzyme activity plot
    if enzyme_acts:
        traces = []
        colors = ["#d946ef", "#f97316", "#eab308", "#22d3ee", "#a855f7"]
        for i, (gene_id, activities) in enumerate(enzyme_acts.items()):
            final_a = activities[-1] if activities else 1.0
            traces.append({
                "x": times[:len(activities)],
                "y": activities,
                "type": "scatter", "mode": "lines",
                "line": {"color": colors[i % len(colors)], "width": 2},
                "name": f"{gene_id} ({final_a:.0%})",
            })

        children.append(dcc.Graph(
            figure={
                "data": traces,
                "layout": {
                    "title": {"text": "Enzyme Activity", "font": {"size": 11, "color": "#888"}},
                    "xaxis": {"title": "Time (s)", "color": "#555", "gridcolor": "#1a1a3e"},
                    "yaxis": {"title": "Activity", "range": [-0.05, 1.05], "color": "#555", "gridcolor": "#1a1a3e"},
                    "paper_bgcolor": "rgba(0,0,0,0)",
                    "plot_bgcolor": "rgba(13,13,36,0.8)",
                    "margin": {"l": 50, "r": 10, "t": 30, "b": 35},
                    "height": 180,
                    "font": {"color": "#888"},
                    "showlegend": True,
                    "legend": {"font": {"size": 9, "color": "#888"}, "x": 0.55, "y": 0.95},
                    "shapes": [{
                        "type": "line", "x0": times[0], "x1": times[-1],
                        "y0": 0.5, "y1": 0.5,
                        "line": {"color": "#ef4444", "width": 1, "dash": "dot"},
                    }],
                },
            },
            config={"displayModeBar": False},
            style={"marginBottom": "8px"},
        ))

    # 3) Growth rate plot
    wt_growth = sim_data.get("wt_growth", 0)
    children.append(dcc.Graph(
        figure={
            "data": [
                {
                    "x": times[:len(growth_rates)],
                    "y": growth_rates,
                    "type": "scatter", "mode": "lines",
                    "line": {"color": "#3b82f6", "width": 2},
                    "name": "Growth Rate",
                    "fill": "tozeroy",
                    "fillcolor": "rgba(59, 130, 246, 0.1)",
                },
                {
                    "x": [times[0], times[-1]],
                    "y": [wt_growth, wt_growth],
                    "type": "scatter", "mode": "lines",
                    "line": {"color": "#44aaff", "width": 1, "dash": "dash"},
                    "name": f"WT ({wt_growth:.2f})",
                },
            ],
            "layout": {
                "title": {"text": "Growth Rate (FBA)", "font": {"size": 11, "color": "#888"}},
                "xaxis": {"title": "Time (s)", "color": "#555", "gridcolor": "#1a1a3e"},
                "yaxis": {"title": "Growth (1/h)", "color": "#555", "gridcolor": "#1a1a3e"},
                "paper_bgcolor": "rgba(0,0,0,0)",
                "plot_bgcolor": "rgba(13,13,36,0.8)",
                "margin": {"l": 50, "r": 10, "t": 30, "b": 35},
                "height": 180,
                "font": {"color": "#888"},
                "showlegend": True,
                "legend": {"font": {"size": 9, "color": "#888"}, "x": 0.55, "y": 0.95},
            },
        },
        config={"displayModeBar": False},
        style={"marginBottom": "8px"},
    ))

    # 4) Simulation summary
    final_drug = drug_intra[-1] if drug_intra else 0
    final_growth = growth_rates[-1] if growth_rates else 0
    children.append(html.Div(style={
        "padding": "8px", "backgroundColor": "#111133",
        "borderRadius": "6px", "border": "1px solid #222255",
        "fontSize": "10px", "color": "#556", "marginTop": "8px",
    }, children=[
        html.Div("SIMULATION DETAILS", style={"fontWeight": "bold", "letterSpacing": "1px", "marginBottom": "4px"}),
        html.Div(f"⏱ Duration: {times[-1]:.0f}s ({len(times)} steps)"),
        html.Div(f"📊 Peak [Drug]ᵢₙ: {max(drug_intra):.2f} µM"),
        html.Div(f"📊 Final [Drug]ᵢₙ: {final_drug:.2f} µM"),
        html.Div(f"📉 Final Growth: {final_growth:.4f} h⁻¹"),
        html.Div(f"📉 WT Growth: {wt_growth:.4f} h⁻¹"),
        html.Div(f"🔬 Mode: {sim_data.get('simulation_mode', 'unknown')}"),
    ]))

    return html.Div(children)


def _metric_card(title, value, subtitle=None, accent_color="#aaa"):
    children = [
        html.Div(title, style={
            "color": "#556", "fontSize": "9px", "textTransform": "uppercase",
            "letterSpacing": "0.5px",
        }),
        html.Div(value, style={
            "color": accent_color, "fontSize": "18px", "fontWeight": "bold",
            "marginTop": "1px",
        }),
    ]
    if subtitle:
        children.append(
            html.Div(subtitle, style={"color": "#445", "fontSize": "10px", "marginTop": "1px"})
        )
    return html.Div(style={
        "backgroundColor": "#0d0d28", "borderRadius": "6px",
        "padding": "10px 12px", "border": f"1px solid {accent_color}15",
    }, children=children)


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
            return html.Div(
                dcc.Markdown(hover.replace("<b>", "**").replace("</b>", "**").replace("<br>", "\n\n")),
                style={"color": "#ddd"},
            )
        return html.Span("Click a node for details", style={"color": "#444"})

    nid = custom[0]
    ntype = custom[1]

    wt_fluxes = sim_data.get("wt_fluxes", {})
    drug_fluxes = sim_data.get("drug_fluxes", {})
    enzyme_acts = sim_data.get("enzyme_activities", {})

    if ntype == "gene":
        rxn_ids = _gene_rxn_map.get(nid, [])
        pathways = set()
        for rid in rxn_ids:
            pathways.add(_rxn_pathway_map.get(rid, "Unknown"))

        enzyme_info = []
        if nid in enzyme_acts:
            acts = enzyme_acts[nid]
            final_act = acts[-1] if acts else 1.0
            enzyme_info = [
                html.Hr(style={"borderColor": "#333", "margin": "6px 0"}),
                html.Div("🔬 Enzyme Kinetics:", style={"color": "#a855f7", "fontSize": "11px", "fontWeight": "bold"}),
                html.Div(f"Final activity: {final_act:.1%}", style={
                    "color": "#ef4444" if final_act < 0.3 else ("#eab308" if final_act < 0.7 else "#22c55e"),
                    "fontSize": "12px", "fontWeight": "bold",
                }),
                html.Div(f"Inhibited by {(1 - final_act) * 100:.1f}%", style={"color": "#888", "fontSize": "11px"}),
            ]

        return html.Div([
            html.B(f"🧬 Gene: {nid}", style={"fontSize": "14px", "color": "#ff6688"}),
            html.P(f"Reactions: {len(rxn_ids)}", style={"margin": "4px 0", "color": "#aaa"}),
            html.P(f"Pathways: {', '.join(sorted(pathways))}", style={"color": "#888", "fontSize": "11px"}),
            *enzyme_info,
        ])

    elif ntype == "reaction":
        try:
            rxn_obj = model.reactions.get_by_id(nid)
            full_name = rxn_obj.name
        except KeyError:
            full_name = nid
        pathway = _rxn_pathway_map.get(nid, "Unknown")
        flux_wt = wt_fluxes.get(nid, 0)
        flux_drug = drug_fluxes.get(nid, 0)
        pct = ((flux_drug - flux_wt) / flux_wt * 100) if abs(flux_wt) > 1e-9 else 0
        fc = "#ff4444" if pct < -10 else ("#ffaa00" if pct < -1 else "#44ff44")

        gene_ids = _rxn_gene_map.get(nid, [])
        constraint_info = []
        for gid in gene_ids:
            if gid in enzyme_acts:
                acts = enzyme_acts[gid]
                final_act = acts[-1] if acts else 1.0
                constraint_info.append(
                    html.Div(f"  ⊘ {gid}: {final_act:.1%} active", style={
                        "color": "#a855f7", "fontSize": "11px",
                    })
                )

        return html.Div([
            html.B(f"⚗️ {full_name}", style={"fontSize": "13px", "color": "#ffaa44"}),
            html.P(f"Pathway: {pathway}", style={"color": "#888", "margin": "3px 0", "fontSize": "11px"}),
            html.Div(style={"display": "flex", "gap": "15px", "marginTop": "4px", "fontSize": "12px"}, children=[
                html.Span(f"WT: {flux_wt:.4f}", style={"color": "#44aaff"}),
                html.Span(f"Drug: {flux_drug:.4f}", style={"color": fc}),
                html.Span(f"({pct:+.1f}%)", style={"color": fc, "fontWeight": "bold"}),
            ]),
            *constraint_info,
        ])

    elif ntype == "metabolite":
        met_name = get_met_name(nid)
        connected = len(_met_rxn_map.get(nid, []))
        return html.Div([
            html.B(f"🧪 {met_name}", style={"fontSize": "13px", "color": "#ffdd55"}),
            html.P(f"ID: {nid}", style={"color": "#666", "margin": "3px 0", "fontSize": "11px"}),
            html.P(f"Connected reactions: {connected}", style={"color": "#888", "fontSize": "11px"}),
        ])

    return html.Span(f"Node: {nid}", style={"color": "#888"})


if __name__ == "__main__":
    app.run(debug=True, port=8050)