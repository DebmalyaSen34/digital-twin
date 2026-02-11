import dash
from dash import html, dcc, Input, Output, State, ctx, no_update, callback
import dash_cytoscape as cyto
import cobra
import pandas as pd
import json
import numpy as np
import sys
import os

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


def build_cytoscape_elements():
    """Build the organism graph with pathway grouping as parent nodes."""
    elements = []
    seen_mets = set()
    node_count = 0
    edge_count = 0

    pathway_ids = {}
    for pw in _pathway_rxns:
        pw_id = f"pw_{pw.replace(' ', '_').replace('&', 'and')}"
        pathway_ids[pw] = pw_id
        elements.append({
            "data": {"id": pw_id, "label": pw, "node_type": "pathway"},
            "classes": "pathway",
        })

    for gene in model.genes:
        rxn_ids = _gene_rxn_map.get(gene.id, [])
        parent = None
        if rxn_ids:
            pw = _rxn_pathway_map.get(rxn_ids[0], "Unknown")
            parent = pathway_ids.get(pw)
        node_data = {"id": gene.id, "label": gene.id, "node_type": "gene"}
        if parent:
            node_data["parent"] = parent
        elements.append({"data": node_data, "classes": "gene"})
        node_count += 1

    for rxn in model.reactions:
        if not rxn.genes:
            continue
        pathway = _rxn_pathway_map.get(rxn.id, "Unknown")
        parent = pathway_ids.get(pathway)
        node_data = {
            "id": rxn.id, "label": rxn.name[:30], "node_type": "reaction",
            "pathway": pathway, "full_name": rxn.name,
        }
        if parent:
            node_data["parent"] = parent
        elements.append({"data": node_data, "classes": "reaction"})
        node_count += 1

        for gene in rxn.genes:
            edge_id = f"eg_{gene.id}_{rxn.id}"
            elements.append({
                "data": {"id": edge_id, "source": gene.id, "target": rxn.id, "edge_type": "gene_rxn"},
                "classes": "edge-default",
            })
            edge_count += 1

        for met in rxn.metabolites:
            mid = met.id
            if mid in HUB_METABOLITES:
                continue
            if mid not in seen_mets:
                met_name = get_met_name(mid)
                first_rxn = _met_rxn_map.get(mid, [None])[0]
                met_parent = None
                if first_rxn:
                    mp = _rxn_pathway_map.get(first_rxn, "Unknown")
                    met_parent = pathway_ids.get(mp)
                met_data = {"id": mid, "label": met_name[:25], "node_type": "metabolite"}
                if met_parent:
                    met_data["parent"] = met_parent
                elements.append({"data": met_data, "classes": "metabolite"})
                seen_mets.add(mid)
                node_count += 1

            edge_id = f"em_{rxn.id}_{mid}"
            elements.append({
                "data": {"id": edge_id, "source": rxn.id, "target": mid, "edge_type": "rxn_met"},
                "classes": "edge-default",
            })
            edge_count += 1

    print(f"Graph built: {node_count} nodes, {edge_count} edges, {len(pathway_ids)} pathway groups")
    return elements


print("Building base graph (one-time)...")
BASE_ELEMENTS = build_cytoscape_elements()


# --- STYLESHEET ---
STYLESHEET = [
    # === PATHWAY COMPOUND NODES ===
    {"selector": "node.pathway", "style": {
        "background-color": "#1e1b4b", "background-opacity": 0.3,
        "border-width": 1, "border-color": "#4f46e5", "border-style": "dashed",
        "shape": "round-rectangle", "label": "data(label)",
        "font-size": 14, "color": "#818cf8", "text-valign": "top", "text-halign": "center",
        "text-margin-y": -10, "text-opacity": 0.8, "padding": 40,
        "compound-sizing-wrt-labels": "include",
    }},
    # === DEFAULT DIMMED ===
    {"selector": "node.gene", "style": {
        "background-color": "#0ea5e9", "width": 12, "height": 12, "shape": "diamond",
        "label": "", "opacity": 0.3,
        "transition-property": "background-color, width, height, opacity, border-width, border-color",
        "transition-duration": "0.5s",
    }},
    {"selector": "node.reaction", "style": {
        "background-color": "#64748b", "width": 8, "height": 8, "shape": "rectangle",
        "label": "", "opacity": 0.25,
        "transition-property": "background-color, width, height, opacity",
        "transition-duration": "0.5s",
    }},
    {"selector": "node.metabolite", "style": {
        "background-color": "#9ca3af", "width": 5, "height": 5, "shape": "ellipse",
        "label": "", "opacity": 0.2,
        "transition-property": "background-color, width, height, opacity",
        "transition-duration": "0.5s",
    }},
    {"selector": "edge.edge-default", "style": {
        "line-color": "#1e293b", "width": 0.5, "opacity": 0.15, "curve-style": "haystack",
    }},
    # === DRUG TARGET GENE ===
    {"selector": "node.gene-target", "style": {
        "background-color": "#d946ef", "width": 60, "height": 60,
        "border-width": 6, "border-color": "#fdf4ff",
        "label": "data(label)", "font-size": 20, "color": "#f0abfc",
        "font-weight": "bold", "text-valign": "top", "text-margin-y": -12,
        "opacity": 1.0, "z-index": 999,
        "transition-property": "all", "transition-duration": "0.8s",
    }},
    # === DIRECT REACTION ===
    {"selector": "node.reaction-hit", "style": {
        "background-color": "#f97316", "width": 30, "height": 30,
        "border-width": 3, "border-color": "#ffedd5",
        "label": "data(label)", "font-size": 12, "color": "#fdba74",
        "text-wrap": "ellipsis", "text-max-width": "100px",
        "opacity": 1.0, "z-index": 100,
        "transition-property": "all", "transition-duration": "0.6s", "transition-delay": "0.1s",
    }},
    # === DIRECT METABOLITE ===
    {"selector": "node.metabolite-hit", "style": {
        "background-color": "#eab308", "width": 20, "height": 20,
        "border-width": 2, "border-color": "#fef08a",
        "label": "data(label)", "font-size": 10, "color": "#fde047",
        "text-wrap": "ellipsis", "text-max-width": "80px",
        "opacity": 1.0, "z-index": 90,
        "transition-property": "all", "transition-duration": "0.6s", "transition-delay": "0.2s",
    }},
    # === SECONDARY CASCADE ===
    {"selector": "node.reaction-secondary", "style": {
        "background-color": "#c2410c", "width": 16, "height": 16,
        "label": "data(label)", "font-size": 8, "color": "#fdba74", "opacity": 0.8,
        "transition-property": "all", "transition-duration": "0.5s", "transition-delay": "0.4s",
    }},
    {"selector": "node.metabolite-secondary", "style": {
        "background-color": "#a16207", "width": 10, "height": 10, "opacity": 0.6,
        "transition-property": "background-color, width, height, opacity",
        "transition-duration": "0.5s", "transition-delay": "0.5s",
    }},
    # === PATHWAY HIT ===
    {"selector": "node.pathway-hit", "style": {
        "background-color": "#312e81", "background-opacity": 0.4,
        "border-width": 4, "border-color": "#818cf8", "border-opacity": 0.8,
        "label": "data(label)", "font-size": 16, "color": "#c7d2fe", "font-weight": "bold",
        "transition-property": "border-color, background-opacity, border-opacity, border-width",
        "transition-duration": "0.8s",
    }},
    # === HIGHLIGHTED EDGES ===
    {"selector": "edge.edge-hit", "style": {
        "line-color": "#fb923c", "width": 4, "opacity": 0.9, "z-index": 95,
        "transition-property": "line-color, width, opacity", "transition-duration": "0.4s",
    }},
    {"selector": "edge.edge-ripple", "style": {
        "line-color": "#ca8a04", "width": 1.5, "opacity": 0.5,
        "transition-property": "line-color, width, opacity",
        "transition-duration": "0.4s", "transition-delay": "0.3s",
    }},
    # === INTERACTION ===
    {"selector": "node:active", "style": {"overlay-opacity": 0.2, "overlay-color": "#ffffff"}},
    {"selector": "node:selected", "style": {"border-width": 4, "border-color": "#22d3ee"}},
]


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
                ]),
            ]),

            # --- CENTER: Network ---
            html.Div(style={"flex": 1, "position": "relative"}, children=[
                cyto.Cytoscape(
                    id="network-graph",
                    elements=BASE_ELEMENTS,
                    stylesheet=STYLESHEET,
                    layout={
                        "name": "cola",
                        "animate": True,
                        "refresh": 1,
                        "maxSimulationTime": 4000,
                        "ungrabifyWhileSimulating": False,
                        "fit": True,
                        "padding": 50,
                        "randomize": True,
                        "avoidOverlap": True,
                        "handleDisconnected": True,
                        "nodeSpacing": 45,
                        "edgeLength": 220,
                        "infinite": False,
                    },
                    style={"width": "100%", "height": "100%", "backgroundColor": CLR_BG},
                    responsive=True, minZoom=0.05, maxZoom=4.0,
                    autoRefreshLayout=False, boxSelectionEnabled=False,
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

                html.Div("Scroll to zoom · Drag to pan · Click nodes for details", style={
                    "position": "absolute", "top": "10px", "right": "15px",
                    "color": "#333", "fontSize": "10px",
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
    # The hybrid simulation returns "fluxes" (final snapshot), not "drug_fluxes"
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
    Output("network-graph", "elements"),
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

    hit_pathway_ids = set()
    for pw in hit_pathways:
        pw_id = f"pw_{pw.replace(' ', '_').replace('&', 'and')}"
        hit_pathway_ids.add(pw_id)

    # --- Update classes ---
    updated = []
    for el in BASE_ELEMENTS:
        data = el["data"]
        nid = data.get("id")
        old_cls = el.get("classes", "")

        if "source" not in data and nid:
            nt = data.get("node_type", "")
            if nt == "pathway":
                cls = "pathway-hit" if nid in hit_pathway_ids else "pathway"
            elif nt == "gene":
                cls = "gene-target" if nid in highlighted_genes else "gene"
            elif nt == "reaction":
                if nid in highlighted_rxns:
                    cls = "reaction-hit"
                elif nid in secondary_rxns:
                    cls = "reaction-secondary"
                else:
                    cls = "reaction"
            elif nt == "metabolite":
                if nid in highlighted_mets:
                    cls = "metabolite-hit"
                elif nid in secondary_mets:
                    cls = "metabolite-secondary"
                else:
                    cls = "metabolite"
            else:
                cls = old_cls
        else:
            src = data.get("source", "")
            tgt = data.get("target", "")
            if (src in highlighted_genes and tgt in highlighted_rxns) or \
               (src in highlighted_rxns and tgt in highlighted_mets):
                cls = "edge-hit"
            elif src in secondary_rxns or tgt in secondary_rxns:
                cls = "edge-ripple"
            else:
                cls = "edge-default"

        if cls != old_cls:
            updated.append({**el, "classes": cls})
        else:
            updated.append(el)

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
    enzyme_activities = sim_data.get("enzyme_activities", {})
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

    return updated, info_children, metrics, header, kinetics_children, timecourse


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
                    # Add 50% threshold line
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
    Input("network-graph", "tapNodeData"),
    State("sim-store", "data"),
)
def show_node_info(data, sim_data):
    if not data:
        return html.Span("Click any node for details", style={"color": "#444"})

    nid = data.get("id", "")
    ntype = data.get("node_type", "")
    label = data.get("label", nid)

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
        full_name = data.get("full_name", label)
        pathway = data.get("pathway", "Unknown")
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
        connected = len(_met_rxn_map.get(nid, []))
        return html.Div([
            html.B(f"🧪 {label}", style={"fontSize": "13px", "color": "#ffdd55"}),
            html.P(f"ID: {nid}", style={"color": "#666", "margin": "3px 0", "fontSize": "11px"}),
            html.P(f"Connected reactions: {connected}", style={"color": "#888", "fontSize": "11px"}),
        ])

    elif ntype == "pathway":
        rxn_count = len(_pathway_rxns.get(label, []))
        return html.Div([
            html.B(f"📂 Pathway: {label}", style={"fontSize": "13px", "color": "#8888cc"}),
            html.P(f"Reactions in pathway: {rxn_count}", style={"color": "#888", "fontSize": "11px"}),
        ])

    return html.Span(f"Node: {nid}", style={"color": "#888"})


if __name__ == "__main__":
    app.run(debug=True, port=8050)