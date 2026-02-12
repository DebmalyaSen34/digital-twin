import dash
from dash import html, dcc, Input, Output, State, ctx, no_update, callback
import dash_cytoscape as cyto
import cobra
import pandas as pd
import json

# --- Load extra layouts for Cytoscape ---
cyto.load_extra_layouts()

# --- DRUG DATABASE ---
DRUG_DB = {
    "Control (No Treatment)": {"targets": [], "desc": "Baseline metabolism."},
    "Trimethoprim": {"targets": ["MG228"], "desc": "Inhibits DHFR (Folate synthesis)."},
    "Methotrexate": {"targets": ["MG228", "MG006"], "desc": "Blocks dTMP production."},
    "Fosmidomycin": {"targets": ["MG066"], "desc": "Inhibits Isoprenoid synthesis."},
    "Cerulenin": {"targets": ["MG212", "MG114"], "desc": "Inhibits Fatty Acid Synthase."},
    "Mupirocin": {"targets": ["MG345"], "desc": "Inhibits Isoleucyl-tRNA synthetase."},
    "Generic Glycolysis Inhibitor": {"targets": ["MG041", "MG429"], "desc": "Blocks Glucose import."},
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

# Exclude metabolites connected to > 15 reactions (cofactors/currency metabolites)
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

# --- Pre-compute pathway → reaction groups (for pathway-based clustering) ---
_pathway_rxns = {}
for rxn_id, pw in _rxn_pathway_map.items():
    if pw not in _pathway_rxns:
        _pathway_rxns[pw] = []
    _pathway_rxns[pw].append(rxn_id)


def run_fba(drug_name, efficacy):
    """Run FBA for wild-type and drug condition."""
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

    return wt_growth, wt_fluxes, drug_growth, drug_fluxes


# --- Cache FBA results to avoid redundant runs ---
_fba_cache = {}


def run_fba_cached(drug_name, efficacy):
    key = (drug_name, round(efficacy, 3))
    if key not in _fba_cache:
        _fba_cache[key] = run_fba(drug_name, efficacy)
        # Keep cache small
        if len(_fba_cache) > 50:
            oldest = next(iter(_fba_cache))
            del _fba_cache[oldest]
    return _fba_cache[key]


def build_cytoscape_elements():
    """Build the organism graph with pathway grouping as parent nodes."""
    elements = []
    seen_mets = set()
    node_count = 0
    edge_count = 0

    # --- Add pathway compound (parent) nodes for clustering ---
    pathway_ids = {}
    for pw in _pathway_rxns:
        pw_id = f"pw_{pw.replace(' ', '_').replace('&', 'and')}"
        pathway_ids[pw] = pw_id
        elements.append({
            "data": {
                "id": pw_id,
                "label": pw,
                "node_type": "pathway",
            },
            "classes": "pathway",
        })

    # --- Gene nodes ---
    for gene in model.genes:
        # Assign gene to the pathway of its first reaction, for clustering
        rxn_ids = _gene_rxn_map.get(gene.id, [])
        parent = None
        if rxn_ids:
            pw = _rxn_pathway_map.get(rxn_ids[0], "Unknown")
            parent = pathway_ids.get(pw)

        node_data = {"id": gene.id, "label": gene.id, "node_type": "gene"}
        if parent:
            node_data["parent"] = parent
        elements.append({
            "data": node_data,
            "classes": "gene",
        })
        node_count += 1

    # --- Reaction nodes ---
    for rxn in model.reactions:
        if not rxn.genes:
            continue
        pathway = _rxn_pathway_map.get(rxn.id, "Unknown")
        parent = pathway_ids.get(pathway)

        node_data = {
            "id": rxn.id,
            "label": rxn.name[:30],
            "node_type": "reaction",
            "pathway": pathway,
            "full_name": rxn.name,
        }
        if parent:
            node_data["parent"] = parent
        elements.append({
            "data": node_data,
            "classes": "reaction",
        })
        node_count += 1

        # gene → reaction edges
        for gene in rxn.genes:
            edge_id = f"eg_{gene.id}_{rxn.id}"
            elements.append({
                "data": {
                    "id": edge_id,
                    "source": gene.id,
                    "target": rxn.id,
                    "edge_type": "gene_rxn",
                },
                "classes": "edge-default",
            })
            edge_count += 1

        # reaction → metabolite edges (no hubs)
        for met in rxn.metabolites:
            mid = met.id
            if mid in HUB_METABOLITES:
                continue
            if mid not in seen_mets:
                met_name = get_met_name(mid)
                # Assign metabolite to pathway of first connected reaction
                first_rxn = _met_rxn_map.get(mid, [None])[0]
                met_parent = None
                if first_rxn:
                    mp = _rxn_pathway_map.get(first_rxn, "Unknown")
                    met_parent = pathway_ids.get(mp)

                met_data = {"id": mid, "label": met_name[:25], "node_type": "metabolite"}
                if met_parent:
                    met_data["parent"] = met_parent
                elements.append({
                    "data": met_data,
                    "classes": "metabolite",
                })
                seen_mets.add(mid)
                node_count += 1

            edge_id = f"em_{rxn.id}_{mid}"
            elements.append({
                "data": {
                    "id": edge_id,
                    "source": rxn.id,
                    "target": mid,
                    "edge_type": "rxn_met",
                },
                "classes": "edge-default",
            })
            edge_count += 1

    print(f"Graph built: {node_count} nodes, {edge_count} edges, {len(pathway_ids)} pathway groups")
    return elements


print("Building base graph (one-time)...")
BASE_ELEMENTS = build_cytoscape_elements()


# --- STYLESHEET ---
STYLESHEET = [
    # === PATHWAY COMPOUND NODES (clusters) ===
    {"selector": "node.pathway", "style": {
        "background-color": "#1e1b4b",  # Very dark indigo
        "background-opacity": 0.3,
        "border-width": 1,
        "border-color": "#4f46e5",      # Indigo 600
        "border-style": "dashed",
        "shape": "round-rectangle",
        "label": "data(label)",
        "font-size": 14,                # Larger font
        "color": "#818cf8",             # Indigo 400
        "text-valign": "top",
        "text-halign": "center",
        "text-margin-y": -10,
        "text-opacity": 0.8,
        "padding": 20,
        "compound-sizing-wrt-labels": "include",
    }},

    # === DEFAULT DIMMED ===
    {"selector": "node.gene", "style": {
        "background-color": "#0ea5e9",  # Sky Blue
        "width": 12, "height": 12,
        "shape": "diamond",
        "label": "",
        "opacity": 0.3,                 # Keep dim
        "transition-property": "background-color, width, height, opacity, border-width, border-color",
        "transition-duration": "0.5s",
    }},
    {"selector": "node.reaction", "style": {
        "background-color": "#64748b",  # Slate
        "width": 8, "height": 8,
        "shape": "rectangle",
        "label": "",
        "opacity": 0.25,
        "transition-property": "background-color, width, height, opacity",
        "transition-duration": "0.5s",
    }},
    {"selector": "node.metabolite", "style": {
        "background-color": "#9ca3af",  # Light Gray
        "width": 5, "height": 5,
        "shape": "ellipse",
        "label": "",
        "opacity": 0.2,
        "transition-property": "background-color, width, height, opacity",
        "transition-duration": "0.5s",
    }},
    {"selector": "edge.edge-default", "style": {
        "line-color": "#1e293b",        # Slate 800
        "width": 0.5,
        "opacity": 0.15,
        "curve-style": "haystack",
    }},

    # === DRUG TARGET GENE (The Source) ===
    {"selector": "node.gene-target", "style": {
        "background-color": "#d946ef",  # Fuchsia 500 (Neon Pink)
        "width": 60, "height": 60,
        "border-width": 6, "border-color": "#fdf4ff", # White glow
        "label": "data(label)",
        "font-size": 20,
        "color": "#f0abfc",             # Light Pink Text
        "font-weight": "bold",
        "text-valign": "top", "text-margin-y": -12,
        "opacity": 1.0,
        "z-index": 999,
        "transition-property": "all",
        "transition-duration": "0.8s",
    }},

    # === DIRECT REACTION (primary impact) ===
    {"selector": "node.reaction-hit", "style": {
        "background-color": "#f97316",  # Orange 500
        "width": 30, "height": 30,
        "border-width": 3, "border-color": "#ffedd5", # Pale Orange
        "label": "data(label)",
        "font-size": 12,
        "color": "#fdba74",             # Light Orange Text
        "text-wrap": "ellipsis", "text-max-width": "100px",
        "opacity": 1.0,
        "z-index": 100,
        "transition-property": "all",
        "transition-duration": "0.6s",
        "transition-delay": "0.1s",
    }},

    # === DIRECT METABOLITE (primary impact) ===
    {"selector": "node.metabolite-hit", "style": {
        "background-color": "#eab308",  # Yellow 500
        "width": 20, "height": 20,
        "border-width": 2, "border-color": "#fef08a", # Pale Yellow
        "label": "data(label)",
        "font-size": 10,
        "color": "#fde047",             # Light Yellow Text
        "text-wrap": "ellipsis", "text-max-width": "80px",
        "opacity": 1.0,
        "z-index": 90,
        "transition-property": "all",
        "transition-duration": "0.6s",
        "transition-delay": "0.2s",
    }},

    # === SECONDARY CASCADE (Ripple) ===
    {"selector": "node.reaction-secondary", "style": {
        "background-color": "#c2410c",  # Dark Orange (700)
        "width": 16, "height": 16,
        "label": "data(label)",
        "font-size": 8, "color": "#fdba74",
        "opacity": 0.8,
        "transition-property": "all",
        "transition-duration": "0.5s",
        "transition-delay": "0.4s",
    }},
    {"selector": "node.metabolite-secondary", "style": {
        "background-color": "#a16207",  # Dark Yellow (700)
        "width": 10, "height": 10,
        "opacity": 0.6,
        "transition-property": "background-color, width, height, opacity",
        "transition-duration": "0.5s",
        "transition-delay": "0.5s",
    }},

    # === PATHWAY HIT (Parent Container Highlight) ===
    {"selector": "node.pathway-hit", "style": {
        "background-color": "#312e81",  # Dark Indigo fill
        "background-opacity": 0.4,
        "border-width": 4,
        "border-color": "#818cf8",      # Bright Indigo border
        "border-opacity": 0.8,
        "label": "data(label)",
        "font-size": 16,
        "color": "#c7d2fe",             # Bright text
        "font-weight": "bold",
        "transition-property": "border-color, background-opacity, border-opacity, border-width",
        "transition-duration": "0.8s",
    }},

    # === HIGHLIGHTED EDGES ===
    {"selector": "edge.edge-hit", "style": {
        "line-color": "#fb923c",        # Orange 400
        "width": 4,
        "opacity": 0.9,
        "z-index": 95,
        "transition-property": "line-color, width, opacity",
        "transition-duration": "0.4s",
    }},
    {"selector": "edge.edge-ripple", "style": {
        "line-color": "#ca8a04",        # Yellow 600
        "width": 1.5,
        "opacity": 0.5,
        "transition-property": "line-color, width, opacity",
        "transition-duration": "0.4s",
        "transition-delay": "0.3s",
    }},

    # === INTERACTION ===
    {"selector": "node:active", "style": {
        "overlay-opacity": 0.2, "overlay-color": "#ffffff",
    }},
    {"selector": "node:selected", "style": {
        "border-width": 4, "border-color": "#22d3ee", # Cyan select
    }},
]


def _legend_row(color, symbol, label):
    return html.Div(style={"display": "flex", "alignItems": "center", "gap": "6px", "marginBottom": "3px"}, children=[
        html.Span(symbol, style={"color": color, "fontSize": "12px", "width": "14px", "textAlign": "center"}),
        html.Span(label, style={"color": "#666"}),
    ])

# --- DASH APP ---
app = dash.Dash(
    __name__,
    title="iPS189 Drug Simulator",
    update_title=None,
)

app.layout = html.Div(
    style={"backgroundColor": CLR_BG, "minHeight": "100vh", "fontFamily": "'Segoe UI', sans-serif",
           "overflow": "hidden"},
    children=[
        # Hidden store for cached FBA results
        dcc.Store(id="fba-store", data={}),

        # --- HEADER ---
        html.Div(style={
            "background": "linear-gradient(135deg, #0a0a1a 0%, #111144 100%)",
            "padding": "15px 30px", "borderBottom": "1px solid #222255",
            "display": "flex", "alignItems": "center", "justifyContent": "space-between",
        }, children=[
            html.Div([
                html.H1("🦠 iPS189 Whole-Cell Simulation", style={
                    "color": "#fff", "margin": 0, "fontSize": "24px",
                }),
                html.P(
                    "System-wide drug response · Mycoplasma genitalium metabolic network",
                    style={"color": "#556", "margin": "2px 0 0 0", "fontSize": "12px"},
                ),
            ]),
            # Quick stats in header
            html.Div(id="header-stats", style={
                "display": "flex", "gap": "30px", "alignItems": "center",
            }),
        ]),

        html.Div(style={"display": "flex", "height": "calc(100vh - 70px)"}, children=[

            # --- LEFT PANEL ---
            html.Div(style={
                "width": "280px", "flexShrink": 0, "padding": "16px",
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

                html.Label("Inhibition Power:", style={"color": "#888", "fontSize": "11px"}),
                html.Div(style={"display": "flex", "alignItems": "baseline", "gap": "8px"}, children=[
                    html.Div(id="efficacy-val", style={
                        "color": "#ff6600", "fontSize": "22px", "fontWeight": "bold",
                    }, children="95%"),
                ]),
                dcc.Slider(
                    id="efficacy-slider", min=0, max=100, step=5, value=95,
                    marks={0: {"label": "0%", "style": {"color": "#444"}},
                           100: {"label": "100%", "style": {"color": "#444"}}},
                ),

                html.Label("Cascade Depth:", style={"color": "#888", "fontSize": "11px", "marginTop": "8px"}),
                dcc.Slider(
                    id="cascade-slider", min=0, max=3, step=1, value=1,
                    marks={i: {"label": str(i), "style": {"color": "#444"}} for i in range(4)},
                ),

                html.Hr(style={"borderColor": "#1a1a3e", "margin": "16px 0"}),

                # Drug info card
                html.Div(id="drug-info", style={
                    "padding": "12px", "backgroundColor": "#111133",
                    "borderRadius": "8px", "border": "1px solid #222255",
                    "fontSize": "12px",
                }),

                html.Hr(style={"borderColor": "#1a1a3e", "margin": "16px 0"}),

                # Metrics
                html.Div(id="metrics-panel"),

                html.Hr(style={"borderColor": "#1a1a3e", "margin": "16px 0"}),

                # Legend
                html.Div(style={"fontSize": "10px", "color": "#555"}, children=[
                    html.Div("LEGEND", style={"color": "#444", "fontWeight": "bold",
                                               "letterSpacing": "1px", "marginBottom": "6px"}),
                    _legend_row("#d946ef", "◆", "Target Gene"),
                    _legend_row("#f97316", "■", "Direct Reaction"),
                    _legend_row("#eab308", "●", "Direct Metabolite"),
                    _legend_row("#c2410c", "■", "Cascade Reaction"),
                    _legend_row("#a16207", "●", "Cascade Metabolite"),
                    _legend_row("#334155", "●", "Unaffected"),
                ]),
            ]),

            # --- RIGHT: Network ---
            html.Div(style={"flex": 1, "position": "relative"}, children=[
                cyto.Cytoscape(
                    id="network-graph",
                    elements=BASE_ELEMENTS,
                    stylesheet=STYLESHEET,
                    layout={
                        "name": "fcose",
                        "animate": True,
                        "animationDuration": 800,
                        "quality": "proof",
                        "randomize": True,
                        "nodeRepulsion": 8000,
                        "idealEdgeLength": 50,
                        "edgeElasticity": 0.45,
                        "gravity": 0.25,
                        "gravityRange": 3.8,
                        "nestingFactor": 0.1,
                        "numIter": 5000,
                        "tile": True,
                        "tilingPaddingVertical": 20,
                        "tilingPaddingHorizontal": 20,
                        "packComponents": True,
                    },
                    style={"width": "100%", "height": "100%", "backgroundColor": CLR_BG},
                    responsive=True,
                    minZoom=0.03,
                    maxZoom=5.0,
                    autoRefreshLayout=False,
                    boxSelectionEnabled=False,
                ),

                # Floating tooltip
                html.Div(id="node-tooltip", style={
                    "position": "absolute", "bottom": "12px", "left": "12px",
                    "color": "#ddd", "padding": "10px 14px",
                    "backgroundColor": "rgba(10, 10, 30, 0.95)",
                    "borderRadius": "8px", "border": "1px solid #333366",
                    "fontSize": "12px", "maxWidth": "320px",
                    "backdropFilter": "blur(10px)",
                    "pointerEvents": "none",
                    "transition": "opacity 0.3s ease",
                }),

                # Zoom controls hint
                html.Div("Scroll to zoom · Drag to pan · Click nodes for details", style={
                    "position": "absolute", "top": "10px", "right": "15px",
                    "color": "#333", "fontSize": "10px",
                }),
            ]),
        ]),
    ],
)


@app.callback(
    Output("efficacy-val", "children"),
    Input("efficacy-slider", "value"),
)
def update_efficacy_label(val):
    return f"{val}%"


@app.callback(
    Output("fba-store", "data"),
    Input("drug-select", "value"),
    Input("efficacy-slider", "value"),
)
def compute_fba(drug_name, efficacy_pct):
    """Run FBA and store results. Separate callback so it only runs when inputs change."""
    efficacy = efficacy_pct / 100.0
    wt_growth, wt_fluxes, drug_growth, drug_fluxes = run_fba_cached(drug_name, efficacy)

    # Only store the fluxes we actually need (affected reactions) to keep Store small
    drug_data = DRUG_DB[drug_name]
    relevant_rxn_ids = set()
    for tid in drug_data["targets"]:
        for rxn_id in _gene_rxn_map.get(tid, []):
            relevant_rxn_ids.add(rxn_id)
            for mid in _rxn_met_map.get(rxn_id, []):
                for r2 in _met_rxn_map.get(mid, []):
                    relevant_rxn_ids.add(r2)

    wt_subset = {rid: float(wt_fluxes.get(rid, 0)) for rid in relevant_rxn_ids}
    drug_subset = {rid: float(drug_fluxes.get(rid, 0)) for rid in relevant_rxn_ids}

    return {
        "wt_growth": float(wt_growth),
        "drug_growth": float(drug_growth),
        "wt_fluxes": wt_subset,
        "drug_fluxes": drug_subset,
    }


@app.callback(
    Output("network-graph", "elements"),
    Output("drug-info", "children"),
    Output("metrics-panel", "children"),
    Output("header-stats", "children"),
    Input("drug-select", "value"),
    Input("cascade-slider", "value"),
    Input("fba-store", "data"),
)
def update_graph(drug_name, cascade_depth, fba_data):
    drug_data = DRUG_DB[drug_name]

    wt_growth = fba_data.get("wt_growth", 0)
    drug_growth = fba_data.get("drug_growth", 0)

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

    # Compute pathway parent IDs that should glow
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
            # Node
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
            # Edge
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
    info_children = [
        html.Div(f"💊 {drug_name}", style={
            "fontWeight": "bold", "color": "#ff6600", "fontSize": "14px", "marginBottom": "6px",
        }),
        html.Div(drug_data['desc'], style={"color": "#999", "fontSize": "11px"}),
    ]
    if drug_data["targets"]:
        info_children.append(
            html.Div(style={"marginTop": "8px", "display": "flex", "flexWrap": "wrap", "gap": "4px"}, children=[
                html.Span("Targets: ", style={"color": "#777", "fontSize": "11px", "marginRight": "2px"}),
                *[html.Span(t, style={
                    "backgroundColor": "#ff0040", "color": "#fff",
                    "padding": "1px 8px", "borderRadius": "10px",
                    "fontSize": "10px", "fontWeight": "bold",
                }) for t in drug_data["targets"]]
            ])
        )
    if hit_pathways - {"Unknown"}:
        info_children.append(
            html.Div(style={"marginTop": "8px"}, children=[
                html.Span("Pathways: ", style={"color": "#777", "fontSize": "11px"}),
                html.Span(", ".join(sorted(hit_pathways - {"Unknown"})),
                          style={"color": "#ff8844", "fontSize": "11px"}),
            ])
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
        ]

    return updated, info_children, metrics, header


def _metric_card(title, value, subtitle=None, accent_color="#aaa"):
    children = [
        html.Div(title, style={"color": "#556", "fontSize": "9px", "textTransform": "uppercase",
                                "letterSpacing": "0.5px"}),
        html.Div(value, style={"color": accent_color, "fontSize": "18px", "fontWeight": "bold",
                                "marginTop": "1px"}),
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
    State("fba-store", "data"),
)
def show_node_info(data, fba_data):
    if not data:
        return html.Span("Click any node for details", style={"color": "#444"})

    nid = data.get("id", "")
    ntype = data.get("node_type", "")
    label = data.get("label", nid)

    wt_fluxes = fba_data.get("wt_fluxes", {})
    drug_fluxes = fba_data.get("drug_fluxes", {})

    if ntype == "gene":
        rxn_ids = _gene_rxn_map.get(nid, [])
        pathways = set()
        for rid in rxn_ids:
            pathways.add(_rxn_pathway_map.get(rid, "Unknown"))
        return html.Div([
            html.B(f"🧬 Gene: {nid}", style={"fontSize": "14px", "color": "#ff6688"}),
            html.P(f"Reactions: {len(rxn_ids)}", style={"margin": "4px 0", "color": "#aaa"}),
            html.P(f"Pathways: {', '.join(sorted(pathways))}", style={"color": "#888", "fontSize": "11px"}),
        ])

    elif ntype == "reaction":
        full_name = data.get("full_name", label)
        pathway = data.get("pathway", "Unknown")
        flux_wt = wt_fluxes.get(nid, 0)
        flux_drug = drug_fluxes.get(nid, 0)
        pct = ((flux_drug - flux_wt) / flux_wt * 100) if abs(flux_wt) > 1e-9 else 0
        fc = "#ff4444" if pct < -10 else ("#ffaa00" if pct < -1 else "#44ff44")
        return html.Div([
            html.B(f"⚗️ {full_name}", style={"fontSize": "13px", "color": "#ffaa44"}),
            html.P(f"Pathway: {pathway}", style={"color": "#888", "margin": "3px 0", "fontSize": "11px"}),
            html.Div(style={"display": "flex", "gap": "15px", "marginTop": "4px", "fontSize": "12px"}, children=[
                html.Span(f"WT: {flux_wt:.4f}", style={"color": "#44aaff"}),
                html.Span(f"Drug: {flux_drug:.4f}", style={"color": fc}),
                html.Span(f"({pct:+.1f}%)", style={"color": fc, "fontWeight": "bold"}),
            ]),
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