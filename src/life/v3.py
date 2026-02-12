import time
import sys
import random
import re
import cobra
import numpy as np
import hashlib
import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from Bio import SeqIO
from Bio.Seq import Seq

# --- CONFIGURATION ---
GENOME_FASTA = "data/ncbi_dataset/data/GCA_000027325.1/cds_from_genomic.fna"
MODEL_FILE = "src/iPS189.xml"

# --- Hub threshold (same as dashboard.py) ---
HUB_THRESHOLD = 15


def typing_effect(text, delay=0.005):
    for char in text:
        sys.stdout.write(char)
        sys.stdout.flush()
        time.sleep(delay)
    print()


def mutate_dna(dna_seq):
    """Introduces a mutation."""
    bases = ['A', 'T', 'C', 'G']
    seq_list = list(dna_seq)
    changes = []

    pos = random.randint(0, len(seq_list) - 1)
    original = seq_list[pos]
    new_base = random.choice([b for b in bases if b != original])
    seq_list[pos] = new_base
    changes.append(f"Position {pos}: {original} -> {new_base}")

    return Seq("".join(seq_list)), changes


def analyze_mutation_impact(original_prot, mutated_prot):
    if original_prot == mutated_prot:
        return "SILENT"
    elif "*" in mutated_prot[:-1]:
        return "NONSENSE"
    else:
        return "MISSENSE"


# =====================================================================
# NETWORK PRECOMPUTATION (mirrors dashboard.py logic)
# =====================================================================

def _stable_hash(s, mod=2**31):
    """Deterministic hash for reproducible layout."""
    return int(hashlib.md5(s.encode()).hexdigest(), 16) % mod


def precompute_mappings(model, hub_metabolites):
    """Build gene→rxn, rxn→met, rxn→pathway maps (same as dashboard.py)."""
    # Load master_map for pathway info
    try:
        rxn_df = pd.read_csv("data/master_map.csv")
    except Exception:
        rxn_df = pd.DataFrame()

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

    gene_rxn_map = {}
    for gene in model.genes:
        gene_rxn_map[gene.id] = [rxn.id for rxn in gene.reactions]

    rxn_met_map = {}
    rxn_gene_map = {}
    rxn_pathway_map = {}
    for rxn in model.reactions:
        rxn_met_map[rxn.id] = [m.id for m in rxn.metabolites if m.id not in hub_metabolites]
        rxn_gene_map[rxn.id] = [g.id for g in rxn.genes]
        rxn_pathway_map[rxn.id] = get_pathway(rxn.name)

    met_rxn_map = {}
    for rxn_id, mets in rxn_met_map.items():
        for mid in mets:
            if mid not in met_rxn_map:
                met_rxn_map[mid] = []
            met_rxn_map[mid].append(rxn_id)

    pathway_rxns = {}
    for rxn_id, pw in rxn_pathway_map.items():
        if pw not in pathway_rxns:
            pathway_rxns[pw] = []
        pathway_rxns[pw].append(rxn_id)

    return {
        "gene_rxn_map": gene_rxn_map,
        "rxn_met_map": rxn_met_map,
        "rxn_gene_map": rxn_gene_map,
        "rxn_pathway_map": rxn_pathway_map,
        "met_rxn_map": met_rxn_map,
        "pathway_rxns": pathway_rxns,
        "get_pathway": get_pathway,
        "get_met_name": get_met_name,
    }


def compute_hub_metabolites(model):
    """Identify hub metabolites (cofactors) connected to many reactions."""
    met_degree = {}
    for rxn in model.reactions:
        for met in rxn.metabolites:
            met_degree[met.id] = met_degree.get(met.id, 0) + 1
    return {mid for mid, deg in met_degree.items() if deg > HUB_THRESHOLD}


# =====================================================================
# 3D LAYOUT (same force-directed approach as dashboard.py)
# =====================================================================

def compute_3d_layout(model, maps, hub_metabolites):
    """Compute 3D positions using pathway-grouped force-directed layout."""
    np.random.seed(42)

    gene_rxn_map = maps["gene_rxn_map"]
    rxn_pathway_map = maps["rxn_pathway_map"]
    met_rxn_map = maps["met_rxn_map"]

    node_ids = []
    node_types = {}
    node_pathway = {}

    for gene in model.genes:
        node_ids.append(gene.id)
        node_types[gene.id] = "gene"
        rxn_ids = gene_rxn_map.get(gene.id, [])
        if rxn_ids:
            node_pathway[gene.id] = rxn_pathway_map.get(rxn_ids[0], "Unknown")
        else:
            node_pathway[gene.id] = "Unknown"

    seen_mets = set()
    for rxn in model.reactions:
        if not rxn.genes:
            continue
        node_ids.append(rxn.id)
        node_types[rxn.id] = "reaction"
        node_pathway[rxn.id] = rxn_pathway_map.get(rxn.id, "Unknown")

        for met in rxn.metabolites:
            if met.id in hub_metabolites:
                continue
            if met.id not in seen_mets:
                node_ids.append(met.id)
                node_types[met.id] = "metabolite"
                first_rxn = met_rxn_map.get(met.id, [None])[0]
                if first_rxn:
                    node_pathway[met.id] = rxn_pathway_map.get(first_rxn, "Unknown")
                else:
                    node_pathway[met.id] = "Unknown"
                seen_mets.add(met.id)

    # Pathway centers on a sphere
    pathways = sorted(set(node_pathway.values()))
    pw_centers = {}
    n_pw = len(pathways)
    golden_angle = np.pi * (3.0 - np.sqrt(5.0))
    sphere_radius = 80.0

    for i, pw in enumerate(pathways):
        y = 1.0 - (2.0 * i / max(n_pw - 1, 1))
        r = np.sqrt(max(0, 1.0 - y * y))
        theta = golden_angle * i
        px = sphere_radius * r * np.cos(theta)
        py = sphere_radius * y
        pz = sphere_radius * r * np.sin(theta)
        pw_centers[pw] = np.array([px, py, pz])

    # Place nodes with jitter
    positions = {}
    for nid in node_ids:
        pw = node_pathway.get(nid, "Unknown")
        center = pw_centers.get(pw, np.array([0, 0, 0]))
        nt = node_types.get(nid, "metabolite")

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
            f = (dist - 8.0) * 0.01
            force = diff / dist * f
            forces[i] += force
            forces[j] -= force

        if n_nodes < 500:
            for i in range(n_nodes):
                for j in range(i + 1, n_nodes):
                    diff = pos_array[j] - pos_array[i]
                    dist_sq = np.dot(diff, diff) + 1e-6
                    if dist_sq < 900:
                        f = -50.0 / dist_sq
                        dist = np.sqrt(dist_sq)
                        force = diff / dist * f
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
                    f = -50.0 / dist_sq
                    dist = np.sqrt(dist_sq)
                    force = diff / dist * f
                    forces[i] += force
                    forces[j] -= force

        for idx, nid in enumerate(node_list):
            pw = node_pathway.get(nid, "Unknown")
            center = pw_centers.get(pw, np.array([0, 0, 0]))
            diff = center - pos_array[idx]
            forces[idx] += diff * 0.005

        damping = 0.8 * (1.0 - iteration / 30.0)
        pos_array += forces * damping

    for idx, nid in enumerate(node_list):
        positions[nid] = pos_array[idx]

    return positions, node_types, node_pathway, pw_centers


# =====================================================================
# COMPUTE MUTATION CASCADE (which genes, rxns, mets are affected)
# =====================================================================

def compute_mutation_cascade(target_gene_id, maps, cascade_depth=2):
    """Compute affected genes, reactions, metabolites from a mutation."""
    gene_rxn_map = maps["gene_rxn_map"]
    rxn_met_map = maps["rxn_met_map"]
    rxn_pathway_map = maps["rxn_pathway_map"]
    met_rxn_map = maps["met_rxn_map"]

    highlighted_genes = {target_gene_id}
    highlighted_rxns = set()
    highlighted_mets = set()
    secondary_rxns = set()
    secondary_mets = set()
    hit_pathways = set()

    frontier_mets = set()
    for rxn_id in gene_rxn_map.get(target_gene_id, []):
        highlighted_rxns.add(rxn_id)
        pw = rxn_pathway_map.get(rxn_id, "Unknown")
        hit_pathways.add(pw)
        for mid in rxn_met_map.get(rxn_id, []):
            highlighted_mets.add(mid)
            frontier_mets.add(mid)

    visited_rxns = set(highlighted_rxns)
    for _ in range(cascade_depth):
        next_frontier = set()
        for mid in frontier_mets:
            for rxn_id in met_rxn_map.get(mid, []):
                if rxn_id in visited_rxns:
                    continue
                secondary_rxns.add(rxn_id)
                visited_rxns.add(rxn_id)
                pw = rxn_pathway_map.get(rxn_id, "Unknown")
                hit_pathways.add(pw)
                for m2 in rxn_met_map.get(rxn_id, []):
                    if m2 not in highlighted_mets:
                        secondary_mets.add(m2)
                    next_frontier.add(m2)
        frontier_mets = next_frontier

    return {
        "highlighted_genes": highlighted_genes,
        "highlighted_rxns": highlighted_rxns,
        "highlighted_mets": highlighted_mets,
        "secondary_rxns": secondary_rxns,
        "secondary_mets": secondary_mets,
        "hit_pathways": hit_pathways,
    }


# =====================================================================
# BUILD 3D FIGURE (mirrors dashboard.py build_3d_figure)
# =====================================================================

def build_3d_figure(model, positions, node_types, node_pathways, pw_centers,
                    maps, hub_metabolites, cascade, mutation_info):
    """Build a Plotly 3D scatter figure showing mutation impact on the network."""
    highlighted_genes = cascade.get("highlighted_genes", set())
    highlighted_rxns = cascade.get("highlighted_rxns", set())
    highlighted_mets = cascade.get("highlighted_mets", set())
    secondary_rxns = cascade.get("secondary_rxns", set())
    secondary_mets = cascade.get("secondary_mets", set())
    hit_pathways = cascade.get("hit_pathways", set())

    gene_rxn_map = maps["gene_rxn_map"]
    rxn_met_map = maps["rxn_met_map"]
    met_rxn_map = maps["met_rxn_map"]
    rxn_pathway_map = maps["rxn_pathway_map"]
    get_met_name = maps["get_met_name"]

    # === EDGES ===
    edge_x, edge_y, edge_z = [], [], []
    edge_hit_x, edge_hit_y, edge_hit_z = [], [], []
    edge_ripple_x, edge_ripple_y, edge_ripple_z = [], [], []

    for rxn in model.reactions:
        if rxn.id not in positions:
            continue
        rx, ry, rz = positions[rxn.id]

        for gene in rxn.genes:
            if gene.id not in positions:
                continue
            gx, gy, gz = positions[gene.id]
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
            if met.id not in positions:
                continue
            mx, my, mz = positions[met.id]
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
    if edge_hit_x:
        traces.append(go.Scatter3d(
            x=edge_hit_x, y=edge_hit_y, z=edge_hit_z,
            mode="lines", name="Direct impact",
            line=dict(color="rgba(251,146,60,0.8)", width=3),
            hoverinfo="skip", showlegend=False,
        ))
    if edge_ripple_x:
        traces.append(go.Scatter3d(
            x=edge_ripple_x, y=edge_ripple_y, z=edge_ripple_z,
            mode="lines", name="Cascade",
            line=dict(color="rgba(202,138,4,0.4)", width=1.5),
            hoverinfo="skip", showlegend=False,
        ))

    # === NODE HELPER ===
    def _add_node_trace(ids, color, size, symbol, name, opacity=1.0,
                        border_color=None, border_width=0):
        if not ids:
            return
        xs = [positions[n][0] for n in ids if n in positions]
        ys = [positions[n][1] for n in ids if n in positions]
        zs = [positions[n][2] for n in ids if n in positions]
        labels = [n for n in ids if n in positions]

        hover_texts = []
        for nid in labels:
            nt = node_types.get(nid, "")
            pw = node_pathways.get(nid, "Unknown")
            if nt == "gene":
                n_rxns = len(gene_rxn_map.get(nid, []))
                hover_texts.append(f"<b>🧬 {nid}</b><br>Type: Gene<br>Reactions: {n_rxns}<br>Pathway: {pw}")
            elif nt == "reaction":
                rxn_obj = model.reactions.get_by_id(nid) if nid in [r.id for r in model.reactions] else None
                rxn_name = rxn_obj.name if rxn_obj else nid
                hover_texts.append(f"<b>⚗️ {rxn_name}</b><br>ID: {nid}<br>Pathway: {pw}")
            elif nt == "metabolite":
                mname = get_met_name(nid)
                n_conn = len(met_rxn_map.get(nid, []))
                hover_texts.append(f"<b>🧪 {mname}</b><br>ID: {nid}<br>Connected rxns: {n_conn}")
            else:
                hover_texts.append(nid)

        marker = dict(size=size, color=color, opacity=opacity, symbol=symbol)
        if border_color and border_width > 0:
            marker["line"] = dict(color=border_color, width=border_width)

        traces.append(go.Scatter3d(
            x=xs, y=ys, z=zs,
            mode="markers+text" if size >= 10 else "markers",
            text=labels if size >= 10 else None,
            textposition="top center" if size >= 10 else None,
            textfont=dict(size=max(7, min(12, size // 2)), color=color) if size >= 10 else None,
            name=name, marker=marker,
            hovertext=hover_texts, hoverinfo="text",
        ))

    # Categorize nodes
    dimmed_genes, dimmed_rxns, dimmed_mets = [], [], []
    target_genes, hit_rxn_list, hit_met_list = [], [], []
    sec_rxn_list, sec_met_list = [], []

    for nid, ntype in node_types.items():
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

    # Dimmed background
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

    # Mutated gene (largest, most prominent)
    _add_node_trace(target_genes, "#ef4444", 18, "diamond", "🧬 Mutated gene",
                    opacity=1.0, border_color="#fecaca", border_width=3)

    # Pathway annotations
    annotations = []
    for pw, center in pw_centers.items():
        if pw == "Unknown":
            continue
        is_hit = pw in hit_pathways
        annotations.append(dict(
            x=center[0], y=center[1], z=center[2],
            text=f"<b>{pw}</b>" if is_hit else pw,
            showarrow=False,
            font=dict(size=13 if is_hit else 10,
                      color="#c7d2fe" if is_hit else "#4a4a7a"),
            bgcolor="rgba(49,46,129,0.6)" if is_hit else "rgba(20,20,50,0.4)",
            borderpad=4,
        ))

    # Layout
    axis_style = dict(
        showbackground=True, backgroundcolor="#0a0a1a",
        gridcolor="#111133", showticklabels=False, title="",
        showspikes=False, zeroline=False, showline=False,
    )

    # Build title from mutation info
    gene_id = mutation_info.get("gene_id", "?")
    impact = mutation_info.get("impact", "?")
    base_growth = mutation_info.get("base_growth", 0)
    new_growth = mutation_info.get("new_growth", 0)
    loss_pct = mutation_info.get("loss_pct", 0)

    if loss_pct < 1:
        status = "ROBUST ✅"
        status_color = "#22c55e"
    elif loss_pct > 99:
        status = "LETHAL ☠️"
        status_color = "#ef4444"
    else:
        status = f"SICK ({loss_pct:.1f}% loss) ⚠️"
        status_color = "#eab308"

    title_text = (
        f"🧬 Mutation Impact: {gene_id} ({impact})<br>"
        f"<span style='font-size:12px; color:{status_color}'>"
        f"Growth: {base_growth:.5f} → {new_growth:.5f} | {status}</span>"
    )

    fig = go.Figure(data=traces)
    fig.update_layout(
        title=dict(text=title_text, font=dict(size=16, color="#ddd"), x=0.5),
        scene=dict(
            xaxis=axis_style, yaxis=axis_style, zaxis=axis_style,
            bgcolor="#0a0a1a", annotations=annotations,
            camera=dict(eye=dict(x=1.6, y=1.6, z=1.0), up=dict(x=0, y=1, z=0)),
            aspectmode="data", dragmode="orbit",
        ),
        paper_bgcolor="#0a0a1a", plot_bgcolor="#0a0a1a",
        margin=dict(l=0, r=0, t=60, b=0),
        showlegend=True,
        legend=dict(x=0.01, y=0.99, bgcolor="rgba(10,10,26,0.8)",
                    font=dict(color="#888", size=10),
                    bordercolor="#222255", borderwidth=1, itemsizing="constant"),
        hovermode="closest",
    )

    return fig


# =====================================================================
# BUILD SUMMARY PANELS (flux comparison, pathway bar chart)
# =====================================================================

def build_summary_figures(model, maps, cascade, wt_fluxes, mut_fluxes,
                          base_growth, new_growth, mutation_info):
    """Build supplementary Plotly figures: flux comparison & pathway impact."""
    rxn_pathway_map = maps["rxn_pathway_map"]
    highlighted_rxns = cascade["highlighted_rxns"]
    secondary_rxns = cascade["secondary_rxns"]
    hit_pathways = cascade["hit_pathways"]

    # --- 1. Flux comparison for directly affected reactions ---
    affected_rxn_ids = sorted(highlighted_rxns | secondary_rxns)
    if len(affected_rxn_ids) > 30:
        # Show most-changed reactions
        changes = []
        for rid in affected_rxn_ids:
            wf = abs(wt_fluxes.get(rid, 0))
            mf = abs(mut_fluxes.get(rid, 0))
            changes.append((rid, abs(wf - mf)))
        changes.sort(key=lambda x: x[1], reverse=True)
        affected_rxn_ids = [c[0] for c in changes[:30]]

    wt_vals = [wt_fluxes.get(rid, 0) for rid in affected_rxn_ids]
    mut_vals = [mut_fluxes.get(rid, 0) for rid in affected_rxn_ids]
    rxn_labels = []
    for rid in affected_rxn_ids:
        try:
            rxn_labels.append(model.reactions.get_by_id(rid).name[:25])
        except KeyError:
            rxn_labels.append(rid[:25])

    flux_fig = go.Figure()
    flux_fig.add_trace(go.Bar(
        name="Wild-Type", x=rxn_labels, y=wt_vals,
        marker_color="#3b82f6", opacity=0.8,
    ))
    flux_fig.add_trace(go.Bar(
        name="Mutant", x=rxn_labels, y=mut_vals,
        marker_color="#ef4444", opacity=0.8,
    ))
    flux_fig.update_layout(
        title=dict(text="Flux Comparison: Wild-Type vs Mutant", font=dict(size=14, color="#ddd")),
        barmode="group",
        paper_bgcolor="#0a0a1a", plot_bgcolor="#111133",
        xaxis=dict(tickangle=-45, color="#888", gridcolor="#1a1a3e"),
        yaxis=dict(title="Flux (mmol/gDW/h)", color="#888", gridcolor="#1a1a3e"),
        legend=dict(font=dict(color="#888")),
        font=dict(color="#888"), margin=dict(l=60, r=20, t=50, b=120),
        height=450,
    )

    # --- 2. Pathway impact summary ---
    pathway_impact = {}
    for rid in highlighted_rxns | secondary_rxns:
        pw = rxn_pathway_map.get(rid, "Unknown")
        if pw == "Unknown":
            continue
        wf = abs(wt_fluxes.get(rid, 0))
        mf = abs(mut_fluxes.get(rid, 0))
        if pw not in pathway_impact:
            pathway_impact[pw] = {"total_wt": 0, "total_mut": 0, "n_rxns": 0}
        pathway_impact[pw]["total_wt"] += wf
        pathway_impact[pw]["total_mut"] += mf
        pathway_impact[pw]["n_rxns"] += 1

    pw_names = sorted(pathway_impact.keys(), key=lambda p: pathway_impact[p]["total_wt"], reverse=True)
    if len(pw_names) > 15:
        pw_names = pw_names[:15]

    pw_wt = [pathway_impact[p]["total_wt"] for p in pw_names]
    pw_mut = [pathway_impact[p]["total_mut"] for p in pw_names]
    pw_labels = [f"{p} ({pathway_impact[p]['n_rxns']}rxn)" for p in pw_names]

    pathway_fig = go.Figure()
    pathway_fig.add_trace(go.Bar(
        name="Wild-Type", y=pw_labels, x=pw_wt,
        orientation="h", marker_color="#3b82f6", opacity=0.8,
    ))
    pathway_fig.add_trace(go.Bar(
        name="Mutant", y=pw_labels, x=pw_mut,
        orientation="h", marker_color="#ef4444", opacity=0.8,
    ))
    pathway_fig.update_layout(
        title=dict(text="Pathway-Level Flux Impact", font=dict(size=14, color="#ddd")),
        barmode="group",
        paper_bgcolor="#0a0a1a", plot_bgcolor="#111133",
        xaxis=dict(title="Total |Flux|", color="#888", gridcolor="#1a1a3e"),
        yaxis=dict(color="#888", autorange="reversed"),
        legend=dict(font=dict(color="#888")),
        font=dict(color="#888"), margin=dict(l=200, r=20, t=50, b=40),
        height=max(300, len(pw_names) * 30 + 100),
    )

    # --- 3. Growth rate gauge ---
    loss_pct = mutation_info.get("loss_pct", 0)
    gauge_fig = go.Figure(go.Indicator(
        mode="gauge+number+delta",
        value=new_growth,
        delta={"reference": base_growth, "decreasing": {"color": "#ef4444"}},
        gauge={
            "axis": {"range": [0, base_growth * 1.2], "tickcolor": "#888"},
            "bar": {"color": "#ef4444" if loss_pct > 50 else ("#eab308" if loss_pct > 10 else "#22c55e")},
            "bgcolor": "#111133",
            "bordercolor": "#222255",
            "steps": [
                {"range": [0, base_growth * 0.5], "color": "rgba(239,68,68,0.2)"},
                {"range": [base_growth * 0.5, base_growth], "color": "rgba(234,179,8,0.2)"},
            ],
            "threshold": {
                "line": {"color": "#3b82f6", "width": 3},
                "thickness": 0.8, "value": base_growth,
            },
        },
        title={"text": f"Growth Rate (Loss: {loss_pct:.1f}%)", "font": {"size": 14, "color": "#ddd"}},
        number={"font": {"size": 28, "color": "#ddd"}, "suffix": " h⁻¹"},
    ))
    gauge_fig.update_layout(
        paper_bgcolor="#0a0a1a", font=dict(color="#888"),
        height=250, margin=dict(l=30, r=30, t=50, b=20),
    )

    return flux_fig, pathway_fig, gauge_fig


# =====================================================================
# MAIN SIMULATION + VISUALIZATION
# =====================================================================

def run_strict_evolution():
    print("\n" + "=" * 60)
    print("   🧬 STRICT EVOLUTIONARY SIMULATOR: FORCED BOTTLENECK")
    print("=" * 60)

    print("[PHASE 1] LOADING METABOLIC MODEL...")
    model = cobra.io.read_sbml_model(MODEL_FILE)
    model.objective = "Biomass"
    for rxn in model.exchanges:
        rxn.lower_bound = -10.0

    # Baseline
    print("   > Calculating Healthy Flux...")
    solution = model.optimize()
    base_growth = solution.objective_value
    print(f"   ✅ Healthy Growth Rate: {base_growth:.5f}")

    fluxes = solution.fluxes.abs()
    wt_fluxes_dict = {rxn.id: float(solution.fluxes[rxn.id]) for rxn in model.reactions}

    # Find active genes
    active_genes = []
    for gene in model.genes:
        is_active = False
        for rxn in gene.reactions:
            if fluxes[rxn.id] > 0.0001:
                is_active = True
                break
        if is_active:
            cleaned_id = gene.id[:2] + "_" + gene.id[2:] if not gene.id.startswith("MG_") else gene.id
            active_genes.append(cleaned_id)

    print(f"   > Active Genes Detected: {len(active_genes)} (Skipping lazy ones)")

    # Pick target
    target_id = random.choice(active_genes)
    print(f"   > Target Selected: {target_id} (Carries significant flux)")

    # Load genome data
    print(f"\n[PHASE 2] EXPOSING {target_id} TO MUTAGEN...")
    original_records = list(SeqIO.parse(GENOME_FASTA, "fasta"))
    candidate_record = None

    for rec in original_records:
        if target_id in rec.description:
            candidate_record = rec
            break

    impact = "SILENT"
    new_growth = base_growth
    mutation_log = ""

    if not candidate_record:
        print(f"   ❌ Error: {target_id} not found in FASTA. Using generic mutation.")
    else:
        original_dna = candidate_record.seq
        mutated_dna, logs = mutate_dna(original_dna)
        mutation_log = logs[0]
        print(f"   ⚠️  MUTATION: {mutation_log}")

        orig_prot = original_dna.translate(table=4, to_stop=True)
        mut_prot = mutated_dna.translate(table=4, to_stop=True)
        impact = analyze_mutation_impact(orig_prot, mut_prot)
        print(f"   Mutation Type: {impact}")

        # Apply penalty
        print("\n[PHASE 3] APPLYING PHYSIOLOGICAL PENALTY...")
        target_gene = model.genes.get_by_id(target_id.replace("_", ""))

        if impact == "SILENT":
            print("   > Effect: None.")
        elif impact == "MISSENSE":
            print("   > Effect: Enzyme damaged. Limiting to 10% capacity.")
            max_flux_needed = 0
            affected_rxns = []
            for rxn in target_gene.reactions:
                current_flux = fluxes[rxn.id]
                if current_flux > max_flux_needed:
                    max_flux_needed = current_flux
                affected_rxns.append(rxn)

            new_limit = max_flux_needed * 0.1
            for rxn in affected_rxns:
                if rxn.upper_bound > 0:
                    rxn.upper_bound = new_limit
                if rxn.lower_bound < 0:
                    rxn.lower_bound = -new_limit

        elif impact == "NONSENSE":
            print("   > Effect: Knockout (Zero capacity).")
            target_gene.knock_out()

        mut_solution = model.optimize()
        new_growth = mut_solution.objective_value if mut_solution.status == "optimal" else 0.0

    # Collect mutant fluxes
    try:
        mut_sol = model.optimize()
        mut_fluxes_dict = {rxn.id: float(mut_sol.fluxes[rxn.id]) for rxn in model.reactions}
    except Exception:
        mut_fluxes_dict = {rxn.id: 0.0 for rxn in model.reactions}

    loss = (base_growth - new_growth) / base_growth * 100 if base_growth > 0 else 0

    print("\n" + "=" * 60)
    print("   EVOLUTIONARY OUTCOME")
    print("=" * 60)
    print(f"   Gene:       {target_id}")
    print(f"   Mutation:   {impact}")
    print(f"   Old Rate:   {base_growth:.5f}")
    print(f"   New Rate:   {new_growth:.5f}")

    if loss < 1:
        print("   Status:     ROBUST (Redundant pathways found)")
    elif loss > 99:
        print("   Status:     DEAD (Lethal Mutation)")
    else:
        print(f"   Status:     SICK ({loss:.1f}% Fitness Loss)")
    print("=" * 60)

    # --- VISUALIZATION ---
    print("\n[PHASE 4] BUILDING 3D NETWORK VISUALIZATION...")
    hub_metabolites = compute_hub_metabolites(model)
    print(f"   Hub metabolites excluded: {len(hub_metabolites)}")

    maps = precompute_mappings(model, hub_metabolites)
    print("   Precomputed gene/rxn/met mappings.")

    # Use the gene ID as it appears in the model (no underscore)
    model_gene_id = target_id.replace("_", "")

    print("   Computing 3D layout (this may take a moment)...")
    positions, node_types, node_pathways, pw_centers = compute_3d_layout(
        model, maps, hub_metabolites
    )
    print(f"   Layout computed for {len(positions)} nodes.")

    cascade = compute_mutation_cascade(model_gene_id, maps, cascade_depth=2)
    print(f"   Cascade: {len(cascade['highlighted_rxns'])} direct rxns, "
          f"{len(cascade['secondary_rxns'])} cascade rxns, "
          f"{len(cascade['hit_pathways'])} pathways")

    mutation_info = {
        "gene_id": target_id,
        "impact": impact,
        "base_growth": base_growth,
        "new_growth": new_growth,
        "loss_pct": loss,
        "mutation_log": mutation_log,
    }

    # Build figures
    network_fig = build_3d_figure(
        model, positions, node_types, node_pathways, pw_centers,
        maps, hub_metabolites, cascade, mutation_info,
    )

    flux_fig, pathway_fig, gauge_fig = build_summary_figures(
        model, maps, cascade, wt_fluxes_dict, mut_fluxes_dict,
        base_growth, new_growth, mutation_info,
    )

    # Save as interactive HTML
    print("\n[PHASE 5] SAVING VISUALIZATIONS...")

    network_fig.write_html("debug/mutation_network_3d.html", include_plotlyjs="cdn")
    print("   ✅ 3D Network → debug/mutation_network_3d.html")

    # Combined summary page
    from plotly.subplots import make_subplots
    import plotly.io as pio

    # Write individual figures
    flux_fig.write_html("debug/mutation_flux_comparison.html", include_plotlyjs="cdn")
    print("   ✅ Flux Comparison → debug/mutation_flux_comparison.html")

    pathway_fig.write_html("debug/mutation_pathway_impact.html", include_plotlyjs="cdn")
    print("   ✅ Pathway Impact → debug/mutation_pathway_impact.html")

    gauge_fig.write_html("debug/mutation_growth_gauge.html", include_plotlyjs="cdn")
    print("   ✅ Growth Gauge → debug/mutation_growth_gauge.html")

    # Build a combined dashboard HTML
    network_html = pio.to_html(network_fig, full_html=False, include_plotlyjs=False)
    flux_html = pio.to_html(flux_fig, full_html=False, include_plotlyjs=False)
    pathway_html = pio.to_html(pathway_fig, full_html=False, include_plotlyjs=False)
    gauge_html = pio.to_html(gauge_fig, full_html=False, include_plotlyjs=False)

    combined_html = f"""<!DOCTYPE html>
<html>
<head>
    <title>🧬 Mutation Impact Dashboard — {target_id}</title>
    <script src="https://cdn.plot.ly/plotly-latest.min.js"></script>
    <style>
        body {{ background: #0a0a1a; color: #ddd; font-family: 'Segoe UI', sans-serif; margin: 0; padding: 20px; }}
        h1 {{ text-align: center; color: #fff; }}
        .subtitle {{ text-align: center; color: #888; font-size: 14px; margin-bottom: 20px; }}
        .grid {{ display: grid; grid-template-columns: 1fr 1fr; gap: 16px; max-width: 1600px; margin: 0 auto; }}
        .full-width {{ grid-column: 1 / -1; }}
        .card {{ background: #0d0d24; border: 1px solid #222255; border-radius: 8px; padding: 12px; overflow: hidden; }}
        .info-bar {{ display: flex; justify-content: center; gap: 40px; margin: 16px 0; }}
        .info-item {{ text-align: center; }}
        .info-label {{ font-size: 10px; color: #556; text-transform: uppercase; letter-spacing: 1px; }}
        .info-value {{ font-size: 22px; font-weight: bold; }}
        .status-robust {{ color: #22c55e; }}
        .status-sick {{ color: #eab308; }}
        .status-dead {{ color: #ef4444; }}
    </style>
</head>
<body>
    <h1>🧬 Evolutionary Mutation Impact Dashboard</h1>
    <div class="subtitle">Gene: <b>{target_id}</b> | Mutation: <b>{impact}</b> | {mutation_log}</div>

    <div class="info-bar">
        <div class="info-item">
            <div class="info-label">Wild-Type Growth</div>
            <div class="info-value" style="color: #3b82f6;">{base_growth:.5f}</div>
        </div>
        <div class="info-item">
            <div class="info-label">Mutant Growth</div>
            <div class="info-value" style="color: #ef4444;">{new_growth:.5f}</div>
        </div>
        <div class="info-item">
            <div class="info-label">Fitness Loss</div>
            <div class="info-value {'status-robust' if loss < 1 else ('status-dead' if loss > 99 else 'status-sick')}">{loss:.1f}%</div>
        </div>
        <div class="info-item">
            <div class="info-label">Affected Pathways</div>
            <div class="info-value" style="color: #f97316;">{len(cascade['hit_pathways'])}</div>
        </div>
        <div class="info-item">
            <div class="info-label">Direct Reactions</div>
            <div class="info-value" style="color: #eab308;">{len(cascade['highlighted_rxns'])}</div>
        </div>
        <div class="info-item">
            <div class="info-label">Cascade Reactions</div>
            <div class="info-value" style="color: #c2410c;">{len(cascade['secondary_rxns'])}</div>
        </div>
    </div>

    <div class="grid">
        <div class="card full-width">{network_html}</div>
        <div class="card">{gauge_html}</div>
        <div class="card">{pathway_html}</div>
        <div class="card full-width">{flux_html}</div>
    </div>
</body>
</html>"""

    with open("debug/mutation_dashboard.html", "w") as f:
        f.write(combined_html)
    print("   ✅ Combined Dashboard → debug/mutation_dashboard.html")

    print("\n" + "=" * 60)
    print("   🎉 VISUALIZATION COMPLETE")
    print("   Open debug/mutation_dashboard.html in your browser")
    print("=" * 60 + "\n")


if __name__ == "__main__":
    run_strict_evolution()