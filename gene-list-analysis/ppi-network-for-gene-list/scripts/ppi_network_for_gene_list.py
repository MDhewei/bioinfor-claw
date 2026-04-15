#!/usr/bin/env python3
"""
PPI Network Analysis for Gene Lists
Uses STRING database API to build and analyze protein-protein interaction networks.
Computes network metrics, identifies hubs and modules, generates visualizations.
"""

import argparse
import json
import sys
import os
from pathlib import Path
from collections import defaultdict
import math
import random

import requests
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
import matplotlib.patches as mpatches


def load_gene_list(genes_input):
    """Load genes from comma-separated string or file."""
    if os.path.isfile(genes_input):
        with open(genes_input, 'r') as f:
            genes = [line.strip() for line in f if line.strip()]
    else:
        genes = [g.strip() for g in genes_input.split(',') if g.strip()]
    return list(set(genes))  # deduplicate


def map_genes_to_string_ids(genes, species_id):
    """Map gene symbols to STRING protein IDs."""
    url = "https://string-db.org/api/json/get_string_ids"
    params = {
        "identifiers": ",".join(genes),
        "species": species_id,
    }
    try:
        resp = requests.get(url, params=params, timeout=10)
        resp.raise_for_status()
        results = resp.json()
        # results is a list of dicts with 'queryIndex', 'queryItem', 'stringId', 'ncbiTaxonId', 'queryItem', 'preferredName'
        mapping = {}
        for item in results:
            mapping[item['queryItem']] = item['stringId']
        return mapping
    except Exception as e:
        print(f"Error mapping genes via STRING API: {e}", file=sys.stderr)
        return {}


def fetch_ppi_network(string_ids, species_id, score_cutoff, network_type, expand_nodes):
    """Fetch PPI network from STRING API."""
    url = "https://string-db.org/api/json/network"
    params = {
        "identifiers": ",".join(string_ids),
        "species": species_id,
        "required_score": score_cutoff,
        "network_type": network_type,
        "add_nodes": expand_nodes,
    }
    try:
        resp = requests.get(url, params=params, timeout=15)
        resp.raise_for_status()
        data = resp.json()
        return data
    except Exception as e:
        print(f"Error fetching network from STRING API: {e}", file=sys.stderr)
        return None


def build_graph_from_api_response(api_response):
    """
    Build adjacency list from STRING API JSON response.
    Returns: dict of {node_id: {'neighbors': [list of neighbor IDs], 'metadata': {...}}}
           edges: list of (source, target, score, interaction_type)
    """
    graph = defaultdict(lambda: {'neighbors': set(), 'metadata': {}})
    edges = []

    if not api_response:
        return graph, edges

    # Extract nodes
    if 'nodes' in api_response:
        for node in api_response['nodes']:
            node_id = node['stringId']
            graph[node_id]['metadata'] = {
                'preferred_name': node.get('preferredName', node_id),
                'ncbi_taxon_id': node.get('ncbiTaxonId', ''),
            }

    # Extract edges
    if 'links' in api_response:
        for edge in api_response['links']:
            source = edge['stringId_a']
            target = edge['stringId_b']
            score = int(edge.get('score', 0)) / 1000.0  # STRING scores are 0-1000; normalize to 0-1
            nscore = score  # For simplicity, use score as weight
            # interaction type is implicit from the data; we use generic 'interaction'
            edges.append((source, target, score, 'interaction'))
            graph[source]['neighbors'].add(target)
            graph[target]['neighbors'].add(source)

    return graph, edges


def compute_degree_centrality(graph):
    """Compute degree for all nodes."""
    degree = {node: len(graph[node]['neighbors']) for node in graph}
    return degree


def compute_clustering_coefficient(graph):
    """Compute clustering coefficient for all nodes (C_v = 2*triangles / (k_v * (k_v - 1)))."""
    clustering = {}
    for node in graph:
        neighbors = list(graph[node]['neighbors'])
        k = len(neighbors)
        if k < 2:
            clustering[node] = 0.0
            continue
        # Count triangles involving this node
        triangles = 0
        for i, u in enumerate(neighbors):
            for v in neighbors[i+1:]:
                if v in graph[u]['neighbors']:
                    triangles += 1
        clustering[node] = (2 * triangles) / (k * (k - 1))
    return clustering


def compute_betweenness_centrality_approx(graph, num_samples=100):
    """Approximate betweenness centrality via random-walk sampling."""
    random.seed(42)
    np.random.seed(42)
    betweenness = {node: 0.0 for node in graph}
    nodes_list = list(graph.keys())

    if len(nodes_list) < 2:
        return betweenness

    for _ in range(num_samples):
        source, target = random.sample(nodes_list, 2)
        # BFS to find shortest path
        path = bfs_shortest_path(graph, source, target)
        if path:
            for node in path[1:-1]:  # Exclude source and target
                betweenness[node] += 1.0 / num_samples

    return betweenness


def bfs_shortest_path(graph, source, target):
    """Find shortest path using BFS."""
    if source == target:
        return [source]
    if source not in graph or target not in graph:
        return None

    visited = {source}
    queue = [(source, [source])]

    while queue:
        current, path = queue.pop(0)
        for neighbor in graph[current]['neighbors']:
            if neighbor == target:
                return path + [neighbor]
            if neighbor not in visited:
                visited.add(neighbor)
                queue.append((neighbor, path + [neighbor]))

    return None


def label_propagation_community_detection(graph, max_iter=100):
    """Detect communities via greedy label propagation."""
    labels = {node: node for node in graph}  # Initialize with unique labels
    converged = False
    iterations = 0

    while not converged and iterations < max_iter:
        converged = True
        nodes_list = list(graph.keys())
        np.random.shuffle(nodes_list)  # Random order for stability

        for node in nodes_list:
            if not graph[node]['neighbors']:
                continue
            # Find majority label among neighbors
            neighbor_labels = [labels[n] for n in graph[node]['neighbors']]
            if neighbor_labels:
                label_counts = defaultdict(int)
                for label in neighbor_labels:
                    label_counts[label] += 1
                new_label = max(label_counts, key=label_counts.get)
                if new_label != labels[node]:
                    labels[node] = new_label
                    converged = False

        iterations += 1

    # Merge small communities (< 3 nodes) with nearest large communities
    community_sizes = defaultdict(int)
    for node, label in labels.items():
        community_sizes[label] += 1

    # Relabel small communities
    node_to_community = labels.copy()
    for node, label in labels.items():
        if community_sizes[label] < 3:
            # Find closest large community
            neighbors = graph[node]['neighbors']
            if neighbors:
                neighbor_labels = [labels[n] for n in neighbors]
                largest_neighbor_label = max(set(neighbor_labels),
                                             key=lambda x: community_sizes[x])
                node_to_community[node] = largest_neighbor_label

    return node_to_community


def compute_network_metrics(graph, edges):
    """Compute all network-level and node-level metrics."""
    n_nodes = len(graph)
    n_edges = len(edges)
    density = 2 * n_edges / (n_nodes * (n_nodes - 1)) if n_nodes > 1 else 0.0

    # Connected components via DFS
    visited = set()
    n_components = 0
    for node in graph:
        if node not in visited:
            dfs(graph, node, visited)
            n_components += 1

    degree = compute_degree_centrality(graph)
    betweenness = compute_betweenness_centrality_approx(graph, num_samples=100)
    clustering = compute_clustering_coefficient(graph)
    communities = label_propagation_community_detection(graph)

    hub_threshold = max(1, int(0.1 * len(graph)))  # Top 10% by degree
    top_degree_nodes = sorted(degree.items(), key=lambda x: x[1], reverse=True)
    hub_nodes = set([n for n, d in top_degree_nodes[:hub_threshold]])

    return {
        'n_nodes': n_nodes,
        'n_edges': n_edges,
        'density': density,
        'n_components': n_components,
        'degree': degree,
        'betweenness': betweenness,
        'clustering': clustering,
        'communities': communities,
        'hub_nodes': hub_nodes,
    }


def dfs(graph, node, visited):
    """Depth-first search to mark connected component."""
    visited.add(node)
    for neighbor in graph[node]['neighbors']:
        if neighbor not in visited:
            dfs(graph, neighbor, visited)


def fruchterman_reingold_layout(graph, edges, k=1.0, iterations=50, seed=42):
    """
    Force-directed layout using Fruchterman-Reingold algorithm.
    Returns dict of {node: (x, y)} positions.
    """
    np.random.seed(seed)
    nodes = list(graph.keys())
    n = len(nodes)

    if n == 0:
        return {}
    if n == 1:
        return {nodes[0]: np.array([0.0, 0.0])}

    # Initialize positions randomly
    pos = {node: np.random.randn(2) for node in nodes}

    # Build edge dict for fast lookup
    edge_set = set()
    for u, v, _, _ in edges:
        edge_set.add((u, v))
        edge_set.add((v, u))

    t_max = 0.5
    t_min = 0.01
    t = t_max

    for iteration in range(iterations):
        # Compute forces
        forces = {node: np.array([0.0, 0.0]) for node in nodes}

        # Repulsive forces (all pairs)
        for i, u in enumerate(nodes):
            for v in nodes[i+1:]:
                delta = pos[u] - pos[v]
                dist = np.linalg.norm(delta)
                if dist > 0.01:  # Avoid division by zero
                    force_magnitude = k * k / dist
                    force = delta / dist * force_magnitude
                    forces[u] += force
                    forces[v] -= force

        # Attractive forces (edges only)
        for u, v, _, _ in edges:
            delta = pos[v] - pos[u]
            dist = np.linalg.norm(delta)
            if dist > 0.01:
                force_magnitude = dist * dist / k
                force = delta / dist * force_magnitude
                forces[u] += force
                forces[v] -= force

        # Apply forces with damping
        for node in nodes:
            force_mag = np.linalg.norm(forces[node])
            if force_mag > 0:
                displacement = forces[node] / force_mag * min(force_mag, t)
                pos[node] += displacement

        # Cool down
        t = t_max * (1.0 - iteration / iterations) + t_min * (iteration / iterations)

    return pos


def visualize_network(graph, edges, metrics, pos, outdir, node_color_by='degree', show_labels=None):
    """Visualize network with matplotlib."""
    fig, ax = plt.subplots(figsize=(14, 12))

    # Extract positions
    nodes = list(graph.keys())
    positions = np.array([pos[n] for n in nodes])

    degree = metrics['degree']
    betweenness = metrics['betweenness']
    hub_nodes = metrics['hub_nodes']

    # Determine node colors based on selection
    if node_color_by == 'degree':
        node_colors = np.array([degree[n] for n in nodes])
        cmap = 'YlOrRd'
        color_label = 'Degree'
    elif node_color_by == 'betweenness':
        node_colors = np.array([betweenness[n] for n in nodes])
        cmap = 'Blues'
        color_label = 'Betweenness Centrality'
    else:
        node_colors = np.ones(len(nodes))
        cmap = 'gray'
        color_label = 'Color'

    # Draw edges
    for u, v, score, _ in edges:
        if u in pos and v in pos:
            x_vals = [pos[u][0], pos[v][0]]
            y_vals = [pos[u][1], pos[v][1]]
            ax.plot(x_vals, y_vals, 'k-', alpha=0.2, linewidth=score * 2)

    # Draw nodes
    node_sizes = np.array([30 + degree[n] * 8 for n in nodes])
    scatter = ax.scatter(positions[:, 0], positions[:, 1],
                        s=node_sizes, c=node_colors, cmap=cmap,
                        alpha=0.8, edgecolors='black', linewidth=1.5)

    # Add colorbar
    cbar = plt.colorbar(scatter, ax=ax, label=color_label)

    # Add labels
    if show_labels is None:
        show_labels = len(nodes) < 50

    if show_labels:
        for node in nodes:
            label = graph[node]['metadata'].get('preferred_name', node)
            ax.text(pos[node][0], pos[node][1], label, fontsize=7,
                   ha='center', va='center', fontweight='bold' if node in hub_nodes else 'normal')
    else:
        # Label only hub nodes
        for node in hub_nodes:
            label = graph[node]['metadata'].get('preferred_name', node)
            ax.text(pos[node][0], pos[node][1], label, fontsize=8,
                   ha='center', va='center', fontweight='bold', color='darkred')

    ax.set_xlabel('X (force-directed layout)', fontsize=11)
    ax.set_ylabel('Y (force-directed layout)', fontsize=11)
    ax.set_title(f'Protein-Protein Interaction Network\n({metrics["n_nodes"]} nodes, {metrics["n_edges"]} edges)',
                fontsize=13, fontweight='bold')
    ax.grid(True, alpha=0.2)

    # Add legend
    hub_patch = mpatches.Patch(color='darkred', label=f'Hub genes (top 10%, n={len(hub_nodes)})')
    ax.legend(handles=[hub_patch], loc='upper right', fontsize=10)

    plt.tight_layout()
    outpath = os.path.join(outdir, 'network_plot.png')
    plt.savefig(outpath, dpi=300, bbox_inches='tight')
    print(f"Saved network plot to {outpath}")
    plt.close()


def save_results(graph, edges, metrics, outdir):
    """Save results to TSV and text files."""
    os.makedirs(outdir, exist_ok=True)

    # Save edges
    edges_df = pd.DataFrame(edges, columns=['Gene1', 'Gene2', 'Score', 'InteractionType'])
    edges_df.to_csv(os.path.join(outdir, 'network_edges.tsv'), sep='\t', index=False)
    print(f"Saved edges to {os.path.join(outdir, 'network_edges.tsv')}")

    # Save node metrics
    nodes_list = list(graph.keys())
    node_data = []
    for node in nodes_list:
        preferred_name = graph[node]['metadata'].get('preferred_name', node)
        node_data.append({
            'GeneID': node,
            'PreferredName': preferred_name,
            'Degree': metrics['degree'][node],
            'BetweennessCentrality': round(metrics['betweenness'][node], 4),
            'ClusteringCoefficient': round(metrics['clustering'][node], 4),
            'ModuleID': metrics['communities'][node],
            'IsHub': 'Yes' if node in metrics['hub_nodes'] else 'No',
        })

    nodes_df = pd.DataFrame(node_data)
    nodes_df = nodes_df.sort_values('Degree', ascending=False)
    nodes_df.to_csv(os.path.join(outdir, 'node_metrics.tsv'), sep='\t', index=False)
    print(f"Saved node metrics to {os.path.join(outdir, 'node_metrics.tsv')}")

    # Save summary
    summary_path = os.path.join(outdir, 'network_summary.txt')
    with open(summary_path, 'w') as f:
        f.write("=== PPI Network Summary ===\n\n")
        f.write(f"Nodes: {metrics['n_nodes']}\n")
        f.write(f"Edges: {metrics['n_edges']}\n")
        f.write(f"Density: {metrics['density']:.4f}\n")
        f.write(f"Connected Components: {metrics['n_components']}\n")
        f.write(f"Hub Genes (top 10%): {len(metrics['hub_nodes'])}\n\n")
        f.write("Hub Genes:\n")
        hub_data = [(node, metrics['degree'][node]) for node in metrics['hub_nodes']]
        hub_data.sort(key=lambda x: x[1], reverse=True)
        for node, deg in hub_data:
            pref_name = graph[node]['metadata'].get('preferred_name', node)
            f.write(f"  {pref_name}: degree={deg}\n")
    print(f"Saved summary to {summary_path}")


def main():
    parser = argparse.ArgumentParser(description='PPI Network Analysis for Gene Lists')
    parser.add_argument('--genes', required=True, help='Comma-separated genes or path to file (one per line)')
    parser.add_argument('--species', type=int, default=9606, help='Species taxid (9606=human, 10090=mouse)')
    parser.add_argument('--score-cutoff', type=int, default=400, help='STRING score cutoff (0-1000)')
    parser.add_argument('--network-type', choices=['functional', 'physical'], default='functional')
    parser.add_argument('--expand', type=int, default=0, help='Expand network by N neighbors')
    parser.add_argument('--outdir', default='./ppi_network_output', help='Output directory')
    parser.add_argument('--layout', choices=['spring', 'circular', 'shell'], default='spring')
    parser.add_argument('--node-color-by', choices=['degree', 'betweenness', 'input_vs_expanded'],
                       default='degree')
    parser.add_argument('--min-degree', type=int, default=1, help='Minimum degree to retain node')
    parser.add_argument('--show-labels', action='store_true', help='Show all gene labels')

    args = parser.parse_args()

    print("Loading gene list...")
    genes = load_gene_list(args.genes)
    print(f"Loaded {len(genes)} genes")

    print("Mapping genes to STRING IDs...")
    mapping = map_genes_to_string_ids(genes, args.species)
    if not mapping:
        print("Error: Failed to map any genes to STRING IDs", file=sys.stderr)
        return 1

    mapped_genes = list(mapping.values())
    print(f"Successfully mapped {len(mapped_genes)} genes")

    print("Fetching PPI network from STRING...")
    api_response = fetch_ppi_network(mapped_genes, args.species, args.score_cutoff,
                                    args.network_type, args.expand)
    if not api_response:
        print("Error: Failed to fetch network", file=sys.stderr)
        return 1

    print("Building graph...")
    graph, edges = build_graph_from_api_response(api_response)
    print(f"Built graph with {len(graph)} nodes and {len(edges)} edges")

    # Filter by min degree
    if args.min_degree > 1:
        degree = compute_degree_centrality(graph)
        keep_nodes = set([n for n, d in degree.items() if d >= args.min_degree])
        filtered_graph = {n: graph[n] for n in keep_nodes}
        filtered_edges = [(u, v, s, t) for u, v, s, t in edges if u in keep_nodes and v in keep_nodes]
        graph = filtered_graph
        edges = filtered_edges
        print(f"After min-degree filter ({args.min_degree}): {len(graph)} nodes, {len(edges)} edges")

    if len(graph) == 0:
        print("Error: No nodes remain after filtering", file=sys.stderr)
        return 1

    print("Computing network metrics...")
    metrics = compute_network_metrics(graph, edges)

    print("Computing layout...")
    if args.layout == 'spring':
        pos = fruchterman_reingold_layout(graph, edges, iterations=50)
    else:
        # Simple circular/shell layout
        nodes = list(graph.keys())
        n = len(nodes)
        if args.layout == 'circular':
            angles = np.linspace(0, 2 * np.pi, n, endpoint=False)
            pos = {nodes[i]: np.array([np.cos(angles[i]), np.sin(angles[i])]) for i in range(n)}
        else:  # shell
            pos = {node: np.array([np.random.randn(), np.random.randn()]) for node in nodes}

    print("Generating visualization...")
    visualize_network(graph, edges, metrics, pos, args.outdir,
                     node_color_by=args.node_color_by, show_labels=args.show_labels)

    print("Saving results...")
    save_results(graph, edges, metrics, args.outdir)

    print(f"\nCompleted! Results saved to {args.outdir}")
    return 0


if __name__ == '__main__':
    sys.exit(main())
