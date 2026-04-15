---
name: ppi-network-for-gene-list
description: Build and analyze a protein-protein interaction (PPI) network for a gene list using the STRING database API. Returns network edges, computes key network metrics, identifies hub genes, detects modules, and visualizes the network.
---

# PPI Network Analysis for Gene Lists

## Purpose

Construct and analyze protein-protein interaction networks from a set of query genes using the STRING database. This skill computes network-level statistics (density, components, degree distribution), identifies hub genes and functional modules, and generates publication-quality network visualizations.

## When to Use

- You have a list of candidate genes (from GWAS, RNA-seq, CRISPR screens, etc.) and want to understand their functional relationships
- You need to prioritize genes by network centrality (hubs are often essential or disease-critical)
- You want to detect cohesive functional modules or pathways within your gene set
- You aim to expand discovery by identifying unsequenced neighbors of your query genes
- You need network metrics to assess whether your gene list represents a tightly connected functional unit

## When NOT to Use

- Your gene list is very small (< 3 genes): network statistics are not meaningful
- You need detailed mechanistic interaction information (e.g., domain interactions, phosphorylation sites) — use curated databases like PhosphoSite instead
- You work exclusively with non-model organisms not in STRING (only ~5000 species supported)

## Expected Inputs & Outputs

**Inputs:**
- Gene list (comma-separated or file with one gene per line; required)
- Species taxid (9606=human, 10090=mouse; default: 9606)
- STRING score cutoff (0–1000; default: 400 = medium confidence)
- Network expansion parameter (add N-th neighbors; default: 0)
- Visualization parameters (layout, coloring, node labeling)

**Outputs:**
- `network_edges.tsv`: tab-separated edges with interaction scores and types
- `node_metrics.tsv`: per-gene degree, betweenness, clustering coefficient, module ID, hub status
- `network_summary.txt`: global metrics (density, components, hub list)
- `network_plot.png`: publication-ready network visualization

## Procedure

### Step 1: Input Parsing
- Load gene list from file or command line
- Validate gene symbols; warn on ambiguous names

### Step 2: STRING API Queries
1. **ID Mapping:** Call `get_string_ids` endpoint to normalize gene symbols to STRING protein IDs
2. **Network Query:** Fetch interactions for the mapped IDs using `network` endpoint with user-specified:
   - `required_score` (score cutoff)
   - `network_type` (functional or physical)
   - `add_nodes` (optional expansion to N-th neighbors)

### Step 3: Graph Construction
- Parse JSON response; build adjacency list (undirected)
- Annotate edges with interaction scores and types

### Step 4: Network Metrics (Pure Python/NumPy)
- **Degree:** count incident edges per node
- **Betweenness Centrality:** approximate via 100 random-walk pairs; full computation is O(V³)
- **Clustering Coefficient:** (triangles incident to node) / C(degree, 2)
- **Global Metrics:**
  - Density = 2*|E| / (|V| * (|V| − 1))
  - Connected components via DFS

### Step 5: Hub Identification & Module Detection
- **Hubs:** top 10% by degree
- **Modules:** greedy community detection via label propagation:
  1. Initialize each node with unique label
  2. Iteratively update each node's label to majority label among neighbors
  3. Converge when labels stabilize
  4. Post-process to merge small modules (< 3 nodes) with nearest large module

### Step 6: Force-Directed Layout
- Fruchterman-Reingold algorithm in pure NumPy (50 iterations):
  - Repulsive force ∝ 1/r²
  - Attractive force ∝ r (only within edges)
  - Gradient descent with cooling

### Step 7: Visualization (Matplotlib)
- **Node size:** proportional to degree
- **Node color:** user-selected metric (degree, betweenness, or input vs. expanded)
- **Edge width:** proportional to interaction score
- **Hub labels:** draw only for hub nodes (or all if < 50 nodes)
- **Legend & colorbar:** include for interpretability

### Step 8: Output
- Write TSV files with numeric precision to 3 decimal places
- Generate summary text report
- Save PNG at 300 dpi for publication

## Key Execution Patterns

### Example 1: Quick network of a 10-gene GWAS hit list
```bash
python ppi_network_for_gene_list.py \
  --genes TP53,BRCA1,PTEN,MYC,E2F1,RB1,MDM2,CDKN1A,GADD45A,BBC3 \
  --species 9606 \
  --score-cutoff 400 \
  --outdir ./gwas_net
```

**Expected output:** compact network with ~40–80 edges; 2–3 connected components; TP53 and BRCA1 as likely hubs.

### Example 2: Expand network to include one layer of neighbors; use physical interactions only
```bash
python ppi_network_for_gene_list.py \
  --genes /path/to/crispr_hits.txt \
  --species 9606 \
  --network-type physical \
  --expand 1 \
  --score-cutoff 500 \
  --outdir ./crispr_network
```

**Expected output:** larger network integrating your query + 1 neighbor layer; node_metrics now includes many non-query genes; useful for pathway discovery.

### Example 3: Mouse data; color nodes by betweenness (centrality-based prioritization); relax stringency
```bash
python ppi_network_for_gene_list.py \
  --genes Tp53,Brca1,Pten,Myc,E2f1 \
  --species 10090 \
  --score-cutoff 300 \
  --node-color-by betweenness \
  --layout spring \
  --show-labels \
  --outdir ./mouse_net
```

**Expected output:** all nodes labeled; betweenness-colored gradient highlights non-obvious bottleneck genes.

### Example 4: Minimal filtering; focus on largest connected component
```bash
python ppi_network_for_gene_list.py \
  --genes /data/disease_genes.txt \
  --score-cutoff 200 \
  --min-degree 0 \
  --outdir ./lenient_net
```

## Parameter Decision Guide

| Scenario | Recommendation |
|----------|---|
| **Protein count < 10** | Use `--score-cutoff 500+` to focus on strongest interactions; avoid expansion |
| **GWAS/rare-variant list (50–100 genes)** | Default settings (score 400, no expansion); enough power to detect hubs |
| **RNA-seq DEG list (1000+ genes)** | Consider `--score-cutoff 600+`; pre-filter to top 500 by adjusted p-value; expansion not recommended (too dense) |
| **Candidate prioritization** | `--node-color-by betweenness`; sort output TSV by betweenness centrality |
| **Pathway annotation** | Save module assignments; cross-check each module against GO/KEGG databases |
| **Publication figure** | `--layout spring` (most aesthetically pleasing); `--show-labels` if N < 60 nodes; set `--min-degree 2` to remove singletons |
| **Exploratory analysis** | `--expand 1` to discover functionally related genes; review neighbor genes for literature mining |

## Failure Modes & Solutions

| Issue | Cause | Solution |
|-------|-------|----------|
| **"No genes found" after ID mapping** | Gene symbols not recognized by STRING (typos, aliases) | Check official gene symbols; use HGNC or MGI lookup; test a few genes in STRING web UI |
| **Empty network (0 edges)** | Score cutoff too high; genes not interacting in STRING | Lower `--score-cutoff` to 300; check if genes are in STRING database; add expansion layer |
| **Very large network (>1000 nodes)** | Too many expansion layers or very loose cutoff | Increase `--score-cutoff`; set `--min-degree 2+`; reduce expansion; pre-filter gene list |
| **Layout fails to converge** | Graph too large or disconnected; FR algorithm stuck | Reduce graph size; use `--layout circular` instead; check for isolated components separately |
| **Betweenness all zeros** | Graph is fully disconnected (many isolated nodes) | Set `--min-degree 2`; check expansion layer; may indicate genes aren't functionally related |
| **Module detection yields one giant module** | Low score cutoff creates dense hairball | Increase `--score-cutoff`; inspect module assignment quality in TSV output |
| **String encoding errors in gene names** | Non-ASCII characters in filenames or gene symbols | Ensure input file is UTF-8 encoded; use only alphanumeric and underscore in gene names |

## Technical Notes

- **API Rate Limiting:** STRING API allows ~1 request/second. For large gene lists (>500 genes), requests may be batched automatically with delays.
- **Network Type:** "functional" includes all interaction evidence (transcriptomic, curated, computational); "physical" includes only direct PPIs (binding, complexes). For discovery, use "functional"; for mechanistic detail, use "physical".
- **Reproducibility:** All calculations are deterministic except random-walk betweenness sampling. For reproducible runs, set seed = 42 (fixed in code).
- **Scalability:** Script handles up to ~5000 input genes on modest hardware (~8 GB RAM). Beyond that, consider subsampling or local tools like igraph.
