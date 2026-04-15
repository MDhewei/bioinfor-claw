---
name: protein-structure-visualizer
description: Visualize and analyse a protein structure from RCSB PDB or AlphaFold. Modules include an interactive HTML 3D viewer (full protein or zoomed subregion), a Cα contact map, a per-residue B-factor/pLDDT plot, secondary structure timeline and composition, SASA-based binding pocket detection with 3D surface viewer, evolutionary conservation coloring (via EBI HMMER), and a STRING protein-protein interaction network. Accepts a PDB ID, UniProt accession (AlphaFold), or local PDB file.
---

# Protein Structure Visualizer

## Purpose

Produce a comprehensive visual and analytical summary of a protein 3D structure.

Supported modules:
- `view` — interactive HTML 3D viewer (opens in any browser; full protein or zoomed subregion)
- `contact_map` — Cα–Cα distance contact map with configurable threshold
- `bfactor` — per-residue B-factor (crystal structures) or pLDDT (AlphaFold) plot
- `secondary` — secondary structure timeline (helix / strand / coil) + composition chart
- `pocket` — SASA-based binding pocket / cavity detection with 3D surface viewer
- `conservation` — evolutionary conservation via EBI HMMER jackhmmer; colours structure variable→conserved (cyan→maroon)
- `ppi` — STRING protein-protein interaction network (requires `--gene`)
- `all` — run all modules (default)

Default: `all`

## Reuse policy

This skill is designed for:
- Structural visualisation of any PDB or AlphaFold protein entry
- Zoomed-in analysis of specific chains, domains, or residue ranges
- Pocket detection as a first-pass for binding site identification
- PPI context for any protein of interest
- Structural validation of results from `protein-structure-for-gene` or DE analysis

This skill does **not** require local data files for PDB/AlphaFold structures —
all structure data is fetched live from RCSB or EBI.

## Data sources

| Module       | Source                        | API / method                                  |
|--------------|-------------------------------|-----------------------------------------------|
| view         | RCSB / EBI                    | PDB file → py3Dmol HTML viewer                |
| contact_map  | PDB file (parsed)             | Biopython + Cα distances                      |
| bfactor      | PDB file (parsed)             | Biopython B-factor / AlphaFold pLDDT          |
| secondary    | DSSP (mkdssp) + PDB HELIX/SHEET records | Biopython DSSP, graceful fallback   |
| pocket       | PDB file (parsed)             | Shrake-Rupley SASA + spatial clustering       |
| conservation | EBI HMMER (jackhmmer) + UniProt | Shannon entropy over UniRef90 MSA           |
| ppi          | STRING REST API               | https://string-db.org/api                     |

## Inputs

### Structure source (provide exactly one)
- `--pdb-id` — RCSB PDB ID (e.g. `4ZJH`, `1TUP`, `6VXX`)
- `--uniprot` — UniProt accession for AlphaFold structure (e.g. `P04637`, `P00533`)
- `--pdb-file` — path to a local PDB file

### Required
- `--outdir` — output directory

### Optional — module selection
- `--modules`

  Default: `all`
  Comma-separated list: `view,contact_map,bfactor,secondary,pocket,ppi`

### Optional — region / zoom
- `--chain`

  Chain ID to focus on (e.g. `A`, `B`). Default: first chain.

- `--zoom-start`

  Residue number start for the zoomed viewer and contact map.

- `--zoom-end`

  Residue number end for the zoomed viewer and contact map.

### Optional — 3D viewer
- `--color-scheme`

  Choices: `spectrum`, `chain`, `bfactor`
  Default: `spectrum` (rainbow N→C terminus)
  Use `bfactor` for AlphaFold pLDDT confidence colouring (applied automatically for `--uniprot`).

### Optional — contact map
- `--contact-threshold`

  Default: `8.0` (Å)
  Cα–Cα distance threshold for contact definition.

### Optional — pocket detection
- `--top-pockets`

  Default: `5`
  Number of top candidate pockets to report.

- `--probe-radius`

  Default: `1.4` (Å, water probe)

### Optional — PPI
- `--gene`

  Gene symbol for STRING PPI search (required for `ppi` module).

- `--ppi-species`

  Default: `9606` (human)
  NCBI taxon ID. Common values: `10090` (mouse), `10116` (rat).

- `--ppi-score`

  Default: `400` (medium confidence)
  Minimum STRING combined score (0–1000). Use `700` for high confidence.

- `--ppi-limit`

  Default: `50`
  Maximum number of interaction partners to retrieve from STRING.

## Outputs

### Always
- `structure_summary.tsv` — one-row summary of structure and analysis metadata
- `summary.txt` — human-readable version

### Structure files (fetched automatically)
- `<PDB_ID>.pdb` or `AF_<ACCESSION>.pdb` — downloaded structure file

### view module
- `<PREFIX>.full_view.html` — interactive 3D viewer for full protein
- `<PREFIX>.zoomed_view.html` — zoomed view (when `--zoom-start` / `--zoom-end` provided)

  Open directly in Chrome, Firefox, or Safari. No installation needed.
  Controls: drag to rotate · scroll to zoom · right-click for options.

### contact_map module
- `<PREFIX>.contact_map.png` / `.pdf` — Cα contact map image
- `<PREFIX>.ca_distances.tsv` — full pairwise Cα distance matrix (TSV)

### bfactor module
- `<PREFIX>.bfactor.png` / `.pdf` — B-factor or pLDDT plot
- `<PREFIX>.bfactor.tsv` — per-residue B-factor / pLDDT values

### secondary module
- `<PREFIX>.secondary_structure.png` / `.pdf` — timeline + composition chart
- `<PREFIX>.secondary_structure.tsv` — per-residue secondary structure assignment

### pocket module
- `<PREFIX>.pockets.png` / `.pdf` — SASA distribution + top pocket ranking
- `<PREFIX>.pockets.tsv` — pocket summary (rank, residue count, centroid, chains)
- `<PREFIX>.pocket_residues.tsv` — per-residue pocket assignments
- `<PREFIX>.sasa_per_residue.tsv` — per-residue SASA values
- `<PREFIX>.pocket_surface.html` — interactive 3D viewer: each pocket shown as a coloured VDW
  surface on the semi-transparent full protein cartoon, with per-residue labels
  (3-letter amino-acid type + sequence number, e.g. "HIS57") anchored to Cα atoms;
  open in any modern browser

### conservation module
- `<PREFIX>.conservation_scores.tsv` — per-residue conservation score (0–1) and ConSurf-style grade (1–9)
- `<PREFIX>.conservation.png` / `.pdf` — bar plot coloured variable (cyan) → conserved (dark red)
- `<PREFIX>.conservation.html` — interactive 3D viewer coloured by conservation; a gradient legend bar is embedded in the footer

### ppi module
- `<PREFIX>.ppi_network.png` / `.pdf` — PPI interaction network
- `<PREFIX>.ppi_interactions.tsv` — raw STRING interaction data

## Parameter decision guide

Choose modules and options based on the specific question:

**Module selection**

| Signal in user request | Modules to activate |
|---|---|
| "show me the structure" / general 3D view | `--modules view` |
| "zoom in on residues 50–200" | `--modules view --zoom-start 50 --zoom-end 200` |
| "contact map" | `--modules contact_map` |
| "B-factor / flexibility" | `--modules bfactor` |
| "secondary structure" | `--modules secondary` |
| "pockets / binding sites" | `--modules pocket` |
| "evolutionary conservation" | `--modules conservation` |
| "protein interactions / STRING" | `--modules ppi --gene <SYMBOL>` |
| "everything" / full analysis | `--modules all --gene <SYMBOL>` |

**Structure source**

| Signal in user request | Parameter to set |
|---|---|
| PDB ID known | `--pdb-id <ID>` |
| No crystal structure, use predicted | `--uniprot <ACC>` (AlphaFold; auto-sets `--color-scheme bfactor` for pLDDT) |
| Local file | `--pdb-file path/to/file.pdb` |

**Module-specific parameters**

| Signal in user request | Parameter to set |
|---|---|
| "colour by pLDDT / confidence" (AlphaFold) | `--color-scheme bfactor` (set automatically for `--uniprot`) |
| "colour by chain" | `--color-scheme chain` |
| "colour by sequence position" | `--color-scheme spectrum` (default) |
| "contact map at 6 Å" | `--contact-threshold 6.0` |
| "top 10 pockets" | `--top-pockets 10` |
| "smaller probe" (strict pocket) | `--probe-radius 1.2` |
| "larger probe" (more surface buried) | `--probe-radius 1.8` |
| "conservation with UniProt sequence" | `--uniprot-for-conservation <ACC>` |
| "more thorough conservation MSA" | `--hmmer-iterations 3` |
| "high-confidence PPI only" | `--ppi-score 700` |
| "medium-confidence PPI" | `--ppi-score 400` (default) |
| "more PPI partners" | `--ppi-limit 100` |
| "focus on one chain" | `--chain A` |

## Execution policy

Provide exactly one structure source (`--pdb-id`, `--uniprot`, or `--pdb-file`).

### Command template

```
python scripts/protein_structure_visualizer.py \
  --pdb-id <PDB_ID> | --uniprot <ACCESSION> | --pdb-file <FILE> \
  [--modules <MODULES>] \
  [--chain <CHAIN>] \
  [--zoom-start <N>] [--zoom-end <N>] \
  [--color-scheme <SCHEME>] \
  [--contact-threshold <A>] \
  [--top-pockets <N>] \
  [--gene <GENE>] [--ppi-score <N>] [--ppi-limit <N>] \
  --outdir <OUTDIR>
```

### Example commands

#### Full analysis of a PDB entry
```
python scripts/protein_structure_visualizer.py \
  --pdb-id 1TUP \
  --gene TP53 \
  --outdir results/
```

#### AlphaFold structure with pLDDT colouring
```
python scripts/protein_structure_visualizer.py \
  --uniprot P04637 \
  --gene TP53 \
  --outdir results/
```

#### Zoomed view of a specific domain (residues 100–200, chain A)
```
python scripts/protein_structure_visualizer.py \
  --pdb-id 4ZJH \
  --modules view,contact_map,bfactor \
  --chain A \
  --zoom-start 100 \
  --zoom-end 200 \
  --outdir results/
```

#### Pocket search only
```
python scripts/protein_structure_visualizer.py \
  --pdb-id 3GFR \
  --modules pocket \
  --top-pockets 10 \
  --outdir results/
```

#### High-confidence PPI network
```
python scripts/protein_structure_visualizer.py \
  --pdb-id 1TUP \
  --modules ppi \
  --gene TP53 \
  --ppi-score 700 \
  --ppi-limit 30 \
  --outdir results/
```

#### Contact map for chain B only
```
python scripts/protein_structure_visualizer.py \
  --pdb-id 6VXX \
  --modules contact_map,secondary \
  --chain B \
  --contact-threshold 6.0 \
  --outdir results/
```

#### Local PDB file
```
python scripts/protein_structure_visualizer.py \
  --pdb-file /path/to/my_model.pdb \
  --modules view,bfactor,pocket \
  --outdir results/
```

## Module behaviour details

### view
- Generates a self-contained HTML file using py3Dmol (JavaScript, no server needed)
- `--color-scheme spectrum` colours residues rainbow from N (blue) to C (red)
- `--color-scheme bfactor` maps B-factor / pLDDT onto a colour gradient (applied automatically for AlphaFold)
- When `--zoom-start` / `--zoom-end` are given, the zoomed region is highlighted in orange with a semi-transparent surface

### contact_map
- Only Cα atoms are used; water and ligands are excluded
- When `--zoom-start` / `--zoom-end` are given, only that residue range is plotted
- The full Cα distance matrix is always saved as a TSV

### bfactor
- For AlphaFold structures, B-factor columns store pLDDT scores
- pLDDT confidence bands: > 90 (very high, blue), 70–90 (confident, cyan), 50–70 (low, yellow), < 50 (very low, orange)

### secondary
- Uses DSSP (mkdssp) when available on the system PATH
- Falls back to parsing HELIX / SHEET records from the PDB file when DSSP is unavailable
- DSSP codes are simplified: H = helix, E = strand (B/G/I/T/S → C), C = coil

### pocket
- Uses the Shrake-Rupley algorithm (Biopython) to compute per-residue SASA
- Residues with SASA < 25 Å² are considered buried
- Buried residues within 8 Å of each other are grouped into pocket clusters
- Pocket-lining residues are visualised as coloured VDW surfaces in `pocket_surface.html`:
  each pocket receives a distinct colour; the full protein is shown as a semi-transparent
  gray cartoon + faint white surface for structural context; an inline colour legend is
  embedded in the page footer
- This is a geometry-based heuristic; results should be validated with fpocket or DoGSiteScorer for publication

### conservation
- Submits the protein sequence to EBI HMMER jackhmmer against UniRef90 (public API, no key needed)
- Parses the resulting Stockholm MSA and computes per-column Shannon entropy
- Conservation score = 1 − (normalised entropy); mapped to a ConSurf-style 1–9 grade
- Structure is coloured from cyan (grade 1, variable) through white (grade 5) to dark red (grade 9, conserved)
- Pass `--uniprot-for-conservation <ACC>` to use the canonical UniProt sequence instead of the PDB chain sequence for more reliable MSA results
- HMMER jobs typically complete in 30–120 seconds; the module polls EBI automatically

### ppi
- Uses the STRING REST API (public, no authentication required)
- `--ppi-score 400` = medium confidence, `700` = high confidence, `900` = very high confidence
- Network layout uses spring layout (networkx). The query protein is shown in red.
- Edge width reflects the STRING combined score

## Failure conditions

Fail clearly if:
- No structure source is provided (`--pdb-id`, `--uniprot`, or `--pdb-file`)
- `--pdb-id` is not found at RCSB
- `--uniprot` has no AlphaFold entry
- `--pdb-file` does not exist
- `--modules ppi` is requested but `--gene` is not provided
- An unknown module name is given in `--modules`
- `--zoom-start` > `--zoom-end`
- Output directory cannot be created or written to

## Agent trigger examples

**Trigger this skill when the user asks:**
- "Visualize the 3D structure of TP53 (PDB: 1TUP)"
- "Show me the AlphaFold structure of EGFR"
- "Generate a contact map for PDB 4ZJH chain A"
- "Zoom in on residues 200–300 of the KRAS structure"
- "What pockets does this protein have?"
- "Find binding cavities in PDB 3GFR"
- "Show me the PPI network for BRCA1"
- "Which proteins interact with MYC according to STRING?"
- "Plot the pLDDT scores for the AlphaFold model of P53"
- "Show secondary structure composition for this PDB"
- "Analyse the B-factor profile of chain B in 6VXX"

**Do NOT trigger this skill when the user asks:**
- "What domains does EGFR have?" → use `protein-structure-for-gene`
- "List PDB structures for BRCA1" → use `protein-structure-for-gene`
- "Perform molecular docking" → requires specialised docking software
- "Run MD simulation" → out of scope
- "Align two protein structures" → structural alignment skill (not yet available)
- "Predict a structure from sequence" → AlphaFold CLI / ESMFold skill

## Notes

- **HTML viewer**: open the `.html` file in any modern browser. It works offline after download.
- **DSSP**: secondary structure assignment requires `mkdssp` on the system PATH.
  Install with: `conda install -c salilab dssp` or `apt install dssp`.
  The skill falls back gracefully to PDB HELIX/SHEET records if DSSP is absent.
- **py3Dmol**: required for the `view` module. Install with `pip install py3Dmol`.
- **Pocket detection**: the built-in method is fast but approximate. For rigorous pocket
  analysis, pipe the PDB file through `fpocket` or the DoGSiteScorer web server.
- **PPI**: STRING API rate limits apply. For large `--ppi-limit` values, allow extra time.
- For AlphaFold structures, `--color-scheme bfactor` is applied automatically to show pLDDT.
- Multi-chain proteins use all chains by default. Use `--chain` to focus on one chain.
