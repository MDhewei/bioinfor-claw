---
name: protein-structure-for-gene
description: Retrieve protein structure and domain information for a single gene. Resolves the gene to a UniProt accession, fetches domain/feature annotations, lists experimental PDB structures from RCSB, fetches the AlphaFold2 predicted structure entry from EBI, and renders a publication-quality linear domain map. No local data files required — all data is fetched live from UniProt, RCSB, and AlphaFold APIs.
---

# Protein Structure Analysis for Gene

## Purpose

For a given **gene symbol**, retrieve comprehensive protein structural information:

- UniProt canonical accession + protein metadata
- Domain and functional feature annotations (UniProt)
- Experimental PDB structures (RCSB) with method and resolution
- AlphaFold2 predicted structure entry (EBI)
- Publication-quality linear domain map figure

Supported modules:
- `uniprot` — protein metadata + domain/feature annotations
- `pdb` — experimental PDB structures from RCSB
- `alphafold` — AlphaFold2 predicted structure entry
- `domain_map` — linear domain figure (requires `uniprot`)
- `all` — run all modules (default)

## Reuse policy

This skill is designed for:
- Structural context for any gene of interest
- Validating domain boundaries before designing experiments
- AlphaFold availability check for a protein
- Downstream workflows that need PDB IDs or UniProt accessions
- Summarising experimental vs predicted structural coverage

This skill does **not** require local data files.
All data is fetched live from public APIs.

## Data sources

| Module      | Source             | API                                          |
|-------------|-------------------|----------------------------------------------|
| uniprot     | UniProt REST       | https://rest.uniprot.org                     |
| pdb         | RCSB PDB           | https://search.rcsb.org, https://data.rcsb.org |
| alphafold   | AlphaFold EBI      | https://alphafold.ebi.ac.uk/api              |
| domain_map  | Derived from UniProt features | —                                 |

## Inputs

### Required
- `--gene` — gene symbol (e.g. `EGFR`, `TP53`, `BRCA1`, `KRAS`)
- `--outdir` — output directory

### Optional
- `--organism`

  Default: `human`
  Also supports: `mouse`, `rat`, or any NCBI taxon ID

- `--modules`

  Default: `all`
  Comma-separated list: `uniprot`, `pdb`, `alphafold`, `domain_map`
  Example: `--modules uniprot,domain_map`

- `--max-pdb`

  Default: `20`
  Maximum number of PDB structures to retrieve.

## Outputs

### Always
- `protein_structure_summary.tsv` — one-row summary of all results
- `summary.txt` — human-readable version of the summary

### UniProt module
- `<GENE>.<ACCESSION>.features.tsv`

  All annotated features with columns:
  `type`, `description`, `start`, `end`, `length`

  Feature types include: Domain, Region, Motif, Binding site, Active site,
  Signal, Transmembrane, Coiled coil, Compositional bias, and more.

### PDB module
- `<GENE>.<ACCESSION>.pdb_structures.tsv`

  Columns: `pdb_id`, `title`, `method`, `resolution_A`, `n_protein_chains`,
  `n_atoms`, `pubmed_id`, `authors`, `deposition_date`, `rcsb_url`

### AlphaFold module
- `<GENE>.<ACCESSION>.alphafold.tsv`

  Columns: `alphafold_id`, `uniprot_accession`, `gene`, `organism`,
  `seq_length`, `model_created`, `latest_version`,
  `pdb_url`, `cif_url`, `pae_image_url`, `alphafold_page`

### Domain map module
- `<GENE>.<ACCESSION>.domain_map.png` (300 dpi)
- `<GENE>.<ACCESSION>.domain_map.pdf`

  Linear schematic with colour-coded feature types, row-stacked to prevent
  overlaps, log-rank annotation style consistent with other skills.

## Parameter decision guide

| Signal in user request | Parameter to set |
|---|---|
| "domains / features of gene X" | `--gene X` — resolves UniProt automatically |
| "what PDB structures exist for X?" | default includes PDB; add `--max-pdb 20` to see more entries |
| "AlphaFold model for X" | default includes AlphaFold; use `--uniprot <ACC>` if auto-resolution fails |
| "high-resolution structures only" | filter output TSV by `resolution_A` after running |
| "only experimental structures (no AlphaFold)" | omit `--uniprot`; script fetches PDB by default |
| "mouse gene" | `--organism mouse` |
| Gene symbol not resolving | supply `--uniprot <ACC>` directly to bypass symbol resolution |
| "show all PDB entries" | `--max-pdb 50` |

## Execution policy

Construct the command from gene, organism, and module selection.

### Command template

```
python scripts/protein_structure_for_gene.py \
  --gene <GENE> \
  [--organism <ORGANISM>] \
  [--modules <MODULES>] \
  [--max-pdb <N>] \
  --outdir <OUTDIR>
```

### Example commands

#### Full analysis (default)
```
python scripts/protein_structure_for_gene.py \
  --gene EGFR \
  --outdir results/
```

#### UniProt domains + domain map only
```
python scripts/protein_structure_for_gene.py \
  --gene TP53 \
  --modules uniprot,domain_map \
  --outdir results/
```

#### PDB structures only (up to 50)
```
python scripts/protein_structure_for_gene.py \
  --gene BRCA1 \
  --modules pdb \
  --max-pdb 50 \
  --outdir results/
```

#### AlphaFold check for mouse gene
```
python scripts/protein_structure_for_gene.py \
  --gene Trp53 \
  --organism mouse \
  --modules alphafold \
  --outdir results/
```

#### All modules for a kinase
```
python scripts/protein_structure_for_gene.py \
  --gene BRAF \
  --modules all \
  --outdir results/
```

## Domain map behaviour

- Feature types drawn: Domain, Region, Motif, Binding site, Active site,
  Signal, Transmembrane, Coiled coil
- Each feature type is assigned a distinct colour (shown in legend)
- Features are row-stacked to prevent overlap
- Feature labels are shown inside boxes when wide enough (> 4% of sequence length)
- Both PNG (300 dpi) and PDF are saved

## Failure conditions

Fail clearly if:
- `--gene` is missing or cannot be resolved in UniProt for the given organism
- `--outdir` is missing or cannot be written to
- `--modules` contains an unrecognised module name
- UniProt API returns no results (try alternate gene name or organism)
- `domain_map` is requested but no features are found in UniProt

## Agent trigger examples

**Trigger this skill when the user asks:**
- "What domains does EGFR have?"
- "Show me the protein structure information for TP53"
- "Is there an AlphaFold structure for KRAS?"
- "How many PDB structures exist for BRCA1?"
- "Draw a domain map for MYC"
- "What's the best resolution crystal structure of CDK2?"
- "Get UniProt annotation and PDB entries for PTEN"
- "What functional regions are annotated in PIK3CA?"

**Do NOT trigger this skill when the user asks:**
- "Predict a protein structure from sequence" → requires local folding tool (e.g. AlphaFold2 CLI)
- "Compare two protein structures" → structural alignment, out of scope
- "Show me TCGA expression for EGFR" → use `tcga-expression-for-gene`
- "Perform docking analysis" → requires specialised docking software

## Notes

- An internet connection is required; all data is fetched live.
- UniProt resolution prefers **reviewed (Swiss-Prot)** entries.
  For unreviewed proteins, results may be less complete.
- PDB search is based on UniProt cross-references in RCSB.
  Structures deposited without UniProt mapping will not appear.
- AlphaFold coverage is nearly complete for human proteins
  and expanding for other organisms. A `404` means no entry exists yet.
- For multi-domain proteins, the domain map may have multiple stacked rows.
  This is expected behaviour, not an error.
- Mouse gene symbols often use title-case (e.g. `Trp53`, not `TP53`).
  Specify `--organism mouse` for correct UniProt resolution.
