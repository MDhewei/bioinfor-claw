# Protein Variant Mapper

Maps amino-acid substitution variants (e.g. A123V, T790M) onto a protein 3D structure. Annotates each variant with solvent exposure (SASA), secondary structure context, B-factor, and proximity to binding pockets. Produces a dual-view output: a coloured lollipop map along the sequence and an interactive 3D HTML viewer with labelled spheres.

## Trigger examples
- "Map these EGFR mutations onto the crystal structure: T790M, L858R, G719S"
- "Show me where the pathogenic variants of TP53 sit in the protein"
- "Which of my variants are buried vs surface-exposed?"
- "Are these mutations near the active site pocket?"
- "Visualize the ClinVar variants for BRCA1 on the AlphaFold model"

## Non-triggers (use a different skill)
- "Find all known variants for a gene" → use UniProt / ClinVar directly
- "Predict whether a variant is pathogenic" → use AlphaMissense or EVE
- "Align two structures" → use protein-structure-alignment

## Command template
```bash
cd protein-structure-analysis/protein-variant-mapper
pip install biopython matplotlib numpy pandas py3Dmol requests

# Map variants onto a PDB structure
python scripts/protein_variant_mapper.py \
  --pdb-id <PDB_ID> \
  --variants "A23V,G45S,K132E" \
  --outdir results/

# AlphaFold structure + auto-fetch UniProt variant classifications
python scripts/protein_variant_mapper.py \
  --uniprot <ACC> \
  --variants "T790M,L858R,G719S" \
  --fetch-uniprot-variants \
  --outdir results/

# ★ Map TCGA mutations from mutation-analysis-for-gene onto 3D structure
# Step 1: Run mutation analysis to get hotspot_details.tsv
cd ../../gene-centered-analysis/mutation-analysis-for-gene
python scripts/mutation_analysis_for_gene.py --gene <GENE> --outdir results/<gene>
# Step 2: Feed hotspots into variant mapper
cd ../../protein-structure-analysis/protein-variant-mapper
python scripts/protein_variant_mapper.py \
  --uniprot <ACC> \
  --hotspot-file ../../gene-centered-analysis/mutation-analysis-for-gene/results/<gene>/hotspot_details.tsv \
  --top-hotspots 20 \
  --outdir results/

# Load variants from file + include pocket proximity
python scripts/protein_variant_mapper.py \
  --pdb-id <PDB_ID> \
  --variants-file variants.txt \
  --pocket-residues results/<PREFIX>.pocket_residues.tsv \
  --outdir results/

# Restrict to chain A
python scripts/protein_variant_mapper.py \
  --pdb-id <PDB_ID> \
  --variants "A23V,G45S" \
  --chain A \
  --outdir results/
```

## Parameter decision guide

| Signal in user request | Parameter to set |
|---|---|
| PDB experimental structure | `--pdb-id <ID>` |
| No crystal structure available / use predicted | `--uniprot <ACC>` (AlphaFold) |
| Local PDB file | `--pdb-file path/to/protein.pdb` |
| Variants as text (few variants) | `--variants "A123V,G45S,R280*"` |
| Variants from a file | `--variants-file variants.txt` |
| **"map TCGA mutations onto structure"** / no specific variants given | **Run mutation-analysis-for-gene first**, then use `--hotspot-file` with the resulting `hotspot_details.tsv`. Add `--top-hotspots 20` to limit to top recurrent mutations |
| Variants with known clinical significance | add `--fetch-uniprot-variants` to auto-classify; also supply `--uniprot-acc` if using `--pdb-id` |
| "are these variants near the binding pocket?" | add `--pocket-residues results/<PREFIX>.pocket_residues.tsv` (run pocket module first) |
| Multi-chain protein, variants on one chain | `--chain A` |
| AlphaFold structure (no chain letter needed) | omit `--chain` |
| Variants classified in a TSV file | use `--variants-file` with `variant` and `impact_class` columns |
| Stop codons / truncations | include `*` in variant string (e.g. `R280*`); automatically classified as `stop` |

## Variant input format
Comma-separated strings: `<WT_1letter><position><MUT_1letter>`
- `A123V` — Ala123 → Val substitution
- `R280*` — Arg280 → stop codon
- `G45=` — synonymous / no change

Alternatively, provide a file (--variants-file):
- One variant per line, OR
- TSV with columns `variant` and optionally `impact_class`

`impact_class` values: `pathogenic`, `likely_pathogenic`, `benign`, `likely_benign`, `uncertain`, `custom`, `stop`

## All CLI arguments
| Argument | Required | Default | Description |
|---|---|---|---|
| `--pdb-id` / `--uniprot` / `--pdb-file` | Yes (one of) | — | Structure source |
| `--variants` / `--variants-file` / `--hotspot-file` | Yes (one of) | — | Variant input |
| `--top-hotspots` | No | `0` (all) | When using --hotspot-file, keep only top N hotspots by mutation count |
| `--chain` | No | all | Restrict to a specific chain |
| `--fetch-uniprot-variants` | No | off | Pull UniProt natural variants and auto-classify |
| `--uniprot-acc` | No | — | UniProt accession for variant fetching (if using `--pdb-id`) |
| `--pocket-residues` | No | — | TSV from visualizer pocket module for proximity annotation |
| `--probe-radius` | No | `1.4` | SASA probe radius in Å |
| `--outdir` | No | `results` | Output directory |
| `--prefix` | No | auto | File prefix |

## Outputs
- `<PREFIX>.variant_summary.tsv` — all variants with chain, resname, SASA, solvent exposure, secondary structure, B-factor, pocket proximity (+ TCGA mutation count when using --hotspot-file)
- `<PREFIX>.variant_linear_map.png/.pdf` — lollipop plot: variant positions along the backbone, coloured by impact class. When using TCGA data, marker sizes scale by mutation count.
- `<PREFIX>.variant_map.html` — interactive 3D viewer: coloured spheres at Cα, stick side-chains, floating labels with mutation count (e.g. "V600E (n=247)"); full protein shown as gray transparent cartoon

## Colour scheme
| Class | Colour |
|---|---|
| pathogenic | Red |
| likely_pathogenic | Orange |
| benign | Green |
| likely_benign | Light green |
| uncertain | Purple |
| custom | Blue |
| stop | Black |

## Chaining with other skills

**TCGA mutations → 3D structure (recommended when user doesn't provide specific variants):**
1. Run `mutation-analysis-for-gene` with `--gene <GENE>` to get `hotspot_details.tsv`
2. Pass the hotspot file to this script via `--hotspot-file`. Marker sizes and labels will automatically reflect TCGA mutation counts.
3. Use `--top-hotspots 20` to focus on the most recurrent mutations.

**Pocket proximity:**
Run `protein-structure-visualizer --modules pocket` first to generate `*.pocket_residues.tsv`, then pass it via `--pocket-residues` to annotate which variants are near binding pockets.

## Failure conditions
- "Residue N not found in structure" → position numbering mismatch between variant string and PDB resseq; check chain with `--chain`
- "wt_matches_structure = False" → the wild-type aa in the variant string doesn't match the PDB residue at that position; verify the PDB is the right isoform
- AlphaFold fetch fails → check UniProt accession at https://alphafold.ebi.ac.uk
