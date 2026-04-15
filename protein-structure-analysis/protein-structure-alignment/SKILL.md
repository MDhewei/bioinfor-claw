# Protein Structure Alignment

Compare two protein structures by superimposing them via Cα atoms (Biopython Superimposer). Computes global RMSD, highlights divergent regions, and produces a dual-colour interactive 3D viewer.

## Trigger examples
- "Align the apo and holo structures of adenylate kinase"
- "Compare 1AKE and 4AKE, what is the RMSD?"
- "Superimpose the AlphaFold model for P00940 against the crystal structure 1TIM"
- "How similar are the WT and mutant protein structures?"
- "Find flexible regions between two conformations of the same protein"

## Non-triggers (use a different skill)
- "Show me the structure of a single protein" → use protein-structure-visualizer
- "What domains does my protein have?" → use protein-structure-for-gene
- "Align protein sequences" → use standard sequence alignment tools

## Command template
```bash
cd protein-structure-analysis/protein-structure-alignment
pip install biopython matplotlib numpy pandas py3Dmol requests

# Two PDB IDs
python scripts/protein_structure_alignment.py \
  --pdb1 <PDB_ID_1> \
  --pdb2 <PDB_ID_2> \
  --outdir results/

# AlphaFold predicted vs experimental
python scripts/protein_structure_alignment.py \
  --pdb1 <PDB_ID> \
  --uniprot2 <UNIPROT_ACC> \
  --outdir results/

# Local files
python scripts/protein_structure_alignment.py \
  --file1 wt.pdb \
  --file2 mutant.pdb \
  --outdir results/

# Restrict to one chain and a residue range
python scripts/protein_structure_alignment.py \
  --pdb1 1AKE --pdb2 4AKE \
  --chain1 A --chain2 A \
  --res-start 1 --res-end 180 \
  --outdir results/
```

## All CLI arguments
| Argument | Required | Default | Description |
|---|---|---|---|
| `--pdb1` / `--uniprot1` / `--file1` | Yes (one of) | — | Structure 1 source |
| `--pdb2` / `--uniprot2` / `--file2` | Yes (one of) | — | Structure 2 source |
| `--chain1` | No | all chains | Chain to use in structure 1 |
| `--chain2` | No | all chains | Chain to use in structure 2 |
| `--res-start` | No | — | First residue number to include |
| `--res-end` | No | — | Last residue number to include |
| `--outdir` | No | `results` | Output directory |
| `--prefix` | No | auto | File prefix (default: `<ID1>_vs_<ID2>`) |

## Parameter decision guide

| Signal in user request | Parameter to set |
|---|---|
| Two PDB IDs to compare | `--pdb1 <ID> --pdb2 <ID>` |
| AlphaFold vs experimental structure | `--uniprot2 <ACC>` (or `--uniprot1`) |
| Local PDB files | `--file1 wt.pdb --file2 mutant.pdb` |
| Multi-chain structure, compare one chain | `--chain1 A --chain2 A` |
| Compare only a domain / specific region | `--res-start 1 --res-end 200` |
| Apo vs holo (ligand-induced conformational change) | use full chain; inspect per-residue RMSD for loop regions |
| Homolog comparison (different residue numbering) | omit `--chain` — let the script pair by resseq; if numbering differs badly, use `--res-start/--res-end` for the conserved core |
| "RMSD of just the kinase domain" | `--res-start 50 --res-end 350 --chain1 A --chain2 A` |

## Outputs
- `<PREFIX>.alignment_summary.tsv` — global RMSD, n Cα pairs, divergent residue counts
- `<PREFIX>.per_residue_rmsd.tsv` — per-residue Cα distance after superimposition
- `<PREFIX>.per_residue_rmsd.png/.pdf` — bar plot coloured green (<1Å) / orange (1–3Å) / red (>3Å)
- `<PREFIX>.superimposed.html` — interactive 3D viewer: structure 1 in blue, structure 2 in red; divergent residues (>3Å) highlighted in orange sticks
- `<PREFIX>.superimposed.pdb` — coordinates of structure 2 in the superimposed frame

## How it works
1. Fetches both structures from RCSB PDB or AlphaFold EBI (or reads local files)
2. Extracts Cα atoms from both; pairs residues by (chain, resseq)
3. Superimposes structure 2 onto structure 1 using Biopython's `Superimposer` (SVD-based rigid-body alignment)
4. Computes global RMSD and per-residue Cα distance
5. Writes all atoms of structure 2 in the new (rotated/translated) frame
6. Visualises both structures in py3Dmol with divergent regions highlighted

## Failure conditions
- "No common Cα pairs found" → chains or residue numbering don't overlap; specify `--chain1/--chain2` or `--res-start/--res-end`
- Very high RMSD on short alignments → try aligning a specific domain with `--res-start/--res-end`
- AlphaFold fetch fails → UniProt accession not found; verify at https://alphafold.ebi.ac.uk
