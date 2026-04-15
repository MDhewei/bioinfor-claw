---
name: design-base-editor-sgrnas
description: Design base editor sgRNAs for a given gene. Fetches exonic sequence from Ensembl, scans for all NGG (or NG/NGA) PAM sites, filters for guides where a cytosine (CBE) or adenine (ABE) falls in the editing window, scores guides for efficiency and bystander risk, and returns a ranked table. Supports CBE (C→T), ABE (A→G), and dual-base editors for human, mouse, and monkey.
---

# Base Editor sgRNA Design

## Purpose

Design sgRNAs for base editing experiments. Given a gene symbol, fetch its coding exon sequences from Ensembl, identify all valid protospacers with a target base (C for CBE, A for ABE) in the editing window, score each guide, and return a ranked table with full annotation of editing outcomes.

## Use when

- the user wants to design base editor sgRNAs for a gene
- the user mentions CBE, ABE, cytosine base editor, adenine base editor, base editing
- the user wants to introduce a specific point mutation (e.g., "install a stop codon", "correct a missense variant")
- the user wants to model a disease mutation via base editing
- the user asks about BE3, BE4max, ABE7.10, ABE8e, AncBE4max, or similar editors
- the user wants guides that avoid bystander edits

## Do not use when

- the user wants standard CRISPR KO guides (use design-sgrnas-by-gene)
- the user wants prime editor pegRNA design (different mechanism)
- the user wants guides for a non-coding region without specifying coordinates
- the user wants CRISPRi or CRISPRa (dCas9-based, no editing window)

## Supported editors and their properties

| Editor | Type | Edit | PAM | Window (protospacer positions) |
|---|---|---|---|---|
| `BE3` | CBE | C→T | NGG | 4–8 |
| `BE4max` / `AncBE4max` | CBE | C→T | NGG | 4–8 |
| `ABE7.10` | ABE | A→G | NGG | 4–7 |
| `ABE8e` | ABE | A→G | NGG | 4–8 |
| `ABE8.20m` | ABE | A→G | NGG | 4–8 |
| `NG-CBE` | CBE | C→T | NG | 4–8 |
| `NG-ABE` | ABE | A→G | NG | 4–8 |
| `SpRY-CBE` | CBE | C→T | NRN/NYN | 4–8 |
| `SpRY-ABE` | ABE | A→G | NRN/NYN | 4–8 |
| `dual` | CBE+ABE | C→T and A→G | NGG | 4–8 |

Positions are counted from the PAM-distal end (position 1 = most 5′ of protospacer).

## Inputs

- `--gene` — gene symbol (e.g., TP53, KRAS, BRCA1). Required.
- `--organism` — `human` (default), `mouse`, or `monkey`
- `--editor` — editor name or type: `CBE`, `ABE`, `dual`, `BE3`, `BE4max`, `ABE7.10`, `ABE8e`, `NG-CBE`, `NG-ABE`, `SpRY-CBE`, `SpRY-ABE` (default: `CBE`)
- `--window-start` — editing window start position (default: editor-specific)
- `--window-end` — editing window end position (default: editor-specific)
- `--pam` — PAM sequence: `NGG` (default), `NG`, or `NRN`
- `--region` — sequence scope: `cds` (coding sequence only, default), `exon` (all exons), `transcript`
- `--allow-bystander` — include guides with bystander C/A in window (default: False)
- `--max-bystander` — max allowed bystander editable bases (default: 0; use 1-2 if --allow-bystander)
- `--top-n` — number of top guides to return (default: 20)
- `--transcript-id` — restrict to a specific Ensembl transcript ID
- `--target-aa-pos` — target a specific amino acid position (e.g., 12 for KRAS G12)
- `--outdir` — output directory (default: `./base_editor_results/`)
- `--prefix` — output file prefix (default: gene symbol)

## Outputs

| File | Description |
|---|---|
| `{prefix}_base_editor_guides.tsv` | Full ranked guide table |
| `{prefix}_top_guides.tsv` | Top N guides (filtered, scored) |
| `{prefix}_editing_summary.png` | Editing window heatmap across protospacer |
| `{prefix}_guide_map.png` | Guide positions mapped onto gene exon structure |
| `{prefix}_bystander_report.tsv` | Bystander edit analysis per guide |

### Guide table columns

| Column | Description |
|---|---|
| `rank` | Final rank (1 = best) |
| `sgrna_seq` | 20 nt protospacer sequence (5′→3′) |
| `pam` | PAM sequence |
| `strand` | `+` or `-` |
| `chrom` | Chromosome |
| `start` / `end` | Genomic coordinates of protospacer |
| `exon_id` | Ensembl exon ID |
| `target_base` | Target base (C or A) |
| `window_positions` | Positions of editable bases in window |
| `edit_outcomes` | Predicted amino acid changes |
| `bystander_count` | Number of bystander editable bases in window |
| `gc_content` | GC% of protospacer |
| `efficiency_score` | Predicted on-target efficiency (0–1) |
| `bystander_penalty` | Penalty for co-editing risk (lower = cleaner) |
| `final_score` | Composite ranking score |

## Execution

```bash
# Basic CBE guide design for TP53
python scripts/design_base_editor_sgrnas.py \
  --gene TP53 \
  --editor CBE \
  --outdir results/

# ABE8e guides for KRAS, targeting codon 12
python scripts/design_base_editor_sgrnas.py \
  --gene KRAS \
  --editor ABE8e \
  --target-aa-pos 12 \
  --outdir results/

# Dual-base editor, no bystander edits allowed
python scripts/design_base_editor_sgrnas.py \
  --gene EGFR \
  --editor dual \
  --allow-bystander \
  --max-bystander 1 \
  --outdir results/

# Mouse gene, NG-CBE with relaxed PAM
python scripts/design_base_editor_sgrnas.py \
  --gene Trp53 \
  --organism mouse \
  --editor NG-CBE \
  --pam NG \
  --outdir results/
```

## Parameter decision guide

| Signal in user request | Parameter to set |
|---|---|
| "CBE", "cytosine base editor", "C to T", "C→T" | `--editor CBE` |
| "ABE", "adenine base editor", "A to G", "A→G" | `--editor ABE` |
| "BE3" | `--editor BE3` (window 4–8, NGG PAM) |
| "BE4max", "AncBE4max" | `--editor BE4max` (window 4–8, NGG PAM) |
| "ABE7.10" | `--editor ABE7.10` (window 4–7) |
| "ABE8e", "ABE8" | `--editor ABE8e` (window 4–8) |
| "NG PAM", "NG editor", "more targets" | `--editor NG-CBE` or `--editor NG-ABE` with `--pam NG` |
| "SpRY", "PAM-flexible" | `--editor SpRY-CBE` or `--editor SpRY-ABE` with `--pam NRN` |
| "both CBE and ABE", "dual base editor" | `--editor dual` |
| "clean edit", "no bystander", "single edit" | `--allow-bystander` NOT set (default off) |
| "allow some bystander", "1 bystander ok" | `--allow-bystander --max-bystander 1` |
| "codon 12", "position 12", "G12D" | `--target-aa-pos 12` |
| "coding sequence only", "CDS" | `--region cds` (default) |
| "all exons" | `--region exon` |
| "mouse", "Mus musculus" | `--organism mouse` |
| "monkey", "cynomolgus" | `--organism monkey` |
| "top 10", "best 10" | `--top-n 10` |
| "specific transcript" | `--transcript-id ENST...` |
