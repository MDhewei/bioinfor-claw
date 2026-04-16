---
name: gene-list-curator
description: Curate literature-backed gene or protein lists for a requested function, pathway, molecular class, or regulator type, then normalize them into a structured table or workbook. Use when the user asks for lists such as transcription factors, pathway genes, receptor families, chromatin regulators, methylation readers, or other function-specific gene sets with identifiers and references.
---

# Gene List Curator

Use this skill to build a curated, reproducible gene or protein list for a requested biological category.

## When to use

Use this skill when the user asks for:
- a list of genes or proteins with a shared function
- pathway members
- transcription factors
- receptor or enzyme families
- chromatin regulators
- reader, writer, or eraser proteins
- a workbook, CSV, or table with identifiers and references

## Default output schema

Use these columns unless the user asks for a different schema:
- `Gene name`
- `Protein ID`
- `Organism`
- `Functional class`
- `Evidence or role`
- `UniProt accession`
- `PMID`
- `References`

Optional columns:
- `Pathway`
- `Molecular function`
- `Complex membership`
- `Marker or substrate`
- `Confidence`
- `Inclusion tier`

## Operating rule

The skill produces a **populated candidate list** directly — it does NOT
require a human-in-the-loop fill-in step. The pipeline:

1. Queries UniProt (reviewed entries) for the requested category in the
   chosen organism — this is the **primary source** of candidate rows.
2. Queries NCBI Gene as a secondary, broader source.
3. Merges the two on gene symbol (case-insensitive), preferring UniProt
   metadata and back-filling from NCBI Gene.
4. Collects PubMed bibliographic evidence (review papers etc.) into a
   side JSONL for curator reference.
5. Writes a fully populated `curated_candidates.csv`. Every row is tagged
   `Inclusion tier = Candidate`. Downstream curation (manual or LLM) can
   re-classify rows into `Core` / `Extended` / `Disputed`.

## Processing steps

1. Define scope.
   - Organism: human by default unless requested otherwise.
   - Inclusion rule: `canonical`, `high-confidence curated`, or `broad`.
   - Boundaries: decide whether to include indirect regulators, cofactors, paralogs, or inferred members.

2. Run the discovery pipeline.
   - `scripts/run_gene_list_pipeline.py` queries UniProt + NCBI Gene + PubMed in one pass.
   - The output `curated_candidates.csv` is **already populated** with all
     candidates from both databases.
   - If you want to widen the search, pass `--include-unreviewed` to add
     UniProt TrEMBL entries (noisier, more hits).

3. Review and re-classify.
   - Open `curated_candidates.csv` and re-tag the `Inclusion tier` column
     for each row: `Core` / `Extended` / `Disputed`.
   - Use the side files for evidence: `uniprot_hits.jsonl`,
     `ncbi_gene_hits.jsonl`, `pubmed_summaries.jsonl`.

4. (Optional) Resolve ambiguous members.
   - Open `pubmed_summaries.jsonl` to see which review papers covered the field.
   - For borderline rows, run targeted PubMed lookups by gene symbol.
   - Move unresolved entries to `Extended` or `Disputed`.

5. (Already done by step 2, but if you imported a CSV from elsewhere)
   - Use `scripts/normalize_with_uniprot.py` to add UniProt accessions
     to a CSV that doesn't have them.

7. Annotate roles.
   - Use short, comparable labels:
     - `sequence-specific transcription factor`
     - `pathway core transducer`
     - `coactivator`
     - `DNA-binding repressor`
     - `receptor tyrosine kinase`

8. Quality control.
   - Collapse aliases.
   - Remove duplicates.
   - Verify that every included row matches the stated inclusion rule.
   - Ensure references support the role claimed.

9. Export.
   - Use `script/build_gene_list_workbook.py` to create an Excel workbook from the curated CSV.
   - Use separate sheets for `Core`, `Extended`, or other logical groups when useful.

## Search templates

Start with general discovery:

```text
[category] review human genes proteins
[category] comprehensive review pathway genes
[category] domain family review
[category] census paper
```

Then confirm specific members:

```text
[gene symbol] [category] PMID
[gene symbol] function in [pathway] PMID
[gene symbol] domain annotation UniProt
[gene symbol] NCBI Gene UniProtKB/Swiss-Prot
```

See:
- [query_templates.md](references/query_templates.md)
- [source_priority.md](references/source_priority.md)

## Minimum quality bar

Before delivering the list, verify:
- every row has a normalized gene symbol
- every row matches the requested category under the stated inclusion rule
- duplicate aliases are collapsed
- ambiguous members are labeled or excluded
- at least one stable source is attached for each row or class

## Parameter decision guide

| Signal in user request | Decision |
|---|---|
| "quick list" / "just give me the genes" | use **Fast mode**: one anchor source, core members only, minimal columns |
| "comprehensive" / "publication-ready" / "with references" | use **Full mode**: multiple sources, per-row references, normalized identifiers |
| "human genes" (default) | `--organism "Homo sapiens"` |
| "mouse genes" | `--organism "Mus musculus"` |
| "transcription factors" | include `sequence-specific transcription factor` as functional class |
| "chromatin regulators / readers / writers / erasers" | include reader/writer/eraser domain annotation in output |
| "receptor family" | include receptor tyrosine kinase or GPCR class labels |
| "disputed members" | label as `Disputed` or `Extended` tier, keep separate from `Core` |
| "Excel workbook output" | run `build_gene_list_workbook.py` as final step |
| "normalize identifiers" | run `normalize_with_uniprot.py` after curation |

## Fast mode versus full mode

### Fast mode

Use when the user wants a quick answer.
- one anchor source
- core members only
- minimal columns
- references can be shared per class

### Full mode

Use when the user wants a comprehensive or publication-ready list.
- multiple anchor sources
- targeted primary-paper checks for ambiguous members
- normalized identifiers
- per-row references where possible
- explicit notes on inclusion boundaries

## Scripts

Use the scripts in `./script`:

- `run_gene_list_pipeline.py`
  - generate standardized PubMed searches and save JSONL/CSV evidence files
- `normalize_with_uniprot.py`
  - add UniProt accessions and canonical protein names to a curated CSV
- `build_gene_list_workbook.py`
  - convert curated CSV into an `.xlsx` workbook

Typical flow:

```bash
# One-shot discovery — produces a populated curated_candidates.csv
python3 scripts/run_gene_list_pipeline.py \
    --category "DNA methylation readers" \
    --organism "Homo sapiens" \
    --email "you@example.com" \
    --output-dir ./out

# Optional: add/refresh UniProt accessions on an existing CSV
python3 scripts/normalize_with_uniprot.py \
    --input ./out/curated_candidates.csv \
    --organism-id 9606 \
    --output ./out/curated_candidates_normalized.csv

# Export to Excel workbook
python3 scripts/build_gene_list_workbook.py \
    --input ./out/curated_candidates_normalized.csv \
    --output ./out/gene_list.xlsx
```

### Tips for getting good results

| Symptom | Try |
|---|---|
| Zero candidates returned | Broaden the category — e.g. `"methyl-CpG binding"` instead of `"DNA methylation readers"` |
| Too few candidates (< 10) | Add `--include-unreviewed` to widen UniProt to TrEMBL |
| Too many candidates (> 500) | Use a more specific category, or filter the CSV downstream |
| Non-human organism | Pass canonical scientific name, e.g. `--organism "Mus musculus"` |
| NCBI rate-limited | Increase `--delay-seconds` (default 0.4 → try 1.0) |
| Want only databases, no PubMed | Pass `--skip-pubmed` for ~10× speedup |

Exit code 2 means zero candidates were found — check the diagnostic
message printed to stderr for the most likely cause.

## References

Load these only when needed:
- [source_priority.md](references/source_priority.md)
- [query_templates.md](references/query_templates.md)
