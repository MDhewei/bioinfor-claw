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

The skill is designed for reproducibility, not blind automation.

- Use one or more anchor sources to define the universe.
- Use scripted PubMed searches to gather repeatable evidence.
- Normalize identifiers with UniProt after the candidate list is assembled.
- Keep disputed members separate from core members.

## Processing steps

1. Define scope.
   - Organism: human by default unless requested otherwise.
   - Inclusion rule: `canonical`, `high-confidence curated`, or `broad`.
   - Boundaries: decide whether to include indirect regulators, cofactors, paralogs, or inferred members.

2. Find anchor sources.
   - Prefer a census paper, authoritative review, or official pathway/gene-family resource.
   - Use the anchor source to generate the first candidate list.

3. Run standardized literature search.
   - Use `script/run_gene_list_pipeline.py` to generate repeatable PubMed queries and store results.
   - Review the resulting JSONL file before final curation.

4. Build candidate list.
   - Extract gene symbols from the anchor source.
   - Expand only when justified by family-level or pathway-step evidence.

5. Resolve ambiguous members.
   - Run targeted searches:
     - `GENE + category + PMID`
     - `GENE + role in pathway`
     - `GENE + domain annotation + UniProt`
   - Move unresolved entries to `Extended` or `Disputed`.

6. Normalize identifiers.
   - Use `script/normalize_with_uniprot.py` on the curated CSV.
   - Prefer canonical reviewed entries for the chosen organism.

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
python3 script/run_gene_list_pipeline.py --category "transcription factors" --organism "Homo sapiens" --email "you@example.com" --output-dir ./out
python3 script/normalize_with_uniprot.py --input ./out/curated_candidates.csv --organism-id 9606 --output ./out/curated_candidates_normalized.csv
python3 script/build_gene_list_workbook.py --input ./out/curated_candidates_normalized.csv --output ./out/gene_list.xlsx
```

## References

Load these only when needed:
- [source_priority.md](references/source_priority.md)
- [query_templates.md](references/query_templates.md)
