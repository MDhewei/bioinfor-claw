---
name: gene-list-curator
description: Curate literature-backed gene or protein lists for a requested function, pathway, molecular class, or regulator type, then normalize them into a structured table or workbook. Use when the user asks for lists such as transcription factors, pathway genes, receptor families, chromatin regulators, methylation readers, or other function-specific gene sets with identifiers and references.
---

# Gene List Curator

Use this skill to build a curated, evidence-backed gene or protein list for a requested biological category.

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

Optional domain-specific columns:
- `Pathway`
- `Molecular function`
- `Reader domain annotation`
- `Marker or substrate`
- `Complex membership`
- `Confidence`

## Workflow

1. Define scope before collecting rows.
   - Organism: human by default unless requested otherwise.
   - Inclusion rule: direct experimental evidence, review-supported membership, database-defined pathway, or a hybrid.
   - Boundary rule: decide whether to include indirect regulators, paralogs with weak evidence, or inferred homologs.

2. Start with review-level discovery.
   - Search for high-quality reviews that define the class.
   - Extract canonical families, synonyms, and inclusion boundaries.
   - Note disagreements in the field before building the table.

3. Expand with targeted confirmation.
   - Search specific genes, domains, complexes, or pathway steps.
   - Use primary papers for disputed or high-value entries.
   - For pathway lists, use authoritative pathway resources if available, then verify literature support when needed.

4. Normalize identifiers.
   - `Gene name`: HGNC-style symbol when working on human genes.
   - `Protein ID`: common protein symbol or protein family label used in the literature.
   - `UniProt accession`: canonical reviewed entry where possible.
   - `PMID`: use one or two representative PMIDs per row or per class.
   - `References`: include stable links, preferably PubMed, PMC, UniProt, or official pathway pages.

5. Separate certainty levels when needed.
   - `Core`: direct, canonical members.
   - `Extended`: likely or context-dependent members.
   - `Excluded or disputed`: only if the user asks for rationale.

6. Export cleanly.
   - For Excel output, use separate sheets when the category naturally splits:
     - core vs extended
     - pathway branches
     - organism-specific sets
     - direct vs indirect evidence

## Search patterns

Adapt the query text to the requested category.

### General discovery

```text
[category] review human genes proteins
[category] comprehensive review pathway genes
[category] key regulators review
[category] domain family review
```

Examples:

```text
transcription factors review human gene list
Wnt signaling pathway human genes review
chromatin regulator families review human
RNA binding proteins review human list
```

### Specific confirmation

```text
[gene symbol] [category] PMID
[gene symbol] function in [pathway] PMID
[gene symbol] domain annotation UniProt
[gene symbol] NCBI Gene UniProtKB/Swiss-Prot
```

Examples:

```text
GATA3 transcription factor PMID
TCF7 Wnt signaling PMID
SMAD4 TGF-beta pathway PMID
TP53 transcription factor cofactor PMID
```

### Identifier normalization

```text
site:uniprot.org/uniprotkb [gene symbol] Homo sapiens UniProt
NCBI Gene [gene symbol] human UniProtKB/Swiss-Prot
[PMCID] PMID
```

## Decision rules

- Prefer a strong review to establish the universe and primary papers to resolve ambiguous members.
- Do not mix pathway members with broad upstream regulators unless the user wants both.
- If the field lacks a universally accepted boundary, state the inclusion rule in the notes.
- For “complete” lists, label them as `curated high-confidence` unless every borderline member has been reviewed.
- If the user wants a quick answer, return a concise core list first and offer an expanded sheet second.

## Output templates

### Minimal table

- `Gene name`
- `Protein ID`
- `Functional class`
- `Evidence or role`
- `References`

### Spreadsheet-ready table

- `Gene name`
- `Protein ID`
- `Organism`
- `Functional class`
- `Evidence or role`
- `UniProt accession`
- `PMID`
- `References`

## Bundled script

Use the template workbook generator when you want a fast starting point:

`python3 scripts/build_gene_list_workbook_template.py`

Edit the row list in the script for the requested category, then rerun it.
