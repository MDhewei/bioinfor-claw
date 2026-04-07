---
name: depmap-download-data
description: Inspect DepMap releases and files, then download only the datasets required for downstream analysis. The agent should choose parameters dynamically based on the user's request instead of relying on one fixed command.
---

# DepMap Download Data

## Purpose

Prepare DepMap datasets for downstream analysis.

This skill supports:
- listing available DepMap releases
- listing files within a release
- downloading selected datasets
- downloading a full recommended bundle
- reusing cached files
- overwriting cached files when needed
- writing a manifest for downstream skills

This skill is for **data preparation only**. It does not perform biological interpretation.

## Supported datasets

Logical dataset types:
- `expression`
- `mutations`
- `copy_number`
- `essentiality`
- `metadata`

These are logical names. Actual DepMap file names may vary by release and should be matched flexibly.

## Core rule for the agent

The agent should **not** use one static command by default.

Instead, it should:
1. determine whether the user wants inspection or download
2. determine which release is needed
3. determine which datasets are needed
4. determine whether cached files are acceptable
5. set parameters accordingly
6. run only the necessary download command

Default principle:

**download only what is needed, reuse cached files whenever possible**

## When to use

Use this skill when the user asks to:
- see available DepMap releases
- inspect files in a release
- download one or more DepMap datasets
- prepare DepMap data for downstream analysis
- refresh cached DepMap files

This skill may also be used internally by another skill, such as:
- `depmap-analysis-for-gene`

## When not to use

Do not use this skill when the user wants:
- interpretation of a gene
- expression ranking results
- mutation summary for a gene
- copy number analysis
- essentiality analysis
- co-expression or co-essentiality outputs

Use a DepMap analysis skill for those tasks.

## Data source

Primary data source:
- DepMap downloads API catalog

Recommended workflow:
1. query the downloads catalog
2. identify the requested release or latest release
3. identify the matching downloadable files
4. normalize download URLs
5. download selected files
6. save files locally
7. optionally write a manifest

## Parameters

### `--release`
DepMap release name.

Use when:
- the user specifies a release
- reproducibility matters
- downstream analysis should use a fixed release

Example:
- `DepMap Public 26Q1`

If omitted:
- use the latest release

### `--outdir`
Local directory for downloaded files.

Use when:
- files need to be cached for reuse
- downstream skills need predictable file locations

Example:
- `data/depmap_26Q1`

Recommended:
- use one directory per release

### `--manifest`
Optional JSON output path describing downloaded or reused files.

Use when:
- downstream analysis needs structured file paths
- reproducibility matters
- multiple files are downloaded
- another skill may use the outputs later

Example:
- `data/depmap_26Q1/manifest.json`

### `--overwrite`
Force re-download of files even if cached copies already exist.

Use when:
- the user explicitly asks to refresh data
- cached files may be outdated or corrupted
- reproducibility requires a clean re-download

Do not use by default.

### `--list-releases`
List available DepMap releases and exit.

Use when:
- the user asks what releases are available
- the release is unclear
- the agent needs to discover current release options

### `--list-files`
List files in the selected release and exit.

Use when:
- the user asks what files are available in a release
- dataset matching needs inspection
- debugging release contents

Usually used with `--release`, but may also use the latest release if `--release` is omitted.

### `--expression`
Download the expression dataset.

Use when downstream analysis needs:
- cell line expression
- lineage expression summary
- co-expression

### `--mutations`
Download the mutations dataset.

Use when downstream analysis needs:
- mutation presence
- mutation type summary
- recurrent protein changes

### `--copy-number`
Download the copy number dataset.

Use when downstream analysis needs:
- amplification / deletion analysis
- copy number distribution
- lineage-level CN summary

### `--essentiality`
Download the essentiality dataset.

Use when downstream analysis needs:
- dependency analysis
- most dependent cell lines
- co-essentiality

### `--metadata`
Download the metadata dataset.

Use when downstream analysis needs:
- cell line names
- lineage
- primary disease
- model annotations

This should usually be included with any biological dataset.

### `--all`
Download the full recommended bundle:
- expression
- mutations
- copy number
- essentiality
- metadata

Use when:
- the user explicitly requests full DepMap preparation
- a downstream skill needs a full profile
- the user wants broad offline exploration

Do not use by default for narrow requests.

## Parameter selection rules for the agent

### If the user wants release inspection
Use:
- `--list-releases`

### If the user wants file inspection
Use:
- `--list-files`
- and `--release` if the user specified a release

### If the user wants expression-related analysis
Use:
- `--expression`
- `--metadata`

### If the user wants mutation-related analysis
Use:
- `--mutations`
- `--metadata`

### If the user wants copy-number-related analysis
Use:
- `--copy-number`
- `--metadata`

### If the user wants essentiality-related analysis
Use:
- `--essentiality`
- `--metadata`

### If the user wants co-expression
Use:
- `--expression`
- `--metadata`

### If the user wants co-essentiality
Use:
- `--essentiality`
- `--metadata`

### If the user wants a full DepMap dataset bundle
Use:
- `--all`

### If the user wants refresh / re-download
Add:
- `--overwrite`

### If downstream analysis will need structured file references
Add:
- `--manifest`

## Examples of agent decisions

### User asks:
“What DepMap releases are available?”
Agent should run:
`python scripts/depmap_download_from_api.py --outdir data/depmap --list-releases`

### User asks:
“What files are in DepMap Public 26Q1?”
Agent should run:
`python scripts/depmap_download_from_api.py --release "DepMap Public 26Q1" --outdir data/depmap_26Q1 --list-files`

### User asks:
“Download expression data for DepMap.”
Agent should run:
`python scripts/depmap_download_from_api.py --expression --metadata --release "DepMap Public 26Q1" --outdir data/depmap_26Q1`

### User asks:
“Prepare data for KRAS essentiality analysis.”
Agent should run:
`python scripts/depmap_download_from_api.py --essentiality --metadata --release "DepMap Public 26Q1" --outdir data/depmap_26Q1`

### User asks:
“Download everything for DepMap 26Q1.”
Agent should run:
`python scripts/depmap_download_from_api.py --all --release "DepMap Public 26Q1" --outdir data/depmap_26Q1 --manifest data/depmap_26Q1/manifest.json`

### User asks:
“Refresh my cached DepMap files.”
Agent should add:
- `--overwrite`

## Execution template

The agent should construct commands from the following template:

`python scripts/depmap_download_from_api.py [DATASET FLAGS] [INSPECTION FLAGS] [RELEASE OPTIONS] [OUTPUT OPTIONS] [CACHE OPTIONS]`

Where:

### Dataset flags
Choose one or more as needed:
- `--expression`
- `--mutations`
- `--copy-number`
- `--essentiality`
- `--metadata`
- `--all`

### Inspection flags
Choose when needed:
- `--list-releases`
- `--list-files`

### Release options
Choose when needed:
- `--release "<release name>"`

### Output options
Usually include:
- `--outdir <path>`

Optional:
- `--manifest <path>`

### Cache options
Choose when needed:
- `--overwrite`

## Outputs

Depending on the mode, the skill should produce:

### Release inspection output
- list of available release names

### File inspection output
- list of file names in the selected release

### Download output
- downloaded files in `outdir`
- optional manifest JSON

### Manifest content
A manifest should contain:
- dataset key
- release name
- selected file name
- download URL
- local path
- status
- note

Typical statuses:
- `downloaded`
- `exists`
- `failed`

## Best practices

- prefer selective download over `--all`
- include `--metadata` with most biological datasets
- reuse cache by default
- use `--overwrite` only when needed
- use `--manifest` when outputs will be consumed by another skill
- keep one release per directory

## Failure conditions

Fail clearly if:
- the downloads catalog cannot be reached
- the requested release does not exist
- the requested dataset cannot be matched to a file
- the output directory cannot be written
- one or more downloads fail

## Agent guidance

When using this skill, the agent should:
- infer what the user actually wants
- choose only the required parameters
- avoid one fixed default command
- prefer the minimum valid dataset bundle
- reuse cached files when possible
- expose release and file inspection when uncertainty exists

The agent should avoid:
- always calling `--all`
- always calling `--overwrite`
- downloading datasets unrelated to the user’s request
