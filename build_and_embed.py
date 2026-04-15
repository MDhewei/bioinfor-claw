#!/usr/bin/env python3
"""
bioinfor-claw  —  Build & Embed Script
========================================
Run from your bioinfor-claw folder:

    python3 build_and_embed.py

It reads all your SKILL.md files and writes a SINGLE self-contained
bioinfor-claw-ready.html with every skill embedded inside.

Open bioinfor-claw-ready.html in any browser — instant, no server needed.
"""

import os, sys, json, re
from pathlib import Path
from datetime import datetime

SKILL_SET_META = {
    "public-datasets-access-and-download":  {"emoji": "🗄️",  "desc": "DepMap, GTEx, TCGA, GEO download & cache"},
    "multiomics-data-analysis":             {"emoji": "🧫",  "desc": "RNA-seq, ATAC-seq, CNV, multi-omics"},
    "crispr-design-and-analysis":           {"emoji": "✂️",  "desc": "sgRNA design, screen QC & hit analysis"},
    "gene-list-analysis":                   {"emoji": "📋",  "desc": "GO, pathway enrichment, annotation"},
    "gene-centered-analysis":               {"emoji": "🔬",  "desc": "Expression, mutation, CNV, dependency"},
    "protein-structure-analysis":           {"emoji": "🧩",  "desc": "Domain mapping, mutation annotation, PDB"},
    "machine-learning-and-deep-learning":   {"emoji": "🤖",  "desc": "Classification, clustering, SHAP"},
    "bioinformatics-plot-generator":        {"emoji": "📊",  "desc": "Publication-quality figures & charts"},
    "paper-search-and-digest":             {"emoji": "📄",  "desc": "PubMed retrieval & summarization"},
    "lab-search-and-track":                {"emoji": "🏛️", "desc": "Find & monitor labs by field"},
}
SKIP = {'.git','__pycache__','Assets','assets','node_modules',
        'web','web_results','.github','docs','tests','examples','.DS_Store'}

def scan(repo):
    tree_ui, flat = {}, {}
    print(f"\n🔍  Scanning: {repo}\n{'─'*58}")
    for sd in sorted(Path(repo).iterdir()):
        if not sd.is_dir() or sd.name in SKIP or sd.name.startswith('.'): continue
        found = {}
        for kd in sorted(sd.iterdir()):
            if not kd.is_dir() or kd.name.startswith('.'): continue
            md = kd / 'SKILL.md'
            if md.exists():
                txt = md.read_text(encoding='utf-8', errors='replace')
                found[kd.name] = txt
                flat[f"{sd.name}/{kd.name}"] = txt
                print(f"  ✅  {sd.name}/{kd.name}  ({len(txt):,} chars)")
            else:
                print(f"  ⬜  {sd.name}/{kd.name}  (no SKILL.md yet)")
        if found:
            m = SKILL_SET_META.get(sd.name, {"emoji":"⚙️","desc":sd.name.replace('-',' ').title()})
            tree_ui[sd.name] = {"emoji":m["emoji"],"desc":m["desc"],"skills":list(found.keys())}
    return tree_ui, flat

def build_html(tree_ui, flat, template_path):
    n  = len(flat)
    kb = sum(len(v) for v in flat.values()) // 1024
    ts = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

    # Serialize safely for JS embedding
    tree_js  = json.dumps(tree_ui, ensure_ascii=False)
    flat_js  = json.dumps(flat,    ensure_ascii=False)

    # The block we inject — goes right before </head>
    inject = f"""
  <!-- ══════════════════════════════════════════════════════════
       bioinfor-claw Embedded Skills  ·  {n} skills  ·  {kb} KB
       Generated: {ts}
       Re-run build_and_embed.py to update after adding skills
  ═══════════════════════════════════════════════════════════ -->
  <script id="embedded-skills">
  // Skills are embedded directly — no network requests needed
  window.BUNDLED_SKILL_TREE    = {tree_js};
  window.BUNDLED_LOADED_SKILLS  = {flat_js};
  window.SKILLS_BUNDLE_META    = {{
    generated   : "{ts}",
    totalSkills : {n},
    totalSets   : {len(tree_ui)},
    totalChars  : {sum(len(v) for v in flat.values())},
    bundled     : true,
  }};
  console.log('[bioinfor-claw] ✅ ' + {n} + ' skills embedded in HTML ({kb} KB) — ready instantly');
  </script>"""

    html = Path(template_path).read_text(encoding='utf-8')

    # Remove any previous embedded block so re-running is idempotent
    html = re.sub(
        r'\n  <!-- ══+.*?embedded-skills.*?</script>',
        '', html, flags=re.DOTALL
    )
    # Remove old external bundle loaders
    html = re.sub(r'\n?<!-- Skills bundle.*?</script>\n?', '', html, flags=re.DOTALL)
    html = re.sub(r'\n?<script>\s*\(function\(\)\s*\{[\s\S]*?tryLoad[\s\S]*?\}\)\(\);[\s\S]*?</script>\n?', '', html, flags=re.DOTALL)

    if '</head>' not in html:
        print("❌  </head> not found in template"); sys.exit(1)

    html = html.replace('</head>', inject + '\n</head>', 1)

    # Also patch the welcome message so it says "Embedded" not "GitHub"
    html = html.replace(
        'Loading skill library from GitHub…',
        f'Loading {n} embedded skills…'
    )
    return html

def main():
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument('--repo',     default='.',                       help='bioinfor-claw folder')
    ap.add_argument('--template', default='bioinfor-claw.html',      help='HTML template file')
    ap.add_argument('--output',   default='bioinfor-claw-ready.html',help='Output file name')
    args = ap.parse_args()

    repo     = Path(args.repo).expanduser().resolve()
    template = (repo / args.template) if not Path(args.template).is_absolute() else Path(args.template)
    out      = repo / args.output

    if not repo.is_dir():
        print(f"❌  Folder not found: {repo}"); sys.exit(1)
    if not template.exists():
        # Try Downloads
        dl = Path.home() / 'Downloads' / 'bioinfor-claw.html'
        if dl.exists():
            template = dl
            print(f"ℹ️   Using template from Downloads: {template}")
        else:
            print(f"❌  Template not found: {template}")
            print(f"    Copy bioinfor-claw.html into {repo} and re-run.")
            sys.exit(1)

    tree_ui, flat = scan(repo)
    if not flat:
        print("❌  No SKILL.md files found. Check the --repo path."); sys.exit(1)

    html = build_html(tree_ui, flat, template)
    out.write_text(html, encoding='utf-8')

    size_kb = out.stat().st_size // 1024
    n = len(flat)
    print(f"""
{'═'*58}
  ✅  Done!

  Output file:  {out.name}
  Size:         {size_kb} KB  ({n} skills embedded)

  ┌─────────────────────────────────────────────────┐
  │  Open this file in your browser:                │
  │                                                 │
  │  open {out.name:<40} │
  │                                                 │
  │  All {n} skills load instantly — no server,      │
  │  no GitHub, no internet required.               │
  │                                                 │
  │  Re-run after adding new skills:                │
  │  python3 build_and_embed.py                     │
  └─────────────────────────────────────────────────┘
{'═'*58}
""")

if __name__ == '__main__':
    main()
