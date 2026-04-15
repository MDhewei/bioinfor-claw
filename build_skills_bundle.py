#!/usr/bin/env python3
"""
bioinfor-claw Skill Bundler
============================
Run this script once from your local bioinfor-claw folder.
It reads every SKILL.md and produces skills-bundle.js.

Usage (from your terminal):
    cd /Users/whe3/Documents/bioinfor-claw
    python build_skills_bundle.py

    # Or specify paths manually:
    python build_skills_bundle.py --repo /Users/whe3/Documents/bioinfor-claw --out /Users/whe3/Desktop

    # Preview without writing:
    python build_skills_bundle.py --dry-run
"""

import os, sys, json, argparse
from pathlib import Path
from datetime import datetime

# ── Display metadata for known skill sets ────────────────────────────────────
SKILL_SET_META = {
    "public-datasets-access-and-download":  {"emoji": "🗄️",  "desc": "DepMap, GTEx, TCGA, GEO download & cache"},
    "multiomics-data-analysis":             {"emoji": "🧫",  "desc": "RNA-seq, ATAC-seq, CNV, multi-omics"},
    "crispr-design-and-analysis":           {"emoji": "✂️",  "desc": "sgRNA design, screen QC & hit analysis"},
    "gene-list-analysis":                   {"emoji": "📋",  "desc": "GO, pathway enrichment, annotation"},
    "gene-centered-analysis":               {"emoji": "🔬",  "desc": "Expression, mutation, CNV, dependency"},
    "protein-structure-analysis":           {"emoji": "🧩",  "desc": "Domain mapping, mutation annotation, PDB"},
    "machine-learning-and-deep-learning":   {"emoji": "🤖",  "desc": "Classification, clustering, SHAP"},
    "bioinformatics-plot-generator":        {"emoji": "📊",  "desc": "Publication-quality figures & charts"},
    "paper-search-and-digest":              {"emoji": "📄",  "desc": "PubMed retrieval & summarization"},
    "lab-search-and-track":                 {"emoji": "🏛️", "desc": "Find & monitor labs by field"},
}

# Folders to always skip
SKIP_DIRS = {
    'Assets', 'assets', '.git', '.github', '__pycache__',
    'node_modules', 'docs', 'tests', 'examples', '.DS_Store',
}

def scan_repo(repo_root: Path):
    skill_tree = {}
    found = missing = 0
    print(f"\n🔍  Scanning: {repo_root}\n{'─'*60}")

    for set_dir in sorted(repo_root.iterdir()):
        if not set_dir.is_dir() or set_dir.name in SKIP_DIRS or set_dir.name.startswith('.'):
            continue

        skills = {}
        for skill_dir in sorted(set_dir.iterdir()):
            if not skill_dir.is_dir() or skill_dir.name.startswith('.'):
                continue
            md_path = skill_dir / 'SKILL.md'
            if md_path.exists():
                content = md_path.read_text(encoding='utf-8', errors='replace')
                skills[skill_dir.name] = content
                found += 1
                print(f"  ✅  {set_dir.name}/{skill_dir.name}  ({len(content):,} chars)")
            else:
                missing += 1
                print(f"  ⬜  {set_dir.name}/{skill_dir.name}  (no SKILL.md yet)")

        if skills:
            meta = SKILL_SET_META.get(set_dir.name, {
                "emoji": "⚙️",
                "desc": set_dir.name.replace('-', ' ').title()
            })
            skill_tree[set_dir.name] = {
                "emoji":  meta["emoji"],
                "desc":   meta["desc"],
                "skills": skills,
            }

    print(f"\n{'─'*60}")
    print(f"  📦  {found} skills with SKILL.md bundled")
    if missing:
        print(f"  ⬜  {missing} skill folders without SKILL.md (skipped)")
    return skill_tree, found

def write_bundle(skill_tree, out_path: Path, dry_run=False):
    # Sidebar structure (names only, no content)
    tree_ui = {
        sk: {"emoji": sv["emoji"], "desc": sv["desc"], "skills": list(sv["skills"].keys())}
        for sk, sv in skill_tree.items()
    }
    # Flat map "set/skill" -> content
    flat = {
        f"{sk}/{sn}": sc
        for sk, sv in skill_tree.items()
        for sn, sc in sv["skills"].items()
    }

    n_skills = len(flat)
    n_sets   = len(skill_tree)
    n_chars  = sum(len(v) for v in flat.values())
    ts = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

    js = f"""// ================================================================
// bioinfor-claw  •  Skills Bundle  (AUTO-GENERATED)
// Generated : {ts}
// Skills    : {n_skills}  across {n_sets} skill sets
// Size      : {n_chars/1024:.1f} KB
//
// DO NOT EDIT MANUALLY — re-run build_skills_bundle.py to update.
// ================================================================
(function(){{

// Sidebar UI structure ──────────────────────────────────────────
window.BUNDLED_SKILL_TREE = {json.dumps(tree_ui, indent=2, ensure_ascii=False)};

// Full SKILL.md contents  "skillSet/skillName" → markdown text ─
window.BUNDLED_LOADED_SKILLS = {json.dumps(flat, ensure_ascii=False)};

// Bundle metadata ───────────────────────────────────────────────
window.SKILLS_BUNDLE_META = {{
  generated   : "{ts}",
  totalSkills : {n_skills},
  totalSets   : {n_sets},
  totalChars  : {n_chars},
  bundled     : true,
}};

console.log(`[bioinfor-claw] ✅ ${{window.SKILLS_BUNDLE_META.totalSkills}} skills loaded from bundle (${{(window.SKILLS_BUNDLE_META.totalChars/1024).toFixed(0)}} KB)`);
}})();
"""
    if dry_run:
        print(f"\n  [dry-run] Would write {len(js):,} chars → {out_path}")
    else:
        out_path.parent.mkdir(parents=True, exist_ok=True)
        out_path.write_text(js, encoding='utf-8')
        kb = out_path.stat().st_size / 1024
        print(f"\n✅  skills-bundle.js written → {out_path}  ({kb:.1f} KB)")
    return n_skills

def main():
    ap = argparse.ArgumentParser(description="Bundle bioinfor-claw SKILL.md files into skills-bundle.js")
    ap.add_argument("--repo",    default=".", help="Path to your bioinfor-claw folder (default: current dir)")
    ap.add_argument("--out",     default=None, help="Output directory for skills-bundle.js (default: same as --repo)")
    ap.add_argument("--dry-run", action="store_true", help="Preview without writing files")
    args = ap.parse_args()

    repo_root = Path(args.repo).expanduser().resolve()
    out_dir   = Path(args.out).expanduser().resolve() if args.out else repo_root
    out_path  = out_dir / "skills-bundle.js"

    if not repo_root.is_dir():
        print(f"❌  Folder not found: {repo_root}")
        sys.exit(1)

    skill_tree, n = scan_repo(repo_root)
    if n == 0:
        print("❌  No SKILL.md files found. Check the --repo path.")
        sys.exit(1)

    write_bundle(skill_tree, out_path, dry_run=args.dry_run)

    print(f"""
┌─────────────────────────────────────────────────────┐
│  Next steps:                                        │
│                                                     │
│  1. Place these two files in the same folder:       │
│       • bioinfor-claw.html                         │
│       • skills-bundle.js                           │
│                                                     │
│  2. Open bioinfor-claw.html in your browser.        │
│     All {n:3d} skills load instantly — no GitHub.   │
│                                                     │
│  3. Re-run this script any time you add/update      │
│     a skill to refresh the bundle.                  │
└─────────────────────────────────────────────────────┘
""")

if __name__ == "__main__":
    main()
