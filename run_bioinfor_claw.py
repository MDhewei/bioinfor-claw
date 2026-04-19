#!/usr/bin/env python3
"""
bioinfor-claw  —  One-Command Launcher
=======================================
Run this from your bioinfor-claw folder:

    cd /Users/whe3/Documents/bioinfor-claw
    python3 run_bioinfor_claw.py

It will:
  1. Scan all your SKILL.md files
  2. Build the complete app in memory
  3. Start a local web server
  4. Open your browser automatically

No files to copy. No paths to get right. Just run and go.
Press Ctrl+C to stop.
"""

import os, sys, json, re, webbrowser, threading, time
from pathlib import Path
from datetime import datetime
from http.server import HTTPServer, BaseHTTPRequestHandler

# ── Skill set metadata ──────────────────────────────────────────────────────
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
    "lab-search-and-track":               {"emoji": "🏛️", "desc": "Find & monitor labs by field"},
}
SKIP = {".git","__pycache__","Assets","assets","node_modules",
        "web","web_results",".github","docs","tests","examples",".DS_Store"}

# ── Scan SKILL.md files ─────────────────────────────────────────────────────
def scan(repo):
    tree_ui, flat = {}, {}
    print(f"\n🔍  Scanning {repo}\n" + "─"*50)
    for sd in sorted(Path(repo).iterdir()):
        if not sd.is_dir() or sd.name in SKIP or sd.name.startswith("."): continue
        found = {}
        for kd in sorted(sd.iterdir()):
            if not kd.is_dir() or kd.name.startswith("."): continue
            md = kd / "SKILL.md"
            if md.exists():
                txt = md.read_text(encoding="utf-8", errors="replace")
                found[kd.name] = txt
                flat[f"{sd.name}/{kd.name}"] = txt
                print(f"  ✅  {sd.name}/{kd.name}")
            else:
                print(f"  ⬜  {sd.name}/{kd.name}  (no SKILL.md)")
        if found:
            m = SKILL_SET_META.get(sd.name, {"emoji":"⚙️","desc":sd.name.replace("-"," ").title()})
            tree_ui[sd.name] = {"emoji":m["emoji"],"desc":m["desc"],"skills":list(found.keys())}
    print(f"\n  📦  {len(flat)} skills found\n")
    return tree_ui, flat

# ── Build the complete HTML in memory ──────────────────────────────────────
def build_html(tree_ui, flat):
    n  = len(flat)
    kb = sum(len(v) for v in flat.values()) // 1024
    ts = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

    skills_js = f"""
<script id="embedded-skills">
// {n} skills embedded at {ts}
window.BUNDLED_SKILL_TREE    = {json.dumps(tree_ui, ensure_ascii=False)};
window.BUNDLED_LOADED_SKILLS  = {json.dumps(flat,    ensure_ascii=False)};
window.SKILLS_BUNDLE_META    = {{
  generated  : "{ts}",
  totalSkills: {n},
  totalSets  : {len(tree_ui)},
  totalChars : {sum(len(v) for v in flat.values())},
  bundled    : true,
}};
console.log("[bioinfor-claw] ✅ " + {n} + " skills embedded ({kb} KB) — ready instantly");
</script>"""

    # Minimal but complete standalone HTML app
    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>bioinfor-claw</title>
<link href="https://fonts.googleapis.com/css2?family=DM+Sans:wght@300;400;500;600;700&family=DM+Mono:wght@400;500&display=swap" rel="stylesheet">
{skills_js}
<style>
*, *::before, *::after {{ box-sizing: border-box; margin: 0; padding: 0; }}
:root {{
  --bg:#f6f9f6; --surf:#fff; --surf2:#eef2ee; --surf3:#e6ece6;
  --border:#d8e4d8; --border2:#c4d4c4;
  --g900:#1b3a1f; --g700:#2e6b35; --g500:#4caf50; --g300:#a5d6a7;
  --g100:#e8f5e9; --g50:#f1faf1;
  --teal:#00796b; --teal-l:#e0f2f1;
  --amber:#e65100; --amber-l:#fff3e0;
  --blue:#1565c0; --blue-l:#e3f2fd;
  --red:#b71c1c; --red-l:#ffebee;
  --text:#1a2e1c; --text2:#3d5c40; --muted:#6b8c6e;
  --sans:'DM Sans',sans-serif; --mono:'DM Mono',monospace;
  --r:10px; --sh:0 1px 4px rgba(0,0,0,.07),0 0 0 1px rgba(0,0,0,.04);
  --sh-md:0 4px 16px rgba(0,0,0,.09);
}}
html,body{{height:100%;overflow:hidden;background:var(--bg);}}
body{{font-family:var(--sans);font-size:13px;color:var(--text);display:flex;flex-direction:column;}}
::-webkit-scrollbar{{width:5px;height:5px;}}
::-webkit-scrollbar-thumb{{background:var(--border2);border-radius:3px;}}

/* HEADER */
header{{display:flex;align-items:center;justify-content:space-between;padding:0 18px;height:52px;background:var(--surf);border-bottom:1px solid var(--border);box-shadow:var(--sh);flex-shrink:0;z-index:20;gap:12px;}}
.logo-area{{display:flex;align-items:center;gap:10px;flex-shrink:0;}}
.logo-img{{height:30px;object-fit:contain;}}
.logo-fb{{font-family:var(--mono);font-weight:500;font-size:15px;color:var(--g700);display:none;}}
.llm-selector-wrap{{display:flex;align-items:center;gap:8px;flex:1;justify-content:center;}}
.llm-selector-label{{font-size:11px;color:var(--muted);white-space:nowrap;}}
.llm-tabs{{display:flex;background:var(--surf2);border:1px solid var(--border);border-radius:8px;padding:3px;gap:2px;}}
.llm-tab{{padding:5px 11px;border-radius:6px;border:none;background:transparent;font-family:var(--sans);font-size:11.5px;font-weight:500;color:var(--muted);cursor:pointer;transition:all .15s;white-space:nowrap;display:flex;align-items:center;gap:5px;}}
.llm-tab:hover{{color:var(--text);background:var(--surf3);}}
.llm-tab.active{{background:var(--surf);color:var(--g700);box-shadow:var(--sh);font-weight:600;}}
.llm-tab .provider-dot{{width:7px;height:7px;border-radius:50%;flex-shrink:0;}}
.header-right{{display:flex;align-items:center;gap:8px;flex-shrink:0;}}
.hbtn{{padding:5px 12px;border-radius:7px;border:1px solid var(--border);background:var(--surf2);font-family:var(--sans);font-size:11.5px;font-weight:500;color:var(--text2);cursor:pointer;transition:all .15s;white-space:nowrap;}}
.hbtn:hover{{border-color:var(--g500);color:var(--g700);background:var(--g100);}}
.hbadge{{padding:4px 10px;border-radius:6px;background:var(--surf2);border:1px solid var(--border);font-size:11px;color:var(--text2);font-family:var(--mono);font-weight:500;}}
.hbadge.ok{{background:var(--g100);border-color:var(--g300);color:var(--g700);}}
.status-dot{{width:8px;height:8px;border-radius:50%;background:var(--g500);box-shadow:0 0 0 2px var(--g300);animation:pulse 2s infinite;}}
@keyframes pulse{{0%,100%{{opacity:1}}50%{{opacity:.4}}}}

/* LAYOUT */
.workspace{{flex:1;display:grid;grid-template-columns:240px 1fr 260px;overflow:hidden;}}

/* SIDEBAR */
.sidebar{{border-right:1px solid var(--border);background:var(--surf);display:flex;flex-direction:column;overflow:hidden;}}
.mode-toggle-bar{{padding:10px 12px;border-bottom:1px solid var(--border);flex-shrink:0;}}
.mode-label{{font-size:10px;font-weight:600;letter-spacing:1px;text-transform:uppercase;color:var(--muted);margin-bottom:7px;}}
.mode-toggle{{display:flex;background:var(--surf2);border:1px solid var(--border);border-radius:8px;padding:3px;gap:2px;}}
.mode-btn{{flex:1;padding:7px 6px;border-radius:6px;border:none;background:transparent;font-family:var(--sans);font-size:11.5px;font-weight:500;color:var(--muted);cursor:pointer;transition:all .15s;text-align:center;}}
.mode-btn:hover{{color:var(--text);background:var(--surf3);}}
.mode-btn.active{{background:var(--g700);color:#fff;font-weight:600;}}
#autoPanel{{flex:1;overflow-y:auto;padding:12px;display:flex;flex-direction:column;gap:14px;}}
.auto-info-card{{background:var(--g50);border:1px solid var(--g300);border-radius:10px;padding:11px 12px;}}
.aic-title{{font-size:12px;font-weight:700;color:var(--g700);display:flex;align-items:center;gap:6px;margin-bottom:5px;}}
.aic-desc{{font-size:11px;color:var(--text2);line-height:1.65;}}
.route-live{{margin-top:9px;padding:7px 10px;background:var(--surf);border:1px solid var(--border);border-radius:7px;font-size:11px;display:flex;align-items:center;gap:7px;color:var(--muted);}}
.rl-dot{{width:7px;height:7px;border-radius:50%;background:var(--border2);flex-shrink:0;transition:background .2s;}}
.rl-dot.routing{{background:var(--amber);animation:pulse .7s infinite;}}
.rl-dot.done{{background:var(--g500);}}
#manualPanel{{flex:1;overflow-y:auto;padding:10px 11px;display:none;flex-direction:column;gap:14px;}}
.sec-lbl{{font-size:10px;font-weight:700;letter-spacing:1.2px;text-transform:uppercase;color:var(--muted);margin-bottom:7px;}}
.skill-set{{margin-bottom:4px;}}
.ss-header{{display:flex;align-items:center;gap:8px;padding:8px 10px;border:1px solid var(--border);border-radius:8px;cursor:pointer;background:transparent;width:100%;text-align:left;transition:all .15s;}}
.ss-header:hover{{background:var(--g50);border-color:var(--g300);}}
.ss-header.active-set{{background:var(--g100);border-color:var(--g700);}}
.ss-icon{{font-size:15px;flex-shrink:0;}}
.ss-info{{flex:1;min-width:0;}}
.ss-name{{font-size:11.5px;font-weight:600;color:var(--text);line-height:1.2;}}
.ss-desc{{font-size:10px;color:var(--muted);line-height:1.3;white-space:nowrap;overflow:hidden;text-overflow:ellipsis;}}
.ss-header.active-set .ss-name{{color:var(--g700);}}
.ss-arrow{{font-size:9px;color:var(--muted);transition:transform .2s;flex-shrink:0;}}
.ss-header.open .ss-arrow{{transform:rotate(90deg);}}
.sub-list{{display:none;margin:3px 0 2px 10px;border-left:2px solid var(--g300);padding-left:8px;}}
.sub-list.open{{display:block;}}
.sub-item{{display:flex;align-items:center;gap:7px;width:100%;padding:6px 9px;background:none;border:1px solid transparent;border-radius:6px;cursor:pointer;text-align:left;font-size:11.5px;color:var(--text2);font-family:var(--sans);transition:all .15s;margin-bottom:2px;}}
.sub-item:hover{{background:var(--g50);color:var(--g700);border-color:var(--g300);}}
.sub-item.selected{{background:var(--g100);border-color:var(--g700);color:var(--g700);font-weight:600;}}
.sub-dot{{width:6px;height:6px;border-radius:50%;background:var(--border2);flex-shrink:0;transition:background .3s;}}
.sub-dot.ok{{background:var(--g500);}}
.sub-dot.err{{background:var(--red);}}
.ex-item{{padding:7px 10px;border-left:2px solid var(--border2);font-size:11px;color:var(--text2);cursor:pointer;margin-bottom:5px;transition:all .15s;border-radius:0 5px 5px 0;line-height:1.5;}}
.ex-item:hover{{border-left-color:var(--teal);color:var(--teal);background:var(--teal-l);}}

/* CHAT */
.chat-area{{display:flex;flex-direction:column;overflow:hidden;background:var(--bg);position:relative;}}
.drop-overlay{{position:absolute;inset:0;z-index:50;background:rgba(46,107,53,.08);border:2px dashed var(--g500);border-radius:12px;display:none;align-items:center;justify-content:center;flex-direction:column;gap:10px;pointer-events:none;}}
.drop-overlay.active{{display:flex;}}
.drop-overlay-icon{{font-size:40px;}}
.drop-overlay-text{{font-size:15px;font-weight:600;color:var(--g700);}}
.drop-overlay-hint{{font-size:12px;color:var(--muted);}}
#messages{{flex:1;overflow-y:auto;padding:28px 30px;display:flex;flex-direction:column;gap:20px;scroll-behavior:smooth;}}
.welcome{{margin:auto;display:flex;flex-direction:column;align-items:center;text-align:center;gap:14px;padding:44px 36px;animation:up .5s ease both;}}
@keyframes up{{from{{opacity:0;transform:translateY(14px)}}to{{opacity:1;transform:none}}}}
.welcome-logo{{height:50px;object-fit:contain;filter:drop-shadow(0 2px 10px rgba(46,107,53,.2));}}
.welcome h2{{font-size:21px;font-weight:700;color:var(--g700);letter-spacing:-.4px;}}
.welcome p{{color:var(--text2);font-size:13px;line-height:1.8;max-width:400px;}}
.welcome-note{{display:flex;align-items:center;gap:8px;padding:8px 14px;background:#fff8e1;border:1px solid #ffe082;border-radius:8px;font-size:11.5px;color:#7a5200;margin-top:2px;}}
.load-bar{{width:280px;height:4px;background:var(--border);border-radius:2px;overflow:hidden;}}
.load-fill{{height:100%;background:var(--g500);border-radius:2px;transition:width .3s;width:0%;}}
.load-txt{{font-size:11px;color:var(--muted);}}
.wpills{{display:none;flex-wrap:wrap;gap:6px;justify-content:center;max-width:440px;margin-top:2px;}}
.wp{{padding:4px 10px;border:1px solid var(--border2);border-radius:20px;font-size:11px;color:var(--text2);background:var(--surf);font-weight:500;}}
.msg{{display:flex;gap:11px;animation:up .2s ease both;}}
.msg.user{{flex-direction:row-reverse;align-self:flex-end;max-width:70%;}}
.msg.assistant{{align-self:flex-start;max-width:86%;}}
.avatar{{width:32px;height:32px;border-radius:9px;display:flex;align-items:center;justify-content:center;font-size:15px;flex-shrink:0;border:1px solid var(--border);overflow:hidden;}}
.msg.user .avatar{{background:var(--g100);border-color:var(--g300);}}
.msg.assistant .avatar{{background:var(--surf);}}
.av-img{{width:100%;height:100%;object-fit:contain;padding:4px;}}
.bubble{{padding:12px 15px;border-radius:11px;line-height:1.75;font-size:13px;box-shadow:var(--sh);}}
.msg.user .bubble{{background:var(--g700);color:#fff;border-radius:11px 3px 11px 11px;}}
.msg.assistant .bubble{{background:var(--surf);border:1px solid var(--border);color:var(--text);border-radius:3px 11px 11px 11px;}}
.bubble p{{margin-bottom:8px;}}.bubble p:last-child{{margin-bottom:0;}}
.msg.assistant .bubble strong{{color:var(--g700);}}
.msg.user .bubble strong{{color:#c8e6c9;}}
.bubble em{{color:var(--teal);}}
.bubble h3{{font-size:14px;font-weight:700;color:var(--g700);margin:12px 0 6px;}}
.bubble ul,.bubble ol{{margin:6px 0 6px 18px;}}.bubble li{{margin-bottom:3px;}}
.msg.assistant .bubble code{{background:var(--surf2);padding:2px 5px;border-radius:4px;font-size:11.5px;color:var(--teal);font-family:var(--mono);border:1px solid var(--border);}}
.msg.user .bubble code{{background:rgba(255,255,255,.2);color:#fff;padding:1px 4px;border-radius:3px;}}
.bubble pre code{{background:none;padding:0;border:none;color:inherit;}}
.route-badge{{display:inline-flex;align-items:center;gap:6px;padding:4px 10px;background:var(--g50);border:1px solid var(--g300);border-radius:20px;font-size:10.5px;color:var(--g700);font-weight:600;margin-bottom:10px;}}
.rb-mode{{font-size:9.5px;color:var(--muted);padding:2px 6px;background:var(--surf2);border:1px solid var(--border);border-radius:10px;margin-left:4px;}}
.routing-row{{display:flex;align-items:center;gap:9px;padding:9px 14px;background:var(--amber-l);border:1px solid #ffcc80;border-radius:9px;font-size:11.5px;color:var(--amber);font-style:italic;animation:up .2s ease both;}}
.spin{{width:14px;height:14px;border:2px solid #ffcc80;border-top-color:var(--amber);border-radius:50%;animation:spin .7s linear infinite;flex-shrink:0;}}
@keyframes spin{{to{{transform:rotate(360deg)}}}}
.code-block{{margin:10px 0;border-radius:9px;overflow:hidden;border:1px solid var(--border);box-shadow:var(--sh);}}
.cb-header{{display:flex;align-items:center;justify-content:space-between;padding:7px 12px;background:var(--surf2);border-bottom:1px solid var(--border);}}
.cb-left{{display:flex;align-items:center;gap:8px;}}
.lang-pill{{padding:2px 8px;border-radius:4px;font-size:10px;font-weight:700;letter-spacing:.4px;text-transform:uppercase;font-family:var(--mono);}}
.lp-py{{background:#dbeafe;color:#1d4ed8;}}.lp-sh{{background:var(--g100);color:var(--g700);}}.lp-r{{background:#fce7f3;color:#9d174d;}}
.cb-skill{{font-size:10px;color:var(--muted);font-family:var(--mono);}}
.cb-copy{{font-size:11px;color:var(--muted);background:none;border:1px solid var(--border);cursor:pointer;padding:2px 8px;border-radius:5px;font-weight:500;transition:all .15s;}}
.cb-copy:hover{{color:var(--g700);border-color:var(--g500);background:var(--g50);}}
.cb-body{{padding:14px 16px;background:#f9fdf9;font-size:12px;color:#1a3a1a;font-family:var(--mono);line-height:1.7;white-space:pre-wrap;word-break:break-word;overflow-x:auto;max-height:400px;overflow-y:auto;}}

/* FILE UPLOAD */
.attached-files{{display:flex;flex-wrap:wrap;gap:6px;padding:0 0 8px 0;}}
.file-chip{{display:flex;align-items:center;gap:6px;padding:4px 9px 4px 8px;background:var(--g50);border:1px solid var(--g300);border-radius:20px;font-size:11px;color:var(--g700);max-width:200px;animation:up .15s ease both;}}
.file-chip-name{{font-weight:500;overflow:hidden;text-overflow:ellipsis;white-space:nowrap;flex:1;min-width:0;}}
.file-chip-size{{font-size:10px;color:var(--muted);flex-shrink:0;}}
.file-chip-remove{{width:16px;height:16px;border-radius:50%;background:var(--border2);border:none;cursor:pointer;display:flex;align-items:center;justify-content:center;font-size:10px;color:var(--muted);flex-shrink:0;transition:all .15s;margin-left:2px;}}
.file-chip-remove:hover{{background:var(--red-l);color:var(--red);}}

/* SETTINGS MODAL */
.modal-overlay{{position:fixed;inset:0;background:rgba(0,0,0,.35);z-index:100;display:none;align-items:center;justify-content:center;backdrop-filter:blur(3px);}}
.modal-overlay.open{{display:flex;}}
.modal{{background:var(--surf);border-radius:14px;box-shadow:var(--sh-md);width:520px;max-width:95vw;overflow:hidden;animation:up .2s ease both;}}
.modal-header{{display:flex;align-items:center;justify-content:space-between;padding:16px 20px;border-bottom:1px solid var(--border);}}
.modal-title{{font-size:15px;font-weight:700;color:var(--text);}}
.modal-close{{width:28px;height:28px;border:1px solid var(--border);border-radius:7px;background:var(--surf2);cursor:pointer;font-size:14px;display:flex;align-items:center;justify-content:center;color:var(--muted);transition:all .15s;}}
.modal-close:hover{{background:var(--red-l);border-color:var(--red);color:var(--red);}}
.modal-body{{padding:20px;display:flex;flex-direction:column;gap:18px;max-height:70vh;overflow-y:auto;}}
.provider-grid{{display:grid;grid-template-columns:repeat(3,1fr);gap:8px;}}
.provider-card{{border:1.5px solid var(--border);border-radius:9px;padding:11px 10px;cursor:pointer;transition:all .15s;text-align:center;background:var(--surf);}}
.provider-card:hover{{border-color:var(--g500);background:var(--g50);}}
.provider-card.selected{{border-color:var(--g700);background:var(--g100);box-shadow:0 0 0 2px var(--g300);}}
.pc-icon{{font-size:22px;margin-bottom:5px;}}.pc-name{{font-size:11.5px;font-weight:700;color:var(--text);}}.pc-models{{font-size:10px;color:var(--muted);margin-top:2px;}}
.field-group{{display:flex;flex-direction:column;gap:6px;}}
.field-label{{font-size:11.5px;font-weight:600;color:var(--text2);}}
.field-input{{width:100%;padding:9px 11px;border:1px solid var(--border2);border-radius:8px;font-family:var(--sans);font-size:13px;color:var(--text);background:var(--surf);outline:none;transition:border .15s;}}
.field-input:focus{{border-color:var(--g500);box-shadow:0 0 0 2px var(--g100);}}
select.field-input{{cursor:pointer;}}
.field-hint{{font-size:10.5px;color:var(--muted);line-height:1.5;}}
.key-row{{display:flex;gap:7px;}}.key-row .field-input{{flex:1;}}
.save-btn{{padding:9px 18px;background:var(--g700);color:#fff;border:none;border-radius:8px;font-family:var(--sans);font-size:13px;font-weight:600;cursor:pointer;transition:background .15s;align-self:flex-end;}}
.save-btn:hover{{background:var(--g900);}}
.key-status{{display:flex;align-items:center;gap:6px;padding:8px 11px;border-radius:7px;font-size:11.5px;border:1px solid var(--border);background:var(--surf2);color:var(--muted);}}
.ks-dot{{width:7px;height:7px;border-radius:50%;background:var(--border2);flex-shrink:0;}}.ks-dot.set{{background:var(--g500);}}

/* INPUT */
.input-zone{{padding:12px 26px 16px;border-top:1px solid var(--border);background:var(--surf);box-shadow:0 -1px 5px rgba(0,0,0,.04);flex-shrink:0;}}
.selected-skill-bar{{display:flex;align-items:center;gap:7px;padding:6px 10px;background:var(--g50);border:1px solid var(--g300);border-radius:8px;font-size:11.5px;color:var(--g700);margin-bottom:9px;font-weight:500;}}
.ssb-change{{margin-left:auto;font-size:10.5px;color:var(--muted);cursor:pointer;padding:2px 7px;border:1px solid var(--border);border-radius:5px;background:var(--surf);transition:all .15s;}}
.ssb-change:hover{{color:var(--g700);border-color:var(--g500);background:var(--g50);}}
.input-row{{display:flex;align-items:flex-end;gap:8px;background:var(--bg);border:1.5px solid var(--border2);border-radius:12px;padding:10px 12px;transition:border .2s,box-shadow .2s;}}
.input-row:focus-within{{border-color:var(--g500);box-shadow:0 0 0 3px var(--g100);}}
#userInput{{flex:1;background:none;border:none;outline:none;color:var(--text);font-family:var(--sans);font-size:13px;resize:none;min-height:20px;max-height:130px;line-height:1.6;}}
#userInput::placeholder{{color:var(--muted);}}
.attach-btn{{width:32px;height:32px;background:none;border:1px solid var(--border2);border-radius:8px;cursor:pointer;display:flex;align-items:center;justify-content:center;font-size:15px;flex-shrink:0;transition:all .15s;color:var(--muted);}}
.attach-btn:hover{{border-color:var(--g500);background:var(--g50);color:var(--g700);}}
.send-btn{{width:36px;height:36px;background:var(--g700);border:none;border-radius:9px;cursor:pointer;display:flex;align-items:center;justify-content:center;color:#fff;font-size:16px;flex-shrink:0;transition:all .15s;box-shadow:0 2px 6px rgba(46,107,53,.3);}}
.send-btn:hover{{background:var(--g900);transform:translateY(-1px);}}
.send-btn:disabled{{opacity:.4;cursor:not-allowed;transform:none;box-shadow:none;}}
.input-meta{{display:flex;align-items:center;justify-content:space-between;font-size:10.5px;color:var(--muted);margin-top:7px;padding:0 2px;}}
.im-right{{display:flex;align-items:center;gap:7px;}}
.mode-pill{{padding:2px 8px;border-radius:10px;font-size:10.5px;font-weight:600;}}
.mode-pill.auto{{background:var(--g100);color:var(--g700);border:1px solid var(--g300);}}
.mode-pill.manual{{background:var(--blue-l);color:var(--blue);border:1px solid #90caf9;}}
.auto-badge{{display:flex;align-items:center;gap:5px;padding:2px 9px;border:1px solid var(--g100);border-radius:10px;font-size:10.5px;color:var(--g700);background:var(--g50);font-weight:600;}}
.model-pill{{padding:2px 8px;border-radius:10px;font-size:10.5px;background:var(--surf2);border:1px solid var(--border);color:var(--text2);}}

/* RIGHT PANEL */
.right-panel{{border-left:1px solid var(--border);background:var(--surf);display:flex;flex-direction:column;overflow:hidden;}}
.panel-tabs{{display:flex;border-bottom:1px solid var(--border);flex-shrink:0;background:var(--surf2);}}
.ptab{{flex:1;padding:11px 6px;font-size:10px;font-weight:600;letter-spacing:.8px;text-transform:uppercase;color:var(--muted);text-align:center;cursor:pointer;border:none;background:none;border-bottom:2px solid transparent;font-family:var(--sans);transition:all .15s;}}
.ptab.active{{color:var(--g700);border-bottom-color:var(--g700);background:var(--surf);}}
.pcontent{{flex:1;overflow-y:auto;padding:14px;}}
.rlog-item{{border:1px solid var(--border);border-radius:9px;margin-bottom:9px;overflow:hidden;box-shadow:var(--sh);}}
.rli-head{{display:flex;align-items:center;justify-content:space-between;padding:7px 11px;background:var(--surf2);border-bottom:1px solid var(--border);font-size:10.5px;}}
.rli-skill{{font-weight:700;color:var(--g700);}}.rli-mode{{padding:2px 6px;border-radius:4px;font-size:9.5px;}}
.rli-auto{{background:var(--g100);color:var(--g700);}}.rli-manual{{background:var(--blue-l);color:var(--blue);}}
.rli-body{{padding:8px 11px;font-size:11px;color:var(--muted);line-height:1.5;}}
.script-card{{border:1px solid var(--border);border-radius:9px;margin-bottom:9px;overflow:hidden;box-shadow:var(--sh);}}
.sc-head{{display:flex;align-items:center;justify-content:space-between;padding:7px 11px;background:var(--surf2);border-bottom:1px solid var(--border);font-size:10.5px;color:var(--text2);}}
.sc-preview{{padding:8px 11px;font-size:11px;color:var(--muted);line-height:1.5;max-height:65px;overflow:hidden;white-space:pre-wrap;word-break:break-word;font-family:var(--mono);}}
.info-row{{display:flex;justify-content:space-between;align-items:center;padding:9px 0;border-bottom:1px solid var(--border);font-size:11.5px;}}
.info-row:last-child{{border-bottom:none;}}.info-label{{color:var(--muted);}}.info-value{{color:var(--g700);font-weight:600;font-family:var(--mono);}}
.clr-btn{{width:100%;padding:9px;background:transparent;border:1px solid #ef9a9a;color:var(--red);border-radius:8px;font-size:12px;font-weight:500;cursor:pointer;margin-top:14px;transition:all .15s;}}
.clr-btn:hover{{background:var(--red-l);}}
.empty-st{{text-align:center;padding:30px 14px;color:var(--muted);font-size:12px;line-height:1.8;}}
.empty-ic{{font-size:26px;margin-bottom:8px;opacity:.45;}}
</style>
</head>
<body>

<header>
  <div class="logo-area">
    <img class="logo-img" src="https://raw.githubusercontent.com/MDhewei/bioinfor-claw/main/Assets/long_logo.png"
      alt="bioinfor-claw" onerror="this.style.display='none';document.getElementById('lf').style.display='block'"/>
    <span id="lf" class="logo-fb">bioinfor-claw</span>
  </div>
  <div class="llm-selector-wrap">
    <span class="llm-selector-label">LLM:</span>
    <div class="llm-tabs" id="llmTabs">
      <button class="llm-tab active" onclick="selectProvider('anthropic',this)" data-p="anthropic">
        <span class="provider-dot" style="background:#d97706"></span>Anthropic
      </button>
      <button class="llm-tab" onclick="selectProvider('openai',this)" data-p="openai">
        <span class="provider-dot" style="background:#10a37f"></span>OpenAI
      </button>
      <button class="llm-tab" onclick="selectProvider('google',this)" data-p="google">
        <span class="provider-dot" style="background:#4285f4"></span>Google
      </button>
      <button class="llm-tab" onclick="selectProvider('mistral',this)" data-p="mistral">
        <span class="provider-dot" style="background:#f97316"></span>Mistral
      </button>
      <button class="llm-tab" onclick="selectProvider('minimax',this)" data-p="minimax">
        <span class="provider-dot" style="background:#e91e8c"></span>MiniMax
      </button>
      <button class="llm-tab" onclick="selectProvider('custom',this)" data-p="custom">
        <span class="provider-dot" style="background:#8b5cf6"></span>Custom
      </button>
    </div>
  </div>
  <div class="header-right">
    <div class="hbadge" id="skillBadge">⏳ Loading…</div>
    <div class="hbadge" id="serverBadge">💬 Chat only</div>
    <div class="status-dot" id="statusDot"></div>
    <button class="hbtn" onclick="openSettings()">⚙ Settings</button>
  </div>
</header>

<!-- Settings Modal -->
<div class="modal-overlay" id="settingsModal">
  <div class="modal">
    <div class="modal-header">
      <span class="modal-title">⚙ LLM Settings</span>
      <button class="modal-close" onclick="closeSettings()">✕</button>
    </div>
    <div class="modal-body">
      <div class="field-group">
        <div class="field-label">Provider</div>
        <div class="provider-grid">
          <div class="provider-card selected" id="pc-anthropic" onclick="selectProviderModal('anthropic')"><div class="pc-icon">🟠</div><div class="pc-name">Anthropic</div><div class="pc-models">Claude 3.5/4</div></div>
          <div class="provider-card" id="pc-openai" onclick="selectProviderModal('openai')"><div class="pc-icon">🟢</div><div class="pc-name">OpenAI</div><div class="pc-models">GPT-4o / o1</div></div>
          <div class="provider-card" id="pc-google" onclick="selectProviderModal('google')"><div class="pc-icon">🔵</div><div class="pc-name">Google</div><div class="pc-models">Gemini 1.5/2</div></div>
          <div class="provider-card" id="pc-mistral" onclick="selectProviderModal('mistral')"><div class="pc-icon">🟡</div><div class="pc-name">Mistral</div><div class="pc-models">Mistral Large</div></div>
          <div class="provider-card" id="pc-minimax" onclick="selectProviderModal('minimax')"><div class="pc-icon">🩷</div><div class="pc-name">MiniMax</div><div class="pc-models">M2.5 · M2 · Text-01</div></div>
          <div class="provider-card" id="pc-custom" onclick="selectProviderModal('custom')"><div class="pc-icon">🟣</div><div class="pc-name">Custom</div><div class="pc-models">OpenAI-compat.</div></div>
        </div>
      </div>
      <div class="field-group">
        <div class="field-label">Model</div>
        <select class="field-input" id="modelSelect" onchange="onModelChange()"></select>
      </div>
      <div class="field-group" id="customEndpointGroup" style="display:none">
        <div class="field-label">API Endpoint</div>
        <input class="field-input" id="customEndpoint" type="url" placeholder="https://your-endpoint.com/v1"/>
        <div class="field-hint">Must be OpenAI-compatible (/v1/chat/completions)</div>
      </div>
      <div class="field-group">
        <div class="field-label" id="keyLabel">Anthropic API Key</div>
        <div class="key-row">
          <input class="field-input" id="apiKeyInput" type="password" placeholder="sk-ant-…"/>
          <button class="save-btn" onclick="saveKey()">Save</button>
        </div>
        <div class="key-status" id="keyStatus">
          <div class="ks-dot" id="ksDot"></div>
          <span id="ksText">No key set</span>
        </div>
        <div class="field-hint" id="keyHint">Stored in sessionStorage only — never sent anywhere except the provider API.</div>
      </div>
    </div>
  </div>
</div>

<div class="workspace">
  <aside class="sidebar">
    <div class="mode-toggle-bar">
      <div class="mode-label">Skill Selection Mode</div>
      <div class="mode-toggle">
        <button class="mode-btn active" id="modeAutoBtn" onclick="setMode('auto')">🤖 Auto-Route</button>
        <button class="mode-btn" id="modeManualBtn" onclick="setMode('manual')">🖱 Manual Pick</button>
      </div>
    </div>
    <div id="autoPanel">
      <div class="auto-info-card">
        <div class="aic-title">🤖 Auto-Routing Active</div>
        <div class="aic-desc">Just describe your analysis — the agent picks the best skill automatically.</div>
        <div class="route-live"><div class="rl-dot" id="rlDot"></div><span id="rlText">Waiting…</span></div>
      </div>
      <div><div class="sec-lbl">Skill Library</div><div id="autoSkillTree"></div></div>
      <div>
        <div class="sec-lbl">Try These</div>
        <div class="ex-item" onclick="useEx(this)">Analyze TP53 across DepMap, GTEx and TCGA</div>
        <div class="ex-item" onclick="useEx(this)">Design sgRNAs for BRCA1 knockout</div>
        <div class="ex-item" onclick="useEx(this)">Run GO enrichment and plot top 20 pathways</div>
        <div class="ex-item" onclick="useEx(this)">Map EGFR mutations onto protein domains</div>
        <div class="ex-item" onclick="useEx(this)">Search recent papers on KRAS inhibitors</div>
        <div class="ex-item" onclick="useEx(this)">Generate a volcano plot from DESeq2 CSV</div>
        <div class="ex-item" onclick="useEx(this)">Train a classifier on DepMap dependency scores</div>
      </div>
    </div>
    <div id="manualPanel">
      <div class="sec-lbl">Browse &amp; Select a Skill</div>
      <div id="manualSkillTree"></div>
    </div>
  </aside>

  <div class="chat-area" id="chatArea">
    <div class="drop-overlay" id="dropOverlay">
      <div class="drop-overlay-icon">📂</div>
      <div class="drop-overlay-text">Drop your files here</div>
      <div class="drop-overlay-hint">CSV, TSV, FASTA, FASTQ, VCF, TXT, JSON, BED</div>
    </div>
    <div id="messages">
      <div class="welcome" id="welcome">
        <img class="welcome-logo" src="https://raw.githubusercontent.com/MDhewei/bioinfor-claw/main/Assets/long_logo.png"
          alt="bioinfor-claw" onerror="this.style.display='none'"/>
        <h2>Your 24/7 bioinformatics copilot</h2>
        <p>Chat to perform your daily bioinformatics analysis with ease</p>
        <div class="welcome-note"><span>🔑</span><span>Click <strong>⚙ Settings</strong> to add your API key (Anthropic, OpenAI, Google, or Mistral)</span></div>
        <div class="load-bar"><div class="load-fill" id="loadFill"></div></div>
        <div class="load-txt" id="loadTxt">Initialising…</div>
        <div class="wpills" id="wPills"></div>
      </div>
    </div>
    <div class="input-zone">
      <div class="selected-skill-bar" id="selectedSkillBar" style="display:none">
        <span id="ssEmoji">📋</span><span id="ssName">No skill selected</span>
        <span class="ssb-change" onclick="setMode('manual')">Change</span>
      </div>
      <div class="attached-files" id="attachedFiles" style="display:none"></div>
      <input type="file" id="fileInput" multiple accept=".csv,.tsv,.txt,.fasta,.fa,.fq,.fastq,.vcf,.json,.bed,.gff,.gtf,.sam" style="display:none" onchange="handleFileSelect(event)"/>
      <div class="input-row">
        <button class="attach-btn" id="attachBtn" onclick="document.getElementById('fileInput').click()" title="Attach files">📎</button>
        <textarea id="userInput" rows="1" placeholder="Describe your analysis — the agent picks the right skill automatically…"></textarea>
        <button class="send-btn" id="sendBtn" onclick="sendMessage()">↑</button>
      </div>
      <div class="input-meta">
        <span>Enter to send · Shift+Enter for new line</span>
        <div class="im-right">
          <span class="mode-pill auto" id="modePill">🤖 Auto</span>
          <span class="model-pill" id="modelPill">claude-sonnet-4</span>
        </div>
      </div>
    </div>
  </div>

  <div class="right-panel">
    <div class="panel-tabs">
      <button class="ptab active" onclick="switchTab('routing',this)">Routing</button>
      <button class="ptab" onclick="switchTab('scripts',this)">Scripts</button>
      <button class="ptab" onclick="switchTab('files',this)">Files</button>
      <button class="ptab" onclick="switchTab('session',this)">Session</button>
    </div>
    <div class="pcontent">
      <div id="tab-routing"><div class="empty-st"><div class="empty-ic">🔀</div>Routing decisions appear here</div></div>
      <div id="tab-scripts" style="display:none"><div class="empty-st"><div class="empty-ic">📝</div>Generated scripts appear here</div></div>
      <div id="tab-files" style="display:none"><div class="empty-st"><div class="empty-ic">📂</div>Uploaded files appear here</div></div>
      <div id="tab-session" style="display:none">
        <div class="info-row"><span class="info-label">Skills loaded</span><span class="info-value" id="infoLoaded">0</span></div>
        <div class="info-row"><span class="info-label">Provider</span><span class="info-value" id="infoProv">Anthropic</span></div>
        <div class="info-row"><span class="info-label">Model</span><span class="info-value" id="infoModel">—</span></div>
        <div class="info-row"><span class="info-label">Mode</span><span class="info-value" id="infoMode">Auto</span></div>
        <div class="info-row"><span class="info-label">Messages</span><span class="info-value" id="infoMsgs">0</span></div>
        <div class="info-row"><span class="info-label">Scripts</span><span class="info-value" id="infoScripts">0</span></div>
        <button class="clr-btn" onclick="clearSession()">🗑 Clear Session</button>
      </div>
    </div>
  </div>
</div>

<script>
// ═══════════════════════════════════════════════════════════
//  SKILL TREE — populated from embedded data at top of page
// ═══════════════════════════════════════════════════════════
const REPO_BASE = 'https://raw.githubusercontent.com/MDhewei/bioinfor-claw/main';
const LOGO_URL  = REPO_BASE + '/Assets/long_logo.png';

// Build SKILL_TREE from embedded bundle
const SKILL_TREE = {{}};
if (window.BUNDLED_SKILL_TREE) {{
  Object.assign(SKILL_TREE, window.BUNDLED_SKILL_TREE);
}}

const ALL_SKILLS = [];
Object.entries(SKILL_TREE).forEach(([sk,m]) =>
  m.skills.forEach(s => ALL_SKILLS.push({{setKey:sk, skillName:s, key:sk+'/'+s, emoji:m.emoji}})));

// ═══ PROVIDERS ═══════════════════════════════════════════════
const PROVIDERS = {{
  anthropic: {{name:'Anthropic',color:'#d97706',keyPlaceholder:'sk-ant-…',keyHint:'console.anthropic.com',
    models:['claude-sonnet-4-20250514','claude-opus-4-6','claude-haiku-4-5-20251001'],
    modelLabels:['Claude Sonnet 4 (recommended)','Claude Opus 4.6','Claude Haiku 4.5'],call:callAnthropic}},
  openai:    {{name:'OpenAI',  color:'#10a37f',keyPlaceholder:'sk-…',keyHint:'platform.openai.com',
    models:['gpt-4o','gpt-4o-mini','o1-mini'],modelLabels:['GPT-4o (recommended)','GPT-4o Mini','o1 Mini'],call:callOpenAI}},
  google:    {{name:'Google',  color:'#4285f4',keyPlaceholder:'AIza…',keyHint:'aistudio.google.com',
    models:['gemini-2.0-flash','gemini-1.5-pro','gemini-1.5-flash'],
    modelLabels:['Gemini 2.0 Flash (recommended)','Gemini 1.5 Pro','Gemini 1.5 Flash'],call:callGoogle}},
  mistral:   {{name:'Mistral', color:'#f97316',keyPlaceholder:'your-key',keyHint:'console.mistral.ai',
    models:['mistral-large-latest','mistral-medium-latest','mistral-small-latest'],
    modelLabels:['Mistral Large','Mistral Medium','Mistral Small'],call:callOpenAI}},
  minimax:   {{name:'MiniMax', color:'#e91e8c',keyPlaceholder:'your-minimax-key',keyHint:'platform.minimax.io → API Keys',
    models:['MiniMax-M2.5','MiniMax-M2.5-highspeed','MiniMax-M2','MiniMax-M2-highspeed','MiniMax-Text-01'],
    modelLabels:['MiniMax M2.5 (recommended)','MiniMax M2.5 High Speed','MiniMax M2','MiniMax M2 High Speed','MiniMax Text-01'],
    endpoint:'https://api.minimax.io/v1', call:callOpenAI}},
  custom:    {{name:'Custom',  color:'#8b5cf6',keyPlaceholder:'your-key',keyHint:'Any OpenAI-compatible endpoint',
    models:['custom-model'],modelLabels:['Custom model'],call:callOpenAI}},
}};

// ═══ STATE ═══════════════════════════════════════════════════
let state = {{
  provider:'anthropic', model:'claude-sonnet-4-20250514',
  keys:{{}}, customEndpoint:'',
  mode:'auto', selectedSkillKey:null,
  conversation:[], loadedSkills:{{}}, systemPrompt:'',
  skillsReady:false, msgCount:0, scriptCount:0,
  attachedFiles:[], serverConnected:false, serverUrl:'http://localhost:8000',
}};

// ═══ INIT ════════════════════════════════════════════════════
window.onload = async () => {{
  const saved = sessionStorage.getItem('bfc_state');
  if (saved) try {{ const s=JSON.parse(saved); state.provider=s.provider||'anthropic'; state.model=s.model||PROVIDERS[state.provider].models[0]; state.keys=s.keys||{{}}; state.customEndpoint=s.customEndpoint||''; }} catch{{}}
  
  buildTrees();
  syncProviderUI();
  initDragDrop();

  const ta = document.getElementById('userInput');
  ta.addEventListener('keydown', e => {{ if(e.key==='Enter'&&!e.shiftKey){{e.preventDefault();sendMessage();}} }});
  ta.addEventListener('input', () => {{ ta.style.height='auto'; ta.style.height=Math.min(ta.scrollHeight,130)+'px'; }});

  document.getElementById('settingsModal').addEventListener('click', e => {{ if(e.target===e.currentTarget)closeSettings(); }});

  // Detect server immediately — retry once after 2s if first attempt fails
  detectServer().then(found => {{
    if (!found) setTimeout(() => detectServer(), 2000);
  }});
  await loadAllSkills();
}};

// ═══ LOAD SKILLS — uses embedded data, no GitHub ═════════════
async function loadAllSkills() {{
  // Skills are embedded directly in this HTML file
  if (window.BUNDLED_LOADED_SKILLS && window.SKILLS_BUNDLE_META) {{
    const meta = window.SKILLS_BUNDLE_META;
    state.loadedSkills = Object.assign({{}}, window.BUNDLED_LOADED_SKILLS);
    
    // Mark all skill dots green
    Object.keys(state.loadedSkills).forEach(key => {{
      const parts = key.split('/');
      const [sk, sn] = [parts[0], parts.slice(1).join('/')];
      ['auto','manual'].forEach(t => {{
        const d = document.getElementById(t+'SkillTree-sd-'+sk+'-'+sn);
        if (d) d.classList.add('ok');
      }});
    }});

    document.getElementById('loadFill').style.width = '100%';
    document.getElementById('loadTxt').textContent = '✅ ' + meta.totalSkills + ' skills ready';
    document.getElementById('infoLoaded').textContent = meta.totalSkills;

    // Show pills
    const pills = document.getElementById('wPills');
    Object.entries(SKILL_TREE).forEach(([sk, sv]) => {{
      const p = document.createElement('span'); p.className='wp';
      p.textContent = sv.emoji+' '+fmt(sk); pills.appendChild(p);
    }});
    pills.style.display = 'flex';

    document.getElementById('skillBadge').textContent = '📦 '+meta.totalSkills+' skills';
    document.getElementById('skillBadge').classList.add('ok');

    buildSystemPrompt();
    state.skillsReady = true;
    return;
  }}

  // Fallback — no embedded data found
  document.getElementById('loadTxt').textContent = '⚠️ No embedded skills found — re-run build_and_embed.py';
  document.getElementById('skillBadge').textContent = '⚠️ No skills';
}}

function buildSystemPrompt() {{
  let docs = ''; let count = 0;
  Object.entries(state.loadedSkills).forEach(([key, cv]) => {{
    if (cv) {{
      // Only include skill name + first 3 lines as index (not full content — that goes in turnSystem)
      const preview = cv.split('\n').slice(0,3).join(' ').replace(/#+/g,'').trim().slice(0,120);
      docs += '\n- ' + key + ': ' + preview;
      count++;
    }}
  }});

  const serverNote = state.serverConnected
    ? 'The backend server is running at ' + state.serverUrl + '. '
      + 'When a skill is selected, the web UI automatically executes the Python script and shows real results. '
      + 'You do NOT need to generate commands or pretend to run anything. '
      + 'Just explain the analysis clearly and concisely.'
    : 'No execution server detected. Provide the exact CLI commands the user needs to run locally.';

  state.systemPrompt =
    'You are bioinfor-claw, an AI assistant for bioinformatics research by MDhewei (MD Anderson Cancer Center). '
    + 'Help users understand and run bioinformatics analyses.\n\n'
    + serverNote + '\n\n'
    + 'CRITICAL RULES:\n'
    + '1. Never write [ROUTED_SKILL], [SKILL:], [INPUT:] or any internal tags in your responses.\n'
    + '2. Never say "I will invoke", "let me run", "I am executing" — you cannot execute code.\n'
    + '3. When a skill is active (injected in system context), explain what it does naturally.\n'
    + '4. Keep replies concise — users see real script output in the chat automatically.\n\n'
    + count + ' skills available:\n' + docs;
}}

// ═══ SERVER DETECTION ════════════════════════════════════════
async function detectServer() {{
  // Try these URLs — use /api/skills which we know works
  const candidates = [
    window.location.origin,
    'http://localhost:7860',
    'http://127.0.0.1:7860',
    'http://localhost:8000',
    'http://127.0.0.1:8000',
  ];
  const seen = new Set();
  const urls = [...new Set(candidates)];

  for (const url of urls) {{
    try {{
      // Use /api/skills — more reliable than /api/health
      const r = await fetch(url+'/api/skills', {{
        signal: AbortSignal.timeout(3000),
        cache: 'no-cache'
      }});
      if (r.ok) {{
        const data = await r.json();
        if (data.skills) {{
          state.serverConnected = true;
          state.serverUrl = url;
          buildSystemPrompt();
          const b = document.getElementById('serverBadge');
          if (b) {{
            b.textContent = '⚡ Server live';
            b.style.cssText += ';background:var(--g100);border-color:var(--g300);color:var(--g700)';
            b.classList.add('ok');
          }}
          console.log('[bioinfor-claw] ✅ Execution server at: ' + url);
          appendMiniLog('✅ Server connected at ' + url);
          console.log('[bioinfor-claw] Skill sets:', Object.keys(data.skills));
          return true;
        }}
      }}
    }} catch(e) {{
      console.log('[bioinfor-claw] tried ' + url + ':', e.message);
    }}
  }}

  state.serverConnected = false;
  const b = document.getElementById('serverBadge');
  if (b) {{ b.textContent = '💬 Chat only'; b.classList.remove('ok'); }}
  console.log('[bioinfor-claw] No server found — chat-only mode');
  return false;
}}

// ═══ PROVIDER / MODEL ════════════════════════════════════════
function selectProvider(p,btn) {{
  state.provider=p; state.model=PROVIDERS[p].models[0];
  document.querySelectorAll('.llm-tab').forEach(t=>t.classList.remove('active')); btn.classList.add('active');
  syncProviderUI(); saveState();
}}
function selectProviderModal(p) {{
  state.provider=p; state.model=PROVIDERS[p].models[0];
  document.querySelectorAll('.provider-card').forEach(c=>c.classList.remove('selected'));
  document.getElementById('pc-'+p).classList.add('selected');
  document.querySelectorAll('.llm-tab').forEach(t=>t.classList.toggle('active',t.dataset.p===p));
  updateModelSelect(); updateKeyStatus();
  document.getElementById('customEndpointGroup').style.display=p==='custom'?'flex':'none';
  saveState();
}}
function onModelChange() {{ state.model=document.getElementById('modelSelect').value; document.getElementById('modelPill').textContent=state.model.split('-').slice(0,3).join('-'); document.getElementById('infoModel').textContent=state.model; saveState(); }}
function syncProviderUI() {{ updateModelSelect(); updateKeyStatus(); document.getElementById('infoProv').textContent=PROVIDERS[state.provider].name; document.querySelectorAll('.provider-card').forEach(c=>c.classList.remove('selected')); document.getElementById('pc-'+state.provider)?.classList.add('selected'); document.querySelectorAll('.llm-tab').forEach(t=>t.classList.toggle('active',t.dataset.p===state.provider)); }}
function updateModelSelect() {{ const sel=document.getElementById('modelSelect'); const p=PROVIDERS[state.provider]; sel.innerHTML=p.models.map((m,i)=>`<option value="${{m}}"${{m===state.model?' selected':''}}>${{p.modelLabels[i]}}</option>`).join(''); document.getElementById('modelPill').textContent=state.model.split('-').slice(0,3).join('-'); document.getElementById('infoModel').textContent=state.model; }}
function updateKeyStatus() {{ const p=PROVIDERS[state.provider]; const has=!!state.keys[state.provider]; document.getElementById('keyLabel').textContent=p.name+' API Key'; document.getElementById('apiKeyInput').placeholder=p.keyPlaceholder; document.getElementById('keyHint').textContent=p.keyHint; document.getElementById('ksDot').className='ks-dot'+(has?' set':''); document.getElementById('ksText').textContent=has?'✓ Key saved for '+p.name:'No key set for '+p.name; }}
function saveKey() {{ const v=document.getElementById('apiKeyInput').value.trim(); const ep=document.getElementById('customEndpoint').value.trim(); if(!v){{alert('Please enter an API key');return;}} state.keys[state.provider]=v; if(ep)state.customEndpoint=ep; saveState(); updateKeyStatus(); document.getElementById('apiKeyInput').value=''; closeSettings(); }}
function openSettings() {{ document.getElementById('settingsModal').classList.add('open'); updateModelSelect(); updateKeyStatus(); document.getElementById('customEndpointGroup').style.display=state.provider==='custom'?'flex':'none'; document.getElementById('customEndpoint').value=state.customEndpoint; }}
function closeSettings() {{ document.getElementById('settingsModal').classList.remove('open'); }}
function saveState() {{ sessionStorage.setItem('bfc_state',JSON.stringify({{provider:state.provider,model:state.model,keys:state.keys,customEndpoint:state.customEndpoint}})); }}

// ═══ MODE ════════════════════════════════════════════════════
function setMode(m) {{
  state.mode=m;
  document.getElementById('modeAutoBtn').classList.toggle('active',m==='auto');
  document.getElementById('modeManualBtn').classList.toggle('active',m==='manual');
  document.getElementById('autoPanel').style.display=m==='auto'?'flex':'none';
  document.getElementById('manualPanel').style.display=m==='manual'?'flex':'none';
  const pill=document.getElementById('modePill'); pill.textContent=m==='auto'?'🤖 Auto':'🖱 Manual'; pill.className='mode-pill '+m;
  document.getElementById('selectedSkillBar').style.display=m==='manual'?'flex':'none';
  document.getElementById('infoMode').textContent=m==='auto'?'Auto-Route':'Manual';
  document.getElementById('userInput').placeholder=m==='auto'?'Describe your analysis — the agent picks the right skill…':'Select a skill, then describe your analysis…';
}}

// ═══ SKILL TREES ═════════════════════════════════════════════
function buildTrees() {{ buildTree('autoSkillTree',false); buildTree('manualSkillTree',true); }}
function buildTree(cid,clickable) {{
  const c=document.getElementById(cid);
  Object.entries(SKILL_TREE).forEach(([sk,meta]) => {{
    const g=document.createElement('div'); g.className='skill-set'; g.id=cid+'-sg-'+sk;
    const h=document.createElement('button'); h.className='ss-header';
    h.innerHTML=`<span class="ss-icon">${{meta.emoji}}</span><div class="ss-info"><div class="ss-name">${{fmt(sk)}}</div><div class="ss-desc">${{meta.desc}}</div></div><span class="ss-arrow">▶</span>`;
    h.onclick=()=>toggleSet(cid,sk,h);
    const sub=document.createElement('div'); sub.className='sub-list'; sub.id=cid+'-sl-'+sk;
    meta.skills.forEach(s=>{{
      const b=document.createElement('button'); b.className='sub-item'; b.id=cid+'-si-'+sk+'-'+s;
      b.innerHTML=`<span class="sub-dot" id="${{cid}}-sd-${{sk}}-${{s}}"></span>${{fmt(s)}}`;
      if(clickable) b.onclick=e=>{{e.stopPropagation();pickSkill(sk,s,meta.emoji);}};
      sub.appendChild(b);
    }});
    g.appendChild(h); g.appendChild(sub); c.appendChild(g);
  }});
}}
function toggleSet(tid,sk,h) {{ const sub=document.getElementById(tid+'-sl-'+sk); const open=sub.classList.contains('open'); document.querySelectorAll('#'+tid+' .sub-list').forEach(s=>s.classList.remove('open')); document.querySelectorAll('#'+tid+' .ss-header').forEach(x=>x.classList.remove('open')); if(!open){{sub.classList.add('open');h.classList.add('open');}} }}
function pickSkill(sk,sn,emoji) {{ state.selectedSkillKey=sk+'/'+sn; document.querySelectorAll('#manualSkillTree .sub-item').forEach(b=>b.classList.remove('selected')); const btn=document.getElementById('manualSkillTree-si-'+sk+'-'+sn); if(btn)btn.classList.add('selected'); document.getElementById('ssEmoji').textContent=emoji; document.getElementById('ssName').textContent=fmt(sn); document.getElementById('selectedSkillBar').style.display='flex'; document.getElementById('userInput').focus(); }}

// ═══ SEND MESSAGE ════════════════════════════════════════════
async function sendMessage() {{
  const ta=document.getElementById('userInput');
  const text=ta.value.trim(); if(!text)return;
  const key=state.keys[state.provider]; if(!key){{openSettings();return;}}
  if(!state.skillsReady){{appendMsg('assistant','⏳ Still loading skills…');return;}}

  // Re-check server connection on every message (in case server started after page load)
  if (!state.serverConnected) await detectServer();

  document.getElementById('welcome')?.remove();
  const filesSnap=[...state.attachedFiles];
  appendMsg('user',text,null,null,filesSnap);
  state.conversation.push({{role:'user',content:text}});
  ta.value=''; ta.style.height='auto';
  state.attachedFiles=[]; document.getElementById('attachedFiles').innerHTML=''; document.getElementById('attachedFiles').style.display='none'; updateAttachBtn();
  state.msgCount++; document.getElementById('sendBtn').disabled=true; setStatusBusy(true);
  document.getElementById('infoMsgs').textContent=state.msgCount;

  try {{
    let routing=null, skillKey=null, skillContent=null;
    if(state.mode==='auto') {{
      const rl=document.getElementById('rlDot'); rl.className='rl-dot routing'; document.getElementById('rlText').textContent='Routing…';
      const rsId='rs-'+Date.now(); showRoutingRow(rsId,'Finding best skill…');
      const skillList=ALL_SKILLS.map(s=>{{const c=state.loadedSkills[s.key]; const p=c?c.split('\\n').slice(0,4).join(' ').replace(/#+/g,'').trim().slice(0,90):s.skillName; return '- '+s.key+': '+p;}}).join('\\n');
      const routerSys='You are a skill router for bioinfor-claw. Given a user message, pick the SINGLE most relevant skill.\\nSkills:\\n'+skillList+'\\nRespond ONLY with JSON: {{"skill":"set/name","reason":"one sentence","confidence":0.0-1.0}}\\nFor general chat: {{"skill":"none","reason":"general","confidence":1.0}}';
      const rr=await callProvider(routerSys,[{{role:'user',content:text+(filesSnap.length?' [Files: '+filesSnap.map(f=>f.name).join(', ')+']':'')}}],'fast');
      document.getElementById(rsId)?.remove();
      try{{routing=JSON.parse(rr.replace(/```json|```/g,'').trim());}}catch{{routing={{skill:'none',reason:'parse error',confidence:0.5}};}}
      if(routing.skill!=='none'){{skillKey=routing.skill;skillContent=state.loadedSkills[skillKey]||null;highlightSkill(skillKey);rl.className='rl-dot done';document.getElementById('rlText').textContent='→ '+routing.skill.split('/').pop();}}
      else{{rl.className='rl-dot';document.getElementById('rlText').textContent='General response';}}
      addRouteLog(routing,text,'auto');
    }} else {{
      if(!state.selectedSkillKey){{appendMsg('assistant','⚠️ Please select a skill from the sidebar.');document.getElementById('sendBtn').disabled=false;setStatusBusy(false);state.conversation.pop();return;}}
      skillKey=state.selectedSkillKey; skillContent=state.loadedSkills[skillKey]||null;
      routing={{skill:skillKey,reason:'manually selected',confidence:1.0}}; addRouteLog(routing,text,'manual');
    }}

    const thId='th-'+Date.now(); appendThinking(thId);
    // Build per-turn system — inject skill context here so LLM never echoes it
    let turnSystem = state.systemPrompt;
    if (skillKey && skillContent) {{
      turnSystem += '\n\n=== SKILL CONTEXT FOR THIS REQUEST ===\n'
        + 'Skill: ' + skillKey + '\n'
        + skillContent
        + '\n=== END SKILL CONTEXT ===\n'
        + 'Use the above SKILL.md to answer. Do NOT echo [ROUTED_SKILL], [SKILL:], [INPUT:] or any internal tags in your reply.';
    }}

    // User message is EXACTLY what they typed — no injected tags
    let userMsg = text;
    if (filesSnap.length) {{
      userMsg += '\n\nUploaded files:\n';
      filesSnap.forEach(f => {{
        userMsg += '--- ' + f.name + ' (' + f.label + ', ' + fmtSize(f.size) + ') ---\n'
          + (f.binary ? '[binary file]' : (f.content || '')) + '\n';
      }});
    }}
    const msgs = [...state.conversation.slice(0,-1), {{role:'user', content:userMsg}}];

    // Start script execution in PARALLEL with LLM call
    let execPromise = null;
    if (state.serverConnected && skillKey) {{
      const fileForExec = filesSnap.find(f => !f.binary && f.content);
      execPromise = runSkillOnServer(skillKey,
        fileForExec ? fileForExec.content : null,
        fileForExec ? fileForExec.name : null,
        skillContent);
    }}

    const reply = await callProvider(turnSystem, msgs, 'main');
    document.getElementById(thId)?.remove();
    state.conversation.push({{role:'assistant', content:reply}});

    // Wait for execution result (may already be done since it ran in parallel)
    let execResult = null;
    if (execPromise) {{
      execResult = await execPromise;
    }}

    appendMsg('assistant', reply, routing, state.mode, [], execResult);
    extractScripts(reply);
    state.msgCount++; document.getElementById('infoMsgs').textContent=state.msgCount;
  }} catch(err) {{
    document.querySelectorAll('[id^="rs-"],[id^="th-"]').forEach(e=>e.remove());
    appendMsg('assistant','**Error:** '+err.message);
  }}
  document.getElementById('sendBtn').disabled=false; setStatusBusy(false);
}}

// ═══ API CALLS ═══════════════════════════════════════════════
async function callProvider(sys,msgs,speed='main') {{
  const p=PROVIDERS[state.provider]; let model=state.model;
  if(speed==='fast'){{const fm={{anthropic:'claude-haiku-4-5-20251001',openai:'gpt-4o-mini',google:'gemini-1.5-flash',mistral:'mistral-small-latest',custom:model}};model=fm[state.provider]||model;}}
  return p.call(sys,msgs,model);
}}
async function callAnthropic(sys,msgs,model) {{
  const res=await fetch('https://api.anthropic.com/v1/messages',{{method:'POST',headers:{{'Content-Type':'application/json','x-api-key':state.keys['anthropic'],'anthropic-version':'2023-06-01'}},body:JSON.stringify({{model,max_tokens:4096,system:sys,messages:msgs}})}});
  const d=await res.json(); if(d.error)throw new Error(d.error.message);
  return d.content.map(b=>b.text||'').join('');
}}
async function callOpenAI(sys,msgs,model) {{
  // Resolve endpoint: use provider's built-in endpoint, custom input, or fallback
  const prov = PROVIDERS[state.provider];
  const baseUrls={{
    openai:'https://api.openai.com/v1/chat/completions',
    mistral:'https://api.mistral.ai/v1/chat/completions',
    minimax:'https://api.minimax.io/v1/chat/completions',
    custom:(state.customEndpoint||'https://api.openai.com/v1')+'/chat/completions',
  }};
  const url = baseUrls[state.provider] || (prov.endpoint ? prov.endpoint+'/chat/completions' : baseUrls.openai);
  const res=await fetch(url,{{method:'POST',headers:{{'Content-Type':'application/json','Authorization':'Bearer '+state.keys[state.provider]}},body:JSON.stringify({{model,messages:[{{role:'system',content:sys}},...msgs],max_tokens:4096}})}});
  const d=await res.json(); if(d.error)throw new Error(d.error.message||JSON.stringify(d.error));
  return d.choices[0].message.content;
}}
async function callGoogle(sys,msgs,model) {{
  const url='https://generativelanguage.googleapis.com/v1beta/models/'+model+':generateContent?key='+state.keys['google'];
  const contents=msgs.map(m=>{{return{{role:m.role==='assistant'?'model':'user',parts:[{{text:m.content}}]}}}});
  const res=await fetch(url,{{method:'POST',headers:{{'Content-Type':'application/json'}},body:JSON.stringify({{systemInstruction:{{parts:[{{text:sys}}]}},contents,generationConfig:{{maxOutputTokens:4096}}}})}});
  const d=await res.json(); if(d.error)throw new Error(d.error.message);
  return d.candidates[0].content.parts.map(p=>p.text||'').join('');
}}

// ═══ ROUTER HELPERS ══════════════════════════════════════════
function highlightSkill(key) {{
  document.querySelectorAll('.ss-header').forEach(h=>h.classList.remove('active-set'));
  document.querySelectorAll('.sub-item').forEach(b=>b.classList.remove('selected'));
  const parts=key.split('/'); const sk=parts[0]; const sn=parts.slice(1).join('/');
  const header=document.querySelector('#autoSkillTree-sg-'+sk+' .ss-header');
  const btn=document.getElementById('autoSkillTree-si-'+sk+'-'+sn);
  const sub=document.getElementById('autoSkillTree-sl-'+sk);
  if(header){{header.classList.add('active-set','open');}} if(sub)sub.classList.add('open'); if(btn){{btn.classList.add('selected');btn.scrollIntoView({{block:'nearest',behavior:'smooth'}});}}
}}
function addRouteLog(routing,query,mode) {{
  const panel=document.getElementById('tab-routing'); panel.querySelector('.empty-st')?.remove();
  const sk=routing.skill==='none'?'general':routing.skill.split('/').pop()||routing.skill;
  const emoji=routing.skill!=='none'?(SKILL_TREE[routing.skill?.split('/')[0]]?.emoji||'⚙️'):'💬';
  const conf=Math.round((routing.confidence||0)*100);
  const item=document.createElement('div'); item.className='rlog-item';
  item.innerHTML='<div class="rli-head"><span>'+emoji+' <span class="rli-skill">'+fmt(sk)+'</span></span><span class="rli-mode rli-'+mode+'">'+mode+'</span></div><div class="rli-body"><em>"'+query.slice(0,70)+(query.length>70?'…':'')+'"</em><br>'+(routing.reason||'')+'<span style="color:var(--muted)"> ('+conf+'%)</span></div>';
  panel.insertBefore(item,panel.firstChild);
}}
function showRoutingRow(id,text) {{ const w=document.getElementById('messages'); const d=document.createElement('div'); d.id=id; d.className='routing-row'; d.innerHTML='<div class="spin"></div><span>'+text+'</span>'; w.appendChild(d); w.scrollTop=w.scrollHeight; }}

// ═══ UI HELPERS ══════════════════════════════════════════════
function setStatusBusy(busy) {{ const d=document.getElementById('statusDot'); d.style.background=busy?'var(--amber)':'var(--g500)'; d.style.boxShadow=busy?'0 0 0 2px #ffcc8088':'0 0 0 2px var(--g300)'; }}
function switchTab(name,btn) {{ document.querySelectorAll('.ptab').forEach(t=>t.classList.remove('active')); btn.classList.add('active'); ['routing','scripts','files','session'].forEach(t=>{{document.getElementById('tab-'+t).style.display=t===name?'block':'none';}}); }}
function useEx(el) {{ const ta=document.getElementById('userInput'); ta.value=el.textContent.trim(); ta.style.height='auto'; ta.style.height=Math.min(ta.scrollHeight,130)+'px'; ta.focus(); }}
function fmt(s){{return s.replace(/-/g,' ').replace(/\\b\\w/g,c=>c.toUpperCase());}}
function esc(s){{return s.replace(/&/g,'&amp;').replace(/</g,'&lt;').replace(/>/g,'&gt;');}}
function fmtSize(b){{if(b<1024)return b+' B';if(b<1024*1024)return(b/1024).toFixed(1)+' KB';return(b/(1024*1024)).toFixed(1)+' MB';}}

// ═══ MESSAGE RENDERING ═══════════════════════════════════════
function appendMsg(role,content,routing=null,mode=null,files=[],execResult=null) {{
  const wrap=document.getElementById('messages');
  const div=document.createElement('div'); div.className='msg '+role;
  const av=document.createElement('div'); av.className='avatar';
  if(role==='user'){{av.textContent='👤';}}
  else{{const img=document.createElement('img');img.className='av-img';img.src=LOGO_URL;img.alt='';img.onerror=()=>{{av.innerHTML='🧬';}};av.appendChild(img);}}
  const bub=document.createElement('div'); bub.className='bubble';
  if(role==='user'&&files&&files.length){{
    const fl=document.createElement('div');fl.style.cssText='display:flex;flex-direction:column;gap:5px;margin-bottom:8px;';
    files.forEach(f=>{{const fi=document.createElement('div');fi.style.cssText='display:flex;align-items:center;gap:7px;padding:6px 10px;border-radius:7px;background:rgba(255,255,255,.15);border:1px solid rgba(255,255,255,.2);font-size:11.5px;';fi.innerHTML='<span style="font-size:14px">'+f.icon+'</span><div><div style="font-weight:600">'+f.name+'</div><div style="font-size:10px;opacity:.7">'+f.label+' · '+fmtSize(f.size)+'</div></div>';fl.appendChild(fi);}});
    bub.appendChild(fl);
  }}
  if(role==='assistant'&&routing&&routing.skill!=='none'){{
    const sk=routing.skill.split('/').pop();const emoji=SKILL_TREE[routing.skill?.split('/')[0]]?.emoji||'⚙️';
    const conf=Math.round((routing.confidence||0)*100);
    const tag=document.createElement('div');tag.className='route-badge';
    tag.innerHTML=emoji+' <strong>'+fmt(sk)+'</strong> <span class="rb-mode">'+(mode==='auto'?'auto':'manual')+' · '+conf+'%</span>';
    // Add "Run Script" button when server is live and script wasn't auto-run
    if(state.serverConnected && !execResult){{
      const runBtn=document.createElement('button');
      runBtn.style.cssText='margin-left:8px;padding:2px 9px;background:var(--g700);color:#fff;border:none;border-radius:5px;font-size:10.5px;font-weight:600;cursor:pointer;';
      runBtn.textContent='▶ Run Script';
      const capturedSkillKey=routing.skill;
      const capturedSkillContent=state.loadedSkills[capturedSkillKey]||null;
      runBtn.onclick=async()=>{{
        runBtn.disabled=true; runBtn.textContent='Running…';
        const r=await runSkillOnServer(capturedSkillKey,null,null,capturedSkillContent);
        const resultDiv=document.createElement('div');
        resultDiv.innerHTML=renderExecResult(r);
        runBtn.closest('.bubble').appendChild(resultDiv);
        runBtn.remove();
      }};
      tag.appendChild(runBtn);
    }}
    bub.appendChild(tag);
  }}
  const cd=document.createElement('div');cd.innerHTML=renderMD(content,routing?.skill||'');
  bub.appendChild(cd);
  if (execResult) {{
    const er=document.createElement('div');
    er.innerHTML=renderExecResult(execResult);
    bub.appendChild(er);
  }}
  div.appendChild(av);div.appendChild(bub);
  wrap.appendChild(div);wrap.scrollTop=wrap.scrollHeight;
}}
function appendThinking(id) {{
  const wrap=document.getElementById('messages');
  const div=document.createElement('div');div.className='msg assistant';div.id=id;
  const av=document.createElement('div');av.className='avatar';
  const img=document.createElement('img');img.className='av-img';img.src=LOGO_URL;img.alt='';img.onerror=()=>{{av.innerHTML='🧬';}};av.appendChild(img);
  const bub=document.createElement('div');bub.className='bubble';
  bub.style.cssText='display:flex;align-items:center;gap:10px;font-size:12px;color:var(--muted);font-style:italic;';
  bub.innerHTML='Generating… <div style="display:flex;gap:4px">'+Array(3).fill(0).map((_,i)=>'<span style="width:5px;height:5px;border-radius:50%;background:var(--g500);display:inline-block;animation:pulse 1.2s '+(i*.2)+'s infinite"></span>').join('')+'</div>';
  div.appendChild(av);div.appendChild(bub);wrap.appendChild(div);wrap.scrollTop=wrap.scrollHeight;
}}
function renderMD(text,skillKey) {{
  text=text.replace(/```(\\w+)?\\n([\\s\\S]*?)```/g,(_,lang,code)=>{{
    const l=(lang||'bash').toLowerCase();
    const cls=l==='python'?'lp-py':l==='r'?'lp-r':'lp-sh';
    const label=l==='python'?'Python':l==='r'?'R':'Bash';
    const id='c'+Math.random().toString(36).slice(2,9);
    const tag=skillKey?skillKey.split('/').pop():'';
    return '<div class="code-block"><div class="cb-header"><div class="cb-left"><span class="lang-pill '+cls+'">'+label+'</span>'+(tag?'<span class="cb-skill">'+tag+'</span>':'')+'</div><button class="cb-copy" onclick="copyCode(\\''+id+'\\')">Copy</button></div><div class="cb-body" id="'+id+'">'+esc(code.trim())+'</div></div>';
  }});
  text=text.replace(/`([^`]+)`/g,'<code>$1</code>');
  text=text.replace(/\\*\\*([^*]+)\\*\\*/g,'<strong>$1</strong>');
  text=text.replace(/\\*([^*]+)\\*/g,'<em>$1</em>');
  text=text.replace(/^###? (.+)$/gm,'<h3>$1</h3>');
  text=text.replace(/^[-*] (.+)$/gm,'<li>$1</li>');
  text=text.split(/\\n{{2,}}/).map(p=>{{p=p.trim();if(!p)return'';if(p.match(/^<(h3|li|div)/))return p;return'<p>'+p.replace(/\\n/g,'<br>')+'</p>';}}).join('');
  text=text.replace(/(<li>[\\s\\S]*?<\\/li>)+/g,m=>'<ul>'+m+'</ul>');
  return text;
}}
function copyCode(id) {{ const el=document.getElementById(id);if(!el)return; navigator.clipboard.writeText(el.textContent).then(()=>{{const btn=el.closest('.code-block').querySelector('.cb-copy');if(btn){{btn.textContent='Copied!';setTimeout(()=>btn.textContent='Copy',1500);}}}})}}
function extractScripts(text) {{
  const m=[...text.matchAll(/```(\\w+)?\\n([\\s\\S]*?)```/g)];if(!m.length)return;
  const p=document.getElementById('tab-scripts');p.querySelector('.empty-st')?.remove();
  m.forEach(x=>{{state.scriptCount++;document.getElementById('infoScripts').textContent=state.scriptCount;const lang=(x[1]||'bash').toLowerCase();const card=document.createElement('div');card.className='script-card';card.innerHTML='<div class="sc-head"><div style="display:flex;align-items:center;gap:5px"><div style="width:5px;height:5px;border-radius:50%;background:var(--g500)"></div><span>'+(lang==='python'?'Python':lang==='r'?'R':'Bash')+' · #'+state.scriptCount+'</span></div><span style="font-size:10px;color:var(--muted)">'+new Date().toLocaleTimeString()+'</span></div><div class="sc-preview">'+esc(x[2].trim().slice(0,160))+'…</div>';p.insertBefore(card,p.firstChild);}});
}}

// ═══ FILE UPLOAD ═════════════════════════════════════════════
const FILE_TYPES={{csv:{{icon:'📊',label:'CSV',text:true}},tsv:{{icon:'📊',label:'TSV',text:true}},txt:{{icon:'📄',label:'TXT',text:true}},fasta:{{icon:'🧬',label:'FASTA',text:true}},fa:{{icon:'🧬',label:'FASTA',text:true}},fastq:{{icon:'🧬',label:'FASTQ',text:true}},fq:{{icon:'🧬',label:'FASTQ',text:true}},vcf:{{icon:'🔬',label:'VCF',text:true}},bed:{{icon:'📐',label:'BED',text:true}},gff:{{icon:'📐',label:'GFF',text:true}},gtf:{{icon:'📐',label:'GTF',text:true}},json:{{icon:'📋',label:'JSON',text:true}},sam:{{icon:'🗂️',label:'SAM',text:true}}}};
function getFileType(n){{const ext=n.split('.').pop().toLowerCase();return FILE_TYPES[ext]||{{icon:'📄',label:ext.toUpperCase(),text:true}};}}
async function readFileContent(file) {{
  const ft=getFileType(file.name);
  if(!ft.text)return{{preview:'[Binary — '+fmtSize(file.size)+']',content:null,binary:true}};
  const MAX=100*1024; const blob=file.size>MAX?file.slice(0,MAX):file;
  return new Promise(res=>{{const r=new FileReader();r.onload=e=>{{const full=e.target.result;const lines=full.split('\\n');const preview=lines.slice(0,100).join('\\n');res({{preview:preview+(file.size>MAX?'\\n[truncated]':''),content:full,lines:lines.length,truncated:file.size>MAX,binary:false}});}};r.onerror=()=>res({{preview:'[Read error]',content:null,binary:true}});r.readAsText(blob);}});
}}
async function handleFileSelect(e) {{ const files=Array.from(e.target.files); e.target.value=''; await addFiles(files); }}
async function addFiles(files) {{
  for(const file of files){{
    if(state.attachedFiles.some(f=>f.name===file.name))continue;
    const ft=getFileType(file.name);
    const {{preview,content,lines,truncated,binary}}=await readFileContent(file);
    const fo={{id:'f'+Date.now()+Math.random().toString(36).slice(2,6),name:file.name,size:file.size,ext:file.name.split('.').pop().toLowerCase(),icon:ft.icon,label:ft.label,preview,content,lines,truncated,binary}};
    state.attachedFiles.push(fo); addFileChip(fo);
    const panel=document.getElementById('tab-files');panel.querySelector('.empty-st')?.remove();
    const item=document.createElement('div');item.style.cssText='border:1px solid var(--border);border-radius:8px;margin-bottom:9px;overflow:hidden;';
    item.innerHTML='<div style="display:flex;align-items:center;gap:8px;padding:8px 11px;background:var(--surf2);border-bottom:1px solid var(--border);font-size:11px;"><span style="font-size:16px">'+fo.icon+'</span><div style="flex:1"><div style="font-weight:600;color:var(--text)">'+fo.name+'</div><div style="font-size:10px;color:var(--muted)">'+fo.label+' · '+fmtSize(fo.size)+(fo.lines?' · '+fo.lines.toLocaleString()+' lines':'')+'</div></div></div>'+(fo.binary?'':('<div style="padding:8px 11px;font-size:11px;color:var(--muted);font-family:var(--mono);max-height:80px;overflow:hidden;white-space:pre-wrap;background:var(--g50)">'+esc(fo.preview.slice(0,250))+'</div>'));
    panel.insertBefore(item,panel.firstChild);
  }}
  updateAttachBtn();
}}
function addFileChip(f) {{ const c=document.getElementById('attachedFiles'); c.style.display='flex'; const chip=document.createElement('div'); chip.className='file-chip'; chip.id='chip-'+f.id; chip.innerHTML='<span>'+f.icon+'</span><span class="file-chip-name" title="'+f.name+'">'+f.name+'</span><span class="file-chip-size">'+fmtSize(f.size)+'</span><button class="file-chip-remove" onclick="removeFile(\\''+f.id+'\\')">✕</button>'; c.appendChild(chip); }}
function removeFile(id) {{ state.attachedFiles=state.attachedFiles.filter(f=>f.id!==id); document.getElementById('chip-'+id)?.remove(); if(!state.attachedFiles.length){{document.getElementById('attachedFiles').style.display='none';}} updateAttachBtn(); }}
function updateAttachBtn() {{ const btn=document.getElementById('attachBtn'); btn.classList.toggle('has-files',state.attachedFiles.length>0); }}
function initDragDrop() {{
  const ca=document.getElementById('chatArea'); const ov=document.getElementById('dropOverlay');
  ca.addEventListener('dragenter',e=>{{if(e.dataTransfer.types.includes('Files')){{e.preventDefault();ov.classList.add('active');}}}});
  ca.addEventListener('dragover',e=>e.preventDefault());
  ca.addEventListener('dragleave',e=>{{if(!ca.contains(e.relatedTarget))ov.classList.remove('active');}});
  ca.addEventListener('drop',async e=>{{e.preventDefault();ov.classList.remove('active');const files=Array.from(e.dataTransfer.files);if(files.length)await addFiles(files);}});
}}


// ═══ REAL SCRIPT EXECUTION ═══════════════════════════════════
async function runSkillOnServer(skillKey, inputText, inputFilename, skillContent) {{
  // Debug: log every step
  console.log('[RUN] serverConnected:', state.serverConnected);
  console.log('[RUN] serverUrl:', state.serverUrl);
  console.log('[RUN] skillKey:', skillKey);
  appendMiniLog('🔧 Calling /api/run for: ' + skillKey + ' on ' + state.serverUrl);

  if (!state.serverConnected) return null;

  // Extract script name from SKILL.md
  let scriptName = null;
  if (skillContent) {{
    const m = skillContent.match(/python scripts\/([\S]+\.py)/);
    if (m) scriptName = m[1];
  }}
  if (!scriptName) {{
    const parts = skillKey.split('/');
    scriptName = parts[parts.length-1].replace(/-/g,'_') + '.py';
  }}

  const execId = 'exec-' + Date.now();
  showExecSpinner(execId, '⚙️ Running ' + scriptName + '…');

  try {{
    const res = await fetch(state.serverUrl + '/api/run', {{
      method: 'POST',
      headers: {{ 'Content-Type': 'application/json' }},
      body: JSON.stringify({{
        skill_path:     skillKey,
        script_name:    scriptName,
        input_text:     inputText || null,
        input_filename: inputFilename || null,
        args: [],
      }})
    }});
    document.getElementById(execId)?.remove();
    if (!res.ok) {{
      const err = await res.json().catch(() => ({{}}));
      return {{ success: false, error: err.detail || 'Server error ' + res.status }};
    }}
    return await res.json();
  }} catch(e) {{
    document.getElementById(execId)?.remove();
    return {{ success: false, error: e.message }};
  }}
}}

function showExecSpinner(id, text) {{
  const w = document.getElementById('messages');
  const d = document.createElement('div'); d.id = id;
  d.style.cssText = 'display:flex;align-items:center;gap:9px;padding:9px 14px;background:var(--g50);border:1px solid var(--g300);border-radius:9px;font-size:11.5px;color:var(--g700);margin:0;';
  d.innerHTML = '<div style="width:14px;height:14px;border:2px solid var(--g300);border-top-color:var(--g700);border-radius:50%;animation:spin .7s linear infinite;flex-shrink:0"></div><span>' + text + '</span>';
  w.appendChild(d); w.scrollTop = w.scrollHeight;
}}

function renderExecResult(r) {{
  if (!r) return '';
  const ok = r.success !== false && r.returncode === 0;
  const files = r.output_files || [];
  let h = '<div style="margin:10px 0;padding:11px 13px;background:' + (ok?'var(--g50)':'var(--red-l)') + ';border:1px solid ' + (ok?'var(--g300)':'#ef9a9a') + ';border-radius:9px;font-size:12px;">';
  h += '<div style="font-weight:700;color:' + (ok?'var(--g700)':'var(--red)') + ';margin-bottom:6px;display:flex;align-items:center;gap:6px;">' + (ok?'✅ Script executed successfully':'❌ Script failed') + (r.script?'<span style="font-size:10px;color:var(--muted);font-family:var(--mono);font-weight:400;margin-left:4px">'+r.script+'</span>':'') + '</div>';
  if (r.error) h += '<div style="color:var(--red);font-size:11.5px;margin-bottom:6px;">'+esc(r.error)+'</div>';
  if (r.stdout) h += '<pre style="font-family:var(--mono);font-size:11px;color:var(--text2);background:var(--surf);border:1px solid var(--border);border-radius:5px;padding:8px;max-height:180px;overflow-y:auto;white-space:pre-wrap;margin:6px 0">'+esc(r.stdout.slice(-2500))+'</pre>';
  if (r.stderr && !ok) h += '<pre style="font-family:var(--mono);font-size:11px;color:var(--red);background:#fff5f5;border:1px solid #ef9a9a;border-radius:5px;padding:8px;max-height:100px;overflow-y:auto;white-space:pre-wrap;margin:6px 0">'+esc(r.stderr.slice(-1200))+'</pre>';
  if (files.length) {{
    h += '<div style="font-weight:600;color:var(--g700);margin:9px 0 5px;">📁 Output files</div>';
    files.forEach(f => {{
      const icon = /png|pdf|svg/i.test(f.ext)?'🖼️':/csv|tsv/i.test(f.ext)?'📊':'📄';
      h += '<div style="display:flex;align-items:center;gap:8px;padding:6px 9px;background:var(--surf);border:1px solid var(--border);border-radius:6px;margin-bottom:4px;">'
         + '<span>' + icon + '</span>'
         + '<span style="flex:1;font-size:11.5px;font-weight:500;">' + esc(f.name) + '</span>'
         + '<span style="font-size:10px;color:var(--muted);">' + (f.size/1024).toFixed(1) + ' KB</span>'
         + '<a href="' + state.serverUrl + f.url + '" download="' + f.name + '" style="padding:3px 9px;background:var(--g700);color:#fff;border-radius:5px;font-size:10.5px;text-decoration:none;font-weight:600;">⬇ Download</a>'
         + '</div>';
    }});
  }}
  return h + '</div>';
}}

// ── Mini log helper ─────────────────────────────────────────────────────────
function appendMiniLog(text) {{
  const wrap = document.getElementById('messages');
  const d = document.createElement('div');
  d.style.cssText = 'font-size:10.5px;color:var(--muted);padding:3px 8px;font-family:var(--mono);';
  d.textContent = text;
  wrap.appendChild(d);
  wrap.scrollTop = wrap.scrollHeight;
}}

// ═══ SESSION ═════════════════════════════════════════════════
function clearSession() {{
  state.conversation=[];state.msgCount=0;state.scriptCount=0;state.attachedFiles=[];
  document.getElementById('messages').innerHTML='';
  ['routing','scripts','files'].forEach(t=>{{document.getElementById('tab-'+t).innerHTML='<div class="empty-st"><div class="empty-ic">'+(t==='routing'?'🔀':t==='scripts'?'📝':'📂')+'</div>Appear here</div>';}});
  document.getElementById('attachedFiles').innerHTML='';document.getElementById('attachedFiles').style.display='none';
  updateAttachBtn();
  ['infoMsgs','infoScripts'].forEach(id=>document.getElementById(id).textContent='0');
  const w=document.createElement('div');w.className='welcome';w.id='welcome';
  w.innerHTML='<img class="welcome-logo" src="'+LOGO_URL+'" alt="" onerror="this.style.display=\\'none\\'"/><h2>Your 24/7 bioinformatics copilot</h2><p>Chat to perform your daily bioinformatics analysis with ease</p>';
  document.getElementById('messages').appendChild(w);
  document.getElementById('rlDot').className='rl-dot';document.getElementById('rlText').textContent='Waiting…';
}}
</script>
</body>
</html>"""
    return html

# ── Web + Execution server ──────────────────────────────────────────────────
import json as _json, subprocess as _subprocess, shutil as _shutil
from urllib.parse import urlparse, parse_qs

class Handler(BaseHTTPRequestHandler):
    html_content = None
    repo_root    = None
    results_dir  = None

    # ── Route GET ────────────────────────────────────────────────────────────
    def do_GET(self):
        path = self.path.split('?')[0]

        if path in ('/', '/index.html'):
            self._send_html(self.html_content)

        elif path in ('/api/health', '/api/health/'):
            self._send_json({'status': 'ok', 'bundled': True,
                             'execution': True,
                             'repo': str(self.repo_root),
                             'time': datetime.now().isoformat()})

        elif path == '/api/skills':
            skills = {}
            skip = {'.git','__pycache__','Assets','assets','node_modules',
                    'web','web_results','.github','docs','tests','examples','.DS_Store',
                    '.ipynb_checkpoints'}
            for sd in sorted(self.repo_root.iterdir()):
                if not sd.is_dir() or sd.name in skip or sd.name.startswith('.'): continue
                sub = {}
                # Is this category itself a skill? (top-level SKILL.md)
                # e.g. bioinformatics-plot-generator/SKILL.md
                if (sd/'SKILL.md').exists():
                    scripts_top = [f.name for f in (sd/'scripts').glob('*.py')] if (sd/'scripts').exists() else []
                    # Represent as "_self" → rendered at top of the category
                    sub['_self'] = {'has_skill_md': True, 'scripts': scripts_top}
                for kd in sorted(sd.iterdir()):
                    if not kd.is_dir() or kd.name.startswith('.') or kd.name in skip:
                        continue
                    has_md = (kd/'SKILL.md').exists()
                    scripts = [f.name for f in (kd/'scripts').glob('*.py')] if (kd/'scripts').exists() else []
                    # Skip subdirs that are clearly NOT skills: no SKILL.md AND
                    # no scripts of their own. Filters out shared scripts/ dirs
                    # and stray work directories.
                    if not has_md and not scripts:
                        continue
                    sub[kd.name] = {'has_skill_md': has_md, 'scripts': scripts}
                if sub: skills[sd.name] = sub
            self._send_json({'skills': skills})

        elif path.startswith('/api/results/'):
            parts = path.split('/')
            if len(parts) >= 5:
                run_id, filename = parts[3], '/'.join(parts[4:])
                fpath = self.results_dir / run_id / filename
                if fpath.exists():
                    self._send_file(fpath)
                    return
            self.send_error(404)

        else:
            self.send_error(404)

    # ── Route POST ───────────────────────────────────────────────────────────
    def do_POST(self):
        path = self.path.split('?')[0]
        length = int(self.headers.get('Content-Length', 0))
        body = self.rfile.read(length)

        if path == '/api/run':
            try:
                req = _json.loads(body)
                result = self._run_skill(req)
                self._send_json(result)
            except Exception as e:
                self._send_json({'error': str(e), 'success': False}, 500)

        elif path.startswith('/api/tools/'):
            # ── Agent-loop tool endpoints ─────────────────────────────────────
            # Each tool takes JSON in, returns JSON out. Used by frontend agent
            # loop (OpenAI-compatible tool-use pattern).
            try:
                req = _json.loads(body) if body else {}
                tool_name = path[len('/api/tools/'):].rstrip('/')
                result = self._dispatch_tool(tool_name, req)
                self._send_json(result)
            except Exception as e:
                import traceback
                traceback.print_exc()
                self._send_json({'error': str(e), 'success': False}, 500)

        elif path == '/api/upload':
            # Save uploaded file to shared uploads dir, return the path
            # Client sends: {"filename": "data.tsv", "content": "...file text..."}
            try:
                req = _json.loads(body)
                fname = req.get('filename', 'upload.txt')
                content = req.get('content', '')
                # Sanitize filename
                safe_name = Path(fname).name  # strip any path components
                upload_dir = self.results_dir / '_uploads'
                upload_dir.mkdir(parents=True, exist_ok=True)
                fpath = upload_dir / safe_name
                fpath.write_text(content, encoding='utf-8')
                print(f"  [upload] Saved {safe_name} ({len(content)} chars)")
                self._send_json({
                    'success': True,
                    'path': str(fpath),
                    'filename': safe_name,
                    'size': len(content),
                })
            except Exception as e:
                self._send_json({'error': str(e), 'success': False}, 500)

        else:
            self.send_error(404)

    # ── OPTIONS (CORS preflight) ──────────────────────────────────────────────
    def do_OPTIONS(self):
        self.send_response(200)
        self._cors_headers()
        self.end_headers()

    # ────────────────────────────────────────────────────────────────────────
    # AGENT-LOOP TOOL ENDPOINTS
    # ────────────────────────────────────────────────────────────────────────
    # These are the primitives an LLM agent-loop calls via tool-use.
    # Each one maps 1:1 to a tool schema exposed to the model.
    # The frontend picks a model → gets a list of tools → calls the LLM →
    # dispatches each tool_use block to /api/tools/{name} → feeds results back.

    def _dispatch_tool(self, name, args):
        """Route tool name to handler."""
        handlers = {
            'list_skills':        self._tool_list_skills,
            'read_skill':         self._tool_read_skill,
            'list_skill_scripts': self._tool_list_skill_scripts,
            'run_script':         self._tool_run_script,
            'list_files':         self._tool_list_files,
            'read_file':          self._tool_read_file,
            'write_text_file':    self._tool_write_text_file,
        }
        h = handlers.get(name)
        if not h:
            return {'success': False, 'error': f'Unknown tool: {name}',
                    'available_tools': list(handlers.keys())}
        try:
            return h(args or {})
        except Exception as e:
            import traceback
            return {'success': False, 'error': str(e),
                    'traceback': traceback.format_exc()[-2000:]}

    # ── Tool: list_skills ──────────────────────────────────────────────────
    def _tool_list_skills(self, args):
        """Return all available skills as a flat list.
        Each item: {id, category, name, has_skill_md, scripts[]}

        A "skill" is any directory that has either a SKILL.md file of its own
        or a scripts/ folder with at least one .py script. This covers:
          - Regular skills:    category/name/SKILL.md
          - Super-skills:      category/SKILL.md       (the category IS a skill)
          - Unnamed skills:    category/name/scripts/*.py (no SKILL.md yet)
        Non-skill subdirs (like shared scripts/, .ipynb_checkpoints, etc.) are
        filtered out."""
        skip = {'.git', '__pycache__', 'Assets', 'assets', 'node_modules',
                'web', 'web_results', '.github', 'docs', 'tests',
                'examples', '.DS_Store', '.ipynb_checkpoints'}

        def describe(skill_md_path):
            """Pull one-line description from SKILL.md."""
            if not skill_md_path.exists():
                return ''
            try:
                t = skill_md_path.read_text(encoding='utf-8', errors='replace')
                for line in t.splitlines()[:30]:
                    line = line.strip()
                    if (line and not line.startswith('#')
                            and not line.startswith('---')
                            and not line.startswith('```')
                            and not line.lower().startswith('name:')
                            and not line.lower().startswith('description:')):
                        return line[:200]
            except Exception:
                pass
            return ''

        skills = []
        for sd in sorted(self.repo_root.iterdir()):
            if not sd.is_dir() or sd.name in skip or sd.name.startswith('.'):
                continue

            # Category-as-skill (top-level SKILL.md)
            sm_top = sd / 'SKILL.md'
            if sm_top.exists():
                scripts_top = []
                if (sd / 'scripts').exists():
                    scripts_top = sorted(f.name for f in (sd / 'scripts').glob('*.py'))
                skills.append({
                    'id':           sd.name,
                    'category':     sd.name,
                    'name':         sd.name,
                    'description':  describe(sm_top),
                    'has_skill_md': True,
                    'scripts':      scripts_top,
                    'is_super':     True,
                })

            # Child skills under this category
            for kd in sorted(sd.iterdir()):
                if not kd.is_dir() or kd.name.startswith('.') or kd.name in skip:
                    continue
                sm = kd / 'SKILL.md'
                scripts = []
                if (kd / 'scripts').exists():
                    scripts = sorted(f.name for f in (kd / 'scripts').glob('*.py'))
                # Must have SKILL.md OR its own scripts/*.py — otherwise it's
                # not a skill (shared resources, cache dirs, etc.)
                if not sm.exists() and not scripts:
                    continue
                skills.append({
                    'id':           f'{sd.name}/{kd.name}',
                    'category':     sd.name,
                    'name':         kd.name,
                    'description':  describe(sm),
                    'has_skill_md': sm.exists(),
                    'scripts':      scripts,
                    'is_super':     False,
                })
        return {'success': True, 'skills': skills, 'count': len(skills)}

    # ── Tool: read_skill ───────────────────────────────────────────────────
    def _tool_read_skill(self, args):
        """Return SKILL.md content. args: {skill_id: 'category/name' or 'category'}."""
        skill_id = (args.get('skill_id') or '').strip().strip('/')
        if not skill_id:
            return {'success': False, 'error': 'skill_id is required'}
        # Support both super-skills (bare category) and child skills (category/name)
        if '/' in skill_id:
            cat, name = skill_id.split('/', 1)
            sm = self.repo_root / cat / name / 'SKILL.md'
        else:
            sm = self.repo_root / skill_id / 'SKILL.md'
        if not sm.exists():
            return {'success': False, 'error': f'SKILL.md not found for {skill_id}'}
        try:
            content = sm.read_text(encoding='utf-8', errors='replace')
            # Cap size to avoid overwhelming LLM context
            truncated = len(content) > 12000
            if truncated:
                content = content[:12000] + '\n\n[... truncated ...]'
            return {'success': True, 'skill_id': skill_id,
                    'content': content, 'truncated': truncated}
        except Exception as e:
            return {'success': False, 'error': str(e)}

    # ── Tool: list_skill_scripts ───────────────────────────────────────────
    def _tool_list_skill_scripts(self, args):
        """Return all scripts for a skill with their --help output.
        args: {skill_id: 'category/name'}."""
        skill_id = (args.get('skill_id') or '').strip().strip('/')
        if not skill_id:
            return {'success': False, 'error': 'skill_id is required'}
        if '/' in skill_id:
            cat, name = skill_id.split('/', 1)
            skill_dir = self.repo_root / cat / name
        else:
            skill_dir = self.repo_root / skill_id
        if not skill_dir.exists():
            return {'success': False, 'error': f'Skill not found: {skill_id}'}

        scripts = []
        script_dir = skill_dir / 'scripts'
        paths = []
        if script_dir.exists():
            paths = sorted(script_dir.glob('*.py'))
        else:
            paths = sorted(skill_dir.glob('*.py'))

        for sp in paths:
            rel = sp.relative_to(skill_dir).as_posix()
            help_text = ''
            # Don't actually exec the script for --help (slow + risky);
            # instead parse argparse directly from source.
            try:
                src = sp.read_text(encoding='utf-8', errors='replace')
                import re as _re
                # Capture add_argument() calls with their help text
                help_lines = []
                for m in _re.finditer(
                    r'add_argument\s*\(\s*([^)]{5,400}?)\)',
                    src, _re.DOTALL
                ):
                    inside = m.group(1)
                    flags = _re.findall(r"['\"](-{1,2}[\w\-]+)['\"]", inside)
                    help_m = _re.search(r"help\s*=\s*['\"]([^'\"]+)['\"]", inside)
                    required = 'required=True' in inside
                    default_m = _re.search(r"default\s*=\s*([^,\n)]+)", inside)
                    if flags:
                        parts = [', '.join(flags)]
                        if required:
                            parts.append('(required)')
                        if default_m:
                            parts.append(f'default={default_m.group(1).strip()}')
                        if help_m:
                            parts.append(help_m.group(1).strip())
                        help_lines.append('  ' + '  '.join(parts))
                help_text = '\n'.join(help_lines[:60])
            except Exception as e:
                help_text = f'(could not parse: {e})'
            scripts.append({
                'script':    rel,
                'full_path': str(sp),
                'help':      help_text,
            })
        return {'success': True, 'skill_id': skill_id,
                'scripts': scripts, 'count': len(scripts)}

    # ── Tool: run_script ───────────────────────────────────────────────────
    def _tool_run_script(self, args):
        """Execute a script with explicit CLI args.
        args: {
          skill_id: 'category/name',
          script:   'scripts/foo.py' or 'foo.py',
          args:     ['--flag', 'value', ...],       # model-provided flags
          input_data:   (optional) text to write to a temp file
          input_file:   (optional) path to pre-uploaded file
          input_flag:   (optional) which flag to use for the input file
          timeout:      (optional) seconds, default 300
        }
        Returns: {success, run_id, returncode, stdout, stderr, output_files, command, duration}
        """
        import subprocess as _sp
        import shutil as _shutil
        import time as _time
        import re as _re
        import uuid as _uuid

        skill_id = (args.get('skill_id') or '').strip().strip('/')
        if not skill_id:
            return {'success': False, 'error': 'skill_id is required'}
        if '/' in skill_id:
            cat, name = skill_id.split('/', 1)
            skill_dir = self.repo_root / cat / name
        else:
            skill_dir = self.repo_root / skill_id
        if not skill_dir.exists():
            return {'success': False, 'error': f'Skill not found: {skill_id}'}

        script_name = (args.get('script') or '').strip().lstrip('/')
        if not script_name:
            return {'success': False, 'error': 'script is required'}
        # Resolve: try scripts/ then skill root
        candidates = [
            skill_dir / script_name,
            skill_dir / 'scripts' / script_name,
        ]
        script_path = next((c for c in candidates if c.exists()), None)
        if not script_path:
            return {'success': False,
                    'error': f'Script not found in {skill_id}: {script_name}',
                    'tried': [str(c) for c in candidates]}

        cli_args = args.get('args') or []
        if not isinstance(cli_args, list):
            return {'success': False, 'error': 'args must be a list of strings'}
        cli_args = [str(a) for a in cli_args]

        # ── Auto-redirect to standalone scripts for co-expression/co-essentiality ──
        if 'depmap-analysis-for-gene' in skill_id:
            _has_coexpr_module = False
            _has_coess_module = False
            for idx, a in enumerate(cli_args):
                if a == '--modules' and idx + 1 < len(cli_args):
                    mods = cli_args[idx + 1].lower()
                    if 'coexpression' in mods and 'coessentiality' not in mods:
                        _has_coexpr_module = True
                    elif 'coessentiality' in mods:
                        _has_coess_module = True
            if _has_coexpr_module and 'coexpression' not in script_path.stem:
                alt = skill_dir / 'scripts' / 'depmap_coexpression.py'
                if alt.exists():
                    print(f"  [tool/run_script] ↗ Redirecting --modules coexpression → standalone depmap_coexpression.py")
                    script_path = alt
                    # Remove --modules from cli_args
                    new_args = []
                    skip = False
                    for a in cli_args:
                        if skip: skip = False; continue
                        if a == '--modules': skip = True; continue
                        new_args.append(a)
                    cli_args = new_args
            elif _has_coess_module and 'coessentiality' not in script_path.stem:
                alt = skill_dir / 'scripts' / 'depmap_coessentiality.py'
                if alt.exists():
                    print(f"  [tool/run_script] ↗ Redirecting --modules coessentiality → standalone depmap_coessentiality.py")
                    script_path = alt
                    new_args = []
                    skip = False
                    for a in cli_args:
                        if skip: skip = False; continue
                        if a == '--modules': skip = True; continue
                        new_args.append(a)
                    cli_args = new_args

        # Unique run dir for outputs
        run_id = datetime.now().strftime('%Y%m%d_%H%M%S_') + _uuid.uuid4().hex[:6]
        out_dir = self.results_dir / run_id
        out_dir.mkdir(parents=True, exist_ok=True)

        # Optional: handle input data/file
        input_path = None
        input_data = args.get('input_data')
        input_file = args.get('input_file')
        input_flag = (args.get('input_flag') or '').lstrip('-') or None
        if input_file and Path(input_file).exists():
            input_path = Path(input_file)
        elif input_data:
            ext = args.get('input_ext') or '.txt'
            input_path = out_dir / f'input{ext}'
            input_path.write_text(str(input_data), encoding='utf-8')

        # Auto-wire --output / --outdir if the script has them and the model
        # didn't include them. This frees the LLM from having to know about
        # our results_dir convention.
        flags = self._script_flags(script_path)
        provided_flags = set()
        for a in cli_args:
            if isinstance(a, str) and a.startswith('--'):
                provided_flags.add(a.split('=', 1)[0].lstrip('-'))

        auto_args = []
        if input_path and input_flag and input_flag in flags and input_flag not in provided_flags:
            auto_args += [f'--{input_flag}', str(input_path)]
        elif input_path and not input_flag:
            # Try common input flag names
            for cand in ('input', 'data', 'file', 'counts', 'input-file',
                         'matrix', 'table', 'gene-file'):
                if cand in flags and cand not in provided_flags:
                    auto_args += [f'--{cand}', str(input_path)]
                    break

        if 'outdir' in flags and 'outdir' not in provided_flags:
            auto_args += ['--outdir', str(out_dir)]
        elif 'output' in flags and 'output' not in provided_flags:
            ext = self._script_output_ext(script_path)
            auto_args += ['--output', str(out_dir / f'{script_path.stem}_result{ext}')]
        elif 'output-prefix' in flags and 'output-prefix' not in provided_flags:
            auto_args += ['--output-prefix', str(out_dir / script_path.stem)]

        # ── Auto-wire DepMap / preflight data files ─────────────────────────
        # If the script needs data files (e.g. --expression-file) and the LLM
        # didn't provide them (or gave a non-existent path), fill from cache.
        preflight = self.SKILL_PREFLIGHT.get(skill_id)
        if preflight:
            file_map = preflight['data_provider']['file_map']
            cache_dir = self.repo_root / 'web_results' / preflight['data_provider']['cache_dir']
            for flag, pattern in file_map.items():
                if flag not in flags:
                    continue
                # Check if LLM already provided a VALID path for this flag
                llm_provided_valid = False
                if flag in provided_flags:
                    for idx, a in enumerate(cli_args):
                        if a == f'--{flag}' and idx + 1 < len(cli_args):
                            if Path(cli_args[idx + 1]).exists():
                                llm_provided_valid = True
                            break
                if not llm_provided_valid:
                    # Search in cache and also in recent web_results subdirs
                    cached = self._find_cached_files(pattern, cache_dir)
                    if not cached:
                        # Broader search: look in web_results subdirectories
                        import fnmatch as _fnm
                        wr = self.repo_root / 'web_results'
                        if wr.exists():
                            for sub in sorted(wr.iterdir(), key=lambda d: d.name, reverse=True):
                                if sub.is_dir() and sub.name != preflight['data_provider']['cache_dir']:
                                    for f in sub.iterdir():
                                        if _fnm.fnmatch(f.name, pattern):
                                            cached = [f]
                                            break
                                if cached:
                                    break
                    if cached:
                        # Remove invalid LLM-provided path from cli_args
                        if flag in provided_flags:
                            new_args = []
                            skip = False
                            for a in cli_args:
                                if skip:
                                    skip = False
                                    continue
                                if a == f'--{flag}':
                                    skip = True
                                    continue
                                new_args.append(a)
                            cli_args = new_args
                            provided_flags.discard(flag)
                        auto_args += [f'--{flag}', str(cached[0])]
                        print(f"  [tool/run_script] Auto-wired --{flag} → {cached[0].name}")

        # Build command. Prevent "-value" from being mistaken for a flag by
        # merging with = syntax.
        cmd = [sys.executable, str(script_path)]
        i = 0
        raw = auto_args + cli_args
        while i < len(raw):
            a = raw[i]
            if (isinstance(a, str) and a.startswith('--')
                    and i + 1 < len(raw)
                    and isinstance(raw[i + 1], str)
                    and raw[i + 1].startswith('-')
                    and not raw[i + 1].startswith('--')):
                cmd.append(f'{a}={raw[i + 1]}')
                i += 2
            else:
                cmd.append(a)
                i += 1

        # Snapshot skill_dir for file rescue
        before_files = set()
        try:
            before_files = {f.name for f in skill_dir.iterdir() if f.is_file()}
        except Exception:
            pass
        run_start = _time.time()
        timeout = int(args.get('timeout') or 300)

        print(f"  [tool/run_script] {skill_id} :: {' '.join(cmd[1:])}")

        try:
            _shared_dir = str(Path(self.repo_root) / '_shared')
            _pypath = os.pathsep.join([str(skill_dir), _shared_dir])
            result = _sp.run(cmd, capture_output=True, text=True,
                             timeout=timeout, cwd=str(skill_dir),
                             env={**os.environ, 'PYTHONPATH': _pypath})
        except _sp.TimeoutExpired:
            return {'success': False, 'run_id': run_id,
                    'error': f'timeout after {timeout}s',
                    'command': ' '.join(cmd)}
        except Exception as e:
            return {'success': False, 'run_id': run_id,
                    'error': str(e), 'command': ' '.join(cmd)}

        duration = _time.time() - run_start

        # Rescue new output files from skill_dir → out_dir
        OUTPUT_EXTS = {'.png', '.svg', '.pdf', '.jpg', '.jpeg', '.tsv',
                       '.csv', '.xlsx', '.xls', '.json', '.html',
                       '.txt', '.fasta', '.fa', '.bed', '.gff', '.gtf',
                       '.bam', '.vcf', '.log'}
        try:
            for f in skill_dir.iterdir():
                if not f.is_file():
                    continue
                try:
                    mtime = f.stat().st_mtime
                except Exception:
                    continue
                is_new = f.name not in before_files
                is_fresh = mtime >= run_start - 1
                if (is_new or is_fresh) and f.suffix.lower() in OUTPUT_EXTS:
                    dest = out_dir / f.name
                    if dest.exists():
                        continue
                    try:
                        _shutil.move(str(f), str(dest))
                        print(f"  [tool/run_script] rescued {f.name}")
                    except Exception:
                        pass
        except Exception:
            pass

        output_files = []
        for f in sorted(out_dir.rglob('*')):
            if not f.is_file():
                continue
            if input_path and f.name == input_path.name:
                continue
            rel = f.relative_to(out_dir)
            ext = f.suffix.lstrip('.').upper()
            output_files.append({
                'name': str(rel),
                'size': f.stat().st_size,
                'url':  f'/api/results/{run_id}/{rel}',
                'ext':  ext,
            })

        return {
            'success':      result.returncode == 0,
            'run_id':       run_id,
            'returncode':   result.returncode,
            'stdout':       (result.stdout or '')[-8000:],
            'stderr':       (result.stderr or '')[-4000:],
            'output_files': output_files,
            'command':      ' '.join(cmd),
            'duration':     round(duration, 2),
            'script':       script_path.name,
            'skill_id':     skill_id,
        }

    # ── Tool: list_files ───────────────────────────────────────────────────
    def _tool_list_files(self, args):
        """List files in a results run directory (recursive). args: {run_id: str}."""
        run_id = (args.get('run_id') or '').strip().strip('/')
        if not run_id or '..' in run_id or '/' in run_id:
            return {'success': False, 'error': 'invalid run_id'}
        d = self.results_dir / run_id
        if not d.exists() or not d.is_dir():
            return {'success': False, 'error': f'run_id not found: {run_id}'}
        files = []
        for f in sorted(d.rglob('*')):
            if f.is_file():
                rel = f.relative_to(d)
                files.append({
                    'name': str(rel),
                    'size': f.stat().st_size,
                    'url':  f'/api/results/{run_id}/{rel}',
                    'ext':  f.suffix.lstrip('.').upper(),
                })
        return {'success': True, 'run_id': run_id, 'files': files,
                'count': len(files)}

    # ── Tool: read_file ────────────────────────────────────────────────────
    def _tool_read_file(self, args):
        """Read a text file from a run dir. args: {run_id, filename, max_bytes?}.
        Returns up to max_bytes (default 8000) of text. Binary files: descriptor only."""
        run_id = (args.get('run_id') or '').strip().strip('/')
        fname = (args.get('filename') or '').strip()
        max_bytes = int(args.get('max_bytes') or 8000)
        if not run_id or not fname or '..' in run_id or '..' in fname:
            return {'success': False, 'error': 'invalid run_id or filename'}
        run_dir = self.results_dir / run_id
        fpath = run_dir / fname
        if not fpath.exists() or not fpath.is_file():
            # Fuzzy fallback 1: search for basename anywhere in the run directory
            basename = Path(fname).name
            candidates = list(run_dir.rglob(basename)) if run_dir.is_dir() else []
            if len(candidates) == 1:
                fpath = candidates[0]
                fname = str(fpath.relative_to(run_dir))
            elif len(candidates) > 1:
                return {'success': False,
                        'error': f'file not found at {run_id}/{fname}; '
                                 f'multiple matches for "{basename}": '
                                 + ', '.join(str(c.relative_to(run_dir)) for c in candidates[:5])}
            else:
                # Fuzzy fallback 2: search ALL recent run directories for this file
                # (handles wrong run_id from LLM hallucination)
                global_candidates = []
                try:
                    for d in sorted(self.results_dir.iterdir(), reverse=True):
                        if not d.is_dir() or d.name.startswith('_'):
                            continue
                        match = d / basename
                        if match.exists() and match.is_file():
                            global_candidates.append((d.name, match))
                        else:
                            # Also check subdirs
                            for m in d.rglob(basename):
                                if m.is_file():
                                    global_candidates.append((d.name, m))
                                    break
                        if len(global_candidates) >= 3:
                            break
                except Exception:
                    pass

                if len(global_candidates) == 1:
                    # Auto-resolve: found in exactly one other run dir
                    found_run_id, fpath = global_candidates[0]
                    run_id = found_run_id
                    run_dir = self.results_dir / run_id
                    fname = str(fpath.relative_to(run_dir))
                    print(f"  [read_file] Auto-resolved {basename} → {run_id}/{fname}")
                elif global_candidates:
                    return {'success': False,
                            'error': f'file not found in {run_id}; found "{basename}" in other runs: '
                                     + ', '.join(f'{rid}/{f.relative_to(self.results_dir / rid)}'
                                                 for rid, f in global_candidates)}
                else:
                    avail = [str(f.relative_to(run_dir))
                             for f in sorted(run_dir.rglob('*')) if f.is_file()] if run_dir.is_dir() else []
                    return {'success': False,
                            'error': f'file not found: {run_id}/{fname}',
                            'available_files': avail[:20]}
        size = fpath.stat().st_size
        TEXT_EXTS = {'.txt', '.csv', '.tsv', '.json', '.md', '.log', '.html',
                     '.xml', '.yaml', '.yml', '.bed', '.gff', '.gtf', '.vcf',
                     '.fasta', '.fa', '.py', '.sh', '.r'}
        if fpath.suffix.lower() not in TEXT_EXTS:
            return {'success': True, 'run_id': run_id, 'filename': fname,
                    'size': size, 'binary': True,
                    'url': f'/api/results/{run_id}/{fname}',
                    'preview': f'[binary file, {size} bytes, {fpath.suffix} — view via URL]'}
        try:
            data = fpath.read_bytes()[:max_bytes]
            text = data.decode('utf-8', errors='replace')
            truncated = size > max_bytes
            return {'success': True, 'run_id': run_id, 'filename': fname,
                    'size': size, 'binary': False,
                    'content': text, 'truncated': truncated,
                    'url': f'/api/results/{run_id}/{fname}'}
        except Exception as e:
            return {'success': False, 'error': str(e)}

    # ── Tool: write_text_file ────────────────────────────────────────────────
    def _tool_write_text_file(self, args):
        """Write text content to a temp file in _uploads/ and return its server
        path.  This is meant for the agent to save data the user pasted in
        the chat (e.g. a gene list, a sample table, a config snippet) so it
        can be fed as input_file to run_script.

        args: {content: str, filename?: str}
        Returns: {success, path, filename, size}
        """
        content = args.get('content')
        if content is None or not isinstance(content, str) or not content.strip():
            return {'success': False, 'error': 'content is required and must be non-empty text'}
        fname = (args.get('filename') or 'agent_input.txt').strip()
        # Sanitize
        safe_name = Path(fname).name
        if not safe_name:
            safe_name = 'agent_input.txt'
        # Deduplicate to avoid overwriting a previous file
        upload_dir = self.results_dir / '_uploads'
        upload_dir.mkdir(parents=True, exist_ok=True)
        stem = Path(safe_name).stem
        suffix = Path(safe_name).suffix or '.txt'
        fpath = upload_dir / safe_name
        counter = 1
        while fpath.exists():
            fpath = upload_dir / f'{stem}_{counter}{suffix}'
            counter += 1
        fpath.write_text(content, encoding='utf-8')
        print(f'  [write_text_file] Saved {fpath.name} ({len(content)} chars)')
        return {
            'success': True,
            'path': str(fpath),
            'filename': fpath.name,
            'size': len(content),
        }

    # ── Inspect a script's accepted CLI flags ────────────────────────────────
    @staticmethod
    def _script_flags(script_path):
        """Return set of --flag names accepted by the script (from add_argument calls).
        Captures ALL flag aliases, e.g. add_argument('-o', '--output', ...) → {'o', 'output'}.
        """
        import re as _re
        flags = set()
        try:
            src = script_path.read_text(encoding='utf-8', errors='replace')
            # Match the full add_argument(...) opening up to the first positional-terminating char
            # Then extract all '-flag' or '--flag' strings within
            for m in _re.finditer(r'add_argument\s*\(([^)]*?)(?:,\s*(?:type|default|help|action|nargs|required|choices|const|dest|metavar)\s*=|\))', src):
                inside = m.group(1)
                # Extract every '-...' or '--...' string literal from inside
                for flag_match in _re.finditer(r"['\"](-{1,2}[\w-]+)['\"]", inside):
                    flags.add(flag_match.group(1).lstrip('-'))
            # Fallback: also do the simpler first-flag pattern in case of odd formatting
            for m in _re.finditer(r"add_argument\s*\(\s*['\"](-{1,2}[\w-]+)['\"]", src):
                flags.add(m.group(1).lstrip('-'))
        except Exception:
            pass
        return flags

    @staticmethod
    def _script_required_flags(script_path):
        """Return set of --flag names that are required by the script."""
        import re as _re
        required = set()
        try:
            src = script_path.read_text(encoding='utf-8', errors='replace')
            # Match add_argument("--flag", ..., required=True)
            for m in _re.finditer(
                r"add_argument\s*\(\s*['\"](-{1,2}[\w-]+)['\"]"
                r"[^)]*required\s*=\s*True", src):
                required.add(m.group(1).lstrip('-'))
        except Exception:
            pass
        return required

    @staticmethod
    def _script_output_ext(script_path):
        """Detect the expected output file extension from multiple signals.

        Signal priority (most reliable first):
          1. Script name: plot_*.py / *_plot.py / *_volcano.py → .png
          2. Imports: matplotlib/seaborn/plotly → .png
          3. Calls: savefig/to_csv/to_excel/to_json/write_html
          4. --output help text keywords
          5. Generic fallback: .txt
        """
        import re as _re
        name_lower = script_path.stem.lower()

        # ── Signal 1: script name heuristics (fastest, usually decisive) ────
        PLOT_NAME_HINTS = ('plot', 'chart', 'figure', 'visualiz', 'heatmap',
                           'volcano', 'scatter', 'boxplot', 'barplot', 'violin',
                           'histogram', 'umap', 'tsne', 'pca', 'survival',
                           'manhattan', 'venn')
        if any(h in name_lower for h in PLOT_NAME_HINTS):
            return '.png'

        try:
            src = script_path.read_text(encoding='utf-8', errors='replace')
        except Exception:
            return '.txt'

        src_lower = src.lower()

        # ── Signal 2: plotting library imports ──────────────────────────────
        PLOT_IMPORTS = ('import matplotlib', 'from matplotlib',
                        'import seaborn', 'import plotly',
                        'import plotnine', 'from plotnine')
        if any(imp in src_lower for imp in PLOT_IMPORTS):
            # Extra confirmation: has a savefig or show() call?
            if 'savefig' in src_lower or 'write_image' in src_lower:
                return '.png'
            # Just the import is strong enough for plotting skills
            return '.png'

        # ── Signal 3: explicit file-writing calls ───────────────────────────
        if _re.search(r'\.to_excel\(', src) or 'openpyxl' in src_lower:
            return '.xlsx'
        if _re.search(r'\.to_html\(|write_html\(', src):
            return '.html'
        if _re.search(r'\.to_json\(|json\.dump\(', src):
            return '.json'
        if _re.search(r'\.to_csv\([^)]*sep\s*=\s*["\']\\t["\']', src):
            return '.tsv'
        if _re.search(r'\.to_csv\(', src):
            return '.csv'

        # ── Signal 4: --output help text keywords ───────────────────────────
        m = _re.search(
            r'add_argument\s*\(\s*["\']--output["\']'
            r'[^)]*help\s*=\s*["\']([^"\']+)["\']', src)
        if m:
            h = m.group(1).lower()
            if 'tsv' in h:   return '.tsv'
            if 'csv' in h:   return '.csv'
            if 'png' in h:   return '.png'
            if 'svg' in h:   return '.svg'
            if 'pdf' in h:   return '.pdf'
            if 'xlsx' in h or 'workbook' in h or 'excel' in h: return '.xlsx'
            if 'json' in h:  return '.json'
            if 'html' in h:  return '.html'
            if 'plot' in h or 'figure' in h or 'image' in h: return '.png'

        return '.txt'  # safe fallback

    # ── Intelligent pre-flight: resolve data dependencies ────────────────────
    # Maps skill paths to their data providers — scripts that can download
    # or generate the required input files automatically
    SKILL_PREFLIGHT = {
        # DepMap analysis needs DepMap data files → run downloader first
        'gene-centered-analysis/depmap-analysis-for-gene': {
            'needs_data_flags': ['expression-file', 'mutation-file', 'copy-number-file',
                                 'essentiality-file', 'metadata-file'],
            'data_provider': {
                'script': 'gene-centered-analysis/depmap-analysis-for-gene/scripts/depmap_downloader.py',
                'alt_script': 'public-datasets-access-and-download/depmap-data-download/scripts/depmap_download_from_api.py',
                'args': ['--all'],
                'cache_dir': 'depmap_cache',
                'file_map': {
                    'expression-file':   '*xpression*.csv',
                    'mutation-file':     '*utation*.csv',
                    'copy-number-file':  '*opy*umber*.csv',
                    'essentiality-file': '*ffect*.csv',
                    'metadata-file':     '*odel*.csv',
                },
            }
        },
    }

    def _find_cached_files(self, pattern, cache_dir):
        """Find files matching a glob pattern in cache directory."""
        import fnmatch
        matches = []
        if cache_dir.exists():
            for f in cache_dir.iterdir():
                if fnmatch.fnmatch(f.name, pattern):
                    matches.append(f)
        return sorted(matches, key=lambda f: f.stat().st_mtime, reverse=True)

    # Maps DepMap modules to which data files they need
    DEPMAP_MODULE_FILES = {
        'expression':     ['expression-file', 'metadata-file'],
        'mutation':       ['mutation-file', 'metadata-file'],
        'copy_number':    ['copy-number-file', 'metadata-file'],
        'essentiality':   ['essentiality-file', 'metadata-file'],
        'coexpression':   ['expression-file'],
        'coessentiality': ['essentiality-file'],
    }

    def _run_preflight(self, skill_path, out_dir, params):
        """Check if a skill needs data that must be downloaded first.
        Returns dict of extra CLI flags to add, or empty dict."""
        preflight = self.SKILL_PREFLIGHT.get(skill_path)
        if not preflight:
            return {}

        extra_flags = {}
        cache_dir = self.repo_root / 'web_results' / preflight['data_provider']['cache_dir']
        file_map = preflight['data_provider']['file_map']

        # Determine which files are actually needed based on detected modules
        modules_str = params.get('modules', 'full')
        if modules_str == 'full':
            needed_flags = set(file_map.keys())
        else:
            needed_flags = set()
            for mod in modules_str.split(','):
                mod = mod.strip()
                needed_flags.update(self.DEPMAP_MODULE_FILES.get(mod, []))
            if not needed_flags:
                needed_flags = set(file_map.keys())  # fallback to all

        print(f"  [preflight] modules={modules_str} → need files: {sorted(needed_flags)}")

        # Check if we already have cached data for needed files
        all_found = True
        for flag in needed_flags:
            pattern = file_map.get(flag)
            if not pattern:
                continue
            matches = self._find_cached_files(pattern, cache_dir)
            if matches:
                extra_flags[flag] = str(matches[0])
            else:
                all_found = False

        if all_found:
            print(f"  [preflight] Found cached data for {skill_path}")
            return extra_flags

        # Need to download — find and run the data provider script
        print(f"  [preflight] Downloading data for {skill_path}...")
        provider = preflight['data_provider']
        script = self.repo_root / provider['script']
        if not script.exists():
            script = self.repo_root / provider['alt_script']
        if not script.exists():
            print(f"  [preflight] WARNING: downloader script not found")
            return extra_flags  # return whatever we found

        cache_dir.mkdir(parents=True, exist_ok=True)
        dl_cmd = [sys.executable, str(script), '--outdir', str(cache_dir)] + provider.get('args', [])
        print(f"  [preflight] Running: {' '.join(dl_cmd)}")

        try:
            _shared_dir = str(Path(__file__).resolve().parent / '_shared')
            _pypath = os.pathsep.join([str(script.parent), _shared_dir])
            dl_result = _subprocess.run(
                dl_cmd, capture_output=True, text=True, timeout=600,
                cwd=str(script.parent),
                env={**os.environ, 'PYTHONPATH': _pypath})
            if dl_result.returncode != 0:
                print(f"  [preflight] Download failed: {dl_result.stderr[:500]}")
            else:
                print(f"  [preflight] Download complete")
        except Exception as e:
            print(f"  [preflight] Download error: {e}")

        # Re-scan for downloaded files
        for flag, pattern in file_map.items():
            if flag not in extra_flags:
                matches = self._find_cached_files(pattern, cache_dir)
                if matches:
                    extra_flags[flag] = str(matches[0])

        return extra_flags

    # ── Smart parameter extraction from conversation text ────────────────────
    @staticmethod
    def _extract_gene_from_text(text):
        """Extract gene name from any text — handles lowercase, mixed case, context."""
        import re as _re
        if not text:
            return None
        NON_GENES = {'TCGA','GTEX','GEO','RNA','DNA','CNV','MAF','API','CSV','TSV',
                     'PDF','PNG','SVG','URL','AND','FOR','THE','WITH','FROM','NOT',
                     'USE','ALL','TOP','LOW','HIGH','SHOW','PLOT','RUN','GET','SET',
                     'CHECK','EXPRESSION','NORMAL','TISSUE','CANCER','GENE','GENES',
                     'ANALYSIS','DATA','DOWNLOAD','SEARCH','FIND','LOOK','QUERY',
                     'ABOUT','PAN','MODE','HELP','WHAT','HOW','TRUE','FALSE','NONE',
                     'IN','OF','TO','BY','IS','AS','AT','BE','IT','IF','OR','ON','AN',
                     'ME','MY','DO','NO','SO','UP','CRISPR','CAS9','PUBMED','SHOW',
                     'PROTEIN','STRUCTURE','DEPMAP','DRUG','PAPER','LAB','FIELD',
                     'CELL','LINE','LINES','SCREEN','DESIGN','GUIDE','TOOL','TOOLS'}

        # Strategy 1: context-aware patterns
        for pat in [
            # "for gene kras", "check tp53", "analyze brca1"
            r'(?:gene|genes|analyze|check)\s+(?:gene\s+)?([A-Za-z][A-Za-z0-9]{1,7})\b',
            # "expression of kras", "mutation of tp53", "survival for brca1"
            r'(?:expression|mutation|dependency|survival|structure)\s+(?:of|for)\s+(?:gene\s+)?([A-Za-z][A-Za-z0-9]{1,7})\b',
            # "for kras" (but not "for the" etc — filter via NON_GENES)
            r'\bfor\s+([A-Za-z][A-Za-z0-9]{1,7})\b',
            # "kras expression", "tp53 in tcga", "brca1 mutation", "prnp protein"
            r'\b([A-Za-z][A-Za-z0-9]{1,7})\s+(?:expression|mutation|dependency|survival|analysis|protein|in\s+(?:tcga|gtex|depmap|normal|cancer))',
        ]:
            for m in _re.finditer(pat, text, _re.IGNORECASE):
                w = m.group(1).upper()
                if w not in NON_GENES and len(w) >= 2:
                    return w

        # Strategy 2: ALL-CAPS words (explicit gene mentions)
        for m in _re.finditer(r'\b([A-Z][A-Z0-9]{1,7})\b', text):
            w = m.group(1)
            if w not in NON_GENES:
                return w
        return None

    @staticmethod
    def _extract_gene_list_from_text(text):
        """Extract a multi-gene list from user text (pasted in chatbox).

        Detects patterns like:
          - One gene per line (newline-separated)
          - Space-separated gene symbols
          - Comma-separated gene symbols
          - Mixed separators

        Returns a list of cleaned gene symbols, or empty list if the text
        doesn't look like a gene list (< 3 gene-like tokens).
        """
        import re as _re
        if not text:
            return []

        # Words that look like genes but aren't
        NON_GENES = {
            'TCGA','GTEX','GEO','RNA','DNA','CNV','MAF','API','CSV','TSV',
            'PDF','PNG','SVG','URL','AND','FOR','THE','WITH','FROM','NOT',
            'USE','ALL','TOP','LOW','HIGH','SHOW','PLOT','RUN','GET','SET',
            'CHECK','EXPRESSION','NORMAL','TISSUE','CANCER','GENE','GENES',
            'ANALYSIS','DATA','DOWNLOAD','SEARCH','FIND','LOOK','QUERY',
            'ABOUT','PAN','MODE','HELP','WHAT','HOW','TRUE','FALSE','NONE',
            'DESIGN','LIBRARY','SGRNA','GUIDE','CRISPR','CAS9','HUMAN',
            'MOUSE','SPECIES','ENRICHMENT','PATHWAY','NETWORK','PLEASE',
            'BUBBLE','LIST','FILE','MAKE','CREATE','GENERATE','HELP','ME',
            'CAN','YOU','THIS','THAT','THESE','THOSE','THEM','HERE',
            'IN','OF','TO','BY','IS','AS','AT','BE','IT','IF','OR','ON','AN',
        }

        # Remove the instruction part (e.g. "help me design sgRNA library for
        # the gene list") — keep only the gene-list-looking block.
        # Split on common instruction phrases
        lines = text.strip().split('\n')

        # Collect candidate gene symbols
        candidates = []
        for line in lines:
            line = line.strip()
            if not line:
                continue
            # Split on whitespace, commas, semicolons, tabs
            tokens = _re.split(r'[\s,;\t]+', line)
            for tok in tokens:
                tok = tok.strip().strip('"').strip("'")
                if not tok:
                    continue
                # Gene symbol pattern: starts with letter, alphanumeric + hyphens/dots
                # Typical: TP53, BRCA1, C11orf88, RP11-382A20.3, 1-Mar (Excel-corrupted)
                if _re.match(r'^[A-Za-z0-9][A-Za-z0-9._-]{0,20}$', tok):
                    upper = tok.upper()
                    if upper not in NON_GENES and len(tok) >= 2:
                        candidates.append(tok)

        # Heuristic: if at least 3 tokens look gene-like, treat as gene list
        if len(candidates) >= 3:
            # Deduplicate while preserving order
            seen = set()
            unique = []
            for g in candidates:
                gu = g.upper()
                if gu not in seen:
                    seen.add(gu)
                    unique.append(g)
            return unique
        return []

    @staticmethod
    def _extract_cancer_type(text):
        """Extract TCGA cancer type code from text (case-insensitive)."""
        import re as _re
        if not text:
            return None
        TCGA_CODES = ['BRCA','LUAD','LUSC','COAD','READ','PRAD','KIRC','KIRP',
            'LIHC','STAD','BLCA','UCEC','GBM','LGG','HNSC','THCA','SKCM','CESC',
            'SARC','LAML','PAAD','ESCA','KICH','DLBC','MESO','UVM','ACC','PCPG',
            'TGCT','THYM','UCS','CHOL']
        # Also map common cancer names to codes
        CANCER_NAMES = {
            'colon': 'COAD', 'colorectal': 'COAD', 'colon cancer': 'COAD',
            'breast': 'BRCA', 'breast cancer': 'BRCA',
            'lung adenocarcinoma': 'LUAD', 'lung squamous': 'LUSC', 'lung cancer': 'LUAD',
            'prostate': 'PRAD', 'prostate cancer': 'PRAD',
            'kidney': 'KIRC', 'renal': 'KIRC',
            'liver': 'LIHC', 'hepatocellular': 'LIHC',
            'stomach': 'STAD', 'gastric': 'STAD',
            'bladder': 'BLCA', 'uterine': 'UCEC', 'endometrial': 'UCEC',
            'glioblastoma': 'GBM', 'glioma': 'LGG',
            'head and neck': 'HNSC', 'thyroid': 'THCA',
            'melanoma': 'SKCM', 'skin': 'SKCM',
            'cervical': 'CESC', 'sarcoma': 'SARC',
            'leukemia': 'LAML', 'pancreatic': 'PAAD', 'pancreas': 'PAAD',
            'esophageal': 'ESCA', 'mesothelioma': 'MESO',
            'ovarian': 'OV', 'ovary': 'OV',
        }
        # Try TCGA codes first (case-insensitive)
        for code in TCGA_CODES:
            if _re.search(r'\b' + code + r'\b', text, _re.IGNORECASE):
                return code
        # Try common cancer names
        text_lower = text.lower()
        for name, code in CANCER_NAMES.items():
            if name in text_lower:
                return code
        return None

    @staticmethod
    def _extract_uniprot_id(text):
        """Extract UniProt accession from text (e.g. P01116, Q9Y6K9, A0A0C5B5G6)."""
        import re as _re
        if not text:
            return None
        # UniProt accession patterns:
        # Classic: [OPQ][0-9][A-Z0-9]{3}[0-9] or [A-NR-Z][0-9][A-Z][A-Z0-9]{2}[0-9]
        # New 10-char: [A-NR-Z][0-9][A-Z][A-Z0-9]{2}[0-9][A-Z][A-Z0-9]{2}[0-9]
        # Also check for explicit "UniProt: X" or "uniprot X" patterns first
        explicit = _re.search(r'(?:uniprot|UniProt)\s*[:=]?\s*([A-Z0-9]{6,10})\b', text, _re.IGNORECASE)
        if explicit:
            return explicit.group(1).upper()
        # Standard 6-char accession
        m = _re.search(r'\b([OPQ][0-9][A-Z0-9]{3}[0-9])\b', text)
        if m:
            return m.group(1)
        m = _re.search(r'\b([A-NR-Z][0-9][A-Z][A-Z0-9]{2}[0-9])\b', text)
        if m:
            return m.group(1)
        # 10-char new format
        m = _re.search(r'\b([A-NR-Z][0-9][A-Z][A-Z0-9]{2}[0-9][A-Z][A-Z0-9]{2}[0-9])\b', text)
        if m:
            return m.group(1)
        return None

    @staticmethod
    def _extract_pdb_id(text):
        """Extract PDB ID from text (4-char alphanumeric starting with digit, e.g. 1A2B, 6LU7)."""
        import re as _re
        if not text:
            return None
        # Explicit "PDB: X" or "pdb-id X"
        explicit = _re.search(r'(?:pdb(?:-?id)?|PDB)\s*[:=]?\s*(\d[A-Za-z0-9]{3})\b', text, _re.IGNORECASE)
        if explicit:
            return explicit.group(1).upper()
        # Standalone PDB ID pattern (digit + 3 alphanumeric)
        m = _re.search(r'\b(\d[A-Za-z0-9]{3})\b', text)
        if m:
            candidate = m.group(1).upper()
            # Avoid matching pure numbers or common 4-char words
            if not candidate.isdigit() and not candidate.isalpha():
                return candidate
        return None

    @staticmethod
    def _extract_genome(text):
        """Extract genome assembly from text (e.g. hg38, hg19, mm10, GRCh38, GRCh37)."""
        import re as _re
        if not text:
            return None
        m = _re.search(r'\b(hg(?:19|38)|mm(?:9|10|39)|GRC[hm]\d{2})\b', text, _re.IGNORECASE)
        if m:
            return m.group(1)
        return None

    @staticmethod
    def _extract_species(text):
        """Extract species/organism from text."""
        import re as _re
        if not text:
            return None
        species_map = {
            'human': 'human', 'homo sapiens': 'human',
            'mouse': 'mouse', 'mus musculus': 'mouse',
            'rat': 'rat', 'rattus': 'rat',
            'zebrafish': 'zebrafish', 'danio rerio': 'zebrafish',
            'drosophila': 'drosophila', 'fruit fly': 'drosophila',
            'yeast': 'yeast', 'saccharomyces': 'yeast',
            'e. coli': 'ecoli', 'escherichia': 'ecoli',
            'arabidopsis': 'arabidopsis',
        }
        text_lower = text.lower()
        for name, val in species_map.items():
            if name in text_lower:
                return val
        return None

    @staticmethod
    def _extract_column_for_flag(flag_name, text):
        """Extract a CSV/table column name for plot-related flags.

        Handles patterns like:
          - "log2FoldChange 作为 X 轴", "padj 作为显著性指标"
          - "X axis: log2FoldChange", "p-value column: padj"
          - "use log2FC for X", "significance column = FDR"
          - Known stat-like column names: log2FoldChange, log2FC, LFC, padj,
            pvalue, p_value, FDR, qvalue, p_kat5_vs_aavs1, LFC_kat5_vs_aavs1
        """
        import re as _re
        if not text or not flag_name:
            return None
        fn = flag_name.lower()

        # Determine axis/role the flag refers to
        axis_keywords = []
        if fn in ('x-col', 'x_col', 'xcol', 'x-column', 'x_column', 'x-axis-col', 'x'):
            axis_keywords = ['x[\\s\\-_]?(?:axis|轴|col(?:umn)?)?', 'horizontal']
        elif fn in ('y-col', 'y_col', 'ycol', 'y-column', 'y_column', 'y-axis-col', 'y'):
            axis_keywords = ['y[\\s\\-_]?(?:axis|轴|col(?:umn)?)?', 'vertical']
        elif fn in ('p-col', 'p_col', 'pcol', 'p-column', 'p_column', 'pval-col', 'pvalue-col',
                    'significance-col', 'sig-col'):
            axis_keywords = [r'(?:p[\s\-_]?(?:value|val)?|padj|p\.?adj|fdr|qvalue|q[\s\-_]?value|'
                             r'significance|显著性(?:指标)?)']
        elif fn in ('feature-col', 'feature_col', 'gene-col', 'label-col', 'name-col'):
            axis_keywords = [r'(?:gene|feature|label|name)[\s\-_]?(?:col(?:umn)?|名)?']
        else:
            return None

        # Pattern 1: "VALUE 作为/as/for X轴/X axis/X column"
        # E.g. "log2FoldChange作为X轴", "padj as the significance"
        for kw in axis_keywords:
            for pat in [
                rf'([A-Za-z_][\w\.]*)\s*(?:作为|当作|用作|as(?:\s+the)?|for(?:\s+the)?|用于|是)\s*(?:the\s+)?{kw}',
                rf'{kw}[:：=\s]*(?:use|using|是|为|用|使用|column[:：=\s]*(?:is|=)?|列[:：=\s]*(?:是|为)?|col[:：=\s]*(?:is|=)?)?\s*([A-Za-z_][\w\.]*)',
                rf'{kw}[^\n\r]{{0,30}}?["\']?([A-Za-z_][\w\.]{{2,}})["\']?',
            ]:
                m = _re.search(pat, text, _re.IGNORECASE)
                if m:
                    val = m.group(1).strip().strip('"\'.,;:')
                    # Filter out generic words accidentally captured
                    bad = {'the', 'a', 'an', 'use', 'using', 'is', 'as', 'for', 'of',
                           'this', 'that', 'with', 'and', 'or', 'column', 'col',
                           'axis', 'value', 'significance', 'sig', 'pvalue', 'pval'}
                    if val.lower() not in bad and len(val) >= 3:
                        return val

        # Pattern 2: detect known stat-like column tokens in the text
        # (useful when user pastes column headers or LLM mentions them directly)
        stat_cols_x = [r'log2?FoldChange', r'log2FC', r'log2?_?FC', r'LFC(?:_\w+)?',
                       r'logfc', r'fold[_\-]?change', r'effect[_\-]?size']
        stat_cols_p = [r'padj', r'p\.?adj', r'p[_\-]?adjusted', r'q[_\-]?value',
                       r'qvalue', r'FDR', r'p[_\-]?value(?:_\w+)?', r'pvalue(?:_\w+)?',
                       r'p_\w+_vs_\w+']

        if 'x' in fn or 'fold' in fn or 'lfc' in fn:
            for sc in stat_cols_x:
                m = _re.search(rf'\b({sc}(?:_\w+)?)\b', text, _re.IGNORECASE)
                if m:
                    return m.group(1)
        if 'p' in fn or 'sig' in fn:
            for sc in stat_cols_p:
                m = _re.search(rf'\b({sc})\b', text, _re.IGNORECASE)
                if m:
                    return m.group(1)

        return None

    @staticmethod
    def _extract_generic_value_for_flag(flag_name, text):
        """Try to extract a value for a given flag name from free text.

        Searches for patterns like:
          --flag VALUE, flag: VALUE, flag = VALUE, flag VALUE
        in the text (from LLM reply or user message).
        """
        import re as _re
        if not text or not flag_name:
            return None
        # Normalize flag name for regex
        fn = flag_name.replace('-', r'[\-_]?')
        # Try: --flag VALUE or --flag=VALUE
        m = _re.search(rf'--{fn}\s*[=\s]\s*(\S+)', text, _re.IGNORECASE)
        if m:
            return m.group(1).strip('"\'')
        # Try: "Flag: VALUE" or "Flag = VALUE"
        fn_human = flag_name.replace('-', r'[\-_ ]?')
        m = _re.search(rf'{fn_human}\s*[:=]\s*(\S+)', text, _re.IGNORECASE)
        if m:
            val = m.group(1).strip('"\'.,;')
            if len(val) >= 2:
                return val
        return None

    @staticmethod
    def _resolve_missing_flags(script_path, flags, params, used_flags, cmd,
                                latest_user, llm_reply, conv_text):
        """Generic resolver: for each required flag not yet in the command,
        try to extract its value from text using smart patterns."""
        import re as _re

        # Known identifier extractors keyed by flag name patterns
        EXTRACTORS = {
            'uniprot':     'uniprot_id',
            'pdb-id':      'pdb_id',
            'pdb':         'pdb_id',
            'genome':      'genome',
            'assembly':    'genome',
            'species':     'species',
            'organism':    'species',
            # Plot column flags (dispatched to _extract_column_for_flag)
            'x-col':       'column',
            'x_col':       'column',
            'y-col':       'column',
            'y_col':       'column',
            'p-col':       'column',
            'p_col':       'column',
            'feature-col': 'column',
            'feature_col': 'column',
            'gene-col':    'column',
            'label-col':   'column',
        }

        required = set()
        try:
            src = script_path.read_text(encoding='utf-8', errors='replace')
            for m in _re.finditer(
                r"add_argument\s*\(\s*['\"](-{1,2}[\w-]+)['\"]"
                r"[^)]*required\s*=\s*True", src):
                required.add(m.group(1).lstrip('-'))
        except Exception:
            pass

        # Also check for mutually_exclusive_group required args
        # e.g. src = parser.add_mutually_exclusive_group(required=True)
        #      src.add_argument("--pdb-id", ...)
        try:
            src = script_path.read_text(encoding='utf-8', errors='replace')
            # Find: varname = xxx.add_mutually_exclusive_group(required=True)
            for mex in _re.finditer(
                r'(\w+)\s*=\s*\w+\.add_mutually_exclusive_group\s*\(\s*required\s*=\s*True\s*\)',
                src):
                var = mex.group(1)
                # Find all var.add_argument("--flag") after this line
                after = src[mex.end():mex.end() + 600]
                for am in _re.finditer(
                    rf'{_re.escape(var)}\.add_argument\s*\(\s*["\'](-{{1,2}}[\w-]+)["\']',
                    after):
                    required.add(am.group(1).lstrip('-'))
        except Exception:
            pass

        # Text sources in priority order
        texts = [latest_user, llm_reply, conv_text]

        added = []
        for rf in required:
            if rf in used_flags or f'--{rf}' in cmd:
                continue
            if rf in ('outdir', 'output', 'output-prefix', 'input', 'data', 'file',
                      'counts', 'gene-file', 'input-file', 'matrix', 'table'):
                continue  # handled elsewhere

            val = None

            # Check known extractors first
            extractor_key = None
            for ek, ext_name in EXTRACTORS.items():
                if ek == rf or rf.endswith(ek) or rf.startswith(ek):
                    extractor_key = ext_name
                    break

            if extractor_key == 'uniprot_id':
                for t in texts:
                    val = Handler._extract_uniprot_id(t)
                    if val:
                        break
            elif extractor_key == 'pdb_id':
                for t in texts:
                    val = Handler._extract_pdb_id(t)
                    if val:
                        break
            elif extractor_key == 'genome':
                for t in texts:
                    val = Handler._extract_genome(t)
                    if val:
                        break
            elif extractor_key == 'species':
                for t in texts:
                    val = Handler._extract_species(t)
                    if val:
                        break
            elif extractor_key == 'column':
                for t in texts:
                    val = Handler._extract_column_for_flag(rf, t)
                    if val:
                        break

            # Generic extraction: search for --flag value or "flag: value" in text
            if not val:
                for t in texts:
                    val = Handler._extract_generic_value_for_flag(rf, t)
                    if val:
                        break

            if val and rf in flags:
                cmd += [f'--{rf}', str(val)]
                used_flags.add(rf)
                added.append(f'--{rf} {val}')
                print(f"  [agent] Generic resolver: added --{rf} {val}")

        return added

    # ── Find the right script for a skill ───────────────────────────────────
    def _resolve_script(self, skill_dir, script_name, skill_path):
        """Find the Python script to run. Tries exact match first, then fuzzy, then auto-detect."""
        if script_name:
            # Try exact match in scripts/ subdir and skill dir
            for candidate in [
                skill_dir / 'scripts' / script_name,
                skill_dir / script_name,
            ]:
                if candidate.exists():
                    return candidate

            # Fuzzy match: SKILL.md might say "for" but file says "by", or vice versa
            # Try matching by stem similarity in both locations
            stem = Path(script_name).stem  # e.g. "normal_tissue_expression_for_gene"
            for search_dir in [skill_dir / 'scripts', skill_dir]:
                if not search_dir.exists():
                    continue
                for f in search_dir.glob('*.py'):
                    if '.ipynb_checkpoints' in str(f):
                        continue
                    # Check if most of the name matches (allow 1-2 word differences)
                    f_parts = set(f.stem.split('_'))
                    s_parts = set(stem.split('_'))
                    overlap = len(f_parts & s_parts)
                    if overlap >= max(len(f_parts), len(s_parts)) - 1 and overlap >= 3:
                        print(f"  [resolve] Fuzzy match: {script_name} → {f.name}")
                        return f
            # Fall through to auto-detect instead of returning None

        # Auto-detect: prefer scripts/ subdir, fall back to skill dir
        for search_dir in [skill_dir / 'scripts', skill_dir]:
            if not search_dir.exists():
                continue
            scripts = [f for f in search_dir.glob('*.py')
                       if '.ipynb_checkpoints' not in str(f)
                       and f.name != '__init__.py'
                       and 'downloader' not in f.name.lower()]
            if scripts:
                # Prefer script matching skill folder name
                folder_stem = skill_path.split('/')[-1].replace('-', '_')
                preferred = [s for s in scripts if folder_stem in s.stem]
                return (preferred or scripts)[0]
        return None

    # ── Run a skill script (intelligent agent) ──────────────────────────────
    def _run_skill(self, req):
        import re as _re
        skill_path  = req.get('skill_path', '')
        script_name = req.get('script_name')
        input_text  = req.get('input_text') or req.get('input_data')
        input_fname = req.get('input_filename', 'input.txt')
        extra_args  = req.get('args', [])
        params      = req.get('params', {})
        llm_args    = req.get('llm_args')  # LLM-generated CLI args from EXEC_ARGS block

        skill_dir = self.repo_root / skill_path
        if not skill_dir.exists():
            return {'error': f'Skill not found: {skill_path}', 'success': False}

        # ── Resolve script ──────────────────────────────────────────────────
        resolved = self._resolve_script(skill_dir, script_name, skill_path)
        if not resolved:
            return {'error': f'No Python scripts found in {skill_path}/', 'success': False}
        script_path = resolved

        # ── Auto-redirect to standalone scripts for co-expression/co-essentiality ──
        # The standalone scripts produce rich output (barplot, network, FDR table)
        # while the main script's --modules coexpression/coessentiality only gives
        # a basic top-N TSV. Always prefer standalone scripts.
        if 'depmap-analysis-for-gene' in skill_path:
            # Gather all text signals to detect co-expression / co-essentiality intent
            _all_text = ' '.join([
                req.get('latest_user_text', ''),
                req.get('llm_reply', ''),
                req.get('conversation_text', '') or '',
                str(llm_args or ''),
                str(params),
            ]).lower()
            _is_coexpr = any(kw in _all_text for kw in [
                'co-expression', 'coexpression', 'co-expressed', 'coexpressed',
                'correlated genes', 'co_expression',
            ])
            _is_coess = any(kw in _all_text for kw in [
                'co-essentiality', 'coessentiality', 'co-essential', 'coessential',
                'co-dependency', 'codependency', 'co_essentiality',
            ])
            # Only redirect if the current script is the main one (not already standalone)
            if _is_coexpr and 'coexpression' not in script_path.stem:
                alt = skill_dir / 'scripts' / 'depmap_coexpression.py'
                if alt.exists():
                    print(f"  [agent] ↗ Redirecting to standalone depmap_coexpression.py (richer output)")
                    script_path = alt
                    # Strip --modules from llm_args if present
                    if isinstance(llm_args, list):
                        new_args = []
                        skip = False
                        for a in llm_args:
                            if skip:
                                skip = False
                                continue
                            if a == '--modules':
                                skip = True
                                continue
                            new_args.append(a)
                        llm_args = new_args
            elif _is_coess and 'coessentiality' not in script_path.stem:
                alt = skill_dir / 'scripts' / 'depmap_coessentiality.py'
                if alt.exists():
                    print(f"  [agent] ↗ Redirecting to standalone depmap_coessentiality.py (richer output)")
                    script_path = alt
                    if isinstance(llm_args, list):
                        new_args = []
                        skip = False
                        for a in llm_args:
                            if skip:
                                skip = False
                                continue
                            if a == '--modules':
                                skip = True
                                continue
                            new_args.append(a)
                        llm_args = new_args

        # ── Determine build strategy ────────────────────────────────────────
        # If the LLM provided args (from reading SKILL.md), use them directly.
        # This is the "LLM-first" approach: the LLM already knows the correct
        # flags from SKILL.md usage examples — no need to re-parse on server.
        use_llm_args = isinstance(llm_args, list) and len(llm_args) > 0
        if use_llm_args:
            print(f"  [agent] ✅ Using LLM-provided args: {llm_args}")
        else:
            print(f"  [agent] ⚠️ No LLM args — falling back to server-side extraction")

        # ── Smart param extraction (4-tier priority) ──────────────────────
        # Priority 1: Latest user message (the user explicitly said it)
        # Priority 2: LLM reply (the LLM confirmed what it's analyzing)
        # Priority 3: Client-side params (from JS extractParams)
        # Priority 4: Full conversation text (older context)
        conv_text = req.get('conversation_text', '') or input_text or ''
        llm_reply = req.get('llm_reply', '')
        latest_user = req.get('latest_user_text', '')
        memory_ctx = req.get('memory_context', '')

        # Append memory context to conversation text for richer extraction
        if memory_ctx:
            conv_text = conv_text + '\n' + memory_ctx

        print(f"  [agent] latest_user_text: {latest_user[:100]!r}")
        print(f"  [agent] llm_reply: {llm_reply[:100]!r}")
        print(f"  [agent] client params: {params}")
        if memory_ctx:
            print(f"  [agent] memory context: {memory_ctx[:150]!r}")

        if not use_llm_args:
            # ── Priority 1: Extract from latest user message ────────────────
            if latest_user:
                if not params.get('gene'):
                    gene = self._extract_gene_from_text(latest_user)
                    if gene:
                        params['gene'] = gene
                        print(f"  [agent] Gene from user msg: {gene}")

                if not params.get('cancer_type'):
                    ct = self._extract_cancer_type(latest_user)
                    if ct:
                        params['cancer_type'] = ct
                        print(f"  [agent] Cancer type from user msg: {ct}")

                if not params.get('uniprot_id'):
                    uid = self._extract_uniprot_id(latest_user)
                    if uid:
                        params['uniprot_id'] = uid
                        print(f"  [agent] UniProt ID from user msg: {uid}")

                if not params.get('pdb_id'):
                    pid = self._extract_pdb_id(latest_user)
                    if pid:
                        params['pdb_id'] = pid
                        print(f"  [agent] PDB ID from user msg: {pid}")

                if not params.get('genome'):
                    gn = self._extract_genome(latest_user)
                    if gn:
                        params['genome'] = gn
                        print(f"  [agent] Genome from user msg: {gn}")

                if not params.get('species'):
                    sp = self._extract_species(latest_user)
                    if sp:
                        params['species'] = sp
                        print(f"  [agent] Species from user msg: {sp}")

            # ── Priority 2: Extract from LLM reply ──────────────────────────
            if llm_reply:
                if not params.get('gene'):
                    gene = self._extract_gene_from_text(llm_reply)
                    if gene:
                        params['gene'] = gene
                        print(f"  [agent] Gene from LLM reply: {gene}")

                if not params.get('cancer_type'):
                    ct = self._extract_cancer_type(llm_reply)
                    if ct:
                        params['cancer_type'] = ct
                        print(f"  [agent] Cancer type from LLM reply: {ct}")

                if not params.get('uniprot_id'):
                    uid = self._extract_uniprot_id(llm_reply)
                    if uid:
                        params['uniprot_id'] = uid
                        print(f"  [agent] UniProt ID from LLM reply: {uid}")

                if not params.get('pdb_id'):
                    pid = self._extract_pdb_id(llm_reply)
                    if pid:
                        params['pdb_id'] = pid
                        print(f"  [agent] PDB ID from LLM reply: {pid}")

                if not params.get('genome'):
                    gn = self._extract_genome(llm_reply)
                    if gn:
                        params['genome'] = gn
                        print(f"  [agent] Genome from LLM reply: {gn}")

                if not params.get('species'):
                    sp = self._extract_species(llm_reply)
                    if sp:
                        params['species'] = sp
                        print(f"  [agent] Species from LLM reply: {sp}")

            # ── Priority 3: Client already sent params (from extractParams) ──
            # (already in params dict from req)

            # ── Priority 4: Full conversation text (older context) ──────────
            if not params.get('gene'):
                gene = self._extract_gene_from_text(conv_text)
                if gene:
                    params['gene'] = gene
                    print(f"  [agent] Gene from conversation: {gene}")

            if not params.get('cancer_type'):
                ct = self._extract_cancer_type(conv_text)
                if ct:
                    params['cancer_type'] = ct
                    print(f"  [agent] Cancer type from conversation: {ct}")

            if not params.get('uniprot_id'):
                uid = self._extract_uniprot_id(conv_text)
                if uid:
                    params['uniprot_id'] = uid
                    print(f"  [agent] UniProt ID from conversation: {uid}")

            if not params.get('pdb_id'):
                pid = self._extract_pdb_id(conv_text)
                if pid:
                    params['pdb_id'] = pid
                    print(f"  [agent] PDB ID from conversation: {pid}")

            # Search query (for paper/lab search skills)
            if not params.get('query'):
                qm = _re.search(r'(?:search|find|look for|query|about)\s+(.{10,80}?)(?:\.|,|$)',
                                llm_reply or conv_text, _re.IGNORECASE)
                if qm:
                    params['query'] = qm.group(1).strip()

        # Inspect script flags (needed for output dir and preflight)
        flags = self._script_flags(script_path)
        print(f"  [agent] script flags: {sorted(flags)}")

        # ── Output directory ────────────────────────────────────────────────
        run_id  = datetime.now().strftime('%Y%m%d_%H%M%S_%f')
        out_dir = self.results_dir / run_id
        out_dir.mkdir(parents=True, exist_ok=True)

        # ── Extract modules from LLM args for smart preflight ───────────────
        # If LLM said --modules expression, we use that for preflight too
        if use_llm_args:
            for i, arg in enumerate(llm_args):
                if arg == '--modules' and i + 1 < len(llm_args):
                    params['modules'] = llm_args[i + 1]
                    print(f"  [agent] LLM specified modules: {params['modules']}")

        # ── Pre-flight: auto-download data if needed ────────────────────────
        preflight_flags = self._run_preflight(skill_path, out_dir, params)

        # ── Resolve input file ──────────────────────────────────────────────
        # Priority: uploaded_path (pre-uploaded to server) > inline input_data
        #         > gene list extracted from conversation text
        input_path = None
        uploaded_path = req.get('uploaded_path')
        if uploaded_path and Path(uploaded_path).exists():
            input_path = Path(uploaded_path)
            print(f"  [agent] Using pre-uploaded file: {input_path}")
        elif input_text:
            input_path = out_dir / input_fname
            input_path.write_text(str(input_text), encoding='utf-8')
            print(f"  [agent] Wrote inline file: {input_path} ({len(input_text)} chars)")
        else:
            # ── Auto-extract gene list from conversation text ──────────────
            # When the user pastes a gene list directly in the chatbox (not as
            # a file attachment), extract gene symbols and save to a temp file.
            # This enables skills like crispr-library-design, GO enrichment,
            # etc. to work with inline gene lists.
            extracted_genes = self._extract_gene_list_from_text(latest_user or conv_text)
            if extracted_genes and len(extracted_genes) >= 3:
                input_path = out_dir / 'gene_list.txt'
                input_path.write_text('\n'.join(extracted_genes), encoding='utf-8')
                print(f"  [agent] Extracted {len(extracted_genes)} genes from chat → {input_path}")

        # Detect which flag the script uses for input files
        INPUT_FLAG_CANDIDATES = ['input', 'genes', 'data', 'file', 'counts',
                                  'gene-file', 'gene-list', 'input-file',
                                  'matrix', 'table']
        input_flag = None
        if input_path:
            for cand in INPUT_FLAG_CANDIDATES:
                if cand in flags:
                    input_flag = cand
                    break
            if not input_flag:
                # Last resort: check if any flag name contains 'input' or 'file'
                for f in flags:
                    if 'input' in f or 'file' in f:
                        input_flag = f
                        break
            if input_flag:
                print(f"  [agent] Input flag: --{input_flag}")
            else:
                print(f"  [agent] WARNING: No input flag found in script, available: {sorted(flags)}")

        # ── Build command ───────────────────────────────────────────────────
        python = sys.executable
        cmd = [python, str(script_path)]

        if use_llm_args:
            # ── LLM-FIRST PATH: use LLM's args directly ────────────────────
            # The LLM read SKILL.md and constructed the correct flags.
            # We only add: output dir, input file, and preflight data paths.

            # Output dir (LLM is told not to include this)
            if 'outdir' in flags:
                cmd += ['--outdir', str(out_dir)]
            elif 'output' in flags:
                ext = self._script_output_ext(script_path)
                cmd += ['--output', str(out_dir / f'{script_path.stem}_result{ext}')]
            elif 'output-prefix' in flags:
                cmd += ['--output-prefix', str(out_dir / script_path.stem)]

            # Input file
            if input_path and input_flag:
                cmd += [f'--{input_flag}', str(input_path)]

            # Add the LLM's args (the core of the command)
            # But skip any --outdir/--output/--input the LLM may have accidentally included
            skip_flags = {'--outdir', '--output', '--output-prefix'}
            if input_flag:
                skip_flags.add(f'--{input_flag}')
            skip_next = False
            i = 0
            while i < len(llm_args):
                arg = llm_args[i]
                if skip_next:
                    skip_next = False
                    i += 1
                    continue
                if arg in skip_flags:
                    skip_next = True  # skip the flag and its value
                    i += 1
                    continue
                # ── Fix: if a value starts with '-', merge with flag using '='
                # to prevent argparse from treating it as another flag.
                # e.g.  --ylabel  -log10(p-value)  →  --ylabel=-log10(p-value)
                if (isinstance(arg, str) and arg.startswith('--')
                        and i + 1 < len(llm_args)
                        and isinstance(llm_args[i+1], str)
                        and llm_args[i+1].startswith('-')
                        and not llm_args[i+1].startswith('--')):
                    val = llm_args[i+1]
                    # Don't merge if the next item is itself a known flag
                    cmd.append(f'{arg}={val}')
                    i += 2
                    continue
                cmd.append(arg)
                i += 1

            # Add preflight data file paths (LLM doesn't know about cached files)
            # Only add if not already specified with a VALID path by LLM args.
            # The LLM often hallucinates paths (e.g. /Users/.../file.csv from
            # user's local machine) that don't exist on the server — override those.
            llm_flag_set = set()
            llm_flag_values = {}
            for idx, arg in enumerate(llm_args):
                if arg.startswith('--'):
                    f = arg.lstrip('-')
                    llm_flag_set.add(f)
                    if idx + 1 < len(llm_args) and not llm_args[idx + 1].startswith('--'):
                        llm_flag_values[f] = llm_args[idx + 1]
            for flag, path in preflight_flags.items():
                if flag not in llm_flag_set:
                    cmd += [f'--{flag}', path]
                else:
                    # LLM specified this flag — check if its value is a valid path
                    llm_val = llm_flag_values.get(flag, '')
                    if llm_val and not Path(llm_val).exists():
                        print(f"  [agent] LLM path for --{flag} does not exist: {llm_val}")
                        print(f"  [agent] Replacing with cached: {path}")
                        # Remove invalid path from cmd and replace
                        new_cmd = []
                        skip_next = False
                        for c in cmd:
                            if skip_next:
                                skip_next = False
                                continue
                            if c == f'--{flag}':
                                skip_next = True
                                continue
                            new_cmd.append(c)
                        cmd = new_cmd
                        cmd += [f'--{flag}', path]

            # ── Preemptive resolver: fill in any required flag the LLM forgot ──
            # e.g. LLM's prose mentions "log2FoldChange 作为 X 轴" but EXEC_ARGS
            # omitted --x-col. Let the resolver harvest values from memory/prose
            # BEFORE we hit argparse "required" errors.
            used_flags_pre = set()
            for i, a in enumerate(cmd):
                if isinstance(a, str) and a.startswith('--'):
                    used_flags_pre.add(a.split('=', 1)[0].lstrip('-'))
            try:
                pre_added = self._resolve_missing_flags(
                    script_path, flags, params, used_flags_pre, cmd,
                    latest_user, llm_reply, conv_text)
                if pre_added:
                    print(f"  [agent] Preemptive resolver added: {pre_added}")
            except Exception as e:
                print(f"  [agent] Preemptive resolver error: {e}")

        else:
            # ── FALLBACK PATH: server-side param mapping (old behavior) ─────
            # Input file
            if input_path and input_flag:
                cmd += [f'--{input_flag}', str(input_path)]

            # Output dir
            if 'outdir' in flags:
                cmd += ['--outdir', str(out_dir)]
            elif 'output' in flags:
                ext = self._script_output_ext(script_path)
                cmd += ['--output', str(out_dir / f'{script_path.stem}_result{ext}')]
            elif 'output-prefix' in flags:
                cmd += ['--output-prefix', str(out_dir / script_path.stem)]

            # Intent-aware module detection from conversation
            if 'modules' in flags and conv_text and not params.get('modules'):
                text_lower = conv_text.lower()
                detected_modules = []
                module_keywords = {
                    'expression':     ['expression', 'gene expression', 'rna', 'mrna', 'transcription'],
                    'mutation':       ['mutation', 'mutations', 'somatic', 'variant', 'variants'],
                    'copy_number':    ['copy number', 'cnv', 'amplification', 'deletion', 'copy-number'],
                    'essentiality':   ['essentiality', 'dependency', 'dependencies', 'crispr', 'gene effect', 'essential'],
                    'coexpression':   ['coexpression', 'co-expression', 'correlated genes'],
                    'coessentiality': ['coessentiality', 'co-essentiality', 'codependency'],
                }
                for mod, keywords in module_keywords.items():
                    if any(kw in text_lower for kw in keywords):
                        detected_modules.append(mod)
                if detected_modules and 'full' not in text_lower:
                    params['modules'] = ','.join(detected_modules)
                    print(f"  [agent] Detected intent → --modules {params['modules']}")

            # Map structured params → CLI flags
            param_map = {
                'gene':        ['gene'],
                'genes':       ['genes', 'gene'],
                'cancer_type': ['cancer-type', 'cancer-types'],
                'cancer':      ['cancer-type', 'cancer-types'],
                'query':       ['query', 'keyword', 'term'],
                'keyword':     ['query', 'keyword', 'term'],
                'modules':     ['modules', 'mode'],
                'uniprot_id':  ['uniprot'],
                'pdb_id':      ['pdb-id'],
                'genome':      ['genome', 'assembly'],
                'species':     ['species', 'organism'],
                'input_file':  [],
            }
            used_flags = set()
            for param_key, candidates in param_map.items():
                val = params.get(param_key)
                if not val:
                    continue
                for cand in candidates:
                    if cand in flags and cand not in used_flags:
                        cmd += [f'--{cand}', str(val)]
                        used_flags.add(cand)
                        break

            # Add pre-flight flags (downloaded data files)
            for flag, path in preflight_flags.items():
                if flag in flags and flag not in used_flags:
                    cmd += [f'--{flag}', path]
                    used_flags.add(flag)

            # Any raw extra args from client
            cmd += extra_args

            # ── Generic resolver: fill in any remaining required flags ──────
            # This catches script-specific flags like --uniprot, --pdb-id,
            # --chain, --genome, etc. that aren't in param_map
            resolved = self._resolve_missing_flags(
                script_path, flags, params, used_flags, cmd,
                latest_user, llm_reply, conv_text)
            if resolved:
                print(f"  [agent] Generic resolver added: {resolved}")

            # Safety check for missing required flags
            required_flags = self._script_required_flags(script_path)
            missing = []
            for rf in required_flags:
                if rf not in used_flags and f'--{rf}' not in cmd:
                    if rf in ('outdir', 'output', 'output-prefix'):
                        continue
                    missing.append(rf)
            if missing:
                print(f"  [agent] WARNING: Missing required flags: {missing}")
                print(f"  [agent] Params available: {params}")

        print(f"  [agent] Running: {' '.join(cmd[:10])}{'...' if len(cmd)>10 else ''}")

        # ── Auto-install dependencies if needed ─────────────────────────────
        self._ensure_deps(skill_dir)

        # ── Snapshot skill_dir files BEFORE running (to detect files saved to cwd)
        import time as _time
        run_start_time = _time.time()
        before_files = set()
        try:
            for f in skill_dir.iterdir():
                if f.is_file():
                    before_files.add(f.name)
        except Exception:
            pass

        # ── Smart Debug-and-Retry Loop ──────────────────────────────────────
        # Try executing; if it fails, analyze the error and apply auto-fixes,
        # then retry. Up to 4 retries with different fix strategies.
        debug_log = []   # human-readable debug history shown to user
        max_retries = 4
        attempt = 0
        result = None

        while attempt <= max_retries:
            attempt += 1
            try:
                _shared_dir = str(Path(__file__).resolve().parent / '_shared')
                _pypath = os.pathsep.join([str(skill_dir), _shared_dir])
                result = _subprocess.run(
                    cmd, capture_output=True, text=True, timeout=300,
                    cwd=str(skill_dir),
                    env={**os.environ, 'PYTHONPATH': _pypath})
            except _subprocess.TimeoutExpired:
                debug_log.append(f"⏱ Attempt {attempt}: timed out (>5 min)")
                # Try with shorter scope or skip retry
                return {'error': 'Script timed out (>5 min)', 'success': False,
                        'run_id': run_id, 'debug_log': debug_log,
                        'command': ' '.join(cmd)}
            except Exception as e:
                debug_log.append(f"💥 Attempt {attempt}: {e}")
                return {'error': str(e), 'success': False, 'run_id': run_id,
                        'debug_log': debug_log, 'command': ' '.join(cmd)}

            # Success → exit loop
            if result.returncode == 0:
                if attempt > 1:
                    debug_log.append(f"✅ Attempt {attempt}: succeeded after debugging")
                break

            # Failed → diagnose and try to fix
            print(f"  [debug] Attempt {attempt}/{max_retries+1} failed (rc={result.returncode})")
            stderr_tail = (result.stderr or '')[-2500:]
            stdout_tail = (result.stdout or '')[-1000:]

            fix = self._diagnose_and_fix(
                error_text=stderr_tail + '\n' + stdout_tail,
                cmd=cmd,
                script_path=script_path,
                skill_dir=skill_dir,
                params=params,
                latest_user=latest_user,
                llm_reply=llm_reply,
                conv_text=conv_text,
                attempt=attempt,
            )

            if not fix:
                debug_log.append(f"❌ Attempt {attempt}: no auto-fix available")
                print(f"  [debug] No auto-fix found, giving up")
                break

            debug_log.append(f"🔧 Attempt {attempt}: {fix['description']}")
            print(f"  [debug] Applying fix: {fix['description']}")

            # Apply the fix to cmd (mutates in place via fix function)
            if fix.get('new_cmd'):
                cmd = fix['new_cmd']
                print(f"  [debug] New command: {' '.join(cmd[:12])}{'...' if len(cmd)>12 else ''}")

            # If fix says "stop" (unfixable), break out
            if fix.get('stop'):
                break

        # ── Rescue files saved to skill_dir (cwd) instead of out_dir ────────
        # Some scripts ignore --output or have a default relative path.
        # Find any file in skill_dir that wasn't there before the run.
        try:
            import shutil as _shutil
            import re as _re_out
            OUTPUT_EXTS = {'.png', '.svg', '.pdf', '.jpg', '.jpeg', '.tsv',
                           '.csv', '.xlsx', '.xls', '.json', '.html',
                           '.txt', '.fasta', '.fa', '.bed', '.gff', '.gtf',
                           '.bam', '.vcf', '.log'}
            # Parse stdout for "Saved plot to: X", "Output: X", "Written to X"
            stdout_str = result.stdout or ''
            mentioned_files = set()
            for m in _re_out.finditer(
                r'(?:saved (?:plot |result |file |to |output )?(?:to\s+)?|'
                r'output(?:\s+(?:file|written)?\s*[:=]\s*)?|'
                r'written (?:to\s+)?|'
                r'generated(?:\s+[:=]?\s*))'
                r'["\']?([A-Za-z0-9_\-./]+\.[A-Za-z]{2,5})["\']?', stdout_str, _re_out.IGNORECASE):
                mentioned_files.add(m.group(1).strip())

            # Find new files in skill_dir created during the run
            for f in skill_dir.iterdir():
                if not f.is_file():
                    continue
                # New file (not present before run) OR modified during run
                try:
                    mtime = f.stat().st_mtime
                except Exception:
                    continue
                is_new = f.name not in before_files
                is_fresh = mtime >= run_start_time - 1
                # Only rescue if new/fresh AND has an output-like extension
                if (is_new or is_fresh) and f.suffix.lower() in OUTPUT_EXTS:
                    dest = out_dir / f.name
                    # Avoid overwriting a legit existing file in out_dir
                    if dest.exists():
                        continue
                    try:
                        _shutil.move(str(f), str(dest))
                        print(f"  [rescue] Moved {f.name} from skill_dir → out_dir")
                    except Exception as _e:
                        print(f"  [rescue] Failed to move {f.name}: {_e}")
        except Exception as _e:
            print(f"  [rescue] Error during file rescue: {_e}")

        # Collect output files recursively (skip the input file we wrote)
        output_files = []
        for f in sorted(out_dir.rglob('*')):
            if not f.is_file():
                continue
            if input_path and f.name == input_path.name:
                continue
            rel = f.relative_to(out_dir)
            ext = f.suffix.lstrip('.').upper()
            output_files.append({
                'name': str(rel),
                'size': f.stat().st_size,
                'url':  f'/api/results/{run_id}/{rel}',
                'ext':  ext,
                'type': ext,
            })

        return {
            'run_id':       run_id,
            'script':       script_path.name,
            'skill':        skill_path,
            'returncode':   result.returncode,
            'success':      result.returncode == 0,
            'stdout':       result.stdout[-8000:],
            'stderr':       result.stderr[-3000:],
            'output_files': output_files,
            'command':      ' '.join(cmd),
            'params_used':  params,
            'llm_args_used': use_llm_args,
            'debug_log':    debug_log,
            'attempts':     attempt,
        }

    # ── Smart Diagnostician: analyze script error and propose a fix ────────
    # Tracks fixes already tried per script so we don't loop forever.
    def _diagnose_and_fix(self, error_text, cmd, script_path, skill_dir,
                          params, latest_user, llm_reply, conv_text, attempt):
        """Analyze script failure and return a fix dict, or None if unfixable.

        Returns: {'description': str, 'new_cmd': list|None, 'stop': bool}
        """
        import re as _re
        import subprocess as _sp

        err = error_text or ''
        err_lower = err.lower()

        # Track what we've already tried this run (avoid loops)
        if not hasattr(self, '_debug_tried'):
            self._debug_tried = set()
        run_key = id(cmd)  # per-execution key (cmd list identity changes per fix)

        # ── Fix 1: Missing Python module → install it ───────────────────────
        if 'modulenotfounderror' in err_lower or 'no module named' in err_lower:
            mod = _re.search(r"No module named ['\"]([^'\"]+)['\"]", err)
            if mod:
                pkg = mod.group(1).split('.')[0]
                # Map common import name → pip package name
                PKG_MAP = {
                    'sklearn': 'scikit-learn',
                    'cv2': 'opencv-python',
                    'PIL': 'Pillow',
                    'yaml': 'pyyaml',
                    'Bio': 'biopython',
                    'pysam': 'pysam',
                    'pyBigWig': 'pyBigWig',
                }
                pkg_install = PKG_MAP.get(pkg, pkg)
                key = f'install:{pkg_install}'
                if key in self._debug_tried:
                    return {'description': f'already tried installing {pkg_install}, giving up', 'stop': True}
                self._debug_tried.add(key)

                print(f"  [debug] Installing missing package: {pkg_install}")
                pip = [sys.executable, '-m', 'pip', 'install', pkg_install]
                try:
                    _sp.run(pip, capture_output=True, text=True, timeout=180)
                except Exception:
                    pass
                return {'description': f'installed missing package "{pkg_install}"',
                        'new_cmd': cmd}

        # ── Fix 2: Missing required argument → extract from context ─────────
        # argparse prints: "error: the following arguments are required: --gene"
        # or: "error: argument --gene is required"
        missing_args = []
        for m in _re.finditer(r'arguments? are required:\s*([\-\w,\s/]+)', err):
            for a in m.group(1).split(','):
                # handle forms like --x-col/-x or -x/--x-col → prefer LONG form
                parts = [p.strip() for p in a.strip().split('/')]
                long_forms = [p for p in parts if p.startswith('--')]
                a = long_forms[0] if long_forms else (parts[-1] if parts else '')
                if a.startswith('-'):
                    missing_args.append(a.lstrip('-'))
        for m in _re.finditer(r"argument\s+(--?[\w-]+)\s+is required", err):
            missing_args.append(m.group(1).lstrip('-'))
        # Mutually exclusive: "one of the arguments --pdb-id --uniprot --pdb-file is required"
        mex = _re.search(r'one of the arguments\s+([\-\w\s]+)\s+is required', err)
        if mex:
            for a in mex.group(1).split():
                if a.startswith('-'):
                    missing_args.append(a.lstrip('-'))

        if missing_args:
            key = f'missing_args:{",".join(sorted(set(missing_args)))}'
            if key in self._debug_tried:
                return {'description': f'already tried fixing missing args {missing_args}, giving up', 'stop': True}
            self._debug_tried.add(key)

            # Try generic resolver to extract values from context
            flags = self._script_flags(script_path)
            used_flags = set()
            for i, a in enumerate(cmd):
                if isinstance(a, str) and a.startswith('--'):
                    used_flags.add(a.lstrip('-'))

            new_cmd = list(cmd)
            added = self._resolve_missing_flags(
                script_path, flags, params, used_flags, new_cmd,
                latest_user, llm_reply, conv_text)

            # If mutually exclusive group: only need ONE to be filled
            if added:
                return {'description': f'extracted missing args from context: {added}',
                        'new_cmd': new_cmd}
            else:
                return {'description': f'could not auto-fill missing args: {missing_args}',
                        'stop': True}

        # ── Fix 2b: "expected one argument" → value swallowed by next flag
        # Happens when the value starts with '-' (e.g. --ylabel -log10(p-value))
        # Fix by merging flag and value with '=': --ylabel=-log10(p-value)
        eoa = _re.search(r"argument\s+(--[\w-]+):\s*expected one argument", err)
        if eoa:
            flag = eoa.group(1)  # e.g. --ylabel
            key = f'eoa:{flag}'
            if key in self._debug_tried:
                return {'description': f'already tried merging {flag} value', 'stop': True}
            self._debug_tried.add(key)

            new_cmd = list(cmd)
            for i, a in enumerate(new_cmd):
                if a == flag and i + 1 < len(new_cmd):
                    nxt = new_cmd[i + 1]
                    # If next item starts with '-' but isn't a known long flag, merge
                    if (isinstance(nxt, str) and nxt.startswith('-')
                            and not nxt.startswith('--')):
                        new_cmd[i] = f'{flag}={nxt}'
                        del new_cmd[i + 1]
                        return {'description': f'merged {flag} value with =-syntax (was being misparsed as a flag)',
                                'new_cmd': new_cmd}
                    # If nothing after flag at all
                    if i + 1 >= len(new_cmd):
                        del new_cmd[i]
                        return {'description': f'removed dangling {flag} (no value)',
                                'new_cmd': new_cmd}
            return {'description': f'could not fix --{flag} expected one argument', 'stop': True}

        # ── Fix 2c: matplotlib/pillow "Format 'X' is not supported" ─────────
        # Script tried to save to an unsupported extension (e.g. .txt). Rewrite
        # the --output arg's extension to the first supported format listed in
        # the error, preferring .png.
        fmt_err = _re.search(
            r"Format\s+['\"]?(\w+)['\"]?\s+is\s+not\s+supported"
            r"\s*\(.*?formats?[:\s]*([^)]+)\)",
            err, _re.IGNORECASE)
        if fmt_err:
            bad_ext = '.' + fmt_err.group(1).lower()
            # Parse supported list: "eps, jpeg, jpg, pdf, pgf, png, ps, raw, rgba, svg, svgz, tif, tiff, webp"
            supported = [s.strip().lower() for s in _re.split(r'[,\s]+', fmt_err.group(2)) if s.strip()]
            # Preferred order
            preferred = ['png', 'svg', 'pdf', 'jpg', 'jpeg', 'tif', 'tiff', 'webp']
            pick = next((p for p in preferred if p in supported),
                        supported[0] if supported else 'png')
            new_ext = '.' + pick

            key = f'fmt:{bad_ext}->{new_ext}'
            if key in self._debug_tried:
                return {'description': f'already tried swapping {bad_ext}→{new_ext}', 'stop': True}
            self._debug_tried.add(key)

            new_cmd = list(cmd)
            replaced = None
            for i, a in enumerate(new_cmd):
                if not isinstance(a, str):
                    continue
                # Match either "--output=path.txt" or "--output" followed by "path.txt"
                if a.startswith('--output=') and a.lower().endswith(bad_ext):
                    new_cmd[i] = a[:-len(bad_ext)] + new_ext
                    replaced = (a, new_cmd[i])
                    break
                if a == '--output' and i + 1 < len(new_cmd):
                    v = new_cmd[i + 1]
                    if isinstance(v, str) and v.lower().endswith(bad_ext):
                        new_cmd[i + 1] = v[:-len(bad_ext)] + new_ext
                        replaced = (v, new_cmd[i + 1])
                        break
            if replaced:
                return {'description': f'rewrote --output extension {bad_ext} → {new_ext} (matplotlib does not support {bad_ext[1:]})',
                        'new_cmd': new_cmd}
            return {'description': f'format {bad_ext} not supported but could not locate --output to rewrite',
                    'stop': True}

        # ── Fix 3: Unrecognized argument → remove it ────────────────────────
        # argparse: "error: unrecognized arguments: --foo bar"
        unrec = _re.search(r'unrecognized arguments?:\s*(.+?)(?:\n|$)', err)
        if unrec:
            bad = unrec.group(1).strip().split()
            key = f'unrec:{":".join(bad)}'
            if key in self._debug_tried:
                return {'description': f'already tried removing {bad}, giving up', 'stop': True}
            self._debug_tried.add(key)

            new_cmd = list(cmd)
            i = 0
            removed = []
            while i < len(new_cmd):
                if new_cmd[i] in bad:
                    removed.append(new_cmd[i])
                    # Remove flag and its value (if next item isn't another flag)
                    del new_cmd[i]
                    if i < len(new_cmd) and not new_cmd[i].startswith('--'):
                        del new_cmd[i]
                else:
                    i += 1
            if removed:
                return {'description': f'removed unrecognized args: {removed}',
                        'new_cmd': new_cmd}

        # ── Fix 4: Invalid choice for argument → try synonyms ───────────────
        # argparse: "error: argument --modules: invalid choice: 'expr' (choose from 'expression', 'mutation', ...)"
        invc = _re.search(
            r"argument\s+(--?[\w-]+):\s*invalid choice:\s*['\"]([^'\"]+)['\"]"
            r"\s*\(choose from\s+([^)]+)\)", err)
        if invc:
            flag = invc.group(1).lstrip('-')
            bad_val = invc.group(2)
            choices_str = invc.group(3)
            choices = [c.strip().strip("'\"") for c in choices_str.split(',')]

            key = f'invchoice:{flag}:{bad_val}'
            if key in self._debug_tried:
                return {'description': f'already tried fixing {flag}, giving up', 'stop': True}
            self._debug_tried.add(key)

            # Try to find best matching choice
            bad_lower = bad_val.lower()
            best = None
            for c in choices:
                if c.lower() == bad_lower:
                    best = c; break
                if c.lower().startswith(bad_lower) or bad_lower.startswith(c.lower()):
                    best = c
            if not best:
                # Use first choice as fallback
                best = choices[0] if choices else None

            if best:
                new_cmd = list(cmd)
                for i, a in enumerate(new_cmd):
                    if a == f'--{flag}' and i + 1 < len(new_cmd):
                        new_cmd[i + 1] = best
                        return {'description': f'corrected --{flag} {bad_val} → {best}',
                                'new_cmd': new_cmd}

        # ── Fix 5: Input file not found → check alternate paths ─────────────
        # "FileNotFoundError: [Errno 2] No such file or directory: 'foo.tsv'"
        # or "Input file not found: foo.tsv"
        fnf = (_re.search(r"No such file or directory:\s*['\"]([^'\"]+)['\"]", err)
               or _re.search(r"Input file not found:\s*([^\s\n]+)", err)
               or _re.search(r"FileNotFoundError.*?['\"]([^'\"]+)['\"]", err))
        if fnf:
            missing_path = fnf.group(1)
            key = f'fnf:{missing_path}'
            if key in self._debug_tried:
                return {'description': f'already tried fixing missing file {missing_path}', 'stop': True}
            self._debug_tried.add(key)

            # Search for the file in common locations
            from pathlib import Path as _P
            search_dirs = [skill_dir, skill_dir.parent, self.results_dir,
                           self.results_dir / '_uploads', _P.cwd()]
            base = _P(missing_path).name
            for d in search_dirs:
                for cand in d.rglob(base):
                    if cand.is_file():
                        new_cmd = [str(cand) if a == missing_path else a for a in cmd]
                        return {'description': f'found file at {cand}',
                                'new_cmd': new_cmd}
            return {'description': f'could not locate missing file {missing_path}', 'stop': True}

        # ── Fix 6: numpy/pandas ABI mismatch → reinstall ────────────────────
        if ('_array_api not found' in err_lower
                or 'numpy.dtype size changed' in err_lower
                or 'incompatible with numpy' in err_lower):
            key = 'numpy_abi'
            if key in self._debug_tried:
                return {'description': 'already tried reinstalling numpy stack', 'stop': True}
            self._debug_tried.add(key)
            print(f"  [debug] numpy/pandas ABI mismatch — reinstalling…")
            pip = [sys.executable, '-m', 'pip', 'install', '--upgrade',
                   '--force-reinstall', '--no-deps', 'numpy', 'pandas',
                   'bottleneck', 'numexpr']
            try:
                _sp.run(pip, capture_output=True, text=True, timeout=300)
            except Exception:
                pass
            return {'description': 'reinstalled numpy/pandas/bottleneck/numexpr',
                    'new_cmd': cmd}

        # ── Fix 7: Permission denied → try chmod ────────────────────────────
        if 'permission denied' in err_lower:
            perm = _re.search(r"Permission denied:?\s*['\"]?([^'\"\n]+)", err)
            if perm:
                key = f'perm:{perm.group(1)}'
                if key in self._debug_tried:
                    return {'description': 'already tried chmod, giving up', 'stop': True}
                self._debug_tried.add(key)
                try:
                    import os as _os, stat as _stat
                    _os.chmod(perm.group(1), 0o755)
                    return {'description': f'fixed permissions on {perm.group(1)}',
                            'new_cmd': cmd}
                except Exception:
                    pass

        # ── Fix 8: Network/download error → wait and retry once ─────────────
        if any(s in err_lower for s in ['connection refused', 'connection reset',
                'timed out', 'temporary failure', 'name or service not known',
                'http error 5', 'urlopen error']):
            key = 'network_retry'
            if key in self._debug_tried:
                return {'description': 'network error already retried', 'stop': True}
            self._debug_tried.add(key)
            import time as _t
            _t.sleep(3)
            return {'description': 'network error — waited 3s and retrying',
                    'new_cmd': cmd}

        # ── Fix 9: Bad value type (int parsing, float parsing) → drop arg ───
        # "argument --top-pockets: invalid int value: 'five'"
        bv = _re.search(r"argument\s+(--?[\w-]+):\s*invalid\s+\w+\s+value:\s*['\"]([^'\"]+)['\"]", err)
        if bv:
            flag = bv.group(1).lstrip('-')
            key = f'badval:{flag}'
            if key in self._debug_tried:
                return {'description': f'already tried fixing --{flag}', 'stop': True}
            self._debug_tried.add(key)
            new_cmd = list(cmd)
            i = 0
            while i < len(new_cmd):
                if new_cmd[i] == f'--{flag}':
                    del new_cmd[i]
                    if i < len(new_cmd) and not new_cmd[i].startswith('--'):
                        del new_cmd[i]
                else:
                    i += 1
            return {'description': f'removed --{flag} (bad value, will use default)',
                    'new_cmd': new_cmd}

        # ── Unknown error → no fix ──────────────────────────────────────────
        # Print short diagnosis for the log
        first_err_line = ''
        for line in err.split('\n'):
            line = line.strip()
            if line and not line.startswith('Traceback'):
                first_err_line = line[:200]
                break
        print(f"  [debug] No matching fix pattern. Error: {first_err_line}")
        return None

    # ── Helpers ──────────────────────────────────────────────────────────────
    # ── Auto-install dependencies from requirements.txt ────────────────────
    _deps_installed = set()   # track which skill dirs we've already checked

    def _ensure_deps(self, skill_dir):
        """Install requirements.txt for a skill if not already done this session."""
        key = str(skill_dir)
        if key in self._deps_installed:
            return
        self._deps_installed.add(key)

        # Walk up to find requirements.txt (could be in skill dir or parent set dir)
        for check_dir in [skill_dir, skill_dir.parent]:
            req_file = check_dir / 'requirements.txt'
            if req_file.exists():
                pkgs = []
                for line in req_file.read_text().splitlines():
                    line = line.strip()
                    if line and not line.startswith('#'):
                        # Extract package name (before any version specifier)
                        pkg_name = line.split('>=')[0].split('<=')[0].split('==')[0].split('>')[0].split('<')[0].strip()
                        if pkg_name:
                            pkgs.append(pkg_name)
                if not pkgs:
                    return

                # Quick check: are any packages missing?
                missing = []
                for pkg in pkgs:
                    # Map common package names to importable names
                    import_name = pkg.replace('-', '_').lower()
                    # Special mappings
                    name_map = {
                        'scikit_learn': 'sklearn',
                        'pillow': 'PIL',
                        'py3dmol': 'py3Dmol',
                        'adjusttext': 'adjustText',
                        'umap_learn': 'umap',
                    }
                    import_name = name_map.get(import_name, import_name)
                    try:
                        __import__(import_name)
                    except ImportError:
                        missing.append(pkg)

                if missing:
                    pip = [sys.executable, '-m', 'pip']  # use same Python's pip
                    print(f"  [deps] Installing {len(missing)} missing packages: {', '.join(missing)}")
                    _subprocess.run(
                        pip + ['install', '--quiet'] + [line for line in req_file.read_text().splitlines()
                                                        if line.strip() and not line.strip().startswith('#')],
                        capture_output=True, text=True, timeout=120)
                    print(f"  [deps] Done")
                return

    def _cors_headers(self):
        self.send_header('Access-Control-Allow-Origin', '*')
        self.send_header('Access-Control-Allow-Methods', 'GET, POST, OPTIONS')
        self.send_header('Access-Control-Allow-Headers', 'Content-Type')

    def _send_html(self, html):
        data = html.encode('utf-8')
        self.send_response(200)
        self.send_header('Content-Type', 'text/html; charset=utf-8')
        self.send_header('Content-Length', str(len(data)))
        self._cors_headers()
        self.end_headers()
        self.wfile.write(data)

    def _send_json(self, obj, code=200):
        data = _json.dumps(obj, ensure_ascii=False).encode('utf-8')
        self.send_response(code)
        self.send_header('Content-Type', 'application/json')
        self.send_header('Content-Length', str(len(data)))
        self._cors_headers()
        self.end_headers()
        self.wfile.write(data)

    def _send_file(self, path):
        import mimetypes
        mime, _ = mimetypes.guess_type(str(path))
        mime = mime or 'application/octet-stream'
        data = path.read_bytes()
        # Serve images/SVG inline so they display in the browser; force download for others
        inline_types = {'image/png', 'image/jpeg', 'image/svg+xml', 'image/gif', 'image/webp', 'text/html'}
        disposition = 'inline' if mime in inline_types else f'attachment; filename="{path.name}"'
        self.send_response(200)
        self.send_header('Content-Type', mime)
        self.send_header('Content-Length', str(len(data)))
        self.send_header('Content-Disposition', disposition)
        self._cors_headers()
        self.end_headers()
        self.wfile.write(data)

    def log_message(self, fmt, *args):
        # Only log API calls, not page loads
        if '/api/' in (args[0] if args else ''):
            print(f"  [{datetime.now().strftime('%H:%M:%S')}] {args[0] if args else ''}")

def open_browser(port):
    time.sleep(1.0)
    url = f"http://localhost:{port}"
    print(f"\n  🌐  Opening {url} in your browser…")
    webbrowser.open(url)

def main():
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument('--repo', default='.', help='bioinfor-claw folder (default: current dir)')
    ap.add_argument('--port', default=7860, type=int, help='Port to serve on (default: 7860)')
    ap.add_argument('--host', default='localhost',
                    help='Bind address. Use 0.0.0.0 to expose on LAN / behind a tunnel. '
                         'Default "localhost" = local-only.')
    ap.add_argument('--no-browser', action='store_true', help='Do not open browser automatically')
    args = ap.parse_args()

    repo = Path(args.repo).expanduser().resolve()
    if not repo.is_dir():
        print(f"❌  Folder not found: {repo}"); sys.exit(1)

    # ── Python health check ─────────────────────────────────────────────────
    print(f"🐍  Python: {sys.executable} ({sys.version.split()[0]})")
    try:
        import pandas, numpy
        print(f"  ✅  pandas {pandas.__version__}, numpy {numpy.__version__}")
    except ImportError as e:
        print(f"  ⚠️   {e} — some skills may fail. Run: {sys.executable} -m pip install pandas numpy")
    except Exception as e:
        print(f"  ⚠️   pandas/numpy issue: {e}")
        print(f"      Fix: {sys.executable} -m pip install --upgrade pandas numpy bottleneck")

    tree_ui, flat = scan(repo)
    if not flat:
        print("❌  No SKILL.md files found."); sys.exit(1)

    print("🔨  Building app…")
    # Prefer serving bioinfor-claw.html (has skills pre-embedded + improved UI)
    static_html = repo / "bioinfor-claw.html"
    if static_html.exists():
        html = static_html.read_text(encoding="utf-8")
        print(f"  📄  Serving existing bioinfor-claw.html ({len(html)//1024} KB)")
    else:
        html = build_html(tree_ui, flat)
        print(f"  🏗️   Built HTML in memory ({len(html)//1024} KB)")

    # ── Inject fresh skill data from disk scan ──────────────────────────────
    # This replaces any stale bundled JSON so new/renamed skills are picked up
    # automatically on every server restart — no manual rebuild needed.
    ts = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    n  = len(flat)
    fresh_bundle = (
        f'\n<script id="live-skill-bundle">\n'
        f'// {n} skills live-injected at {ts}\n'
        f'window.BUNDLED_SKILL_TREE   = {json.dumps(tree_ui, ensure_ascii=False)};\n'
        f'window.BUNDLED_LOADED_SKILLS = {json.dumps(flat, ensure_ascii=False)};\n'
        f'window.SKILLS_BUNDLE_META   = {{\n'
        f'  generated: "{ts}", totalSkills: {n}, totalSets: {len(tree_ui)},\n'
        f'  totalChars: {sum(len(v) for v in flat.values())}, bundled: true,\n'
        f'  source: "live-server-inject"\n'
        f'}};\n'
        f'</script>\n'
    )
    # Inject server origin + fresh bundle before </head>
    inject_tag = (
        f'\n<script id="server-inject">'
        f'window.__SERVER_PORT={args.port};'
        f'window.__SERVER_INJECTED=true;'
        f'</script>\n'
        + fresh_bundle
    )
    html = html.replace('</head>', inject_tag + '</head>', 1)
    print(f"  🔄  Injected {n} skills from live disk scan")
    Handler.html_content = html

    port = args.port
    results_dir = repo / "web_results"
    results_dir.mkdir(exist_ok=True)
    Handler.repo_root   = repo
    Handler.results_dir = results_dir
    # Allow quick restart without "Address already in use" errors
    HTTPServer.allow_reuse_address = True
    server = HTTPServer((args.host, port), Handler)
    if args.host not in ('localhost', '127.0.0.1'):
        print(f"\n⚠️  Server is bound to {args.host}:{port} — reachable beyond this machine.")
        print("    Make sure you have auth (Tailscale / Cloudflare Access / basic auth) in front.\n")

    print(f"""
{'═'*58}
  ✅  bioinfor-claw is ready!

  📦  {len(flat)} skills embedded from:
      {repo}

  🌐  Web UI:   http://localhost:{port}
  ⚡  API:      http://localhost:{port}/api/health
  📋  API docs: http://localhost:{port}/docs (if using FastAPI)

  The web UI and execution server run on the SAME port.
  No need to run server.py separately.

  Press Ctrl+C to stop.
{'═'*58}
""")

    if not args.no_browser:
        threading.Thread(target=open_browser, args=(port,), daemon=True).start()

    try:
        server.serve_forever()
    except KeyboardInterrupt:
        print("\n\n  👋  Server stopped. Goodbye!\n")

if __name__ == '__main__':
    main()
