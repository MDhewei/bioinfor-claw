#!/bin/bash
# ============================================================
# bioinfor-claw Local Setup Script
# Run this once from your bioinfor-claw folder:
#   cd /Users/whe3/Documents/bioinfor-claw
#   bash setup_local.sh
# ============================================================

set -e
REPO_DIR="$(cd "$(dirname "$0")" && pwd)"
echo ""
echo "╔══════════════════════════════════════════════════╗"
echo "║     bioinfor-claw  Local Setup                  ║"
echo "╚══════════════════════════════════════════════════╝"
echo ""
echo "📁 Repo folder: $REPO_DIR"

# ── 1. Download latest files if not present ─────────────────
echo ""
echo "── Step 1: Checking required files ─────────────────────"
for f in bioinfor-claw.html build_skills_bundle.py server.py; do
    if [ -f "$REPO_DIR/$f" ]; then
        echo "  ✅  $f found"
    else
        echo "  ❌  $f NOT found — copy it from Downloads first"
        echo "      cp ~/Downloads/$f $REPO_DIR/"
        MISSING=1
    fi
done
[ -n "$MISSING" ] && echo "" && echo "Please copy missing files and re-run." && exit 1

# ── 2. Build the skills bundle ───────────────────────────────
echo ""
echo "── Step 2: Building skills bundle ──────────────────────"
cd "$REPO_DIR"
python3 build_skills_bundle.py
echo "  ✅  skills-bundle.js generated"

# ── 3. Set up web folder ─────────────────────────────────────
echo ""
echo "── Step 3: Setting up web folder ───────────────────────"
mkdir -p "$REPO_DIR/web"
cp "$REPO_DIR/bioinfor-claw.html"  "$REPO_DIR/web/index.html"
cp "$REPO_DIR/skills-bundle.js"    "$REPO_DIR/web/skills-bundle.js"
echo "  ✅  web/index.html"
echo "  ✅  web/skills-bundle.js"

# ── 4. Install Python server deps ───────────────────────────
echo ""
echo "── Step 4: Installing server dependencies ───────────────"
pip3 install fastapi uvicorn python-multipart 2>&1 | tail -3
echo "  ✅  Server dependencies installed"

# ── 5. Done ──────────────────────────────────────────────────
echo ""
echo "╔══════════════════════════════════════════════════╗"
echo "║  ✅  Setup complete!                             ║"
echo "║                                                  ║"
echo "║  Start the server:                               ║"
echo "║    python3 server.py                             ║"
echo "║                                                  ║"
echo "║  Then open in your browser:                      ║"
echo "║    http://localhost:8000                         ║"
echo "╚══════════════════════════════════════════════════╝"
echo ""
