#!/usr/bin/env bash
# =============================================================================
# bioinfor-claw  ·  One-click setup script
# =============================================================================
# Usage:
#   bash setup.sh              # install core dependencies only
#   bash setup.sh --all        # install ALL skill dependencies
#   bash setup.sh --skill rnaseq-differential-expression
#   bash setup.sh --list       # list available skills
# =============================================================================

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PYTHON=${PYTHON:-python3}
VENV_DIR="$REPO_ROOT/.venv"

# ── Colours ─────────────────────────────────────────────────────────────────
GREEN='\033[0;32m'; YELLOW='\033[1;33m'; CYAN='\033[0;36m'
RED='\033[0;31m'; NC='\033[0m'; BOLD='\033[1m'

info()    { echo -e "${CYAN}[INFO]${NC}  $*"; }
success() { echo -e "${GREEN}[OK]${NC}    $*"; }
warn()    { echo -e "${YELLOW}[WARN]${NC}  $*"; }
error()   { echo -e "${RED}[ERROR]${NC} $*"; exit 1; }

# ── Banner ───────────────────────────────────────────────────────────────────
echo -e "${BOLD}"
echo "  ╔══════════════════════════════════════════╗"
echo "  ║         bioinfor-claw  setup             ║"
echo "  ║  Modular bioinformatics skill library    ║"
echo "  ╚══════════════════════════════════════════╝"
echo -e "${NC}"

# ── Parse arguments ──────────────────────────────────────────────────────────
MODE="core"
TARGET_SKILL=""

while [[ $# -gt 0 ]]; do
    case "$1" in
        --all)    MODE="all"; shift ;;
        --skill)  MODE="skill"; TARGET_SKILL="$2"; shift 2 ;;
        --list)   MODE="list"; shift ;;
        --help|-h)
            echo "Usage: bash setup.sh [--all | --skill <name> | --list]"
            echo ""
            echo "  (no args)          Install core dependencies only"
            echo "  --all              Install dependencies for every skill"
            echo "  --skill <name>     Install dependencies for one skill"
            echo "  --list             List all available skills"
            exit 0 ;;
        *) warn "Unknown argument: $1"; shift ;;
    esac
done

# ── Check Python ─────────────────────────────────────────────────────────────
if ! command -v "$PYTHON" &>/dev/null; then
    error "Python 3 not found. Install from https://www.python.org or via conda."
fi

PY_VERSION=$("$PYTHON" -c "import sys; print(f'{sys.version_info.major}.{sys.version_info.minor}')")
info "Python version: $PY_VERSION"

if [[ "${PY_VERSION%%.*}" -lt 3 ]] || \
   ([[ "${PY_VERSION%%.*}" -eq 3 ]] && [[ "${PY_VERSION##*.}" -lt 9 ]]); then
    error "Python 3.9 or higher is required (found $PY_VERSION)."
fi

# ── List skills mode ─────────────────────────────────────────────────────────
list_skills() {
    echo -e "\n${BOLD}Available skills:${NC}"
    find "$REPO_ROOT" -name "requirements.txt" \
        ! -path "*/.venv/*" \
        ! -path "*/example_results/*" | sort | while read -r req; do
        skill_dir="$(dirname "$req")"
        skill_name="$(basename "$skill_dir")"
        skill_set="$(basename "$(dirname "$skill_dir")")"
        printf "  ${CYAN}%-40s${NC} %s\n" "$skill_name" "$skill_set"
    done
    echo ""
}

if [[ "$MODE" == "list" ]]; then
    list_skills
    exit 0
fi

# ── Create virtual environment ────────────────────────────────────────────────
_create_venv() {
    info "Creating virtual environment at .venv ..."
    if "$PYTHON" -m venv "$VENV_DIR" 2>/dev/null; then
        success "Virtual environment created."
    else
        # Fallback: some systems lack ensurepip (e.g. Debian/Ubuntu without python3-venv)
        warn "Standard venv creation failed; trying --without-pip fallback ..."
        "$PYTHON" -m venv --without-pip "$VENV_DIR" || \
            error "Could not create virtual environment. Try: sudo apt install python3-venv"
        success "Virtual environment created (pip will be bootstrapped)."
    fi
}

if [[ ! -d "$VENV_DIR" ]]; then
    _create_venv
elif [[ ! -f "$VENV_DIR/bin/activate" ]]; then
    warn ".venv exists but looks incomplete — recreating ..."
    rm -rf "$VENV_DIR" 2>/dev/null || true
    _create_venv
else
    info "Virtual environment already exists at .venv"
fi

# Activate
# shellcheck disable=SC1091
source "$VENV_DIR/bin/activate"
PIP="$VENV_DIR/bin/pip"

# Bootstrap pip if absent (happens with --without-pip fallback)
if [[ ! -f "$PIP" ]]; then
    info "Bootstrapping pip ..."
    curl -sSL https://bootstrap.pypa.io/get-pip.py | "$PYTHON" - --quiet 2>/dev/null || \
        "$PYTHON" -m ensurepip --upgrade 2>/dev/null || \
        error "Could not install pip. Please install it manually: https://pip.pypa.io/en/stable/installation/"
    success "pip bootstrapped."
fi

# Upgrade pip silently
"$PIP" install --upgrade pip --quiet

# ── Core dependencies (shared across most skills) ─────────────────────────────
install_core() {
    info "Installing core dependencies ..."
    "$PIP" install --quiet \
        numpy>=1.24.0 \
        pandas>=2.0.0 \
        matplotlib>=3.7.0 \
        requests>=2.28.0 \
        scipy>=1.10.0
    success "Core dependencies installed."
}

# ── Install from a requirements.txt ──────────────────────────────────────────
install_requirements() {
    local req_file="$1"
    local skill_label="$2"
    info "Installing: $skill_label"
    # Filter out comments and blank lines; skip optional markers
    grep -v '^\s*#' "$req_file" \
        | grep -v '^\s*$' \
        | grep -v '^\s*#.*optional' \
        | "$PIP" install --quiet -r /dev/stdin 2>/dev/null || \
        warn "Some packages in $skill_label may not have installed cleanly."
}

# ── Install all skills ────────────────────────────────────────────────────────
install_all() {
    info "Installing dependencies for all skills (this may take a few minutes) ..."
    local failed=()
    while IFS= read -r req; do
        skill_name="$(basename "$(dirname "$req")")"
        install_requirements "$req" "$skill_name" || failed+=("$skill_name")
    done < <(find "$REPO_ROOT" -name "requirements.txt" \
                  ! -path "*/.venv/*" \
                  ! -path "*/example_results/*" | sort)

    if [[ ${#failed[@]} -gt 0 ]]; then
        warn "The following skills had install issues (often optional packages like xgboost / pydeseq2):"
        for s in "${failed[@]}"; do echo "    - $s"; done
    fi
    success "All skill dependencies processed."
}

# ── Install single skill ──────────────────────────────────────────────────────
install_skill() {
    local name="$1"
    local req
    req=$(find "$REPO_ROOT" -name "requirements.txt" \
              ! -path "*/.venv/*" \
              ! -path "*/example_results/*" \
          | grep "/$name/" | head -1 || true)

    if [[ -z "$req" ]]; then
        error "Skill '$name' not found. Run 'bash setup.sh --list' to see available skills."
    fi
    install_requirements "$req" "$name"
    success "Skill '$name' ready."
}

# ── Run ───────────────────────────────────────────────────────────────────────
install_core

case "$MODE" in
    core)  ;;   # core already done
    all)   install_all ;;
    skill) install_skill "$TARGET_SKILL" ;;
esac

# ── Claude Code integration hint ──────────────────────────────────────────────
echo ""
echo -e "${BOLD}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${NC}"
echo -e "${BOLD}  Setup complete!${NC}"
echo ""
echo -e "  ${CYAN}Activate the environment:${NC}"
echo -e "    source .venv/bin/activate"
echo ""
echo -e "  ${CYAN}Use with Claude Code:${NC}"
echo -e "    cd bioinfor-claw"
echo -e "    claude           # Claude Code reads SKILL.md files automatically"
echo ""
echo -e "  ${CYAN}Use with OpenClaw:${NC}"
echo -e "    Set BIOINFOR_CLAW_ROOT=$(pwd)"
echo -e "    Point your agent config to this directory"
echo ""
echo -e "  ${CYAN}Run a skill directly:${NC}"
echo -e "    python gene-centered-analysis/tcge_survival_for_gene/scripts/tcga_survival_for_gene.py \\"
echo -e "      --gene TP53 --cancer-type BRCA --mode os --outdir results/"
echo -e "${BOLD}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${NC}"
