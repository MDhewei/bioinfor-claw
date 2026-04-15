#!/usr/bin/env bash
# =============================================================================
#  install-openclaw.sh — Install bioinfor-claw skills permanently into OpenClaw
#
#  OpenClaw is an open-source AI agent: https://openclaw.ai
#  Skills live in ~/.openclaw/skills/ (global) or <workspace>/skills/ (local).
#  This installer clones the repo and registers all skills via extraDirs in
#  ~/.openclaw/openclaw.json so they are discovered automatically on every run.
#
#  One-liner from GitHub (no clone required):
#    bash <(curl -sSL https://raw.githubusercontent.com/MDhewei/bioinfor-claw/main/install-openclaw.sh)
#
#  Or from inside a cloned repo:
#    bash install-openclaw.sh [OPTIONS]
# =============================================================================
set -euo pipefail

# ── Colours ───────────────────────────────────────────────────────────────────
BOLD="\033[1m"; CYAN="\033[0;36m"; GREEN="\033[0;32m"
YELLOW="\033[0;33m"; RED="\033[0;31m"; RESET="\033[0m"

info()    { echo -e "${CYAN}[INFO]${RESET}  $*"; }
success() { echo -e "${GREEN}[OK]${RESET}    $*"; }
warn()    { echo -e "${YELLOW}[WARN]${RESET}  $*"; }
error()   { echo -e "${RED}[ERROR]${RESET} $*" >&2; exit 1; }
step()    { echo -e "\n${BOLD}▶ $*${RESET}"; }

# ── Banner ────────────────────────────────────────────────────────────────────
echo -e "${BOLD}
  ╔══════════════════════════════════════════════════╗
  ║    bioinfor-claw  →  OpenClaw installer          ║
  ║  Registers all 17 skills permanently             ║
  ╚══════════════════════════════════════════════════╝
${RESET}"

# ── Defaults ──────────────────────────────────────────────────────────────────
REPO_URL="https://github.com/MDhewei/bioinfor-claw.git"
INSTALL_DIR="${HOME}/.bioinfor-claw"          # permanent repo home
OPENCLAW_CONFIG="${HOME}/.openclaw/openclaw.json"  # OpenClaw config file
PYTHON="${PYTHON:-python3}"
MODE="extradirs"   # extradirs | copy  (copy places skills in ~/.openclaw/skills/)

# ── Argument parsing ──────────────────────────────────────────────────────────
while [[ $# -gt 0 ]]; do
    case "$1" in
        --install-dir)  INSTALL_DIR="$2"; shift 2 ;;
        --config)       OPENCLAW_CONFIG="$2"; shift 2 ;;
        --copy)         MODE="copy"; shift ;;
        --help|-h)
            echo "Usage: bash install-openclaw.sh [OPTIONS]"
            echo ""
            echo "Options:"
            echo "  --install-dir PATH   Where to clone the repo (default: ~/.bioinfor-claw)"
            echo "  --config PATH        Path to openclaw.json  (default: ~/.openclaw/openclaw.json)"
            echo "  --copy               Copy skills into ~/.openclaw/skills/ instead of using extraDirs"
            echo "  --help               Show this help"
            echo ""
            echo "Default method (extraDirs):"
            echo "  Adds skill parent directories to skills.load.extraDirs in openclaw.json."
            echo "  Skills stay in the cloned repo and are auto-updated with 'git pull'."
            echo ""
            echo "Copy method (--copy):"
            echo "  Physically copies each skill folder into ~/.openclaw/skills/."
            echo "  Re-run the installer to pull in updates."
            echo ""
            echo "One-liner from GitHub:"
            echo "  bash <(curl -sSL https://raw.githubusercontent.com/MDhewei/bioinfor-claw/main/install-openclaw.sh)"
            exit 0 ;;
        *)  error "Unknown argument: $1. Use --help for usage." ;;
    esac
done

# ── Clone or update bioinfor-claw ─────────────────────────────────────────────
step "Setting up bioinfor-claw repository at ${INSTALL_DIR}"
if [[ -d "$INSTALL_DIR/.git" ]]; then
    info "Repository already exists — pulling latest changes ..."
    if git -C "$INSTALL_DIR" pull --quiet --ff-only 2>/dev/null; then
        success "Repository up to date."
    else
        warn "Could not pull (offline or network error). Using existing version."
    fi
else
    info "Cloning ${REPO_URL} ..."
    git clone --quiet "$REPO_URL" "$INSTALL_DIR" || \
        error "git clone failed. Check your internet connection and that git is installed."
    success "Repository cloned to ${INSTALL_DIR}"
fi

# ── Python virtual environment ────────────────────────────────────────────────
step "Setting up Python environment"
VENV_DIR="${INSTALL_DIR}/.venv"

if [[ ! -f "${VENV_DIR}/bin/activate" ]]; then
    info "Creating virtual environment ..."
    if "$PYTHON" -m venv "$VENV_DIR" 2>/dev/null; then
        success "Virtual environment created."
    else
        warn "Standard venv failed — trying --without-pip fallback ..."
        "$PYTHON" -m venv --without-pip "$VENV_DIR" || \
            error "Could not create venv. Try: sudo apt install python3-venv"
    fi
else
    info "Virtual environment already exists."
fi

# shellcheck disable=SC1091
source "${VENV_DIR}/bin/activate"
PIP="${VENV_DIR}/bin/pip"

# Bootstrap pip if missing (happens with --without-pip fallback)
if [[ ! -f "$PIP" ]]; then
    info "Bootstrapping pip ..."
    curl -sSL https://bootstrap.pypa.io/get-pip.py | "$PYTHON" - --quiet 2>/dev/null || \
        error "Could not install pip. See https://pip.pypa.io/en/stable/installation/"
fi
"$PIP" install --upgrade pip --quiet

info "Installing core dependencies ..."
if "$PIP" install --quiet numpy pandas matplotlib requests scipy biopython 2>/dev/null; then
    success "Core dependencies installed."
else
    warn "Could not install core dependencies (network error?). Run manually later:"
    warn "  source ${VENV_DIR}/bin/activate && pip install numpy pandas matplotlib requests scipy biopython"
fi

# ── Discover all skills and their parent directories ──────────────────────────
step "Discovering skills"
declare -a SKILL_DIRS=()
declare -A PARENT_DIRS=()   # unique parent dirs for extraDirs config

while IFS= read -r skill_md; do
    skill_dir=$(dirname "$skill_md")
    [[ "$skill_dir" == *"example_results"* ]] && continue
    SKILL_DIRS+=("$skill_dir")
    parent_dir=$(dirname "$skill_dir")
    PARENT_DIRS["$parent_dir"]=1
done < <(find "$INSTALL_DIR" -name "SKILL.md" | sort)

info "Found ${#SKILL_DIRS[@]} skills across ${#PARENT_DIRS[@]} skill-set directories."

# Install per-skill Python dependencies
for skill_dir in "${SKILL_DIRS[@]}"; do
    req="${skill_dir}/requirements.txt"
    if [[ -f "$req" ]]; then
        if "$PIP" install --quiet -r "$req" 2>/dev/null; then
            info "  Dependencies installed: $(basename "$skill_dir")"
        else
            warn "  Could not install deps for $(basename "$skill_dir")"
        fi
    fi
done

# ── Register skills in OpenClaw ───────────────────────────────────────────────
step "Registering skills in OpenClaw (mode: ${MODE})"

if [[ "$MODE" == "extradirs" ]]; then
    # ── Method 1: extraDirs in openclaw.json ─────────────────────────────────
    # Adds each skill-set parent directory to skills.load.extraDirs.
    # OpenClaw discovers all SKILL.md folders inside those dirs automatically.
    # Lowest precedence — won't override bundled or user skills.

    OPENCLAW_DIR=$(dirname "$OPENCLAW_CONFIG")
    mkdir -p "$OPENCLAW_DIR"

    # Read existing config or start with empty object
    if [[ -f "$OPENCLAW_CONFIG" ]]; then
        existing_json=$(cat "$OPENCLAW_CONFIG")
    else
        existing_json="{}"
    fi

    # Build the list of directories to add
    declare -a NEW_DIRS=()
    for dir in "${!PARENT_DIRS[@]}"; do
        NEW_DIRS+=("$dir")
    done

    # Use Python to merge the new dirs into the existing JSON safely
    "$PYTHON" - <<PYEOF
import json, sys

config_path = "$OPENCLAW_CONFIG"
new_dirs    = [$(printf '"%s",' "${NEW_DIRS[@]}" | sed 's/,$//')]

try:
    with open(config_path) as f:
        cfg = json.load(f)
except (FileNotFoundError, json.JSONDecodeError):
    cfg = {}

cfg.setdefault("skills", {}).setdefault("load", {}).setdefault("extraDirs", [])
existing = set(cfg["skills"]["load"]["extraDirs"])
added = []
for d in new_dirs:
    if d not in existing:
        cfg["skills"]["load"]["extraDirs"].append(d)
        added.append(d)

with open(config_path, "w") as f:
    json.dump(cfg, f, indent=2)

if added:
    print(f"  Added {len(added)} director{'y' if len(added)==1 else 'ies'} to extraDirs:")
    for d in added:
        print(f"    {d}")
else:
    print("  All directories already present in extraDirs.")
PYEOF

    success "openclaw.json updated: ${OPENCLAW_CONFIG}"
    info "  OpenClaw will discover skills from these directories at startup."

else
    # ── Method 2: Copy into ~/.openclaw/skills/ ───────────────────────────────
    SKILLS_TARGET="${HOME}/.openclaw/skills"
    mkdir -p "$SKILLS_TARGET"

    INSTALLED=0; UPDATED=0; SKIPPED=0
    for skill_dir in "${SKILL_DIRS[@]}"; do
        skill_name=$(basename "$skill_dir")
        target="${SKILLS_TARGET}/${skill_name}"

        if [[ -d "$target" ]]; then
            info "  Updating: ${skill_name}"
            rm -rf "$target"
            cp -r "$skill_dir" "$target"
            (( UPDATED++ )) || true
        else
            cp -r "$skill_dir" "$target"
            success "  Installed: ${skill_name}"
            (( INSTALLED++ )) || true
        fi
    done

    success "Skills copied to ${SKILLS_TARGET}"
    info "  Installed: ${INSTALLED}  Updated: ${UPDATED}  Skipped: ${SKIPPED}"
fi

# ── Write env helper ──────────────────────────────────────────────────────────
VENV_PYTHON="${VENV_DIR}/bin/python"
ENV_FILE="${INSTALL_DIR}/.env"
cat > "$ENV_FILE" <<EOF
# bioinfor-claw runtime environment
# Source before running skill scripts manually: source ${ENV_FILE}
export BIOINFOR_CLAW_ROOT="${INSTALL_DIR}"
export BIOINFOR_CLAW_PYTHON="${VENV_PYTHON}"
export PATH="${VENV_DIR}/bin:\$PATH"
EOF

# ── Shell profile ─────────────────────────────────────────────────────────────
step "Shell profile integration"
PROFILE_LINE="source \"${ENV_FILE}\"  # bioinfor-claw"
PROFILE=""
for f in "${HOME}/.zshrc" "${HOME}/.bashrc" "${HOME}/.bash_profile"; do
    if [[ -f "$f" ]]; then PROFILE="$f"; break; fi
done

if [[ -n "$PROFILE" ]]; then
    if grep -q "bioinfor-claw" "$PROFILE" 2>/dev/null; then
        info "Shell profile already configured — skipping."
    else
        printf "\n%s\n" "$PROFILE_LINE" >> "$PROFILE"
        success "Added environment to ${PROFILE}"
    fi
else
    warn "No shell profile detected. Add this line to your .zshrc or .bashrc manually:"
    echo "  ${PROFILE_LINE}"
fi

# ── Verify (if openclaw CLI is available) ─────────────────────────────────────
step "Verifying installation"
if command -v openclaw &>/dev/null; then
    info "Running: openclaw skills list --eligible"
    openclaw skills list --eligible 2>/dev/null | grep -E "$(IFS=\|; echo "${SKILL_DIRS[*]##*/}")" \
        && success "Skills confirmed visible to OpenClaw." \
        || warn "Skills not yet visible — try restarting OpenClaw."
else
    warn "openclaw CLI not found in PATH. Install OpenClaw from https://openclaw.ai"
    warn "Once installed, run: openclaw skills list  (to verify skills are registered)"
fi

# ── Summary ───────────────────────────────────────────────────────────────────
echo ""
echo -e "${BOLD}${GREEN}════════════════════════════════════════════════${RESET}"
echo -e "${BOLD}  Installation complete!${RESET}"
echo -e "${GREEN}════════════════════════════════════════════════${RESET}"
echo ""
if [[ "$MODE" == "extradirs" ]]; then
    echo -e "  Skills registered via ${BOLD}skills.load.extraDirs${RESET} in:"
    echo -e "    ${OPENCLAW_CONFIG}"
    echo ""
    echo -e "  ${BOLD}To update skills:${RESET}"
    echo -e "    git -C ${INSTALL_DIR} pull"
    echo -e "  (No reinstall needed — extraDirs picks up changes automatically)"
else
    echo -e "  Skills copied to ${BOLD}~/.openclaw/skills/${RESET}"
    echo ""
    echo -e "  ${BOLD}To update skills:${RESET}"
    echo -e "    bash <(curl -sSL https://raw.githubusercontent.com/MDhewei/bioinfor-claw/main/install-openclaw.sh) --copy"
fi
echo ""
echo -e "  Restart OpenClaw if it is already running, then ask anything"
echo -e "  bioinformatics-related to invoke a skill automatically."
echo ""
echo -e "  ${BOLD}Repo:${RESET}   ${INSTALL_DIR}"
echo -e "  ${BOLD}Config:${RESET} ${OPENCLAW_CONFIG}"
echo ""
