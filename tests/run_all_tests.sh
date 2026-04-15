#!/bin/bash
#
# Bioinfor-Claw Test Runner
# Runs all skill test scripts and produces a pass/fail summary
#
# Usage:
#   ./tests/run_all_tests.sh                      # Run all tests
#   ./tests/run_all_tests.sh --skill-set <name>  # Run single skill set
#   ./tests/run_all_tests.sh --verbose            # Verbose output
#   ./tests/run_all_tests.sh --dry-run            # Show what would be tested

set -e

# Get the repo root (parent of tests directory)
REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
TESTS_DIR="$REPO_ROOT/tests"
RESULTS_DIR="$TESTS_DIR/results"

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
BOLD='\033[1m'
NC='\033[0m' # No Color

# Default options
SKILL_SET=""
VERBOSE=0
DRY_RUN=0
RUN_API=0

# Parse arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        --skill-set)
            SKILL_SET="$2"
            shift 2
            ;;
        --verbose)
            VERBOSE=1
            shift
            ;;
        --dry-run)
            DRY_RUN=1
            shift
            ;;
        --run-api)
            RUN_API=1
            shift
            ;;
        --help)
            echo "Bioinfor-Claw Test Runner"
            echo ""
            echo "Usage: $0 [OPTIONS]"
            echo ""
            echo "Options:"
            echo "  --skill-set <name>  Run tests for a specific skill set"
            echo "  --verbose           Enable verbose output"
            echo "  --dry-run           Show what would be tested without running"
            echo "  --run-api           Include API-requiring skills"
            echo "  --help              Show this help message"
            exit 0
            ;;
        *)
            echo "Unknown option: $1"
            exit 1
            ;;
    esac
done

# Create results directory
mkdir -p "$RESULTS_DIR"

# Build Python command
PYTHON_CMD="python3 $TESTS_DIR/test_skills.py"

if [ -n "$SKILL_SET" ]; then
    PYTHON_CMD="$PYTHON_CMD --skill-set $SKILL_SET"
fi

if [ $VERBOSE -eq 1 ]; then
    PYTHON_CMD="$PYTHON_CMD --verbose"
fi

if [ $DRY_RUN -eq 1 ]; then
    PYTHON_CMD="$PYTHON_CMD --dry-run"
fi

if [ $RUN_API -eq 1 ]; then
    PYTHON_CMD="$PYTHON_CMD --run-api"
fi

# Run the test suite
echo -e "${BOLD}Bioinfor-Claw Test Suite${NC}"
echo "===================================="
echo ""

cd "$REPO_ROOT"
$PYTHON_CMD
TEST_EXIT_CODE=$?

echo ""
echo "===================================="
if [ $TEST_EXIT_CODE -eq 0 ]; then
    echo -e "${GREEN}All tests passed!${NC}"
    exit 0
else
    echo -e "${RED}Some tests failed!${NC}"
    exit 1
fi
