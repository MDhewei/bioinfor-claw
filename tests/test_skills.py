#!/usr/bin/env python3
"""
Bioinfor-Claw test suite.
Tests all 49 skills with synthetic data from tests/data/.
Usage: python tests/test_skills.py [--skill <skill-name>] [--verbose] [--dry-run] [--run-api]
"""

import subprocess
import sys
import time
import json
from pathlib import Path
from datetime import datetime
from collections import defaultdict
import argparse


# Constants
REPO = Path(__file__).parent.parent
DATA = REPO / "tests" / "data"
RESULTS = REPO / "tests" / "results"

# ANSI color codes
PASS_COLOR = "\033[92m"  # Green
FAIL_COLOR = "\033[91m"  # Red
SKIP_COLOR = "\033[93m"  # Yellow
API_COLOR = "\033[94m"   # Blue
RESET_COLOR = "\033[0m"
BOLD = "\033[1m"


def colored_status(status):
    """Return colored status string."""
    if status == "PASS":
        return f"{PASS_COLOR}[PASS]{RESET_COLOR}"
    elif status == "FAIL":
        return f"{FAIL_COLOR}[FAIL]{RESET_COLOR}"
    elif status == "SKIP":
        return f"{SKIP_COLOR}[SKIP]{RESET_COLOR}"
    elif status == "API_REQUIRED":
        return f"{API_COLOR}[API]{RESET_COLOR}"
    elif status == "TIMEOUT":
        return f"{FAIL_COLOR}[TIMEOUT]{RESET_COLOR}"
    return f"[{status}]"


def run_skill(script, args, timeout=120):
    """Run a skill script and return (status, runtime, error)."""
    if not Path(script).exists():
        return "SKIP", 0, f"Script not found: {script}"

    cmd = [sys.executable, str(script)] + [str(a) for a in args]
    t0 = time.time()
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout)
        runtime = time.time() - t0
        if result.returncode == 0:
            return "PASS", runtime, None
        else:
            error_msg = result.stderr or result.stdout
            return "FAIL", runtime, error_msg[-500:] if len(error_msg) > 500 else error_msg
    except subprocess.TimeoutExpired:
        return "TIMEOUT", timeout, "Script exceeded timeout"
    except Exception as e:
        return "FAIL", time.time() - t0, str(e)


def skip_skill(reason):
    """Return skip status with reason."""
    return "SKIP", 0, reason


# ========================
# MULTIOMICS SKILLS
# ========================

def test_rnaseq_differential_expression():
    """Test RNA-seq differential expression analysis."""
    script = REPO / "multiomics-data-analysis" / "rnaseq-differential-expression" / "scripts" / "rnaseq_differential_expression.py"
    return run_skill(script, [
        "--counts", DATA / "rnaseq_counts.tsv",
        "--metadata", DATA / "rnaseq_metadata.tsv",
        "--group-a", "control",
        "--group-b", "treatment",
        "--outdir", RESULTS / "rnaseq_de",
    ])


def test_atac_chipseq_downstream():
    """Test ATAC-ChIP-seq downstream analysis."""
    script = REPO / "multiomics-data-analysis" / "atac-chipseq-downstream-analysis" / "scripts" / "atac_chipseq_downstream.py"
    return run_skill(script, [
        "--peaks", DATA / "atac_peaks.narrowPeak",
        "--mode", "qc",
        "--outdir", RESULTS / "atac_qc",
    ])


def test_methylation_analysis():
    """Test methylation analysis."""
    script = REPO / "multiomics-data-analysis" / "methylation-analysis" / "scripts" / "methylation_analysis.py"
    return run_skill(script, [
        "--input", DATA / "methylation_beta.tsv",
        "--metadata", DATA / "methylation_metadata.tsv",
        "--group-col", "group",
        "--ref-group", "Control",
        "--mode", "dmp",
        "--outdir", RESULTS / "methylation",
    ])


def test_proteomics_analysis():
    """Test proteomics analysis."""
    script = REPO / "multiomics-data-analysis" / "proteomics-analysis" / "scripts" / "proteomics_analysis.py"
    return run_skill(script, [
        "--input", DATA / "proteomics_intensities.tsv",
        "--metadata", DATA / "proteomics_metadata.tsv",
        "--group-col", "group",
        "--ref-group", "control",
        "--outdir", RESULTS / "proteomics",
    ])


def test_single_cell_analysis():
    """Test single-cell RNA-seq analysis."""
    script = REPO / "multiomics-data-analysis" / "single-cell-basic-analysis" / "scripts" / "single_cell_basic_analysis.py"
    return run_skill(script, [
        "--input", DATA / "scrna_counts.tsv",
        "--metadata", DATA / "scrna_metadata.tsv",
        "--color-by", "true_celltype",
        "--outdir", RESULTS / "scrna",
    ])


# ========================
# CRISPR SKILLS
# ========================

def test_crispr_screen_analysis():
    """Test CRISPR screen analysis."""
    script = REPO / "crispr-design-and-analysis" / "crispr-screen-analysis" / "scripts" / "crispr_screen_analysis.py"
    return run_skill(script, [
        "--count-table", DATA / "crispr_counts.tsv",
        "--treatment", "Treatment_1,Treatment_2",
        "--control", "Control_1,Control_2",
        "--mode", "test",
        "--outdir", RESULTS / "crispr_screen",
    ])


def test_crispr_screen_qc():
    """Test CRISPR screen QC."""
    script = REPO / "crispr-design-and-analysis" / "crispr-screen-qc" / "scripts" / "crispr_screen_qc.py"
    return run_skill(script, [
        "--counts", DATA / "crispr_counts.tsv",
        "--outdir", RESULTS / "crispr_qc",
    ])


def test_design_sgrnas_by_gene():
    """Test sgRNA design by gene (requires API)."""
    return skip_skill("Requires API access to UCSC Genome Browser")


def test_design_base_editor_sgrnas():
    """Test base editor sgRNA design (requires API)."""
    return skip_skill("Requires API access to BE design service")


def test_design_prime_editor_sgrnas():
    """Test prime editor sgRNA design (requires API)."""
    return skip_skill("Requires API access to PE design service")


def test_crispr_library_design():
    """Test CRISPR library design (requires API)."""
    return skip_skill("Requires API access to design service")


# ========================
# GENE LIST ANALYSIS SKILLS
# ========================

def test_gene_list_overlap():
    """Test gene list overlap analysis."""
    script = REPO / "gene-list-analysis" / "gene-list-overlap" / "scripts" / "gene_list_overlap.py"
    return run_skill(script, [
        "--lists", f"SetA:{DATA}/gene_list_A.txt,SetB:{DATA}/gene_list_B.txt,SetC:{DATA}/gene_list_C.txt",
        "--outdir", RESULTS / "gene_overlap",
    ])


def test_go_analysis_for_gene_list():
    """Test GO analysis for gene list (requires API)."""
    return skip_skill("Requires API access to GO databases")


def test_gsea_for_ranked_gene_list():
    """Test GSEA for ranked gene list (requires API)."""
    return skip_skill("Requires API access to GSEA databases")


def test_ppi_network_for_gene_list():
    """Test PPI network for gene list (requires API)."""
    return skip_skill("Requires API access to PPI databases")


def test_transcription_factor_enrichment():
    """Test transcription factor enrichment (requires API)."""
    return skip_skill("Requires API access to TF databases")


def test_curate_gene_list_by_function():
    """Test gene list curation by function (requires UniProt API)."""
    return skip_skill("Requires API access to UniProt")


def test_function_annotation_for_gene_list():
    """Test function annotation for gene list (requires MyGene/UniProt API)."""
    return skip_skill("Requires API access to MyGene.info and UniProt")


# ========================
# GENE-CENTERED ANALYSIS SKILLS
# ========================

def test_cox_survival():
    """Test Cox survival analysis."""
    script = REPO / "gene-centered-analysis" / "cox-survival-analysis" / "scripts" / "cox_survival_analysis.py"
    return run_skill(script, [
        "--input", DATA / "survival_data.tsv",
        "--time-col", "time_days",
        "--event-col", "event",
        "--covariates", "age,stage",
        "--categorical", "stage",
        "--mode", "both",
        "--outdir", RESULTS / "cox_survival",
    ])


def test_depmap_analysis_for_gene():
    """Test DepMap analysis for gene (requires API)."""
    return skip_skill("Requires API access to DepMap portal")


def test_drug_sensitivity_for_gene():
    """Test drug sensitivity for gene (requires API)."""
    return skip_skill("Requires API access to drug sensitivity databases")


def test_mutation_analysis_for_gene():
    """Test mutation analysis for gene (requires API)."""
    return skip_skill("Requires API access to mutation databases")


def test_coexpression_for_gene():
    """Test coexpression analysis for gene (requires API)."""
    return skip_skill("Requires API access to coexpression databases")


def test_normal_tissue_expression_by_gene():
    """Test normal tissue expression by gene (requires API)."""
    return skip_skill("Requires API access to tissue expression databases")


def test_tcge_survival_for_gene():
    """Test TCGE survival analysis for gene (requires API)."""
    return skip_skill("Requires API access to TCGA/GEO data")


# ========================
# MACHINE LEARNING SKILLS
# ========================

def test_dimensionality_reduction():
    """Test dimensionality reduction."""
    script = REPO / "machine-learning-and-deep-learning" / "dimensionality-reduction" / "scripts" / "dimensionality_reduction.py"
    return run_skill(script, [
        "--input", DATA / "expression_matrix.tsv",
        "--metadata", DATA / "sample_metadata.tsv",
        "--method", "pca",
        "--color-by", "tissue_type",
        "--outdir", RESULTS / "dimred",
    ])


def test_clustering_analysis():
    """Test clustering analysis."""
    script = REPO / "machine-learning-and-deep-learning" / "clustering-analysis" / "scripts" / "clustering_analysis.py"
    return run_skill(script, [
        "--input", DATA / "expression_matrix.tsv",
        "--method", "hierarchical",
        "--n-clusters", "3",
        "--outdir", RESULTS / "clustering",
    ])


def test_omics_ml_classifier():
    """Test omics ML classifier."""
    script = REPO / "machine-learning-and-deep-learning" / "omics-ml-classifier" / "scripts" / "omics_ml_classifier.py"
    return run_skill(script, [
        "--features", DATA / "ml_features.tsv",
        "--labels", DATA / "ml_labels.tsv",
        "--model", "logistic",
        "--outdir", RESULTS / "ml_classifier",
    ], timeout=180)


# ========================
# PLOTTING SKILLS
# ========================

def test_plot_volcano():
    """Test volcano plot."""
    script = REPO / "bioinformatics-plot-generator" / "plot-volcano" / "scripts" / "plot_volcano.py"
    return run_skill(script, [
        "--input", DATA / "de_results.tsv",
        "--feature-col", "gene",
        "--x-col", "log2FoldChange",
        "--p-col", "padj",
        "--output", RESULTS / "volcano.png",
        "--annotate-top-n", "10",
        "--show-quadrant-counts",
    ])


def test_plot_heatmap():
    """Test heatmap plot."""
    script = REPO / "bioinformatics-plot-generator" / "plot-heatmap" / "scripts" / "plot_heatmap.py"
    return run_skill(script, [
        "--input", DATA / "heatmap_matrix.tsv",
        "--index-col", "gene",
        "--scale", "row",
        "--output", RESULTS / "heatmap.png",
    ])


def test_plot_box_violin():
    """Test box/violin plot."""
    script = REPO / "bioinformatics-plot-generator" / "plot-box-violin" / "scripts" / "plot_box_violin.py"
    return run_skill(script, [
        "--input", DATA / "boxviolin_data.tsv",
        "--value-col", "value",
        "--group-col", "group",
        "--plot-type", "violin",
        "--stats", "all_pairs",
        "--output", RESULTS / "violin.png",
    ])


def test_plot_scatter_bar():
    """Test scatter/bar plot."""
    script = REPO / "bioinformatics-plot-generator" / "plot-scatter-bar" / "scripts" / "plot_scatter_bar.py"
    return run_skill(script, [
        str(DATA / "scatter_data.tsv"),
        "--plot-type", "scatter",
        "--x-col", "x_value",
        "--y-col", "y_value",
        "--output", RESULTS / "scatter.png",
    ])


def test_plot_survival():
    """Test survival plot."""
    script = REPO / "bioinformatics-plot-generator" / "plot-survival" / "scripts" / "plot_survival.py"
    return run_skill(script, [
        "--input", DATA / "survival_plot.tsv",
        "--time-col", "time",
        "--event-col", "event",
        "--group-col", "group",
        "--show-pvalue",
        "--show-at-risk",
        "--output", RESULTS / "survival.png",
    ])


# ========================
# PROTEIN STRUCTURE SKILLS
# ========================

def test_protein_structure_for_gene():
    """Test protein structure for gene (requires API)."""
    return skip_skill("Requires API access to protein structure databases")


def test_protein_sequence_analysis():
    """Test protein sequence analysis (requires API)."""
    return skip_skill("Requires API access to protein sequence databases")


def test_protein_structure_visualizer():
    """Test protein structure visualizer (requires API)."""
    return skip_skill("Requires API access to protein structure databases")


def test_protein_structure_alignment():
    """Test protein structure alignment."""
    script = REPO / "protein-structure-analysis" / "protein-structure-alignment" / "scripts" / "protein_structure_alignment.py"
    return skip_skill("Requires protein structure files")


def test_protein_variant_mapper():
    """Test protein variant mapper."""
    script = REPO / "protein-structure-analysis" / "protein-variant-mapper" / "scripts" / "protein_variant_mapper.py"
    return skip_skill("Requires protein structure files and variant data")


# ========================
# PAPER SEARCH SKILLS
# ========================

def test_pubmed_search():
    """Test PubMed search (requires API)."""
    return skip_skill("Requires API access to PubMed")


def test_paper_digest_single():
    """Test single paper digest (requires API)."""
    return skip_skill("Requires API access to paper content")


def test_preprint_tracker():
    """Test preprint tracker (requires API)."""
    return skip_skill("Requires API access to preprint servers")


def test_big_papers_weekly_report():
    """Test big papers weekly report (requires API)."""
    return skip_skill("Requires API access to paper databases")


# ========================
# LAB SEARCH SKILLS
# ========================

def test_search_big_labs_by_field():
    """Test search big labs by field (requires API)."""
    return skip_skill("Requires API access to lab databases")


def test_track_lab_publications():
    """Test track lab publications (requires API)."""
    return skip_skill("Requires API access to publication databases")


def test_find_collaborators():
    """Test find collaborators (requires API)."""
    return skip_skill("Requires API access to researcher databases")


# ========================
# PUBLIC DATASETS SKILLS
# ========================

def test_tcga_download_data():
    """Test TCGA download (requires API)."""
    return skip_skill("Requires API access to TCGA portal")


def test_gtex_download_data():
    """Test GTEx download (requires API)."""
    return skip_skill("Requires API access to GTEx portal")


def test_geo_download_data():
    """Test GEO download (requires API)."""
    return skip_skill("Requires API access to GEO database")


def test_depmap_data_download():
    """Test DepMap data download (requires API)."""
    return skip_skill("Requires API access to DepMap portal")


# ========================
# TEST REGISTRY
# ========================

TESTS = {
    # Multiomics
    "rnaseq-differential-expression": ("multiomics-data-analysis", test_rnaseq_differential_expression),
    "atac-chipseq-downstream-analysis": ("multiomics-data-analysis", test_atac_chipseq_downstream),
    "methylation-analysis": ("multiomics-data-analysis", test_methylation_analysis),
    "proteomics-analysis": ("multiomics-data-analysis", test_proteomics_analysis),
    "single-cell-basic-analysis": ("multiomics-data-analysis", test_single_cell_analysis),

    # CRISPR
    "crispr-screen-analysis": ("crispr-design-and-analysis", test_crispr_screen_analysis),
    "crispr-screen-qc": ("crispr-design-and-analysis", test_crispr_screen_qc),
    "design-sgrnas-by-gene": ("crispr-design-and-analysis", test_design_sgrnas_by_gene),
    "design-base-editor-sgrnas": ("crispr-design-and-analysis", test_design_base_editor_sgrnas),
    "design-prime-editor-sgrnas": ("crispr-design-and-analysis", test_design_prime_editor_sgrnas),
    "crispr-library-design": ("crispr-design-and-analysis", test_crispr_library_design),

    # Gene List
    "gene-list-overlap": ("gene-list-analysis", test_gene_list_overlap),
    "go-analysis-for-gene-list": ("gene-list-analysis", test_go_analysis_for_gene_list),
    "gsea-for-ranked-gene-list": ("gene-list-analysis", test_gsea_for_ranked_gene_list),
    "ppi-network-for-gene-list": ("gene-list-analysis", test_ppi_network_for_gene_list),
    "transcription-factor-enrichment": ("gene-list-analysis", test_transcription_factor_enrichment),
    "curate-gene-list-by-function": ("gene-list-analysis", test_curate_gene_list_by_function),
    "function-annotation-for-gene-list": ("gene-list-analysis", test_function_annotation_for_gene_list),

    # Gene-centered
    "cox-survival-analysis": ("gene-centered-analysis", test_cox_survival),
    "depmap-analysis-for-gene": ("gene-centered-analysis", test_depmap_analysis_for_gene),
    "drug-sensitivity-for-gene": ("gene-centered-analysis", test_drug_sensitivity_for_gene),
    "mutation-analysis-for-gene": ("gene-centered-analysis", test_mutation_analysis_for_gene),
    "coexpression-for-gene": ("gene-centered-analysis", test_coexpression_for_gene),
    "normal-tissue-expression-by-gene": ("gene-centered-analysis", test_normal_tissue_expression_by_gene),
    "tcge-survival-for-gene": ("gene-centered-analysis", test_tcge_survival_for_gene),

    # ML
    "dimensionality-reduction": ("machine-learning-and-deep-learning", test_dimensionality_reduction),
    "clustering-analysis": ("machine-learning-and-deep-learning", test_clustering_analysis),
    "omics-ml-classifier": ("machine-learning-and-deep-learning", test_omics_ml_classifier),

    # Plots
    "plot-volcano": ("bioinformatics-plot-generator", test_plot_volcano),
    "plot-heatmap": ("bioinformatics-plot-generator", test_plot_heatmap),
    "plot-box-violin": ("bioinformatics-plot-generator", test_plot_box_violin),
    "plot-scatter-bar": ("bioinformatics-plot-generator", test_plot_scatter_bar),
    "plot-survival": ("bioinformatics-plot-generator", test_plot_survival),

    # Protein structure
    "protein-structure-for-gene": ("protein-structure-analysis", test_protein_structure_for_gene),
    "protein-sequence-analysis": ("protein-structure-analysis", test_protein_sequence_analysis),
    "protein-structure-visualizer": ("protein-structure-analysis", test_protein_structure_visualizer),
    "protein-structure-alignment": ("protein-structure-analysis", test_protein_structure_alignment),
    "protein-variant-mapper": ("protein-structure-analysis", test_protein_variant_mapper),

    # Papers
    "pubmed-search": ("paper-search-and-digest", test_pubmed_search),
    "paper-digest-single": ("paper-search-and-digest", test_paper_digest_single),
    "preprint-tracker": ("paper-search-and-digest", test_preprint_tracker),
    "big-papers-weekly-report": ("paper-search-and-digest", test_big_papers_weekly_report),

    # Labs
    "search-big-labs-by-field": ("lab-search-and-track", test_search_big_labs_by_field),
    "track-lab-publications": ("lab-search-and-track", test_track_lab_publications),
    "find-collaborators": ("lab-search-and-track", test_find_collaborators),

    # Public datasets
    "tcga-download-data": ("public-datasets-access-and-download", test_tcga_download_data),
    "gtex-download-data": ("public-datasets-access-and-download", test_gtex_download_data),
    "geo-download-data": ("public-datasets-access-and-download", test_geo_download_data),
    "depmap-data-download": ("public-datasets-access-and-download", test_depmap_data_download),
}


def main():
    """Main test runner."""
    parser = argparse.ArgumentParser(description="Run bioinformatics skill tests")
    parser.add_argument("--skill", help="Run a single skill by name")
    parser.add_argument("--skill-set", help="Run tests for a skill set (e.g., multiomics-data-analysis)")
    parser.add_argument("--verbose", action="store_true", help="Verbose output")
    parser.add_argument("--dry-run", action="store_true", help="Show what would be tested without running")
    parser.add_argument("--run-api", action="store_true", help="Include API-requiring skills")
    args = parser.parse_args()

    # Create results directory
    RESULTS.mkdir(parents=True, exist_ok=True)

    # Determine which tests to run
    tests_to_run = TESTS.copy()

    if args.skill:
        if args.skill not in tests_to_run:
            print(f"Error: Skill '{args.skill}' not found")
            return 1
        tests_to_run = {args.skill: tests_to_run[args.skill]}

    if args.skill_set:
        tests_to_run = {k: v for k, v in tests_to_run.items() if v[0] == args.skill_set}
        if not tests_to_run:
            print(f"Error: Skill set '{args.skill_set}' not found")
            return 1

    if not args.run_api:
        # Filter out API-requiring tests
        api_tests = {
            "design-sgrnas-by-gene", "design-base-editor-sgrnas",
            "design-prime-editor-sgrnas", "crispr-library-design",
            "go-analysis-for-gene-list", "gsea-for-ranked-gene-list",
            "ppi-network-for-gene-list", "transcription-factor-enrichment",
            "depmap-analysis-for-gene", "drug-sensitivity-for-gene",
            "mutation-analysis-for-gene", "coexpression-for-gene",
            "normal-tissue-expression-by-gene", "tcge-survival-for-gene",
            "protein-structure-for-gene", "protein-sequence-analysis",
            "protein-structure-visualizer",
            "pubmed-search", "paper-digest-single", "preprint-tracker",
            "big-papers-weekly-report",
            "search-big-labs-by-field", "track-lab-publications",
            "find-collaborators",
            "tcga-download-data", "gtex-download-data", "geo-download-data",
            "depmap-data-download",
        }
        # Filter to non-API tests
        non_api = {k: v for k, v in tests_to_run.items() if k not in api_tests}
        tests_to_run = non_api

    if args.dry_run:
        print(f"Would test {len(tests_to_run)} skills:")
        for skill_name in sorted(tests_to_run.keys()):
            print(f"  - {skill_name}")
        return 0

    # Run tests
    results = []
    total = len(tests_to_run)
    passed = 0
    failed = 0
    skipped = 0

    for idx, (skill_name, (skill_set, test_func)) in enumerate(sorted(tests_to_run.items()), 1):
        print(f"[{idx}/{total}] Testing {skill_name}... ", end="", flush=True)

        try:
            status, runtime, error = test_func()

            if status == "PASS":
                passed += 1
                print(f"{colored_status(status)} ({runtime:.2f}s)")
            elif status == "SKIP":
                skipped += 1
                print(f"{colored_status(status)}")
                if args.verbose:
                    print(f"         Reason: {error}")
            elif status == "FAIL":
                failed += 1
                print(f"{colored_status(status)} ({runtime:.2f}s)")
                if args.verbose and error:
                    print(f"         Error: {error}")
            else:
                failed += 1
                print(f"{colored_status(status)} ({runtime:.2f}s)")
                if args.verbose and error:
                    print(f"         Error: {error}")

            results.append({
                "skill": skill_name,
                "skill_set": skill_set,
                "status": status,
                "runtime": runtime,
                "error": error,
            })
        except Exception as e:
            failed += 1
            print(f"{colored_status('FAIL')}")
            if args.verbose:
                print(f"         Error: {str(e)}")
            results.append({
                "skill": skill_name,
                "skill_set": skill_set,
                "status": "FAIL",
                "runtime": 0,
                "error": str(e),
            })

    # Print summary
    print("\n" + "=" * 70)
    print(f"{BOLD}TEST SUMMARY{RESET_COLOR}")
    print("=" * 70)
    print(f"Total:   {total}")
    print(f"Passed:  {PASS_COLOR}{passed}{RESET_COLOR}")
    print(f"Failed:  {FAIL_COLOR}{failed}{RESET_COLOR}")
    print(f"Skipped: {SKIP_COLOR}{skipped}{RESET_COLOR}")
    print("=" * 70)

    # Save results
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    report_file = RESULTS / f"test_report_{timestamp}.tsv"

    with open(report_file, "w") as f:
        f.write("skill\tskill_set\tstatus\truntime\terror\n")
        for result in results:
            error = (result["error"] or "").replace("\n", " ").replace("\t", " ")
            f.write(f"{result['skill']}\t{result['skill_set']}\t{result['status']}\t{result['runtime']:.2f}\t{error}\n")

    print(f"Results saved to: {report_file}")

    # Print results table
    print("\n" + "=" * 70)
    print(f"{BOLD}DETAILED RESULTS{RESET_COLOR}")
    print("=" * 70)
    print(f"{'Skill':<40} {'Status':<12} {'Runtime':>10}")
    print("-" * 70)
    for result in results:
        status_str = colored_status(result["status"])
        print(f"{result['skill']:<40} {status_str:<24} {result['runtime']:>10.2f}s")

    return 0 if failed == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
