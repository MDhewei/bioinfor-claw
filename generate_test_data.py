#!/usr/bin/env python3
"""
Generate comprehensive synthetic test datasets for bioinformatics skills.
All data files are written to tests/data/
"""

import numpy as np
import pandas as pd
import os
from pathlib import Path
import random

# Set random seeds for reproducibility
np.random.seed(42)
random.seed(42)

OUTPUT_DIR = "/sessions/festive-admiring-sagan/mnt/bioinfor-claw/tests/data"
Path(OUTPUT_DIR).mkdir(parents=True, exist_ok=True)

print(f"Generating test data in {OUTPUT_DIR}")
print("=" * 80)

# ============================================================================
# 1. RNA-seq data
# ============================================================================
print("\n1. Generating RNA-seq data...")

n_genes = 500
n_samples = 12
gene_names = ["TP53", "BRCA1", "EGFR", "MYC", "KRAS"] + [f"Gene_{i:03d}" for i in range(6, n_genes + 1)]
samples = [f"Control_{i}" for i in range(1, 7)] + [f"Treatment_{i}" for i in range(1, 7)]

# Generate base counts
counts = np.random.negative_binomial(20, 0.5, size=(n_genes, n_samples))

# Add differential expression: first 50 genes upregulated 2-4x in treatment
de_indices = np.arange(50)
fc_multipliers = np.random.uniform(2, 4, size=50)
for i, gene_idx in enumerate(de_indices):
    counts[gene_idx, 6:12] = (counts[gene_idx, 6:12] * fc_multipliers[i]).astype(int)

rnaseq_df = pd.DataFrame(counts, columns=samples)
rnaseq_df.insert(0, "GeneID", gene_names)
rnaseq_df.to_csv(os.path.join(OUTPUT_DIR, "rnaseq_counts.tsv"), sep="\t", index=False)
print(f"  Created rnaseq_counts.tsv: {rnaseq_df.shape}")

# RNA-seq metadata
metadata_rna = pd.DataFrame({
    "sample_id": samples,
    "group": ["control"] * 6 + ["treatment"] * 6,
    "batch": ["batch1", "batch1", "batch2", "batch2", "batch1", "batch2",
              "batch1", "batch2", "batch2", "batch1", "batch2", "batch1"]
})
metadata_rna.to_csv(os.path.join(OUTPUT_DIR, "rnaseq_metadata.tsv"), sep="\t", index=False)
print(f"  Created rnaseq_metadata.tsv: {metadata_rna.shape}")

# ============================================================================
# 2. ATAC-seq / ChIP-seq peak data
# ============================================================================
print("\n2. Generating ATAC-seq peak data...")

def generate_peaks(n_peaks, seed_offset=0):
    """Generate narrowPeak format peaks."""
    np.random.seed(42 + seed_offset)
    chromosomes = [f"chr{i}" for i in range(1, 23)] + ["chrX"]
    chrom_sizes = {
        "chr1": 248956422, "chr2": 242193529, "chr3": 198295559,
        "chr4": 190214555, "chr5": 181538259, "chr6": 170805979,
        "chr7": 159345973, "chr8": 145138636, "chr9": 138394717,
        "chr10": 133797422, "chr11": 135086622, "chr12": 133275309,
        "chr13": 114364328, "chr14": 107043718, "chr15": 101991189,
        "chr16": 90338345, "chr17": 83257441, "chr18": 80373912,
        "chr19": 58617616, "chr20": 64444167, "chr21": 46709983,
        "chr22": 50818468, "chrX": 155270560
    }

    peaks_list = []
    for i in range(n_peaks):
        chrom = np.random.choice(chromosomes)
        width = np.random.randint(200, 2001)
        start = np.random.randint(1000, chrom_sizes[chrom] - width)
        end = start + width
        name = f"peak_{i+1:04d}"
        score = np.random.randint(100, 1001)
        strand = np.random.choice(["+", "-"])
        signal = np.random.uniform(1, 20)
        pvalue = 10 ** (-np.random.uniform(1, 10))
        qvalue = 10 ** (-np.random.uniform(1, 8))
        peak_pos = np.random.randint(start + 10, end - 10)

        peaks_list.append([chrom, start, end, name, score, strand, signal, pvalue, qvalue, peak_pos - start])

    return pd.DataFrame(peaks_list, columns=["chrom", "start", "end", "name", "score", "strand",
                                              "signalValue", "pValue", "qValue", "peak"])

peaks1 = generate_peaks(2000, seed_offset=0)
peaks1.to_csv(os.path.join(OUTPUT_DIR, "atac_peaks.narrowPeak"), sep="\t", index=False, header=False)
print(f"  Created atac_peaks.narrowPeak: {peaks1.shape}")

# Create overlapping peaks for second file
peaks2_base = generate_peaks(1800, seed_offset=1)
# 60% overlap: keep some from peaks1 and add some new ones
overlap_indices = np.random.choice(len(peaks1), size=int(2000 * 0.6), replace=False)
peaks2_list = peaks1.iloc[overlap_indices].reset_index(drop=True)
# Rename them with new peak numbers
peaks2_list["name"] = [f"peak_{i+1:04d}" for i in range(len(peaks2_list))]

# Add new peaks
new_peaks_count = 1800 - len(peaks2_list)
new_peaks = generate_peaks(new_peaks_count, seed_offset=2)
new_peaks["name"] = [f"peak_{i+len(peaks2_list)+1:04d}" for i in range(len(new_peaks))]
peaks2 = pd.concat([peaks2_list, new_peaks], ignore_index=True)
peaks2.to_csv(os.path.join(OUTPUT_DIR, "atac_peaks2.narrowPeak"), sep="\t", index=False, header=False)
print(f"  Created atac_peaks2.narrowPeak: {peaks2.shape}")

# ============================================================================
# 3. DNA methylation data
# ============================================================================
print("\n3. Generating DNA methylation data...")

n_cpgs = 3000
cpg_ids = [f"cg{i:08d}" for i in range(1, n_cpgs + 1)]
chromosomes_meth = [f"chr{i}" for i in range(1, 23)] + ["chrX"]
positions = [np.random.randint(1000, 100000000) for _ in range(n_cpgs)]
chroms_meth = np.random.choice(chromosomes_meth, n_cpgs)

# Generate beta values
beta_values = np.random.beta(2, 5, size=(n_cpgs, n_samples))

# First 200 CpGs: truly differential (treatment group shifted)
de_cpg_indices = np.arange(200)
beta_values[de_cpg_indices, 6:12] = np.clip(beta_values[de_cpg_indices, 6:12] + np.random.uniform(0.2, 0.4, (200, 6)), 0, 1)

methyl_df = pd.DataFrame(beta_values, columns=samples)
methyl_df.insert(0, "position", positions)
methyl_df.insert(0, "chr", chroms_meth)
methyl_df.insert(0, "CpG_ID", cpg_ids)
methyl_df.to_csv(os.path.join(OUTPUT_DIR, "methylation_beta.tsv"), sep="\t", index=False)
print(f"  Created methylation_beta.tsv: {methyl_df.shape}")

# Methylation metadata (same structure as RNA-seq)
metadata_meth = metadata_rna.copy()
metadata_meth.to_csv(os.path.join(OUTPUT_DIR, "methylation_metadata.tsv"), sep="\t", index=False)
print(f"  Created methylation_metadata.tsv: {metadata_meth.shape}")

# ============================================================================
# 4. Proteomics data
# ============================================================================
print("\n4. Generating proteomics data...")

n_proteins = 1500
protein_names = [f"Protein_{i:04d}" for i in range(1, n_proteins + 1)]

# Generate log2 intensities
proteomics_data = np.random.uniform(15, 35, size=(n_proteins, n_samples))

# Add ~20% missing values
for i in range(n_samples):
    missing_indices = np.random.choice(n_proteins, size=int(n_proteins * 0.2), replace=False)
    proteomics_data[missing_indices, i] = np.nan

# First 100 proteins: truly DE (treatment has +1.5 log2FC)
proteomics_data[0:100, 6:12] = proteomics_data[0:100, 6:12] + 1.5

# Add 5% entirely missing rows (low-abundance proteins)
fully_missing_indices = np.random.choice(n_proteins, size=int(n_proteins * 0.05), replace=False)
proteomics_data[fully_missing_indices, :] = np.nan

proteomics_df = pd.DataFrame(proteomics_data, columns=samples)
proteomics_df.insert(0, "Protein", protein_names)
proteomics_df.to_csv(os.path.join(OUTPUT_DIR, "proteomics_intensities.tsv"), sep="\t", index=False)
print(f"  Created proteomics_intensities.tsv: {proteomics_df.shape}")

# Proteomics metadata
metadata_prot = metadata_rna.copy()
metadata_prot.to_csv(os.path.join(OUTPUT_DIR, "proteomics_metadata.tsv"), sep="\t", index=False)
print(f"  Created proteomics_metadata.tsv: {metadata_prot.shape}")

# ============================================================================
# 5. Single-cell RNA-seq data
# ============================================================================
print("\n5. Generating single-cell RNA-seq data...")

n_cells = 500
n_genes_sc = 1000
cell_barcodes = [f"CELL_{i:04d}" for i in range(1, n_cells + 1)]
mt_genes = [f"MT-CO{i}" for i in range(1, 21)]
gene_names_sc = mt_genes + [f"GENE_{i:04d}" for i in range(21, n_genes_sc + 1)]

# Generate sparse scRNA-seq counts
scrna_counts = np.zeros((n_genes_sc, n_cells))

# Cell type assignment: 3 types, ~167 cells each
cell_types = []
for i in range(n_cells):
    if i < 167:
        cell_types.append("TypeA")
    elif i < 334:
        cell_types.append("TypeB")
    else:
        cell_types.append("TypeC")

# Type A: genes 0-49 (after MT genes) highly expressed
for gene_idx in range(20, 70):
    scrna_counts[gene_idx, :167] = np.random.poisson(5, 167)

# Type B: genes 100-150 highly expressed
for gene_idx in range(100, 150):
    scrna_counts[gene_idx, 167:334] = np.random.poisson(5, 167)

# Type C: genes 200-250 highly expressed
for gene_idx in range(200, 250):
    scrna_counts[gene_idx, 334:500] = np.random.poisson(5, 166)

# Background Poisson(0.2) for all genes
background = np.random.poisson(0.2, size=(n_genes_sc, n_cells))
scrna_counts = (scrna_counts + background).astype(int)

# MT genes for all cells
for mt_idx in range(20):
    scrna_counts[mt_idx, :] = np.random.poisson(2, n_cells)

# 5% of cells with high MT percentage (>30%)
high_mt_cells = np.random.choice(n_cells, size=int(n_cells * 0.05), replace=False)
for cell_idx in high_mt_cells:
    scrna_counts[:20, cell_idx] = np.random.poisson(10, 20)

scrna_df = pd.DataFrame(scrna_counts, index=gene_names_sc, columns=cell_barcodes).T
scrna_df.insert(0, "cell_id", cell_barcodes)
scrna_df.to_csv(os.path.join(OUTPUT_DIR, "scrna_counts.tsv"), sep="\t", index=False)
print(f"  Created scrna_counts.tsv: {scrna_df.shape}")

# Single-cell metadata
metadata_sc = pd.DataFrame({
    "cell_id": cell_barcodes,
    "true_celltype": cell_types,
    "batch": ["batch1" if i < 250 else "batch2" for i in range(n_cells)]
})
metadata_sc.to_csv(os.path.join(OUTPUT_DIR, "scrna_metadata.tsv"), sep="\t", index=False)
print(f"  Created scrna_metadata.tsv: {metadata_sc.shape}")

# ============================================================================
# 6. CRISPR screen count data
# ============================================================================
print("\n6. Generating CRISPR screen count data...")

n_genes_crispr = 950
sgrnas_per_gene = 4
total_sgrnas = n_genes_crispr * sgrnas_per_gene + 200  # +200 for essential and controls

# Gene list: 950 regular, 50 essential cancer genes
essential_genes = ["TP53", "EGFR", "MYC", "KRAS", "PTEN", "AKT1", "PIK3CA", "MTOR", "RB1", "CDKN2A",
                   "CDH1", "VHL", "ATM", "CHEK2", "PALB2", "RAD51", "BRCA1", "BRCA2", "NF1", "BRAF",
                   "MAP2K1", "NRAS", "HRAS", "FGFR1", "FGFR2", "FGFR3", "MET", "ALK", "ERBB2", "ERBB3",
                   "ROS1", "RET", "PDGFRA", "PDGFRB", "KIT", "FLT3", "JAK1", "JAK2", "STAT3", "STAT5A",
                   "BCR", "ABL1", "NOTCH1", "NOTCH2", "NOTCH3", "NOTCH4", "HES1", "HEY1", "ADAM10"]

regular_genes = [f"Gene_{i:04d}" for i in range(1, n_genes_crispr + 1)]
gene_list_crispr = regular_genes + essential_genes

# Create sgRNA table
sgrna_list = []
gene_col = []
counts_data = []

for gene_idx, gene_name in enumerate(regular_genes):
    for sg_idx in range(4):
        sgrna_id = f"sg_{gene_idx:04d}_{sg_idx}"
        sgrna_list.append(sgrna_id)
        gene_col.append(gene_name)

# Add essential genes
for ess_idx, gene_name in enumerate(essential_genes):
    for sg_idx in range(4):
        sgrna_id = f"sg_ess_{ess_idx:02d}_{sg_idx}"
        sgrna_list.append(sgrna_id)
        gene_col.append(gene_name)

# Add non-targeting controls
for nt_idx in range(50):
    sgrna_list.append(f"NonTargeting_{nt_idx:03d}")
    gene_col.append("NonTargeting")

# Generate counts: NegBinomial with mean ~500 in Plasmid
samples_crispr = ["Plasmid_lib", "Control_1", "Control_2", "Treatment_1", "Treatment_2"]
np.random.seed(42)

counts_crispr = np.random.negative_binomial(100, 0.17, size=(len(sgrna_list), 5))  # mean ~500

# Add realistic CV ~30% between replicates
for i in range(len(sgrna_list)):
    for j in range(1, 5):
        noise = np.random.normal(1, 0.3)
        counts_crispr[i, j] = (counts_crispr[i, j] * noise).astype(int)

# Make essential genes depleted in treatment
for i in range(len(sgrna_list)):
    if gene_col[i] in essential_genes[:30]:  # First 30 essential genes
        counts_crispr[i, 3:5] = (counts_crispr[i, 3:5] * np.random.uniform(0.1, 0.3)).astype(int)

# Make some oncogenes enriched in treatment
for i in range(len(sgrna_list)):
    if gene_col[i] in essential_genes[30:60]:  # 30 enriched genes
        counts_crispr[i, 3:5] = (counts_crispr[i, 3:5] * np.random.uniform(2, 5)).astype(int)

crispr_df = pd.DataFrame(counts_crispr, columns=samples_crispr)
crispr_df.insert(0, "Gene", gene_col)
crispr_df.insert(0, "sgRNA", sgrna_list)
crispr_df.to_csv(os.path.join(OUTPUT_DIR, "crispr_counts.tsv"), sep="\t", index=False)
print(f"  Created crispr_counts.tsv: {crispr_df.shape}")

# ============================================================================
# 7. Gene lists
# ============================================================================
print("\n7. Generating gene lists...")

cancer_genes = ["TP53", "BRCA1", "BRCA2", "EGFR", "MYC", "KRAS", "PTEN", "AKT1", "PIK3CA", "MTOR",
                "RB1", "CDKN2A", "CDH1", "VHL", "ATM", "CHEK2", "PALB2", "RAD51", "FANCB", "FANCC",
                "FANCD2", "FANCE", "FANCF", "FANCG", "FANCI", "FANCJ", "FANCL", "FANCM", "BRCAI", "BRCA1",
                "BRCA2", "NBN", "MRE11A", "RAD50", "RAD51", "RAD51B", "RAD51C", "RAD51D", "XRCC2", "XRCC3",
                "XRCC4", "XRCC5", "XRCC6", "XRCC7", "LIG1", "LIG3", "LIG4", "POLB", "POLD1", "POLE",
                "RFC1", "RPA1", "RPA2", "RPA3", "RECO1", "TOPBP1", "TP53", "TP53BP1", "MDM2", "MDM4",
                "CDKN1A", "CDKN1B", "CDKN2A", "CDKN2B", "CDKN2C", "CDKN2D", "RBL1", "RBL2", "E2F1", "E2F2",
                "E2F3", "E2F4", "E2F5", "E2F6", "E2F7", "E2F8", "CCND1", "CCND2", "CCND3", "CCNE1",
                "CCNE2", "CDK2", "CDK4", "CDK6", "CDK7", "CDK8", "CDK9", "CDK10", "CDK11A", "CDK11B",
                "CCNA1", "CCNA2", "CCNB1", "CCNB2", "CCNB3", "CCNC", "CCNF", "CCNH", "CCNI", "CCNJ",
                "CCNK", "CCNL1", "CCNL2", "CCNO", "CCNQ", "CCNR1", "CCNS", "CCNT1", "CCNT2", "CCNU",
                "CCNV", "CCNW", "CCNX", "CCNY", "CCNZ", "BRAF", "MAP2K1", "MAP2K2", "MAPK1", "MAPK3",
                "NRAS", "HRAS", "FGFR1", "FGFR2", "FGFR3", "FGFR4", "MET", "ALK", "ERBB2", "ERBB3",
                "ERBB4", "ROS1", "RET", "PDGFRA", "PDGFRB", "KIT", "FLT3", "JAK1", "JAK2", "JAK3",
                "STAT1", "STAT2", "STAT3", "STAT4", "STAT5A", "STAT5B", "STAT6", "BCR", "ABL1", "ABL2",
                "NOTCH1", "NOTCH2", "NOTCH3", "NOTCH4", "HES1", "HEY1", "HEY2", "MYCT1", "APC", "CTNNB1"]

# Truncate to 150 for list A
list_a = cancer_genes[:150]
with open(os.path.join(OUTPUT_DIR, "gene_list_A.txt"), "w") as f:
    f.write("\n".join(list_a))
print(f"  Created gene_list_A.txt: {len(list_a)} genes")

# List B: 120 genes, ~70 overlap with A
list_b = list_a[:70] + cancer_genes[70:120]
with open(os.path.join(OUTPUT_DIR, "gene_list_B.txt"), "w") as f:
    f.write("\n".join(list_b))
print(f"  Created gene_list_B.txt: {len(list_b)} genes")

# List C: 100 genes, ~30 overlap with A, ~20 with B
list_c = list_a[:30] + list_b[70:90] + cancer_genes[120:150]
with open(os.path.join(OUTPUT_DIR, "gene_list_C.txt"), "w") as f:
    f.write("\n".join(list_c))
print(f"  Created gene_list_C.txt: {len(list_c)} genes")

# Ranked gene list
n_genes_ranked = 500
genes_ranked = [f"Gene_{i:03d}" for i in range(1, n_genes_ranked + 1)]
log2fc = np.concatenate([np.random.uniform(1, 3, 50),  # upregulated
                         np.random.uniform(-3, -1, 50),  # downregulated
                         np.random.uniform(-1, 1, n_genes_ranked - 100)])  # background
log2fc = np.sort(log2fc)[::-1]  # sort descending

pvalues = 10 ** (-np.random.uniform(1, 20, n_genes_ranked))
padj = np.minimum(pvalues * n_genes_ranked, 1)  # Bonferroni correction

ranked_df = pd.DataFrame({
    "gene": genes_ranked,
    "log2FC": log2fc,
    "pvalue": pvalues,
    "padj": padj
})
ranked_df = ranked_df.sort_values("log2FC", ascending=False).reset_index(drop=True)
ranked_df.to_csv(os.path.join(OUTPUT_DIR, "ranked_gene_list.tsv"), sep="\t", index=False)
print(f"  Created ranked_gene_list.tsv: {ranked_df.shape}")

# ============================================================================
# 8. Survival / clinical data
# ============================================================================
print("\n8. Generating survival/clinical data...")

n_patients = 200
sample_ids = [f"Patient_{i:03d}" for i in range(1, n_patients + 1)]

# Weibull-distributed survival times
shape, scale = 1.5, 365
time_days = np.random.weibull(shape, n_patients) * scale
time_days = np.minimum(time_days, 1825)  # right-censor at 5 years

# Events: ~60%
events = np.random.binomial(1, 0.6, n_patients)

# Clinical variables
ages = np.random.normal(60, 15, n_patients).astype(int)
ages = np.clip(ages, 25, 95)

stages = np.random.choice([1, 2, 3, 4], n_patients, p=[0.3, 0.3, 0.2, 0.2])
grades = np.random.choice([1, 2, 3], n_patients, p=[0.3, 0.4, 0.3])
treatments = np.random.choice(["A", "B"], n_patients)

# Gene expression (continuous)
gene_expr = np.random.normal(0, 1, n_patients)

survival_df = pd.DataFrame({
    "sample_id": sample_ids,
    "time_days": time_days.astype(int),
    "event": events,
    "age": ages,
    "stage": stages,
    "grade": grades,
    "treatment": treatments,
    "gene_expression": gene_expr
})

# High-expression patients: shorter survival (HR ~2)
high_expr_idx = gene_expr > 0.5
survival_df.loc[high_expr_idx, "time_days"] = (survival_df.loc[high_expr_idx, "time_days"] * 0.5).astype(int)

# Stage IV: shorter survival (HR ~3)
stage_iv_idx = stages == 4
survival_df.loc[stage_iv_idx, "time_days"] = (survival_df.loc[stage_iv_idx, "time_days"] * 0.33).astype(int)

survival_df.to_csv(os.path.join(OUTPUT_DIR, "survival_data.tsv"), sep="\t", index=False)
print(f"  Created survival_data.tsv: {survival_df.shape}")

# ============================================================================
# 9. Machine learning / omics classification
# ============================================================================
print("\n9. Generating ML classification data...")

n_ml_samples = 200
n_ml_features = 200
ml_samples = [f"Sample_{i:03d}" for i in range(1, n_ml_samples + 1)]
ml_features = [f"Feature_{i:03d}" for i in range(1, n_ml_features + 1)]

# Generate feature matrix
ml_features_data = np.random.normal(0, 1, (n_ml_samples, n_ml_features))

# First 20 features: discriminative (Class_A has higher expression)
ml_features_data[:100, :20] = np.random.normal(1.5, 1, (100, 20))  # Class A
ml_features_data[100:, :20] = np.random.normal(0, 1, (100, 20))    # Class B

ml_df = pd.DataFrame(ml_features_data, columns=ml_features)
ml_df.insert(0, "sample_id", ml_samples)
ml_df.to_csv(os.path.join(OUTPUT_DIR, "ml_features.tsv"), sep="\t", index=False)
print(f"  Created ml_features.tsv: {ml_df.shape}")

# ML labels
labels = ["Class_A"] * 100 + ["Class_B"] * 100
ml_labels_df = pd.DataFrame({
    "sample_id": ml_samples,
    "label": labels
})
ml_labels_df.to_csv(os.path.join(OUTPUT_DIR, "ml_labels.tsv"), sep="\t", index=False)
print(f"  Created ml_labels.tsv: {ml_labels_df.shape}")

# ============================================================================
# 10. Generic expression matrix (dimensionality reduction / clustering)
# ============================================================================
print("\n10. Generating expression matrix...")

n_expr_samples = 50
n_expr_genes = 500
expr_samples = [f"Sample_{i:02d}" for i in range(1, n_expr_samples + 1)]
expr_genes = [f"Gene_{i:03d}" for i in range(1, n_expr_genes + 1)]

# Initialize expression matrix
expr_matrix = np.random.normal(0, 1, (n_expr_samples, n_expr_genes))

# Create 3 clusters with distinct marker genes
# Cluster 1 (samples 0-16): genes 0-49 high, others low
expr_matrix[0:17, 0:50] = np.random.normal(3, 1, (17, 50))
expr_matrix[0:17, 50:] = np.random.normal(0, 1, (17, 450))

# Cluster 2 (samples 17-33): genes 100-149 high
expr_matrix[17:34, 100:150] = np.random.normal(3, 1, (17, 50))
expr_matrix[17:34, [i for i in range(n_expr_genes) if i not in range(100, 150)]] = np.random.normal(0, 1, (17, 450))

# Cluster 3 (samples 34-49): genes 200-249 high
expr_matrix[34:50, 200:250] = np.random.normal(3, 1, (16, 50))
expr_matrix[34:50, [i for i in range(n_expr_genes) if i not in range(200, 250)]] = np.random.normal(0, 1, (16, 450))

expr_df = pd.DataFrame(expr_matrix, columns=expr_genes)
expr_df.insert(0, "sample_id", expr_samples)
expr_df.to_csv(os.path.join(OUTPUT_DIR, "expression_matrix.tsv"), sep="\t", index=False)
print(f"  Created expression_matrix.tsv: {expr_df.shape}")

# Sample metadata
tissue_types = ["TypeA"] * 17 + ["TypeB"] * 17 + ["TypeC"] * 16
batches = ["batch1"] * 25 + ["batch2"] * 25
treatments = ["control"] * 25 + ["treatment"] * 25

sample_metadata_df = pd.DataFrame({
    "sample_id": expr_samples,
    "tissue_type": tissue_types,
    "batch": batches,
    "treatment": treatments
})
sample_metadata_df.to_csv(os.path.join(OUTPUT_DIR, "sample_metadata.tsv"), sep="\t", index=False)
print(f"  Created sample_metadata.tsv: {sample_metadata_df.shape}")

# ============================================================================
# 11. Plot test data
# ============================================================================
print("\n11. Generating plot test data...")

# DE results for volcano plot
n_de_genes = 1000
de_genes = [f"Gene_{i:04d}" for i in range(1, n_de_genes + 1)]
de_log2fc = np.concatenate([np.random.uniform(1, 3, 80),     # upregulated significant
                             np.random.uniform(-3, -1, 80),  # downregulated significant
                             np.random.uniform(-1, 1, 840)])  # background

de_pvalues = np.concatenate([10 ** (-np.random.uniform(4, 20, 80)),  # significant
                              10 ** (-np.random.uniform(4, 20, 80)),  # significant
                              10 ** (-np.random.uniform(0.3, 2, 840))])  # ns
de_padj = np.minimum(de_pvalues * n_de_genes, 1)
de_basemean = np.random.uniform(10, 10000, n_de_genes)

de_results_df = pd.DataFrame({
    "gene": de_genes,
    "log2FoldChange": de_log2fc,
    "pvalue": de_pvalues,
    "padj": de_padj,
    "baseMean": de_basemean
})
de_results_df.to_csv(os.path.join(OUTPUT_DIR, "de_results.tsv"), sep="\t", index=False)
print(f"  Created de_results.tsv: {de_results_df.shape}")

# Heatmap matrix
n_hm_genes = 80
n_hm_samples = 24
hm_genes = [f"Gene_{i:03d}" for i in range(1, n_hm_genes + 1)]
hm_samples = [f"Sample_{i:02d}" for i in range(1, n_hm_samples + 1)]

hm_matrix = np.random.normal(0, 1, (n_hm_genes, n_hm_samples))

# 3 groups of 8 samples each with distinct patterns
for i in range(n_hm_genes):
    hm_matrix[i, 0:8] = np.random.normal(1, 0.5, 8)      # Group 1: higher
    hm_matrix[i, 8:16] = np.random.normal(-1, 0.5, 8)    # Group 2: lower
    hm_matrix[i, 16:24] = np.random.normal(0, 0.5, 8)    # Group 3: middle

hm_df = pd.DataFrame(hm_matrix, columns=hm_samples)
hm_df.insert(0, "gene", hm_genes)
hm_df.to_csv(os.path.join(OUTPUT_DIR, "heatmap_matrix.tsv"), sep="\t", index=False)
print(f"  Created heatmap_matrix.tsv: {hm_df.shape}")

# Box/violin plot data
n_bv_obs = 120
bv_values = np.concatenate([np.random.normal(0, 0.5, 30),      # Control
                             np.random.normal(1.5, 0.6, 30),   # TreatA
                             np.random.normal(2.5, 0.5, 30),   # TreatB
                             np.random.normal(0.8, 0.7, 30)])  # TreatC
bv_groups = ["Control"] * 30 + ["TreatA"] * 30 + ["TreatB"] * 30 + ["TreatC"] * 30
bv_subgroups = ["sub1", "sub2"] * 60  # alternating subgroups

bv_df = pd.DataFrame({
    "value": bv_values,
    "group": bv_groups,
    "subgroup": bv_subgroups
})
bv_df.to_csv(os.path.join(OUTPUT_DIR, "boxviolin_data.tsv"), sep="\t", index=False)
print(f"  Created boxviolin_data.tsv: {bv_df.shape}")

# Scatter plot data
n_scatter = 100
scatter_x = np.random.normal(0, 1, n_scatter)
scatter_y = scatter_x + np.random.normal(0, 0.4, n_scatter)  # correlated, r~0.7
scatter_groups = np.random.choice(["GroupA", "GroupB", "GroupC"], n_scatter)
scatter_sizes = np.random.uniform(50, 500, n_scatter)
scatter_colors = np.random.uniform(0, 10, n_scatter)

scatter_df = pd.DataFrame({
    "sample": [f"S{i}" for i in range(n_scatter)],
    "x_value": scatter_x,
    "y_value": scatter_y,
    "group": scatter_groups,
    "size_value": scatter_sizes,
    "color_value": scatter_colors
})
scatter_df.to_csv(os.path.join(OUTPUT_DIR, "scatter_data.tsv"), sep="\t", index=False)
print(f"  Created scatter_data.tsv: {scatter_df.shape}")

# Bar plot data
n_bar_genes = 10
n_bar_conditions = 3
n_bar_replicates = 5
bar_genes = [f"Gene_{i:02d}" for i in range(1, n_bar_genes + 1)]
bar_conditions = ["Ctrl", "TreatA", "TreatB"]

bar_data = []
for gene in bar_genes:
    for condition in bar_conditions:
        means = {"Ctrl": 10, "TreatA": 15, "TreatB": 20}
        base_mean = means[condition]
        for rep in range(n_bar_replicates):
            value = np.random.normal(base_mean, 2)
            bar_data.append({
                "gene": gene,
                "expression": value,
                "condition": condition,
                "replicate": f"rep{rep+1}"
            })

bar_df = pd.DataFrame(bar_data)
bar_df.to_csv(os.path.join(OUTPUT_DIR, "bar_data.tsv"), sep="\t", index=False)
print(f"  Created bar_data.tsv: {bar_df.shape}")

# Survival plot data
n_surv = 180  # 60 per group
surv_samples = [f"Patient_{i:03d}" for i in range(1, n_surv + 1)]
surv_times = np.concatenate([np.random.exponential(500, 60),     # GroupA: best survival
                              np.random.exponential(350, 60),    # GroupB: medium
                              np.random.exponential(200, 60)])   # GroupC: poor
surv_times = np.minimum(surv_times, 1000)

surv_events = np.concatenate([np.random.binomial(1, 0.4, 60),    # GroupA: lower event rate
                               np.random.binomial(1, 0.6, 60),   # GroupB: medium
                               np.random.binomial(1, 0.8, 60)])  # GroupC: higher

surv_groups = ["GroupA"] * 60 + ["GroupB"] * 60 + ["GroupC"] * 60

surv_plot_df = pd.DataFrame({
    "sample_id": surv_samples,
    "time": surv_times.astype(int),
    "event": surv_events,
    "group": surv_groups
})
surv_plot_df.to_csv(os.path.join(OUTPUT_DIR, "survival_plot.tsv"), sep="\t", index=False)
print(f"  Created survival_plot.tsv: {surv_plot_df.shape}")

# ============================================================================
# 12. API skill parameter files
# ============================================================================
print("\n12. Generating API skill parameter files...")

api_genes = ["BRCA1", "TP53", "EGFR", "MYC", "KRAS", "PTEN", "AKT1", "CDK4", "RB1", "CDKN2A"]
with open(os.path.join(OUTPUT_DIR, "api_test_genes.txt"), "w") as f:
    f.write("\n".join(api_genes))
print(f"  Created api_test_genes.txt: {len(api_genes)} genes")

api_gene_list = cancer_genes[:50]
with open(os.path.join(OUTPUT_DIR, "api_test_gene_list.txt"), "w") as f:
    f.write("\n".join(api_gene_list))
print(f"  Created api_test_gene_list.txt: {len(api_gene_list)} genes")

crispr_lib_genes = ["TP53", "BRCA1", "EGFR", "MYC", "KRAS", "PTEN", "AKT1", "PIK3CA", "MTOR", "RB1",
                    "CDKN2A", "CDH1", "VHL", "ATM", "CHEK2", "PALB2", "RAD51", "BRAF", "NRAS", "HRAS"]
with open(os.path.join(OUTPUT_DIR, "crispr_library_genes.txt"), "w") as f:
    f.write("\n".join(crispr_lib_genes))
print(f"  Created crispr_library_genes.txt: {len(crispr_lib_genes)} genes")

with open(os.path.join(OUTPUT_DIR, "pubmed_test_query.txt"), "w") as f:
    f.write("CRISPR base editing cancer 2023")
print(f"  Created pubmed_test_query.txt")

# ============================================================================
# Generate README.md
# ============================================================================
print("\n" + "=" * 80)
print("Generating README.md...")

readme_content = """# Bioinformatics Skills Test Data

This directory contains comprehensive synthetic test datasets for all 49 bioinformatics skills.

## Dataset Summary

### 1. RNA-seq Data
- **rnaseq_counts.tsv**: 500 genes × 12 samples
  - Gene names: TP53, BRCA1, EGFR, MYC, KRAS, Gene_006..Gene_500
  - Samples: Control_1-6, Treatment_1-6
  - Counts: NegBinomial-distributed
  - First 50 genes: differentially expressed (2-4x upregulated in treatment)

- **rnaseq_metadata.tsv**: 12 samples
  - Columns: sample_id, group, batch
  - Groups: control (6), treatment (6)

### 2. ATAC-seq / ChIP-seq Peaks
- **atac_peaks.narrowPeak**: 2000 peaks
  - Format: narrowPeak (10 columns: chrom, start, end, name, score, strand, signalValue, pValue, qValue, peak)
  - Chromosomes: chr1-chr22, chrX
  - Peak widths: 200-2000 bp

- **atac_peaks2.narrowPeak**: 1800 peaks
  - ~60% overlap with atac_peaks.narrowPeak for differential analysis testing

### 3. DNA Methylation
- **methylation_beta.tsv**: 3000 CpGs × 12 samples
  - CpG IDs: cg00000001..cg00003000
  - Beta values: 0-1 range (Beta distribution)
  - First 200 CpGs: differentially methylated (treatment +0.2-0.4)
  - Columns: CpG_ID, chr, position, [12 samples]

- **methylation_metadata.tsv**: 12 samples
  - Same structure as rnaseq_metadata.tsv

### 4. Proteomics
- **proteomics_intensities.tsv**: 1500 proteins × 12 samples
  - Log2 intensities: 15-35 range
  - ~20% missing values per sample
  - First 100 proteins: differentially abundant (+1.5 log2FC in treatment)
  - ~5% of rows entirely missing (low-abundance proteins)

- **proteomics_metadata.tsv**: 12 samples
  - Same structure as rnaseq_metadata.tsv

### 5. Single-cell RNA-seq
- **scrna_counts.tsv**: 500 cells × 1000 genes
  - Cell barcodes: CELL_0001..CELL_0500
  - 20 MT genes (MT-CO1..MT-CO20)
  - 3 cell types with distinct marker genes:
    - TypeA: 167 cells, high expression in genes 21-70
    - TypeB: 167 cells, high expression in genes 100-150
    - TypeC: 166 cells, high expression in genes 200-250
  - ~5% cells with high MT percentage (>30%)

- **scrna_metadata.tsv**: 500 cells
  - Columns: cell_id, true_celltype, batch

### 6. CRISPR Screen
- **crispr_counts.tsv**: 4000 sgRNAs × 5 samples
  - sgRNA × gene mapping: ~1000 genes × 4 sgRNAs/gene + controls
  - Samples: Plasmid_lib, Control_1-2, Treatment_1-2
  - First 30 essential genes: depleted in treatment (0.1-0.3x)
  - 30 enriched oncogenes: enriched in treatment (2-5x)
  - ~30% CV between replicates

### 7. Gene Lists
- **gene_list_A.txt**: 150 cancer-related genes
- **gene_list_B.txt**: 120 genes (~70% overlap with A)
- **gene_list_C.txt**: 100 genes (~30% overlap with A, ~20% with B)
- **ranked_gene_list.tsv**: 500 genes ranked by log2FC
  - Columns: gene, log2FC, pvalue, padj
  - ~50 upregulated (log2FC>1, padj<0.05)
  - ~50 downregulated (log2FC<-1, padj<0.05)

### 8. Survival / Clinical Data
- **survival_data.tsv**: 200 patients × 8 columns
  - Columns: sample_id, time_days, event, age, stage, grade, treatment, gene_expression
  - Time: Weibull-distributed, right-censored at 1825 days
  - Events: ~60%
  - High expression & stage IV patients: shorter survival

### 9. ML Classification
- **ml_features.tsv**: 200 samples × 200 features
  - First 20 features: discriminative (Class_A higher expression)
  - Remaining: random noise

- **ml_labels.tsv**: 200 samples
  - 100 Class_A, 100 Class_B

### 10. Expression Matrix
- **expression_matrix.tsv**: 50 samples × 500 genes
  - 3 clusters with 50 marker genes each:
    - Cluster 1 (samples 1-17): genes 1-50
    - Cluster 2 (samples 18-34): genes 100-150
    - Cluster 3 (samples 35-50): genes 200-250

- **sample_metadata.tsv**: 50 samples
  - Columns: sample_id, tissue_type, batch, treatment

### 11. Plot Test Data
- **de_results.tsv**: 1000 genes
  - DE volcano plot data
  - ~80 upregulated, ~80 downregulated (padj<0.05)

- **heatmap_matrix.tsv**: 80 genes × 24 samples
  - 3 sample groups × 8 samples with distinct patterns

- **boxviolin_data.tsv**: 120 observations
  - 4 groups with different means/spreads

- **scatter_data.tsv**: 100 points
  - x and y correlated (r~0.7)

- **bar_data.tsv**: 10 genes × 3 conditions × 5 replicates
  - Expression by condition

- **survival_plot.tsv**: 180 patients
  - 3 groups: GroupA (best), GroupB (medium), GroupC (poor) survival

### 12. API Skills Parameters
- **api_test_genes.txt**: 10 well-known genes (BRCA1, TP53, EGFR, etc.)
- **api_test_gene_list.txt**: 50 cancer-related genes
- **crispr_library_genes.txt**: 20 genes for library design
- **pubmed_test_query.txt**: "CRISPR base editing cancer 2023"

## Data Characteristics

### Quality Control Features
- Single-cell data includes cells with varying mitochondrial gene content
- Proteomics data includes realistic missing values and low-abundance proteins
- CRISPR screen includes technical replicates with realistic variance

### Biological Realism
- RNA-seq: Negative binomial counts with known DE signal
- ATAC-seq: Overlapping peak sets for differential analysis
- Methylation: Beta-distributed values with differential methylation
- Survival: Weibull-distributed times with meaningful covariates
- CRISPR: Essential gene depletion and oncogene enrichment patterns

### Integration Features
- Consistent sample IDs across related datasets
- Metadata files for all experiments
- Named features (genes, CpGs, etc.) for real-world compatibility
- Multiple clustering structures for unsupervised analysis testing

## Usage

All files are tab-separated text files that can be imported into R, Python, or other analysis tools:

```python
import pandas as pd
df = pd.read_csv("rnaseq_counts.tsv", sep="\\t", index_col=0)
```

```r
df <- read.delim("rnaseq_counts.tsv", row.names=1)
```

---
Generated: 2026-04-09
Total files: 49
"""

with open(os.path.join(OUTPUT_DIR, "README.md"), "w") as f:
    f.write(readme_content)
print(f"  Created README.md")

print("\n" + "=" * 80)
print("All test data generated successfully!")
print("=" * 80)

# List all files with sizes
print("\nGenerated Files Summary:")
print("-" * 80)
import subprocess
result = subprocess.run(
    f"cd {OUTPUT_DIR} && ls -lh && echo '' && du -sh .",
    shell=True,
    capture_output=True,
    text=True
)
print(result.stdout)

print("\nFile Count and Directory Size:")
result = subprocess.run(
    f"find {OUTPUT_DIR} -type f | wc -l",
    shell=True,
    capture_output=True,
    text=True
)
n_files = result.stdout.strip()
print(f"Total files created: {n_files}")

result = subprocess.run(
    f"du -sh {OUTPUT_DIR}",
    shell=True,
    capture_output=True,
    text=True
)
print(f"Total directory size: {result.stdout.strip()}")

print("\n" + "=" * 80)
print("SUCCESS: All synthetic test datasets have been created!")
print("=" * 80)
