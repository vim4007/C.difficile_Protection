# Figure 5 — Metabolic Niche Composition (BIOLOG PM1)

## Overview
Analysis of binary growth profiles for 21 ST1 *C. difficile* strains across 95 carbon 
sources on the PM1 BIOLOG plate, and their association with the protection phenotype.

## Script
`biolog_analysis.m`

## Required Input Files
- `Biolog_growth_matrix.xlsx` — binary growth matrix (95 metabolites × 21 strains)
- `strain_groups.xlsx` — protection estimates and group assignments (High/Medium/Low)
- `moas.xlsx` — metabolite chemical class annotations

## Analyses
1. **Wilcoxon rank-sum test** — per-metabolite association with protection score; 
   Bonferroni correction applied (threshold = 0.05/95)
2. **Spearman correlation by chemical class** — growth count per class vs protection estimate
3. **Random Forest classifier** — LOO-CV prediction of protection group (High/Medium/Low); 
   200 trees, accuracy = 57%
4. **Lasso regression** — feature selection on continuous protection score; 
   identifies glyoxylic acid as top discriminating metabolite (coef = 3.33)

## Figures Generated
- Chemical class associations with protection (waterfall, Spearman ρ)
- Stacked strain contributions per chemical class
- PCA of growth profiles colored by protection estimate
- Per-strain metabolic niche composition sorted by protection + linear fit

## Key Result
 Protection is positively correlated 
with growth on amino acids (ρ = 0.6, p = 0.0034).