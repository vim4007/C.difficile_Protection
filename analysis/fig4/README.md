# Figure 4 — Genomic contributions to the protection score

## Overview
Pan-genome and phylogenetic analysis of 21 ST1 C. difficile strains
showing that neither core genome phylogeny nor accessory gene content
explains the protection phenotype. Identifies clcA_1 as nominally
enriched in highly protective strains.

## Scripts
- `genomic_analysis.m` — pan-genome pie chart, PCA of accessory
  genes, clcA_1 presence/absence, waterfall plot of gene correlations
  with protection score
- `phylogeny.m` — core genome phylogeny rerooted at VPI
  outgroup and colored by protection score

## Data
Place in `data/raw/genomics/`:
- `binary_gene_pre_abs.csv` — binary gene presence/absence matrix
  (21 strains × 5109 genes) with virulence and protection scores
- `strain_names_mapping.csv` — maps numeric genome IDs to strain names
- `VPI Outgroup_tree.nwk` — core genome newick tree with VPI outgroup

## How to run
Set MATLAB working directory to `data/raw/genomics/` then run scripts

## Dependencies
MATLAB R2025b, Statistics and Machine Learning Toolbox,
Bioinformatics Toolbox (required for phytreeread, reroot, prune)