# Figure 2 — Protection is independent of adaptive immunity

## Overview
Flow cytometry and mouse weight loss analysis showing that protection
conferred by ST1.75 colonization is independent of adaptive immunity.
Includes innate and adaptive immune cell composition analysis and RAG1
KO experiment.

## Scripts
- `fig2_innate_analysis.m` — innate immune panel (DCs, macrophages, monocytes, neutrophils)
- `fig2_adaptive_analysis.m` — adaptive immune panel (B cells, T cells, NK cells)
- `fig2_innate_QC.m` — innate data quality control
- `fig2_adaptive_QC.m` — adaptive data quality control
- `fig2_rag1ko_KS11.m` — RAG1 KO secondary challenge weight trajectories

## Data
- Flow cytometry CSVs: `data/raw/flow_cytometry/`
- RAG1 KO rechallenge weights: `data/raw/mouse/rag1ko/KS11_rechallenge_relweight.xlsx`

## How to run
- Flow scripts: set MATLAB working directory to `data/raw/flow_cytometry/`
- RAG1 KO script: set MATLAB working directory to `data/raw/mouse/rag1ko/`

## Dependencies
MATLAB R2025b, Statistics and Machine Learning Toolbox