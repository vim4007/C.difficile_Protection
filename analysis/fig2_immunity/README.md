# Figure 2 — Protection is independent of adaptive immunity

## Overview
Flow cytometry analysis comparing innate and adaptive immune cell 
composition between uninfected (UI) and ST1.75-colonized mice (n=3 
per group). No significant differences were detected in any cell type 
(Wilcoxon rank-sum, all p >= 0.1), ruling out host immunity as the 
mechanism of protection.

## Scripts
- `fig2_innate_analysis.m` — innate panel (DCs, macrophages, monocytes, neutrophils)
- `fig2_adaptive_analysis.m` — adaptive panel (B cells, T cells, NK cells)
- `fig2_innate_QC.m` — innate data quality control
- `fig2_adaptive_QC.m` — adaptive data quality control

## Data
Place CSV files exported from FlowJo in `data/raw/flow_cytometry/`

## How to run
Set MATLAB working directory to `data/raw/flow_cytometry/` then run scripts

## Dependencies
MATLAB R2025b, Statistics and Machine Learning Toolbox
