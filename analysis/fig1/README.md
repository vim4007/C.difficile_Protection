# Figure 1 — Virulence and Protection Screening of *C. difficile* Strains

## Overview

This folder contains the MATLAB analysis scripts and figure file for Figure 1, which characterizes a panel of *C. difficile* strains based on their **virulence** and **protection** estimates in a mouse co-infection model. 
---

## Contents

| File | Description |
|------|-------------|
| `fig1_analysis.m` | Primary analysis script: computes virulence and protection scores, generates scatter plot, bar charts, weight trajectories, and survival curves |
| `secondary_challnege_analysis.m` | Secondary challenge analysis: plots relative weight trajectories and survival curves for ST1.75, ST1.68, and VPI groups over a re-challenge experiment |
| `Figure_1.ai` | Final assembled figure (Adobe Illustrator) |

---

## Data Dependencies

These scripts read from the following files (relative to `/Users/vmishra/C.difficile_Protection/`):

| Data File | Used In |
|-----------|---------|
| `data/mouse/Scores/ProtectionScreen_CDI_mouse.csv` | `fig1_analysis.m` |
| `data/mouse/Scores/Virulence_screen_clean_table.csv` | `fig1_analysis.m` |
| `data/mouse/Scores/weights.xlsx` | `secondary_challnege_analysis.m` |

---

## Requirements

- MATLAB (tested with R2025b or later)
- Statistics and Machine Learning Toolbox (for `fitlme`, `ecdf`)
