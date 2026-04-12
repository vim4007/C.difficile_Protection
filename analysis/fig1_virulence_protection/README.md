# Figure 1 — Virulence and protection phenotypes across ST1 C. difficile isolates

## What this script does
Fits linear mixed-effects models to quantify virulence and protection 
scores across 21 clinical MLST1 C. difficile isolates, generates the 
virulence vs. protection correlation scatter plot, weight loss 
trajectories, bar plots of ranked scores, and Kaplan-Meier survival 
curves.

## Input files
Place these in `data/raw/mouse/` before running:
- `Virulence_screen_clean_table.csv` — weight loss data from primary 
  infection screen (21 ST1 strains)
- `ProtectionScreen_CDI_mouse.csv` — weight loss data from co-infection 
  and secondary rechallenge experiments

## Output figures
- Scatter plot: virulence vs. protection effect sizes with linear fit 
  (R² = 0.51, p = 0.00019)
- Scatter plot with bidirectional 95% CIs
- Bar plots: virulence and protection scores ranked across all strains
- Weight loss trajectory: mean ± SD for surviving mice, individual 
  traces for mice that died
- Kaplan-Meier survival curves

## How to run
1. Open MATLAB R2025b
2. Set working directory to `data/raw/mouse/`
3. Open and run `fig1_analysis.m`

## Dependencies
MATLAB R2025b  
Statistics and Machine Learning Toolbox (required for fitlme, fitlm, ecdf)

## Notes
- The trajectory plot is hardcoded to one strain at a time via the 
  `strain_name` variable (line 169). Change this variable to plot a 
  different strain.
- The mixed-effects model uses PBS sham-infected mice as the reference 
  category for both virulence and protection models.


## Secondary rechallenge trajectories

**Script:** `secondary_challenge_analysis.m`  
**Input:** `data/raw/mouse/weights.xlsx`  
**What it does:** Plots weight loss trajectories for ST1-75, ST1-68, 
and VPI10463 during the secondary rechallenge experiment. Shows mean 
± SD with individual traces and red X markers for mice that died.  
**How to run:** Set MATLAB working directory to `data/raw/mouse/` 
then run the script.
