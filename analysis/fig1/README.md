# Figure 1 — Virulence and Protection Screen Analysis

Code to reproduce all panels of Figure 1, which characterizes the virulence and protective capacity of 21 clinical ST1 isolates of *C. difficile* against VPI10463 rechallenge in an antibiotic-treated mouse model.

---

## Figure overview

| Panel | Description | Script |
|-------|-------------|--------|
| B | Representative relative weight trajectories during primary infection (sham, ST1-75, VPI) | `fig1_analysis.m` |
| C | Virulence scores for all 21 ST1 isolates (bar plot with 95% CI) | `fig1_analysis.m` |
| D | Survival curves during secondary VPI rechallenge (ST1-75, ST1-68, VPI alone) | `secondary_challenge_analysis.m` |
| F | Representative relative weight trajectories during co-infection | `secondary_challenge_analysis.m` |
| G | Protection scores for all 21 ST1 isolates (bar plot with 95% CI) | `fig1_analysis.m` |
| H | Virulence vs protection scatter plot with linear regression | `fig1_analysis.m` |

---

## Scripts

### `fig1_analysis.m`

Main analysis script. Loads the primary infection (virulence) and co-infection (protection) datasets, fits linear mixed-effects models for both screens, combines the scores, and produces the bar plots, scatter plot, representative weight trajectories, and survival curves for selected strains.


### `secondary_challenge_analysis.m`

Dedicated script for the secondary VPI rechallenge experiment (panels D and F). Loads weight data from `weights.xlsx`, computes relative weight from a day-0 baseline per mouse, assigns group labels (ST1-75, ST1-68, VPI), and plots:


## Data files required

Place the following files at the paths specified in the scripts, or update the paths at the top of each script:

```
data/mouse/Scores/Virulence_screen_clean_table.csv   → used by fig1_analysis.m
data/mouse/Scores/ProtectionScreen_CDI_mouse.csv     → used by fig1_analysis.m
data/mouse/Scores/weights.xlsx                       → used by secondary_challenge_analysis.m
```


## Dependencies

- MATLAB R2025b or later
- Statistics and Machine Learning Toolbox 

---

## How to run

Run the two scripts in MATLAB's cell mode (Ctrl+Enter per section) or as full scripts. They are independent of each other and can be run in either order.