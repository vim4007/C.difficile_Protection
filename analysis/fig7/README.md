# Figure 7 — Consumer-Resource Model of *C. difficile* Strain Competition

Figure 7 presents a consumer-resource (CR) model in which two *C. difficile* strains compete for shared and private metabolic carbon sources. The model predicts competitive outcomes (exclusion vs. coexistence) based on Biolog PM1 niche breadth data, and is validated against qPCR colonization.

---

## Contents

| File | Description |
|------|-------------|
| `MCT_Trajectories.m` | Simulates population and resource dynamics over time for ST1.75 vs VPI10463 and ST1.68 vs VPI10463 (panels A–D) |
| `Phase_portrait.m` | Generates phase-plane portraits with nullclines, vector fields, steady states, and trajectories for each strain pair (panel E) |
| `MCT_all_strains.m` | Loads Biolog data, computes fitness difference (κ) and stabilizing difference (S) for all ST1 strains vs VPI, and plots the MCT coexistence phase plane colored by protection group (panel I) |
| `qPCR.m` | Plots relative weight trajectories (panel H) and stacked bar plots of ST1-75 vs VPI colonization proportions by qPCR over time (panel G) |

---

## Data

```
data/
├── Biolog/
│   ├── Biolog_growth_matrix.xlsx    ← used by MCT_all_strains.m
│   └── strain_groups.xlsx           ← used by MCT_all_strains.m
└── qPCR/
    ├── weights.xlsx                 ← used by qPCR.m
    └── qPCR section.xlsx            ← used by qPCR.m
```

`MCT_Trajectories.m` and `Phase_portrait.m` are pure mathematical models with no data dependencies.

