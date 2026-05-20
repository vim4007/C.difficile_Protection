# Figure 6 — GC-MS Metabolic Profiling of *C. difficile* Strains

Figure 6 compares intracellular and extracellular (secretome) metabolite profiles of ST1.75, ST1.68, and VPI10463 strains using GC-MS, identifying metabolites that are differentially accumulated or secreted across strains.

---

## Contents

| File | Description |
|------|-------------|
| `gcms_figures.m` | Loads GC-MS data, normalizes and log2-transforms metabolite abundances, runs PCA on intracellular and extracellular profiles, and generates dot plots of log2 fold change vs. media for significantly different metabolites (Welch's t-test) |

---

## Data

```
data/
└── GCMS/
    ├── secretome.xlsx    ← extracellular metabolites (used by gcms_figures.m)
    └── intra.xlsx        ← intracellular metabolites (used by gcms_figures.m)
```
