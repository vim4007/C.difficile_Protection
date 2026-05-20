# Figure 2 — Protection by ST1-75 is Independent of Adaptive Immunity


## Contents

| File | Description |
|------|-------------|
| `fig2_innate_analysis.m` | Loads innate flow cytometry data, runs PCA and tSNE, computes cell type fractions per mouse, and generates bar plots comparing UI vs ST1-75 colonized mice (DCs, CD11b+, Lymphoid DCs, Macrophages, Monocytes, Myeloid DCs, Neutrophils) |
| `fig2_adaptive_analysis.m` | Same pipeline as innate script applied to adaptive immune populations (B cells, CD4+, CD8+, T cells, NK cells, CD19-TCRβ-) |
| `fig2_rag1ko_KS11.m` | Plots relative weight trajectories following VPI10463 secondary challenge in PBS/B6, ST1-75/B6, and ST1-75/RAG1-KO groups |

---

## Data

```
data/
├── flow_cytometry/
│   ├── ks10_innate_csv_files/       ← used by fig2_innate_analysis.m
│   └── ks10_adaptive_csv_files/     ← used by fig2_adaptive_analysis.m
└── mouse/
    └── rag1ko/
        └── KS11_rechallenge_relweight.xlsx   ← used by fig2_rag1ko_KS11.m
```

---

## Requirements

- MATLAB (tested with R2025b or later)