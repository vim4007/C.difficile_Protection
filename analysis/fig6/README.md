# Figure 6

MATLAB scripts that generate the panels for Figure 6 — mechanistic comparison of the protective strain **ST1-75** against **ST1-68** and the virulent **VPI 10463** reference.

| Script | Panel(s) | Data source |
| --- | --- | --- |
| `comp_manus.m` | Cerillo OD600 competition curves (mean ± SD) for ST1-75 vs VPI and ST1-68 vs VPI on alanine and glycine as sole carbon sources. | `data/cerillo/competition.xlsx`  |
| `gcms_figures.m` | Secretome and intracellular GC-MS dotplots: log₂ fold-change vs. media for each strain, with pairwise Welch's *t*-tests and significance stars. | `data/GCMS/secretome.xlsx`, `data/GCMS/intra.xlsx` |
| `qPCR.m` | Weight-trajectory plot + stacked bar plots of ST1-75 vs VPI relative abundance (mono- and co-colonized mice) at days 0.5, 1, and 3. | `data/qPCR/weights.xlsx`, `data/qPCR/qPCR section.xlsx` |

Run each script from MATLAB with `data/` as described above
