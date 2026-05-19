# Figure 7

MATLAB scripts that generate the Figure 7 panels — a **Modern Coexistence Theory (MCT)** / consumer–resource model of ST1 strains vs. VPI 10463, parameterised from the Biolog substrate-utilisation data in [`../../data/Biolog/`](../../data/Biolog/).

| Script | Panel(s) | Inputs |
| --- | --- | --- |
| `MCT_Trajectories.m` | Two-species consumer–resource ODE trajectories (population density + shared/private resource concentrations over 24 h) for **ST1-75 vs VPI** and **ST1-68 vs VPI**. Substrate counts `n1`, `n2`, `ns` are from BIOLOG phenotype assyays. | *(values taken from `Biolog_growth_matrix.xlsx`)* |
| `MCT_all_strains.m` | MCT phase-plane scatter of every ST1 strain vs. VPI: **stabilizing difference** *S* = 2·n₁·n₂/(n₁+n₂) on the x-axis vs. **fitness difference** κ = (n₁ − n₂) + U·nₛ on the y-axis, with coexistence / competitive-exclusion regions shaded, and points coloured by the protection group (High / Medium / Low). | `Biolog_growth_matrix.xlsx`, `strain_groups.xlsx` |


Run from MATLAB with `data/Biolog/` 
