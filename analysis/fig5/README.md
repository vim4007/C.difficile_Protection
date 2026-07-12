# Figure 5 — Metabolic Niche Composition (BIOLOG PM1)

## Overview

Analysis of binary growth profiles for 21 ST1 *C. difficile* strains across 95 carbon sources on the PM1 BIOLOG plate, and their association with the protection phenotype. The workflow has two stages: first, raw kinetic absorbance reads from each plate are processed into a binary growth call per carbon source (grows / does not grow); second, the assembled growth matrix is analysed against the protection phenotype.

All scripts locate their data relative to the repository, so no paths need editing. A user who clones the repository can run the full pipeline unchanged.

## Pipeline and scripts

The analysis runs in three steps, each with its own script. Steps 1 and 2 build the growth matrix from raw plate data; step 3 is the phenotype analysis that produces the figures.

### Step 1 — Per-plate growth scoring: `analyze_biolog_plate.m`

Processes a single BIOLOG PM1 plate file and returns a binary growth call for each of the 95 carbon sources. For each plate it subtracts the first (T0) reading from every well to remove baseline offset, fits a linear mixed-effects model of absorbance over time (`AbsorbanceminusT0 ~ Time*Names + (Time|Well)`), and flags a carbon source as supporting growth when its time-by-metabolite interaction is significantly positive after Bonferroni correction (threshold = 0.05/95). Two quality-control checks are printed automatically: the positive control (a-D-Glucose) is expected to score 1 and the negative control 0. This function is normally driven by step 2 rather than called directly.

### Step 2 — Batch processing and replicate handling: `run_biolog_folders.m`

Drives `analyze_biolog_plate.m` over every strain in `data/Biolog/Biolog_PM1_files/`, where each strain has its own subfolder. For a strain run once, it produces a single binary growth file. For a strain run in replicate (multiple PM1 files), it scores each replicate independently and then combines them by logical OR: a carbon source is called positive if growth was detected in any replicate, and negative only if all replicates agree on no growth. Each strain's result is written back into its own subfolder as `<Strain>_binary_growth.xlsx`, with per-replicate files named `<Strain>_binary_growth_<n>.xlsx`.

### Assembling the matrix: `combine_binary_growth.m`

Gathers the strain-level `<Strain>_binary_growth.xlsx` files (ignoring the per-replicate files) into a single metabolites-by-strains matrix, matched by metabolite name, and writes `data/Biolog/Biolog_growth_matrix.xlsx`. This is the input to step 3.

### Helper: `biologDataDir.m`

A small utility used by the scripts above to locate `data/Biolog` relative to the code, so the pipeline is portable across machines with no hardcoded paths.

### Step 3 — Phenotype analysis and figures: `biolog_analysis.m`

Reads the assembled growth matrix and the phenotype annotations and produces the Figure 5 panels and statistics.

## Running the pipeline

From MATLAB, with the code folder on the path:

```matlab
run_biolog_folders();        % Step 2: score every plate, write per-strain growth files
combine_binary_growth();     % assemble data/Biolog/Biolog_growth_matrix.xlsx
biolog_analysis;             % Step 3: phenotype analysis and figures
```

Steps 1 and 2 fit a mixed-effects model per plate and are relatively slow. Because the growth matrix is committed to the repository, this regeneration is only needed to rebuild it from raw plates; to reproduce the figures alone, run `biolog_analysis` directly.

## Input files

Produced by the pipeline:

* `Biolog_growth_matrix.xlsx` — binary growth matrix (95 metabolites × 21 strains), written by `combine_binary_growth.m`

Provided in the repository:

* `Biolog_PM1_files/` — raw PM1 plate reads, one subfolder per strain, plus `Biolog_names_sorted.xlsx` mapping plate wells to carbon-source names
* `strain_groups.xlsx` — protection estimates and group assignments (High / Medium / Low)
* `moas.xlsx` — metabolite chemical class annotations

