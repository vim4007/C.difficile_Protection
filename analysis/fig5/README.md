# Figure 5 — Metabolic Niche Composition (BIOLOG PM1)

## Overview

Analysis of binary growth profiles for 21 ST1 *C. difficile* strains across 95 carbon sources on the PM1 BIOLOG plate, and their association with the protection phenotype. The workflow has two stages: first, raw kinetic absorbance reads from each plate are processed into a binary growth call per carbon source (grows / does not grow); second, the assembled growth matrix is analysed against the protection phenotype.

## Pipeline and scripts

### Step 1 — Per-plate growth scoring: `analyze_biolog_plate.m`

Processes a single BIOLOG PM1 plate file and returns a binary growth call for each of the 95 carbon sources. For each plate it subtracts the first (T0) reading from every well to remove baseline offset, fits a linear mixed-effects model of absorbance over time (`AbsorbanceminusT0 ~ Time*Names + (Time|Well)`), and flags a carbon source as supporting growth when the coefficient associated with metabolite and the coefficient with time-by-metabolite interaction is significantly positive after Bonferroni correction (threshold = 0.05/95). Two quality-control checks are printed automatically: the positive control (a-D-Glucose) is expected to score 1 and the negative control 0. 

### Step 2 — Batch processing and replicate handling: `run_biolog_folders.m`

Drives `analyze_biolog_plate.m` over every strain in `data/Biolog/Biolog_PM1_files/`, where each strain has its own subfolder. For a strain run once, it produces a single binary growth file. For a strain run in replicate (multiple PM1 files), it scores each replicate independently and then combines them by logical OR: a carbon source is called positive if growth was detected in any replicate, and negative only if all replicates agree on no growth. 

### Assembling the matrix: `combine_binary_growth.m`

Gathers the strain-level `<Strain>_binary_growth.xlsx` files (ignoring the per-replicate files) into a single metabolites-by-strains matrix, matched by metabolite name, and writes `data/Biolog/Biolog_growth_matrix.xlsx`. This is the input to step 3.


### Step 3 — Phenotype analysis and figures: `biolog_analysis.m`

Reads the assembled growth matrix and the phenotype annotations and produces the Figure 5 panels and statistics.

## Running the pipeline

From MATLAB, with the code folder on the path:

```matlab
run_biolog_folders();        % Step 2: score every plate, write per-strain growth files
combine_binary_growth();     % assemble data/Biolog/Biolog_growth_matrix.xlsx
biolog_analysis;             % Step 3: phenotype analysis and figures
```
## Input files

Produced by the pipeline:

* `Biolog_growth_matrix.xlsx` — binary growth matrix (95 metabolites × 21 strains), written by `combine_binary_growth.m`

Provided in the repository:

* `Biolog_PM1_files/` — raw PM1 plate reads, one subfolder per strain,
* `Biolog_names_sorted.xlsx` mapping plate wells to carbon-source names,
* `strain_groups.xlsx` — protection estimates 
* `moas.xlsx` — metabolite chemical class annotations

