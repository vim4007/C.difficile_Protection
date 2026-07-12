# Data

This directory contains the raw and processed datasets for the *C. difficile* protection project.

## Directory layout

```
data/
├── 16s_sequencing/                 # 16S rRNA amplicon abundance + matched host weights
├── Biolog/                         # Biolog PM substrate utilization (C-source) panels
│   ├── Biolog_PM1_files/           # raw PM1 plate reads, one subfolder per strain
│   │   ├── <Strain>/               # e.g. ST1-6/ (single run), ST1-75/ and VPI/ (replicate runs _1/_2/_3)
│   │   └── Biolog_names_sorted.xlsx  # maps plate wells (A1..H12) to carbon-source names
│   ├── Biolog_growth_matrix.xlsx   # binary growth matrix (metabolites × strains), generated from the raw reads
│   ├── strain_groups.xlsx          # protection estimates and group assignments (High/Medium/Low)
│   └── moas.xlsx                   # metabolite chemical class (mechanism-of-action) annotations
├── Genomics/                       # Pan-genome presence/absence + phylogenetic trees
├── flow_cytometry/                 # Gated flow cytometry populations (innate + adaptive)
├── mouse/                          # in vivo CDI screens and host readouts
└── qPCR/                           # Strain-specific qPCR quantification + mouse weights
```

Within `Biolog/Biolog_PM1_files/`, each strain has its own subfolder of kinetic absorbance reads: strains run once contain a single PM1 file, strains run in replicate contain numbered files (`_1`, `_2`, `_3`). The processing pipeline writes a per-strain binary growth call back into each subfolder as `<Strain>_binary_growth.xlsx`, combines replicates by logical OR, and assembles the strain-level calls into `Biolog_growth_matrix.xlsx`. That matrix is a generated file and can be rebuilt from the raw reads at any time; `strain_groups.xlsx`, `moas.xlsx`, and `Biolog_names_sorted.xlsx` are provided inputs.
