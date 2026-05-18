# Data Directory

Only small benchmark resources are tracked in git.

Tracked:

- `resources/human_tfs.tsv`: predefined transcription-factor gene set used by
  the TF feature-selection baseline.

Not tracked:

- raw h5ad files
- processed h5ad files
- intermediate matrices
- downloaded spatial transcriptomics archives

The benchmark configs expect local dataset paths under `data/raw/`. Recreate
those files from the original sources before running the full benchmark.

Known dataset sources used during this project include:

- DLPFC: spatialLIBD / Visium DLPFC slices
- MouseBrainSerialSections: 10x Genomics mouse brain serial-section Visium data
- STOmics0212: STOmicsDB STDS0000212 processed h5ad
- STOmics0218: STOmicsDB STDS0000218 processed h5ad
- STOmics0224: STOmicsDB STDS0000224 processed h5ad
- STOmicsVisium5Samples: local multi-sample Visium/STOmics h5ad bundle
- E8p5Embryo and E9p5Embryo: Slide-seq embryo datasets

Dataset-level metadata used for figure generation is available in
`results/current_rank/data/datasets-metadata.tsv`.
