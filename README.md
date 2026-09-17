# Multi-Omics Approaches Applied to Parkinson's Disease

Exploratory integration of Parkinson's Progression Markers Initiative (PPMI) molecular and clinical data to identify patient subgroups and characterize their molecular and clinical differences.

## Scope

The scripts in this repository analyze Parkinson's disease (PD) overall and three genetically defined strata:

- sporadic PD (`sPD`; participants without the mutation label used by the source phenotype table),
- GBA1-associated PD, and
- LRRK2-associated PD.

The represented molecular layers are CSF proteomics, urine proteomics, whole-blood transcriptomics, and plasma metabolomics. Depending on the cohort and the active code path, two or three layers are integrated with `iClusterBayes`. The repository also contains cluster-level differential analyses, clinical comparisons, longitudinal mixed-effects modeling, pathway enrichment, heatmaps, and UPDRS item plots.

This repository does **not** currently contain the controlled-access PPMI data or a fully portable pipeline. 

## Analysis overview

1. Select a baseline clinical cohort (all PD, sporadic PD, GBA1, or LRRK2).
2. Load and harmonize participant identifiers across available omics layers.
3. For urine proteomics and metabolomics, remove features with more than 20% missing values and mean-impute remaining missing values.
4. Log-transform metabolite abundances and select the most variable CSF/urine protein features by coefficient of variation where used.
5. Retain participants shared by the cohort and the active omics layers.
6. Fit `iClusterPlus::tune.iClusterBayes()` over candidate values of `K` and inspect BIC and deviance ratio.
7. Extract cluster assignments and features passing posterior-probability thresholds.
8. Characterize clusters using heatmaps, clinical/demographic comparisons, differential molecular analyses, pathway enrichment, UPDRS item plots, and longitudinal mixed-effects models.

## Cohort-specific integration scripts

| Analysis | Cohort selection | Active integration inputs | Main products |
|---|---|---|---|
| All PD | `COHORT == "Parkinson's Disease"` | top-2,000 variable CSF proteins + top-2,000 coding transcripts | fitted iClusterBayes workspace, BIC/deviance plot, annotated heatmap, clinical plots |
| Sporadic PD | PD with `type_mutation == "nomut"` | metabolomics + top-2,000 variable CSF proteins + top-2,000 coding transcripts | fitted workspace, model-selection plot, feature signatures, Venn diagram, heatmap, exploratory DE/clinical plots |
| GBA1 | baseline records with `type_mutation == "GBA1"` | top-2,000 variable CSF proteins + top-2,000 variable urine proteins + top-2,000 coding transcripts | fitted workspace, model-selection/posterior plots, annotated heatmap, enrichment output |
| LRRK2 | baseline records with `type_mutation == "LRRK2"` | top-2,000 variable CSF proteins + top-2,000 coding transcripts (urine is prepared but not passed to the active model) | fitted workspace, model-selection/posterior plots, annotated heatmap, enrichment output |

The scripts also contain commented alternatives for other layer combinations. These are retained as exploratory history and are not counted as executed analyses above.

## Script map

| Proposed name | Purpose |
|---|---|
| `01_integrate_pd.R` | Prepare PD multi-omics data, fit/select iClusterBayes models, annotate clusters, and create heatmaps and clinical plots. |
| `02_integrate_sporadic_pd.R` | Integrate metabolomics, CSF proteomics, and transcriptomics in sporadic PD; extract signatures and perform exploratory follow-up analyses. |
| `03_integrate_gba1.R` | Integrate CSF/urine proteomics and transcriptomics in the GBA1 stratum; visualize signatures and run enrichment/follow-up analyses. |
| `04_integrate_lrrk2.R` | Integrate CSF proteomics and transcriptomics in the LRRK2 stratum; visualize signatures and run enrichment/follow-up analyses. |
| `05_pd_differential_analysis.R` | Test cluster differences in RNA counts (DESeq2) and metabolomics/urine/CSF proteomics (limma helper). |
| `06_sporadic_pd_clinical_analysis.R` | Reconstruct sporadic-PD/healthy-control phenotypes, summarize clinical/genetic variables, and test proteomic/metabolomic contrasts. |
| `07_pd_longitudinal_models.R` | Fit participant-random-intercept longitudinal models with a time-by-cluster interaction and plot predicted trajectories. |
| `08_plot_updrs_items.R`| Generate baseline violin/box plots and pairwise comparisons for individual UPDRS-related fields. |


## Data requirements

The code expects PPMI-derived files containing:

- baseline and longitudinal phenotype data;
- curated mutation/cohort labels;
- CSF SomaScan measurements and SomaScan annotation;
- urine proteomics measurements and technical principal components;
- metabolomics measurements in long format;
- normalized expression matrices and raw RNA counts;
- levodopa equivalent daily dose (LEDD);
- CSF seed-amplification assay (SAA) status;
- selected auxiliary phenotype files, including caffeine history in the exploratory clinical script.

Participant identifiers are represented both as numeric values and as `PP-<PATNO>`. Harmonize this once during ingestion and validate row order before constructing matrices; the current scripts use a mixture of membership filters, sorting, and merges.

PPMI data are not redistributed here. Obtain data under the applicable PPMI access and data-use terms.

## Software

The scripts use R and the following packages:

- analysis: `iClusterPlus`, `DESeq2`, `limma`, `BiocParallel`, `lme4`, `lmerTest`, `ggeffects`, `enrichR`, `MASS`;
- data handling: `tidyverse`, `data.table`, `GenomicRanges`;
- visualization: `ggplot2`, `ggpubr`, `ComplexHeatmap`, `circlize`, `gplots`, `RColorBrewer`, `VennDiagram`, `ggrepel`, `lattice`, `scales`.

## Intended execution workflow

The original scripts are interactive and cannot yet be run unattended. After refactoring paths and explicit data contracts, the intended order is:

```text
cohort-specific integration (01-04)
        |
        +--> cluster assignment tables + model objects
                    |
                    +--> PD differential analysis (05)
                    +--> sporadic-PD clinical analysis (06)
                    +--> PD longitudinal models (07)
                    +--> baseline UPDRS item plots (08)
```

Each integration script should write a small cluster-assignment table (at minimum `PATNO`, cohort, and cluster), feature-signature tables, and a serialized model object. 

## Data privacy and citation

Do not upload PPMI participant-level data, credentials, or local paths. Add the PPMI acknowledgment and citation language required by the data-use agreement, along with citations for `iClusterPlus`/iClusterBayes and the statistical packages used. Add an explicit open-source license only after confirming that the code and any adapted helpers may be redistributed.
