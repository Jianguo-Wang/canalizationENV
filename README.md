# canalizationENV

### A self-similar equation for the morphogenetic reproducibility of *C. elegans*

> Waddington's canalization landscape, made quantitative: mother–daughter positional noise is suppressed by negative feedback that acts through cell migration velocity, and the resulting velocity oscillations follow a one-parameter self-similar equation whose parameter is the golden ratio conjugate, **φ = 0.618**.

This repository holds the single-cell spatiotemporal data and the complete R analysis code for the study *"A self-similar equation for the morphogenetic reproducibility of C. elegans"*. It tracks **2,039 embryos** — wild type, environmentally perturbed and single-gene RNAi knockdown — from the 4-cell to the 350-cell stage, and contains everything needed to reproduce the parameter estimates, the figures and the gene-discovery analysis.

**Preprint:** [*A kinematic equation for the morphogenetic reproducibility of an animal*](https://doi.org/10.65215/ltspreprints.2026.04.04.000174) — 浪淘沙预印本平台 (Langtaosha Preprint Platform), 4 April 2026. Registered DOI: `10.65215/ltspreprints.2026.04.04.000174`.

---

## Contents

- [Key findings](#key-findings)
- [Study at a glance](#study-at-a-glance)
- [Repository layout](#repository-layout)
- [Code](#code)
- [Data](#data)
- [Getting started](#getting-started)
- [Reproducing the analysis](#reproducing-the-analysis)
- [Notes and conventions](#notes-and-conventions)
- [Data availability](#data-availability)
- [Citation and contact](#citation-and-contact)
- [Funding and declarations](#funding-and-declarations)

---

## Key findings

- **Ubiquitous mother–daughter negative feedback.** Across 2,039 embryos and three coordinates, the noise-transmission slope `k` is almost always negative (median `kx = −0.61`, `ky = −0.66`, `kz = −0.74`; 706 mother–daughter pairs per embryo). Daughter cells inherit positional noise only partially.
- **Feedback is implemented through migration velocity, not timing.** Daughter migration displacement is tuned by velocity, while cell cycle length is tightly conserved — an *isochronous* regulatory design.
- **Velocity oscillates across generations.** Cell velocity alternates in sign along every lineage, and negative correlation with ancestors decays geometrically with generation distance.
- **A one-parameter self-similar equation describes the oscillations**, obtained by adapting the geometric distributed lag model:

  ```math
  V_t = -\varphi \sum_{i=1}^{t-1} \varphi^{\,i}\, V_{t-i} + \varepsilon_t
  ```

  with a single velocity-adjustment parameter φ. Estimated φ = **0.619 ± 0.001** (control), **0.614 ± 0.0003** (RNAi) and **0.603 ± 0.001** (environmental), i.e. the golden ratio conjugate **0.618**.
- **φ = 0.618 is a theoretical optimum.** When α = β = φ the cumulative feedback strength `−αβ/(1−β)` converges to `−1`, meaning exactly 100 % correction of any ancestral deviation and zero net noise accumulation. Smaller values under-correct, larger values over-correct.
- **Violating the equation predicts lethality and names genes.** Per-embryo deviation (mean absolute error of predicted grand-daughter velocities) predicts hatching failure with **AUC = 0.89** (0.55 when outcomes are shuffled), and identifies **24 candidate genes** whose knockdown destabilises the equation — significantly enriched for cytoskeleton and kinesin motor functions.

## Study at a glance

| | |
|---|---|
| Embryos analysed | **2,039** — RNAi 1,823, environmental 111, control 105 |
| Developmental window | 4-cell → 350-cell stage (~9 of 10 embryonic divisions) |
| Lineage coverage | 714 cells, 706 mother–daughter pairs per embryo |
| Perturbations | 40 environmental conditions in 8 stressor classes (temperature, heat shock, bacterial food, stearic acid, glucose, NaCl, FeCl₃, paraquat); ~700 genes by RNAi |
| Imaging | spinning-disk confocal, 75 s intervals, 30 z-planes at 1 µm, 4.75–6 h at 21 °C |
| Lineage reconstruction | StarryNite + AceTree manual curation; tracing accuracy 99.9 % |
| Central estimate | φ = 0.619 ± 0.001 (control) |
| Predictive performance | hatching-lethality AUC = 0.89 |

## Repository layout

```text
canalizationENV/
├── README.md
├── Code/                            53 R scripts
│   ├── raw-1 … raw-3                data preprocessing
│   ├── 1 … 12                       derived variables
│   ├── plot-1 … plot-26             figures
│   └── Mechanism-1 … Mechanism-12   kinematic mechanism
└── Data/                            2,313 files
    ├── Raw_data_Control/            108 files
    │   ├── spatial/                 105 per-embryo tracking files
    │   └── temporal/                3 cross-embryo matrices
    ├── Raw_data_Env/                124 files
    │   ├── spatial/                 120 per-embryo tracking files
    │   ├── temporal/                3 cross-embryo matrices
    │   └── 条件说明.docx            description of the environmental conditions
    └── Raw_data_RNAi/               2,081 files
        ├── spatial/                 2,075 per-embryo tracking files, 698 RNAi conditions
        │   └── A/ … Z/              grouped by the first character of the gene name
        └── temporal/                6 summary tables and 1 R script
```

## Code

All analysis code is in `Code/` and is written in R. Scripts are numbered in the order they should be run.

| Group | Scripts | Purpose |
|---|---|---|
| **A. Data preprocessing** | `raw-1`, `raw-2`, `raw-3` | Build per-embryo 3D coordinates from the StarryNite lineage output, align cells across embryos (translation–rotation–scaling), and assemble cell cycle length (CCL) and cell division timing (CDT) tables |
| **B. Basic variables** | `1` … `12` | Generate the derived quantities used throughout: the 714-cell reference set, mother–daughter noise and feedback coefficients, velocity/time/displacement kinematics, canal depth, lethality prediction and gene identification |
| **C. Visualization** | `plot-1` … `plot-26` | Reproduce the main and supplementary figures — data structure, alignment quality, spatial and temporal negative feedback, variance accumulation, scaling factor and dilution-effect controls, prediction ROC curves, gene enrichment |
| **D. Kinematic mechanism** | `Mechanism-1` … `Mechanism-12` | The self-similar velocity equation and its analysis: theoretical curves, φ estimation, self-similarity analysis, damped-system theory, natural-selection simulation, robustness to sample size, the mother–sister paradox and gene identification |

<details>
<summary><b>Full script list</b></summary>

**A. Data preprocessing**
`raw-1-embryo_coordinates_from_raw.R`, `raw-2-embryo_coordinates_alignment.R`, `raw-3-embryo_CCL_CDT_from_raw_&_alignment.R`

**B. Basic variables**
`1-lineage_info_short.R`, `2-mother_daughter1_daughter2.R`, `3-cell_fate.R`, `4-lineage_time_series_model_spatial.R`, `5-real-time_negative_feedback.R`, `6-spatial_subset_cell_association.R`, `7-embryonic_lethal_prediction.R`, `8-depth_of_canals.R`, `9-lineage_time_series_model_temporal.R`, `10-gene_identification_CCL.R`, `11-embryonic_lethal_prediction_wormbase.R`, `12-gene_identification_xyz.R`

**C. Visualization**
`plot-1-data_structure.R`, `plot-2-negative_feedback_concept.R`, `plot-3-spatial_negative_feedback.R`, `plot-4-effect_of_generation.R`, `plot-5-effect_of_spatial_distribution.R`, `plot-6-effect_of_differentiation.R`, `plot-7-effect_of_scaling_factor.R`, `plot-8-effect_of_dilution_effect.R`, `plot-9-effect_of_tracking_error.R`, `plot-10-effect_of_alignment_process.R`, `plot-11-real-time_negative_feedback.R`, `plot-12-variance_accumulation.R`, `plot-13-embryo_develop_failure_withoutNF.R`, `plot-14-spatial_subset_cell_association.R`, `plot-15-embryonic_lethal_prediction.R`, `plot-16-spatial_high_dimensional_stability.R`, `plot-17-temporal_negative_feedback_high_dimensional_stability.R`, `plot-18-effect_mo_da_cv_sd.R`, `plot-19-common_expecation_for_xyz_CCL.R`, `plot-20-effect_of_tracking_error_on_CCL.R`, `plot-21-spatial_alignment_process_and_quality.R`, `plot-22-temporal_alignment_process_and_quality.R`, `plot-23-gene_identification.R`, `plot-24-phenotype_prediction_wormbase.R`, `plot-25-phenotype_prediction_CCL.R`, `plot-26-simulate_if_noise_accumulation_rather_than_negative_feedback.R`

**D. Kinematic mechanism**
`Mechanism-1-theoretical_curve_displacement_velocity_time-final.R`, `Mechanism-2-velocity_gene_identification_PPI.R`, `Mechanism-3-self_similarity_analysis.R`, `Mechanism-4-simulation_for_phi_based_self_similarity_model.R`, `Mechanism-5-approximation_of_three_generation_model.R`, `Mechanism-6-simulation_for_natural_selectin_of_phi.R`, `Mechanism-7-phi_estimation_robustness_for_varying_sample_size.R`, `Mechanism-8-damped_system.R`, `Mechanism-9-working_model.R`, `Mechanism-10-mother_sister_paradox.R`, `Mechanism-11-gene_identification.R`, `Mechanism-12-gene_identification_shuffle.R`

</details>

## Data

Raw data are organised by perturbation type. Each condition folder has the same two subfolders.

### `spatial/` — per-embryo tracking files

One tab-separated file per embryo, produced by the lineage-tracking pipeline. File naming:

| Folder | Naming | Files |
|---|---|---|
| `Raw_data_Control/spatial/` | `ctr_emb<N>_raw_data.txt` | 105 |
| `Raw_data_Env/spatial/` | `XL_RW<slide>_<temp>C_<date>_<session>_emb<N>_raw_data.txt` | 120 |
| `Raw_data_RNAi/spatial/<letter>/` | `<gene>_emb<N>_raw_data.txt` | 2,075 |

Columns of each file:

| Column | Content |
|---|---|
| `time` | imaging time index of the frame |
| `cell_name` | cell identity in Sulston nomenclature (e.g. `ABal`, `EMS`) |
| `X`, `Y`, `Z` | 3D position of the nucleus, in raw tracking coordinates |
| `size` | nucleus size |
| `raw-expression` | histone::mCherry signal intensity |
| `blot-correction` | photobleaching-corrected intensity |

Coordinates are raw outputs of the tracker; alignment, rotation and scaling are performed by `raw-2` and `raw-3`.

### `temporal/` — cross-embryo summary tables

| File | Content |
|---|---|
| `CCL_control.csv`, `CCL_Env.csv`, `CCL_RNAi.csv` | Cell cycle length matrix — rows are cell identities, columns are embryo names, values in minutes |
| `cdt_control.csv`, `cdt_Env.csv`, `cdt_RNAi.csv` | Cell division timing matrix, same orientation |
| `cell_division_time_wt.csv` | Division times of the wild-type reference embryos |
| `Cell_division_timing_matrix_environment.txt` | Division timing matrix for the environmental series, used for cross-checking |
| `Cell_division_timing_matrix_E03H4.8_and_NED-8_refine.csv` | Refined division timing for the E03H4.8 and NED-8 RNAi series |
| `cell_division_time_mt_revise.csv` | Revised division-time table for the RNAi series |
| `cdt_mt_revise_code.R` | Script that produced the revised table |
| `条件说明.docx` | Description of the 40 environmental conditions (Chinese) |

## Getting started

```bash
git clone https://github.com/Jianguo-Wang/canalizationENV.git
cd canalizationENV
```

The code requires **R ≥ 4.0** (it uses raw string literals) and a set of CRAN and Bioconductor packages.

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c(
  "clusterProfiler", "org.Ce.eg.db", "org.Hs.eg.db", "AnnotationDbi",
  "enrichplot", "ReactomePA", "pathview", "biomaRt", "limma", "ggtree", "treeio"
))
install.packages(c(
  "ggplot2", "dplyr", "tidyr", "readr", "reshape", "reshape2", "gplots", "pheatmap",
  "glmnet", "pROC", "caret", "yardstick", "plotROC", "cowplot", "gridExtra",
  "patchwork", "scales", "stringr", "openxlsx",
  "igraph", "ggraph", "tidygraph", "ggVennDiagram",
  "lavaan", "semPlot", "bootnet", "corpcor", "mnormt", "apcluster",
  "concaveman", "doParallel", "foreach", "rgl", "scatterplot3d", "extrafont",
  "bayesmeta", "MASS", "cluster"
))
```

### Configure the project path

Most scripts derive every input and output path from a single variable declared at the top of the file:

```r
focal_projectFolderPath <- r"(E:\canalizationENV)"
```

Change it to the location of your clone before running anything — 44 of the 53 scripts call `setwd()` relative to this root. All other paths are expressed relative to it.

## Reproducing the analysis

Run the scripts in the following order; each stage writes the intermediate products the next stage expects under `data/`.

1. **Preprocessing** — `raw-1` builds per-embryo coordinate tables from the tracker output in `Data/*/spatial/`; `raw-2` performs the translation–rotation–scaling alignment across embryos (reference embryo `ctr-emb58`); `raw-3` assembles CCL and CDT matrices.
2. **Derived variables** — run `1-…` through `12-…` in numerical order to obtain the 714-cell reference set, mother–daughter noise and feedback coefficients, kinematic quantities and predictions.
3. **Figures and mechanism** — run the `plot-*` scripts for the main and supplementary figures, and the `Mechanism-*` scripts for the self-similar equation, φ estimation, simulations and gene identification.

Intermediate products are read from and written to subfolders of `data/`:

| Folder | References in code |
|---|---|
| `data/basic/` | 176 |
| `data/spatial/` | 166 |
| `data/alignment/` | 31 |
| `data/subsetDiff/` | 30 |
| `data/prediction/` | 27 |
| `data/realTime/` | 24 |
| `data/temporal/` | 23 |
| `data/interpolation/` | 15 |
| `data/velocity/` | 5 |
| `data/depth/` | 3 |
| `data/annotation/` | 1 |

> **These intermediates are not included in this repository** — they run to several gigabytes. They are regenerated by running the pipeline in order, and are available from the authors on request.

## Notes and conventions

- **File counts vs. analysed embryos.** The manuscript analyses 2,039 embryos (1,823 RNAi, 111 environmental, 105 control). `Data/Raw_data_RNAi/spatial/` contains 2,075 per-embryo files covering 698 conditions, i.e. the full tracking set, which includes embryos that did not enter the final analysis. Likewise `Raw_data_Env/spatial/` holds all 120 imaged environmental embryos, of which 111 were lineage-traced to the 350-cell stage.
- **RNAi subfolders.** GitHub's file browser truncates directory listings at 1,000 entries, and the flat RNAi folder holds 2,075 files. Files are therefore grouped into `spatial/A/` … `spatial/Z/` by the first character of the gene name. This is a presentation choice only: the file contents are unchanged, and the largest subfolder holds 234 files.
- **Line endings.** The `.txt` and `.csv` files are stored in the repository with LF line endings (git normalised the original Windows CRLF endings on commit). Numeric content is unaffected; a checkout on Windows restores CRLF automatically.
- **Repository size.** The data expand to about 2.4 GB on disk but the repository itself is only ~64 MB compressed, because the tracked files are text.
- **Spatial folders are flat apart from the RNAi split.** `Raw_data_Control/spatial/` and `Raw_data_Env/spatial/` have fewer than 1,000 files each and are not subdivided.

## Data availability

The raw data, processed data and code for this study are available in this repository and accompany the preprint:

> Wang, J., Xiao, L., Yang, M., Yao, Z., Deng, S., Du, Z., & He, X. (2026). *A kinematic equation for the morphogenetic reproducibility of an animal* [Preprint]. 浪淘沙预印本平台 (Langtaosha Preprint Platform). https://doi.org/10.65215/ltspreprints.2026.04.04.000174

The spatiotemporal data of the RNAi and control embryos were originally generated for our previous study and are described there:

> L. Xiao, D. Fan, H. Qi, Y. Cong, Z. Du, Defect-buffering cellular plasticity increases robustness of metazoan embryogenesis. *Cell Systems* **13**, 615–630.e9 (2022). https://doi.org/10.1016/j.cels.2022.07.001

## Citation and contact

If you use these data or code, please cite the preprint:

> Wang, J., Xiao, L., Yang, M., Yao, Z., Deng, S., Du, Z., & He, X. (2026). *A kinematic equation for the morphogenetic reproducibility of an animal* [Preprint]. 浪淘沙预印本平台 (Langtaosha Preprint Platform). https://doi.org/10.65215/ltspreprints.2026.04.04.000174

```bibtex
@misc{wang2026kinematic,
  title        = {A kinematic equation for the morphogenetic reproducibility of an animal},
  author       = {Wang, Jianguo and Xiao, Long and Yang, Miaoling and Yao, Zeqi and
                  Deng, Shanjun and Du, Zhuo and He, Xionglei},
  year         = {2026},
  month        = apr,
  howpublished = {Preprint, 浪淘沙预印本平台 (Langtaosha Preprint Platform)},
  doi          = {10.65215/ltspreprints.2026.04.04.000174},
  url          = {https://doi.org/10.65215/ltspreprints.2026.04.04.000174}
}
```

A revised manuscript, *A self-similar equation for the morphogenetic reproducibility of C. elegans* (August 2026), is in preparation; please cite the preprint until the journal version appears.

**Corresponding authors**

- Xionglei He — hexiongl@mail.sysu.edu.cn
- Zhuo Du — zdu@genetics.ac.cn

**Affiliations**

1. Hongkong Institute for Advanced Studies, School of Life Sciences, Sun Yat-sen University, Guangzhou 510275, China
2. State Key Laboratory of Molecular Developmental Biology, Institute of Genetics and Developmental Biology, Chinese Academy of Sciences, Beijing 100101, China

Jianguo Wang and Long Xiao contributed equally to this work.

**Licence.** The preprint is released under [CC BY-NC-ND 4.0](https://creativecommons.org/licenses/by-nc-nd/4.0/). No separate licence has been assigned to the data and code in this repository yet; please contact the corresponding authors before redistributing or reusing them.

## Funding and declarations

This work was supported by the National Natural Science Foundation of China (grant nos. 32293191 and 32400497) and the GuangDong Basic and Applied Basic Research Foundation (grant no. 2023A1515110669). The authors declare no competing interests.
