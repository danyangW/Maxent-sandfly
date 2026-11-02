# Maxent species distribution models with spatial cross-validation and forward selection

R code to fit **Maxent** species distribution models per species, screen predictor combinations for multicollinearity, evaluate them with **spatial cross-validation**, and quantify the effect size of each predictor (combination).

> **Status / citation:** *Wang, D., Hof, A.R., Matson, K.D. et al. Understanding and predicting the geographic distributions of phlebotomine sand flies in and around Europe. Climatic Change 178, 205 (2025). https://doi-org.library.ezproxy.wur.nl/10.1007/s10584-025-04009-z*

## What the script does

For every species listed in `sp_config`, the script:

1. **Fits single-variable Maxent models** for a set of candidate regularization multipliers, so you can inspect the response curves (`figures/reg/`) and choose a multiplier per species.
2. **Computes VIF** for each predictor combination of size `nVar` and keeps only combinations with VIF ≤ 10.
3. **Runs spatial cross-validation** (`blockCV::cv_cluster`, k = 5, fixed seed) and records **AUC** and **TSS** on each test fold.
4. **Fits the final models** on all data (jackknife and response curves enabled).
5. **Computes an effect size** per model, based on the Maxent lambda coefficients (linear, quadratic, forward/reverse hinge and product features) evaluated at the presence points.

The "forwRC" in the output folder name stands for *spatial cross-validation, forward selection of all random combinations*. Set `nVar_min` and `nVar_max` in section 3.3 to step from 1 predictor up to larger combinations (the script currently runs `nVar = 1` only).

## Repository structure

```
.
├── cv_forwRC.R              # main script (rename to match your file)
├── data/
│   └── input/
│       └── train.csv        # input data (see "Input data")
├── output/
│   └── cv_forwRC/           # created automatically
│       ├── mVIF_forwRC.csv  # VIF per predictor combination
│       ├── mSize_forwRC.csv # effect size per model
│       ├── e_forwRC.csv     # AUC / TSS per fold, all species
│       └── <species>/       # per-species models, data, evaluation
└── figures/
    └── reg/                 # response curves per species and multiplier
```

Paths are built with the [`here`](https://here.r-lib.org/) package, so run the script from within the project folder (or open the `.Rproj`).

## Requirements

- **R** ≥ 4.3
- **Java** (JDK/JRE) installed and matching your R architecture, required by `rJava`
- **`maxent.jar`** placed in the `java` folder of the `dismo` package, i.e. the path returned by `system.file("java", package = "dismo")` (download from the [Maxent website](https://biodiversityinformatics.amnh.org/open_source/maxent/))
- R packages: `tidyverse`, `dismo`, `usdm`, `rmaxent`, `blockCV` (≥ 3.0, for `cv_cluster`), `sf`, `rJava`, `here`

```r
install.packages(c("tidyverse", "dismo", "usdm", "blockCV", "sf", "rJava", "here"))

# rmaxent is not on CRAN
# install.packages("remotes")
remotes::install_github("johnbaums/rmaxent")
```

For reproducibility, consider running `renv::init()` and committing the resulting `renv.lock`.

## Input data

The script expects `data/input/train.csv` (file name set by `inputfile`) with one row per presence or background point:

| Column | Description |
|---|---|
| `species` | Species name. **`NA` for background points**, which are shared across all species |
| `preab` | 1 = presence, 0 = background/absence |
| `X3035`, `Y3035` | Coordinates in ETRS89-LAEA (EPSG:3035) |
| `year` | Year of the record |
| `*_std` | Standardized environmental predictors (see below) |

Predictor groups defined in the script:

- **Bioclimatic** (`bio1_std` … `bio19_std`, growing degree days, growing season length/temperature, etc.)
- **Land use** (`c3ann_std`, `urban_std`, `primf_std`, … from LUH2)
- **Monthly climate** (`pr_`, `tas_`, `tasmax_`, `tasmin_` for months 01–12)
- **Static** (`lith_std`)
- **Biodiversity** (vertebrate, mammal and bird richness and change metrics)

Commented-out predictors in the script were excluded (e.g. because of 10–30 % missing values). *Sources of predictor variables are documented in the original article (https://doi-org.library.ezproxy.wur.nl/10.1007/s10584-025-04009-z).*

## Usage

1. Place `train.csv` in `data/input/`.
2. Edit the configuration section (2.1 to 2.4): `inputfile`, `k`, `output_format` (default `cloglog`), and `sp_config`, which sets the regularization multiplier per species.
3. Run the script. Section 3.1 writes response curves to `figures/reg/`; inspect them and adjust `regMul` in `sp_config` before running section 3.3.
4. Results appear in `output/cv_forwRC/`.

## Outputs

| File | Content |
|---|---|
| `mVIF_forwRC.csv` | Maximum VIF per predictor combination |
| `e_forwRC.csv`, `<species>/e_<species>.csv` | AUC and TSS per test fold |
| `mSize_forwRC.csv` | Mean absolute effect size per model |
| `<species>/maList_*.rds` | Fitted Maxent models |
| `<species>/d_*.csv` | Data subset used for each species |

## Notes

- Fold assignment uses `set.seed(2024)` for reproducibility.
- Model fitting over many predictor combinations is computationally heavy. A parallel version (`parLapply`) is sketched but commented out in section 3.3.
- Result files are appended to on each run, so delete `output/cv_forwRC/` before re-running from scratch.

## License

*MIT*

## Contact

*danyang.wang@wur.nl*
