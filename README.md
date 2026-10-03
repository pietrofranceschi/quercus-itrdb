# Quercus ring-width correlation vs distance (ITRDB)

Analysis of how the correlation between detrended tree-ring width series of
*Quercus petraea* (QUPE) and *Quercus robur* (QURO) decays with the geographic
distance between sites, and how that relationship can be inverted to estimate
distance from an observed correlation. Data come from the International Tree
Ring Data Bank (ITRDB).

**[Read the rendered analysis, with all figures and results](quercus_analysis_final.md)**

## Contents

| Path | Description |
|---|---|
| `quercus_analysis_final.qmd` | The analysis (Quarto source) |
| `quercus_analysis_final.md` | The rendered analysis (GitHub Markdown, displays directly on GitHub) |
| `quercus_analysis_final_files/` | Figures embedded in the rendered analysis |
| `plots/` | The main figures saved as standalone PNGs |
| `data/itrdb_study.parquet` | One row per ITRDB site/study: code, name, coordinates, elevation, species, year span |
| `data/itrdb_measurements.parquet` | One row per ring: `itrdb_code`, `sample_id`, `global_id`, `year`, `raw_width`, `species` |
| `data/itrdb_sample.parquet` | One row per series, with summary statistics (not read by the analysis; included for convenience) |
| `data/sites_and_contributors.csv` | Human-readable list of the sites used, with contributor names |

The subset has 182 sites (57 QUPE, 125 QURO), 5,028 series and 630,676 ring
measurements.

## Reproducing the analysis

Requirements: R (≥ 4.x) and Quarto, plus these R packages:

```r
install.packages(c("arrow", "tidyverse", "data.table", "dplR", "sf", "tmap",
                   "rnaturalearth", "rnaturalearthdata", "patchwork",
                   "gridExtra", "knitr", "rmarkdown"))
```

The exact package versions used are listed in `sessionInfo.txt`.

Then, from this folder:

```bash
quarto render quercus_analysis_final.qmd
```

This regenerates `quercus_analysis_final.md` together with its figure folder and `plots/`.

## Data provenance

All data come from the International Tree Ring Data Bank (ITRDB), distributed
by NOAA NCEI (<https://www.ncei.noaa.gov/products/paleoclimatology/tree-ring>).
The ring-width files were downloaded on 11 March 2026. The subset contains
the ITRDB ring-width sites recorded as *Quercus petraea* (QUPE) or
*Quercus robur* (QURO). Density, multivariate and malformed files were
excluded.

About the measurements:

- `raw_width` values are stored exactly as they appear in the source files.
  Different contributors used different units (mm, 0.01 mm or 0.001 mm), and
  no unit correction is applied. The analysis explains why the log-spline
  detrending makes this irrelevant.
- Non-positive widths are dropped inside the analysis.
- Some site names in the source files were latin1-encoded; they have been
  re-encoded to UTF-8.

## Licence and attribution

- **Code** (`*.qmd`, `*.R`) is released under the MIT licence (see `LICENSE`).
- **Data** in `data/` is a subset of the ITRDB, distributed by NOAA NCEI. The
  MIT licence does not cover it. ITRDB data are freely available, but the
  original contributors should be credited. The contributor of each site is
  listed in `data/sites_and_contributors.csv` (the first part of
  `study_name`). If you use these data, please cite the ITRDB and the
  relevant contributors, following the NOAA data citation guidance.
