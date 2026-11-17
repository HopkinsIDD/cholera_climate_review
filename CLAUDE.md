# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

# Instructions

- DO NOT TAKE INITIATIVES.
- Only perform the exact actions explicitly requested by the user.
- Do not make changes, edits, or updates unless directly asked.
- __THE PROTECTED FOLDER `~/secure/` IS GATED__: every read there asks the user first, and only the `data-guardian` agent reads it, returning structure and summaries, never raw records. Never copy, link or route around the gate.
- New decisions and dated change history go to the knowledge vault (a `decision` note or a `record`), never into this file.

## Project Overview

VIMC (Vaccine Impact Modelling Consortium) **WP3 deliverable 1**: a systematic review and Bayesian meta-analysis of the association between weather / climate exposures and cholera, led by Qulu Zheng (manuscript "Zheng et al. 2025, cholera weather systematic review"). It is the literature counterpart of the modelling repo `../cholera_climate_outbreaks_africa` (vault [[cholera-climate-outbreaks-africa]]), which covers WP3 deliverable 3.

**Question**: across published studies, how often and how strongly are precipitation, temperature, flood and drought associated with cholera occurrence and incidence, and how does that vary by setting (continent, climatic zone, endemicity, HDI, coast, urbanicity, lag, rainy season, temporal scale)?

- **Exposures** (`exp_class`): `precipitation`, `temperature`, `flood`, `drought` (`utils.R::get_all_exp_class()`).
- **Outcomes** (`out_class_metric`): `occurrence`, `incidence` (`get_all_out_class()`).
- **Effect types** (`get_all_effects()`): direction of effect (`direction_positive` / `direction_negative`, vote counting) and effect sizes (`corr` Fisher-z correlations, `diff` risk differences, `ratio` log10 ratios).
- Estimates before the 7th pandemic (`year(eff_tl) <= 1961`) are excluded; only `national` and `sub-national` estimates enter the meta-regression.

## Status

Analysis is essentially finished; the work is manuscript revision / resubmission. Submission history lives in `manuscript/` (01 Lancet Planetary Health → 02 Nature Climate Change → 03 PLOS Medicine → 04 Nature Health → PLOS Pathogens, most recent). Co-author comments: `manuscript/coauthors/`. `manuscript/` is binary (docx/pdf); do not edit it.

## Pipeline

Scripts in `Analysis/`, run from the repo root (paths are relative to it or built with `here::here()`). All sourced helpers are in `Analysis/utils.R` (single module, ~1.6k lines).

| Script | What it does |
|---|---|
| `00_prepare_shapefiles.R` | one-time: fuzzy-match (`stringdist` JW) extracted study locations to the HopkinsIDD location taxonomy and build location shapefiles. Reads files not in the repo (`location.csv`, `location_periods.csv`, `filtered_best_match.csv`). |
| `001_link_external_data.R` | one-time: link each location to Köppen-Geiger climatic zone, sub-national HDI, distance to coast/inland water (`Data/climatic_zone_raster`, `hdi_raster`, `inland_coast`). |
| `002_validate_extracted_data.R` | sanity checks on the extraction sheet (study period vs publication year, etc.). |
| `01_process_extracted_data.R` | process the extraction CSV → `generated_data/processed_estimate_data.rds`, `processed_effect_size_data_by_type.rdata`, quality scores (`paper_scores.rds`, `sum_scores_long.rds`). Option `-p` for plots. |
| `02_prepare_meta_analysis.R` | build Stan data for one (effect, exposure, outcome, stratum) → `generated_data/stan_data_*.rds`. |
| `03_run_meta_analysis.R` | fit with cmdstanr (4 chains, 250 warm-up + 500 sampling) and run generated quantities → `stan_fit_*.rds`, `stan_genquant_*.rds`. |
| `04_make_figures.R` | manuscript figures and tables (direction-of-effect `gt` tables, effect-size and post-stratified plots, maps) → `figures/`. |
| `submit_models.sh` | loops 02 → 03 over effects × exposures × outcomes × 10 strata (all; continent Asia/Africa/Other; climatic zone arid/temperate/tropical; % years with cholera <50 / 50–95 / >95). |

CLI options of 02/03 (`optparse`): `-e` effect, `-c` exposure class, `-o` outcome class, `-s` stratification column (`none`, `continent`, `new_climatic_zone`, `pct_years_cholera_cat`), `-l` level, `-r` redo.

```bash
Rscript Analysis/01_process_extracted_data.R
Rscript Analysis/02_prepare_meta_analysis.R -e direction_positive -c precipitation -o incidence -s continent -l Africa -r TRUE
Rscript Analysis/03_run_meta_analysis.R    -e direction_positive -c precipitation -o incidence -s continent -l Africa -r TRUE
bash Analysis/submit_models.sh       # full grid, long-running
Rscript Analysis/04_make_figures.R
```

File names for Stan artefacts come from `utils.R::make_stan_data_filename()` / `make_stan_fit_filename()` / `make_stan_genquant_filename()` (pattern `stan_<kind>_effect-<e>_exp-<c>_out-<o>_strat-<s>-<l>.rds`); never build them by hand.

## Models (`Analysis/stan/`)

- `direction_meta_analysis_marginalize_weighted.stan` — Bernoulli-logit meta-regression of P(positive direction) with study random effect (`sigma_re`) and estimate-level overdispersion (`sigma_epsilon`); missing covariates are **marginalised** over their possible values weighted by `p_missing_covar`; post-stratified probabilities (`poststrat`, `poststrat_pred`) and I² in generated quantities. Priors: `beta ~ N(0, 0.5)`, both sds `N(0, 1)`.
- `effect_meta_analysis_marginalize_weighted.stan` — same structure for continuous effect sizes (correlations on the Fisher-z scale, back-transformed with `inv_z_transform`). `effect_meta_analysis[_marginalize].stan` are earlier versions kept side by side.
- Meta-regression covariates (`02_prepare_meta_analysis.R`): `continent + rainy_season + urban_rural_class + new_climatic_zone + hdi_cat + inland_coast + lag_class + out_temp_scale + pct_years_cholera_cat`, minus the stratification column.

## Data

- `Data/extracted_data/` — the extraction workbook (`cholera_climate_systematic_review_extract_into_sheets.xlsx`, sheet `estimates_0329`) and its processed CSVs (`..._manuscript_version.csv` is the frozen manuscript version). The `.xlsx` is the source of truth; co-authors edit it.
- `Data/endemicity/` — WHO annual reports (`who_annual_report.csv`) used for `pct_years_cholera_cat`.
- `Data/climatic_zone_raster/`, `Data/hdi_raster/`, `Data/inland_coast/`, `Data/shapefiles/` — covariate layers.
- `generated_data/` — all intermediate RDS / Stan outputs (largely untracked; regenerate rather than edit).

## Reports

- `notebooks/systematic_review_summary_report.Rmd` — main summary report (incl. meta-regression results).
- `notebooks/Systematic-Review-Supplement-tableS2.Rmd`, `...-table-S5.Rmd` — supplement tables.
- `notebooks/country_report.Rmd` — per-country / continent data reports, rendered by `Analysis/scratchpads/s02_generate_reports.R`.
- `Analysis/scratchpads/` — exploration only (`s01` VIMC presentation, `s02` report rendering, `s03` direction-of-effect plot restyle); nothing in the pipeline calls it.

## Conventions

- Numbered scripts run; `utils.R` is sourced. No script-header or section convention is enforced here (unlike the outbreaks repo).
- Covariate names, exposure / outcome / effect vocabularies and colour palettes are defined once in `utils.R` (`get_*()`, `colors_*()`, `get_covar_dict()`); reuse them.
- Repo root junk to ignore: `tmp.rdata`, `nohup.out`, `Rplots.pdf`, `.Rhistory`.
- Remote: `HopkinsIDD/cholera_climate_review`; Qulu Zheng commits most data / notebook changes directly on `main`.

## Key R packages

`tidyverse`, `cmdstanr`, `optparse`, `here`, `furrr`, `sf`, `terra`, `stringdist`, `gt`, `patchwork`, `cowplot`, `RColorBrewer`, `rmarkdown`.

## Session Management (Knowledge Vault)

Project index note: [[cholera-climate-review]] in the knowledge vault at `../claude_code_vault` (absolute path `/home/perez/projects/claude_code_vault`). See `claude_code_vault/CLAUDE.md` for `/save` and `/resume`, and the note-type taxonomy ([[workflow/note-type-taxonomy]]).
