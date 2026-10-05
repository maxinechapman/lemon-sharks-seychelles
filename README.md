# lemon-sharks-seychelles
This code was written and published as part of an undergraduate project by Maxine Chapman and Giorgos Papageorgiou. The data tracks the movements of 43 lemon sharks over almost a decade in the Seychelles.

## DEPENDENCIES

The annual space-use analysis uses package "tidytracks". Install it once with:

install.packages("pak")
pak::pak("EvolEcolGroup/tidytracks")


## HELPER FILES — KEEP IN THE CODE FOLDER, BUT DO NOT RUN DIRECTLY

03_data_processing_functions_MC.R
Defines shared functions used for shark size estimation and annual spatial analyses.

plot_lmer_diag_loess.R
Defines mixed-model diagnostic and random-effect plotting functions used by the modelling script.

## MAIN RUN ORDER
Please note that the raw data used in 1, 2 and our seasonality analysis cannot be made publicly available because it contains spatially explicit data for an endangered species. Therefore we have provided the intermediate results tables results_annual_kde_coa.csv for the main analysis and results_monthly.csv for the seasonality analysis. Further analysis from 3 onwards or within the seasonality document can be continued in full.

1. 01_COA.Rmd
Reads combined_sharks.csv, projects receiver locations to the local LAEA coordinate system, compares candidate time bins and calculates detection-weighted 60-minute centres of activity (COAs). It validates the calculations and writes data/intermediate/combined_sharks_coa.csv.

2. 02_annual_data_kde.Rmd
Creates one row per shark-year, including estimated size, detection coverage, projected COA-based KDE50/KDE95 areas and COA-based MCP areas. It writes results/results_annual_kde_coa.csv

3. 03_metadata.R
Combines the annual results with capture metadata to make one summary row per shark, including tagging/final sizes, tracking duration, detections and maturity at tagging and final detection. It writes results/metadata/shark_metadata_full.csv and shark_metadata_full.rds, which are required by the seasonality script.

4. 04_filtering.Rmd
Applies the final annual-analysis criteria: valid KDE50 and KDE95, at least 50 COAs and at least six detection months for KUDs; MCP100 needs temporal coverage and a valid positive estimate. It saves separate KDE and MCP model datasets, and COAs joined to the annual inclusion flags (for plotting polygons) in data/intermediate.

5. 05_modelling.Rmd
Reads the filtered RDS files and fits the KDE50, KDE95 and MCP100 mixed models against size, sex and detection months, with shark ID as a random intercept. It also runs diagnostics, plots the spatial relationships and KDE ratio, examines selected large-range cases and tests alternative maturity-threshold effects.

6. 06_figures.Rmd
Creates figures for publication.

## ADDITIONAL ANALYSES

growth_morphometrics.Rmd — RUN AFTER STEP 2
Fits and checks the TL–PCL relationship, evaluates the Brown–Gruber growth curve against Stevens data and converts published life-stage thresholds into PCL.

seasonality.Rmd — RUN AFTER STEP 3
Uses raw detections plus shark_metadata_full.rds to calculate monthly residency, including genuine zero-detection months, for either the full array or the atoll receiver subset. It writes results_monthly.csv/.rds and residency figures, then fits a beta-binomial mixed model testing whether seasonal residency changes with shark size.



