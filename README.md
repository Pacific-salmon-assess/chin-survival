# Chinook Acoustic Telemetry Surival Analysis

Data and code associated with acoustic telemetry-based estimates of marine survival rates for adult Chinook salmon. Manuscript (**Revisiting Marine Mortality Assumptions in Adult Pacific Salmon: High Survival and Minimal Variability Among Multiple Stocks of Chinook Salmon**) submitted to **Ecological Applications**.

Scripts
1. General

  + data_clean.R - import individual detection histories and prep for models
  + map_receivers - map figures
  
2. Cormack-Jolly-Seber

  + cjs_priors_check - prior predictive checks for CJS models
  + fit_cjs - fit primary CJS models presented in main text and generate figures
  + fit_cjs_sensitivity - fit supplementary CJS models presented in online appendix and generate figures

3. Hierarchical Terminal Survival

  + prep_cyer_dat - clean and export calendar year exploitation rate data used as proxy for harvest pressure in HTS model
  + supp_figs - box plots included in online supplement
  + fit_hts - fit primary HTS models presented in main text and generate figures
  + fit_hts_sensitivity - fit supplementary HTS models presented in online appendix and generate figures
  
