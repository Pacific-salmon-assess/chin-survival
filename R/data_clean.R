### Data Cleaning
## Import individual detections data and reform to pass to models

library(tidyverse)


dat_tbl <- readRDS(here::here("data", "det_history_tbl.RDS"))

chin <- readRDS(here::here("data", "cleanTagData_GSI.RDS")) %>%
  filter(
    !charter == "msf",
    !is.na(acoustic),
    !is.na(agg_name),
    !fish == "CK113", # magnet not removed
    (agg_prob > 90 | genetic_source != "GSI")
  )

indicator_key <- read.csv(
  here::here("data", "ctc_stock_decoder.csv")
) %>% 
  mutate(
    ctc_indicator = ifelse(ctc_indicator == "SRH/ELK", "SRH", ctc_indicator)
  )

# generate alternative CYER datasets
cyer_dat <- readRDS(
  here::here("data", "harvest", "cleaned_cyer_dat_adj.rds")
) %>%
  rename(ctc_indicator = stock, focal_er_adj = focal_er) %>%
  filter(year > 2018)


stage_dat <- readRDS(
  here::here("data", "agg_lifestage_df.RDS")) %>% 
  select(
    fish, known_stage = stage, stage_1 = stage_predicted, 
    stage_2 = stage_predicted_sens
  )


## CLEAN AND EXPORT FOR MODEL FITTING ------------------------------------------

fr_sum_yr <- indicator_key %>% 
  filter(agg_name == "Fraser Summer Year.") %>% 
  pull(stock)
fr_spr_yr <- indicator_key %>% 
  filter(agg_name == "Fraser Spring Year.") %>% 
  pull(stock)

chin2 <- left_join(
  chin, 
  indicator_key %>% 
    select(stock, ctc_indicator, ctc_name) %>% 
    distinct(),
  by = "stock") %>% 
  mutate(
    agg_name = case_when(
      stock %in% fr_sum_yr ~ "Fraser Sum. 1.3",
      stock %in% fr_spr_yr ~ "Fraser Spr. 1.x",
      agg_name == "Fraser Sum. 4.1" ~ "Fraser Spr. 0.3",
      # define subyearlings based on run bubble_ts plots
      agg_name == "Fraser 4.1" ~ "Fraser Sum. 0.3",
      acoustic_year %in% c("7719_2019", "7701_2019", "7692_2019", "7691_2019",
                           "7692_2019", "7690_2019") ~ "Fraser Sum. 0.3",
      acoustic_year %in% c("7696_2019", "5353_2022") ~ "Fraser Fall",
      grepl("CAPILANO", stock) ~ "ECVI",
      TRUE ~ agg_name
    ),
    # define injury scores as per SJ's analyses
    adj_inj = case_when(
      hook_loc == "eye" ~ 3,
      injury == "3" | fin_dam == "3" | scale_loss == "3" ~ 2,
      injury == "2" | fin_dam == "2" | scale_loss == "2" ~ 1,
      TRUE ~ 0
    )
  ) %>% 
  left_join(., cyer_dat, by = c("year", "ctc_indicator", "clip")) %>% 
  left_join(., stage_dat, by = "fish")


# define aggregate names for Fraser and add to tbl
agg_names <- chin2 %>%
  filter(grepl("Fraser", agg_name)) %>%
  select(vemco_code = acoustic_year, agg = agg_name)


# export supplementary table summarizing stock breakdown
stock_supp_table <- chin2 %>% 
  mutate(
    agg_name = fct_recode(agg_name,
                          "Spring 1.x" = "Fraser Spr. 1.x",
                          "Summer 1.3" = "Fraser Sum. 1.3", 
                          "Summer 0.3" = "Fraser Sum. 0.3", 
                          "Fall 0.3" = "Fraser Fall")
  ) %>% 
  filter(agg_prob > 90) %>% 
  group_by(
    Stock = agg_name, Population = stock, CTC_Indicator = ctc_name
  ) %>% 
  tally() 
# write.csv(
#   stock_supp_table,
#   here::here(
#     "data", "stock_supp_table.csv"
#   ),
#   row.names = FALSE
# )


## Sample sizes by maturity stage
chin2 %>% 
  group_by(stage_2) %>% 
  tally()


## Logistic Regression Dataset
det_dat1 <- dat_tbl %>% 
  dplyr::select(stock_group, agg_det) %>% 
  unnest(cols = c(agg_det)) %>%
  left_join(., agg_names, by = "vemco_code") %>% 
  mutate(
    stock_group = ifelse(stock_group == "Fraser", agg, stock_group) %>% 
      as.factor()
  ) %>% 
  # remove fish with low stock assignment
  filter(vemco_code %in% chin2$acoustic_year) %>%
  left_join(., 
            chin2 %>% 
              mutate(month = lubridate::month(date)) %>% 
              select(vemco_code = acoustic_year, month, year, acoustic_type, 
                     known_stage, stage_1, stage_2, lat, lon, year_day, 
                     fl, wt, lipid, adj_inj, ctc_indicator, clip,
                     focal_er_adj, comment),
            by = "vemco_code") %>%
  mutate(
    # redeploy = ifelse(acoustic_type %in% c("V13P", "V13"), "no", "yes"),
    terminal_p = case_when(
      stock_group %in% c("South Puget", "Up Col.", "Low Col.") ~ 1,
      grepl("Fraser", stock_group) ~ 1,
      # specify years when extensive in-river detections available
      stock_group == "WCVI" & year %in% c("2021", "2022") ~ 1,
      TRUE ~ 0
    )
  ) %>% 
  arrange(
    year, stock_group
  )


# small number of tags (<2% missing lipid data; impute)
bio_dat <- det_dat1 %>% 
  select(vemco_code, fl, wt, lipid, stock_group, year) %>%
  distinct() 
interp_lipid <- bio_dat %>%
  select(-vemco_code) %>% 
  VIM::kNN(., k = 5) %>% 
  select(-ends_with("imp")) 
det_dat1$lipid <- interp_lipid$lipid

 
# export 
saveRDS(det_dat1, here::here("data", "surv_hts_data.rds"))



## CJS Model Data
chin_no_severe <- chin2 %>% 
  filter(!adj_inj == "3")

# add updated Fraser groupings and imputed lipid content
dat_tbl_trim <- dat_tbl %>% 
  filter(!stock_group %in% c("ECVI", "North Puget", "WA_OR", "WCVI")) %>% 
  mutate(
    bio_dat = purrr::map(bio_dat, function (x) {
      x %>%
        filter(vemco_code %in% chin_no_severe$acoustic_year) %>%
        select(-c(agg, lipid)) %>% 
        mutate(tag_date_z = scale(year_day) %>% as.numeric()) %>% 
        left_join(
          ., 
          det_dat1 %>%
            select(vemco_code, agg = stock_group, lipid,
                   known_stage, stage_1, stage_2)
        )
    }),
    wide_array_dat = purrr::map(wide_array_dat, function (x) {
      x %>%
        filter(vemco_code %in% chin_no_severe$acoustic_year) %>% 
        select(-agg)
    })
  ) %>%
  select(stock_group, bio_dat, wide_array_dat)

saveRDS(dat_tbl_trim, here::here("data", "surv_cjs_data.rds"))
