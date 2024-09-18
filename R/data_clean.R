## Data Cleaning, Exploration and Summary Figures 
# Sep 13, 2024

library(tidyverse)


dat_tbl <- readRDS(here::here("data", "det_history_tbl.RDS"))

chin <- readRDS(here::here("data", "cleanTagData_GSI.RDS")) %>% 
  filter(
    !charter == "msf",
    !is.na(acoustic)
  )

indicator_key <- read.csv(
  here::here("data", "ctc_decoder", "ctc_stock_decoder.csv")
) 


## CLEAN AND EXPORT FOR MODEL FITTING ------------------------------------------


## add exploitation rate indicator and redefine Fraser by SMU
# TODO: add exploitation rate estimates by year
# TODO: eventually replace stocks without stock IDs, but with aggregate IDs with
# mean annual estimates for the aggregate
# TODO: decide whether to split upriver Columbia, lower Columbia and Puget Sound
# similarly to FR
fr_sum_yr <- indicator_key %>% 
  filter(agg_name == "Fraser Summer Year.") %>% 
  pull(stock)
fr_spr_yr <- indicator_key %>% 
  filter(agg_name == "Fraser Spring Year.") %>% 
  pull(stock)

chin2 <- left_join(
  chin, 
  indicator_key %>% 
    select(stock, ctc_indicator) %>% 
    distinct(),
  by = "stock") %>% 
  mutate(
    agg_name = case_when(
      stock %in% fr_sum_yr ~ "Fraser Sum. Yr.",
      stock %in% fr_spr_yr ~ "Fraser Spr. Yr.",
      # define subyearlings based on run bubble_ts plots
      acoustic_year %in% c("7719_2019", "7701_2019", "7692_2019", "7691_2019",
                           "7692_2019", "7690_2019") ~ "Fraser 4.1",
      acoustic_year %in% c("7696_2019", "5353_2022") ~ "Fraser Fall",
      TRUE ~ agg_name
    )
  ) 


# define aggregate names for Fraser and add to tbl
agg_names <- chin2 %>%
  filter(grepl("Fraser", agg_name)) %>%
  select(vemco_code = acoustic_year, agg = agg_name)



## Logistic Regression Dataset
det_dat1 <- dat_tbl %>% 
  dplyr::select(stock_group, agg_det) %>% 
  unnest(cols = c(agg_det)) %>%
  left_join(., agg_names, by = "vemco_code") %>% 
  mutate(
    term_det = ifelse(final_det + river_det > 0, 1, 0),
    stock_group = ifelse(stock_group == "Fraser", agg, stock_group) %>% 
      as.factor()
  ) %>% 
  left_join(., 
            chin %>% 
              mutate(month = lubridate::month(date)) %>% 
              select(vemco_code = acoustic_year, month, year, acoustic_type, 
                     fl, lipid, year_day, hook_loc, fin_dam, injury, scale_loss,
                     comment),
            by = "vemco_code") %>%
  mutate(
    comp_inj = (injury + scale_loss + fin_dam),
    redeploy = ifelse(acoustic_type %in% c("V13P", "V13"), "no", "yes")
  ) %>% 
  arrange(
    year, stock_group
  )

# small number of tags (<2% missing lipid data; impute)
bio_dat <- det_dat1 %>% 
  select(vemco_code, fl, lipid, stock_group, year) %>%
  distinct() 
interp_lipid <- bio_dat %>%
  select(-vemco_code) %>% 
  VIM::kNN(., k = 5) %>% 
  select(-ends_with("imp")) 
det_dat1$lipid <- interp_lipid$lipid

# export 
saveRDS(det_dat1, here::here("data", "surv_log_reg_data.rds"))



## CJS Model Data
# ID tag codes to retain
kept_tags <- dat_tbl %>%
  # unnest and filter out stage 3
  unnest(cols = bio_dat) %>%
  filter(redeploy == "no") %>%
  pull(vemco_code)

# add updated Fraser groupings and imputed lipid content
dat_tbl_trim <- dat_tbl %>% 
  filter(!stock_group %in% c("ECVI", "North Puget")) %>% 
  mutate(
    bio_dat = purrr::map(bio_dat, function (x) {
      x %>%
        filter(vemco_code %in% kept_tags) %>% 
        select(-c(agg, lipid)) %>% 
        left_join(
          ., 
          det_dat1 %>% select(vemco_code, agg = stock_group, lipid)
        )
    }),
    wide_array_dat = purrr::map(wide_array_dat, function (x) {
      # remove aggregate vector if only one level
      dd <- x %>%
        filter(vemco_code %in% kept_tags)
      # if (length(unique(dd$agg)) == "1") {
      dd <- dd %>% select(-agg)
      # }
      return(dd)
    })
  ) %>%
  select(stock_group, bio_dat, wide_array_dat)

saveRDS(dat_tbl_trim, here::here("data", "surv_cjs_data.rds"))



## DATA SHARE ------------------------------------------------------------------

# export summary of stocks, stock groups and indicators to share
indicator_dat <- chin %>% 
  select(stock, cu_name, agg_name) %>% 
  distinct() %>% 
  arrange(agg_name)

write.csv(
  indicator_dat,
  here::here(
    "data", "ctc_decoder", "ctc_stock_decoder_template.csv"
  ),
  row.names = FALSE
)
