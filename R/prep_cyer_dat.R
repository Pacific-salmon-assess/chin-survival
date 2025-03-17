## Prepare CTC harvest rate estimates
# Unlike prep_cyer_dat_old.R use standard CTC mortality tables
# March 17, 2025

library(tidyverse)
library(readxl)

stock_key <- read.csv(
  here::here("data", "ctc_decoder", "ctc_stock_decoder.csv"),
  stringsAsFactors = FALSE
) 


## CLEAN CYER DATA -------------------------------------------------------------

sheet_names <- excel_sheets(
  here::here(
    "data", "harvest", "ctc_esc_data", 
    "TCCHINOOK-25-01-Appendix-C-Mortality-Distribution-Tables-Detailed.xlsx"
  ))

stocks <- stock_key %>% 
  filter(!is.na(ctc_indicator)) %>% 
  pull(ctc_indicator) %>% 
  unique()
matching_sheets <- sheet_names[
  sapply(sheet_names, function(x) any(grepl(paste(stocks, collapse = "|"), x)))
]
matching_sheets2 <- matching_sheets[
  sapply(matching_sheets, function(x) any(grepl("TM", x)))
]
sheet_ids <- which(sheet_names %in% matching_sheets2)


# CWT based CYERs from Laura Tessier
# identify sheets w/ relevant data and associated stock name
new_col_names <- c(
  "year", "cwt_n", "ages", "aabm_seak_t", "aabm_seak_n", "aabm_seak_s", 
  "aabm_nbc_t", "aabm_nbc_s", "aabm_wcvi_t", "aabm_wcvi_s", "isbm_nbc_t", 
  "isbm_nbc_n", "isbm_nbc_s", "isbm_sbc_t", "isbm_sbc_n", "isbm_sbc_s", 
  "isbm_n_falcon_t", "isbm_n_falcon_s", "isbm_s_falcon_t", "isbm_s_falcon_s", 
  "isbm_wac_n", "isbm_puget_n", "isbm_puget_s", "term_seak_t", "term_seak_n",
  "term_seak_s", "term_can_n", "term_can_s", "term_sus_t", "term_sus_n", 
  "term_sus_s", "stray", "esc", "comment"
) 

cwt_dat <- purrr::map2(
  sheet_ids, matching_sheets2, 
  function(x, y) {
    dum <- read_xlsx(
      here::here(
        "data", "harvest", "ctc_esc_data", 
        "TCCHINOOK-25-01-Appendix-C-Mortality-Distribution-Tables-Detailed.xlsx"
      ),
      sheet = x,
      skip = 6,
      col_names = FALSE
    )
    colnames(dum) <- new_col_names
    dum %>% 
      mutate(
        indicator = str_split(y, " ") %>% unlist() %>% .[1],
        mark = str_split(y, " ") %>% unlist() %>% .[2]
      ) %>% 
      # remove five year averages at bottom of table
      filter(!grepl("-", year))
  }
) %>% 
  bind_rows()


# southern US harvest not available for 2023; calculate mean values 2016-22 for
# each fishery, add to original dataset for 23 then rescale
cwt_dat_long <- cwt_dat %>% 
  filter(comment == "ok") %>% 
  mutate(indicator = paste(indicator, mark, sep = "_")) %>% 
  pivot_longer(cols = c(starts_with("aabm"), starts_with("isbm"),
                        starts_with("term"), stray, esc),
               names_to = "strata", values_to = "percent_run") %>% 
  mutate(
    year = as.numeric(year),
    southern_us = ifelse(
      (grepl("falcon", strata) | grepl("_sus_", strata) | 
         grepl("puget", strata) | grepl("wac", strata)),
      TRUE,
      FALSE
    ),
    canadian_er = ifelse(
      (grepl("nbc", strata) | grepl("sbc", strata) | 
         grepl("wcvi", strata) | grepl("term_can", strata)),
      TRUE,
      FALSE
    ),
    missing_from_fmi = ifelse(
      strata %in% c("isbm_nbc_t", "isbm_nbc_n", "isbm_sbc_t", "isbm_sbc_n",
                    "term_can_n"),
      TRUE,
      FALSE
    )
  ) 

# calculate mean southern US exploitation rate to use since 2022 values 
# unavailable
mean_sus <- cwt_dat_long %>% 
  filter(year > 2015 & year < 2023) %>% 
  group_by(strata, indicator) %>% 
  summarize(mean_percent_run = mean(percent_run))


cwt_dat_long2 <- left_join(cwt_dat_long, mean_sus, 
                           by = c("indicator", "strata")) %>% 
  mutate(
    percent_run = ifelse(year == "2023" & southern_us == TRUE,
                         mean_percent_run,
                         percent_run)
  ) %>% 
  group_by(
    indicator, year
  ) %>% 
  mutate(
    total_percent = sum(percent_run)
  ) %>% 
  ungroup() %>% 
  mutate(
    scaled_percent = percent_run / total_percent
  )


#combine stray and escapement, then use to calculate total exploitation
cyer_dat <- cwt_dat_long2 %>% 
  filter(strata %in% c("esc", "stray")) %>% 
  group_by(year, indicator) %>% 
  summarize(
    percent_escaped = sum(scaled_percent)
  ) %>% 
  ungroup() %>% 
  mutate(
    total_er = 1 - percent_escaped
  ) 

ggplot(cyer_dat) + 
  # geom_point(aes(x = year, y= can_er)) + 
  geom_point(aes(x = year, y= total_er), color = "red") + 
  facet_wrap(~indicator) +
  ggsidekick::theme_sleek()


# export exploitation rate for focal domain
cwt_dat_out <- cwt_dat_long2 %>% 
  filter(grepl("isbm", strata) | grepl("aabm_wcvi", strata),
         !grepl("isbm_nbc", strata)) %>% 
  group_by(year, indicator) %>% 
  summarize(
    focal_er = sum(scaled_percent)
  ) %>% 
  ungroup() %>% 
  mutate(
    clip = ifelse(grepl("unmarked", indicator), "N", "Y"),
    stock = str_split(indicator, "_")  %>%
      purrr::map(., head, n = 1) %>%
      unlist() 
  ) 

saveRDS(cwt_dat_out,
        here::here("data", "harvest", "cleaned_cyer_dat.rds"))
