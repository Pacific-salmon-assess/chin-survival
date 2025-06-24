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

stocks <- stock_key %>% 
  filter(!is.na(ctc_indicator)) %>% 
  pull(ctc_indicator) %>% 
  unique()

## import marked data
sheet_names <- excel_sheets(
  here::here(
    "data", "harvest", "ctc_esc_data", 
    "TCCHINOOK-25-XX-Appendix-C-Mortality-Distribution-Tables-Detailed-marked_noTBRSEAK.xlsx"
  ))

matching_sheets <- sheet_names[
  sapply(sheet_names, function(x) any(grepl(paste(stocks, collapse = "|"), x)))
]
matching_sheets2 <- matching_sheets[
  sapply(matching_sheets, function(x) any(grepl("total mort", x)))
]
sheet_ids <- which(sheet_names %in% matching_sheets2)

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

cwt_dat_marked <- purrr::map2(
  sheet_ids, matching_sheets2, 
  function(x, y) {
    dum <- read_xlsx(
      here::here(
        "data", "harvest", "ctc_esc_data", 
        "TCCHINOOK-25-XX-Appendix-C-Mortality-Distribution-Tables-Detailed-marked_noTBRSEAK.xlsx"
      ),
      sheet = x,
      skip = 6,
      col_names = FALSE
    )
    colnames(dum) <- new_col_names
    dum %>% 
      mutate(
        indicator = str_split(y, " ") %>% unlist() %>% .[1],
        mark = "marked"
      ) %>% 
      # remove five year averages at bottom of table
      filter(!grepl("-", year),
             ! year == "2024") 
  }
) %>% 
  bind_rows()

cwt_dat_unmarked <- purrr::map2(
  sheet_ids, matching_sheets2, 
  function(x, y) {
    dum <- read_xlsx(
      here::here(
        "data", "harvest", "ctc_esc_data", 
        "TCCHINOOK-25-XX-Appendix-C-Mortality-Distribution-Tables-Detailed-unmarked_noTBRSEAK.xlsx"
      ),
      sheet = x,
      skip = 6,
      col_names = FALSE
    )
    colnames(dum) <- new_col_names
    dum %>% 
      mutate(
        indicator = str_split(y, " ") %>% unlist() %>% .[1],
        mark = "unmarked"
      ) %>% 
      # remove five year averages at bottom of table
      filter(!grepl("-", year),
             ! year == "2024") 
  }
) %>% 
  bind_rows()


cwt_dat_long <- rbind(cwt_dat_unmarked, cwt_dat_marked) %>%
  filter(comment == "ok") %>%
  mutate(indicator = paste(indicator, mark, sep = "_")) %>%
  pivot_longer(cols = c(starts_with("aabm"), starts_with("isbm"),
                        starts_with("term"), stray, esc),
               names_to = "strata", values_to = "percent_run") %>%
  mutate(
    year = as.numeric(year)#,
    # southern_us = ifelse(
    #   (grepl("falcon", strata) | grepl("_sus_", strata) |
    #      grepl("puget", strata) | grepl("wac", strata)),
    #   TRUE,
    #   FALSE
    # )
  )

cwt_dat_long2 <- cwt_dat_long %>% 
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


# look at relative impact of WCVI fisheries on different indicators
cwt_dat_long2 %>% 
  filter(strata %in% c("aabm_wcvi_s", "aabm_wcvi_s"), 
         year > 2018,
         mark == "unmarked") %>%
  group_by(indicator, year) %>% 
  summarize(total_er = sum(percent_run)) %>% 
  group_by(indicator) %>% 
  summarize(mean_er = mean(total_er)) %>% 
  arrange(mean_er) %>% 
  print(n = Inf)

cwt_dat_long2 %>% 
  filter(strata %in% c("isbm_puget_n", "isbm_puget_s"), 
         year > 2018,
         mark == "unmarked") %>% 
  group_by(indicator, year) %>% 
  summarize(total_er = sum(percent_run)) %>% 
  group_by(indicator) %>% 
  summarize(mean_er = mean(total_er)) %>% 
  arrange(mean_er) %>% 
  print(n = Inf)


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
# 
# ggplot(cyer_dat) + 
#   # geom_point(aes(x = year, y= can_er)) + 
#   geom_point(aes(x = year, y= total_er), color = "red") + 
#   facet_wrap(~indicator) +
#   ggsidekick::theme_sleek()


# export exploitation rate for focal domain
cwt_dat_out <- cwt_dat_long2 %>% 
  filter(grepl("isbm", strata) | grepl("aabm_wcvi", strata),
         !grepl("isbm_nbc", strata),
         !grepl("isbm_puget", strata)
         ) %>% 
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
        here::here("data", "harvest", "cleaned_cyer_dat_no_puget.rds"))



# as above but includes puget
cwt_dat_out <- cwt_dat_long2 %>% 
  filter(grepl("isbm", strata) | grepl("aabm_wcvi", strata),
         !grepl("isbm_nbc", strata)
  ) %>% 
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
