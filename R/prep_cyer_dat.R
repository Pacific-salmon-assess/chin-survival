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
    year = as.numeric(year)
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
# cwt_dat_long2 %>% 
#   filter(strata %in% c("aabm_wcvi_s", "aabm_wcvi_s"), 
#          year > 2018,
#          mark == "unmarked") %>%
#   group_by(indicator, year) %>% 
#   summarize(total_er = sum(percent_run)) %>% 
#   group_by(indicator) %>% 
#   summarize(mean_er = mean(total_er)) %>% 
#   arrange(mean_er) %>% 
#   print(n = Inf)
# 
# cwt_dat_long2 %>% 
#   filter(strata %in% c("isbm_puget_n", "isbm_puget_s"), 
#          year > 2018,
#          mark == "unmarked") %>% 
#   group_by(indicator, year) %>% 
#   summarize(total_er = sum(percent_run)) %>% 
#   group_by(indicator) %>% 
#   summarize(mean_er = mean(total_er)) %>% 
#   arrange(mean_er) %>% 
#   print(n = Inf)


#combine stray and escapement, then use to calculate total exploitation
# cyer_dat <- cwt_dat_long2 %>% 
#   filter(strata %in% c("esc", "stray")) %>% 
#   group_by(year, indicator) %>% 
#   summarize(
#     percent_escaped = sum(scaled_percent)
#   ) %>% 
#   ungroup() %>% 
#   mutate(
#     total_er = 1 - percent_escaped
#   ) 
#  
# ggplot(cyer_dat) + 
#   # geom_point(aes(x = year, y= can_er)) + 
#   geom_point(aes(x = year, y= total_er), color = "red") + 
#   facet_wrap(~indicator) +
#   ggsidekick::theme_sleek()


# function to make proxy Chilko data assumed 50% harvest of SHU
chi_foo <- function(x) {
  x %>% 
    filter(
      stock == "SHU"
    ) %>% 
    group_by(
      year, clip
    ) %>% 
    mutate(
      focal_er = focal_er / 2,
      stock = "SHU_adjusted",
      indicator = ifelse(clip == "Y", "SHU_adjusted_marked", "SHU_adjusted_unmarked")
    ) 
}

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
chi_dat <- chi_foo(cwt_dat_out)
saveRDS(rbind(cwt_dat_out, chi_dat),
        here::here("data", "harvest", "cleaned_cyer_dat_no_puget.rds"))


# as above but includes puget
cwt_dat_out2 <- cwt_dat_long2 %>% 
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
chi_dat2 <- chi_foo(cwt_dat_out2)
saveRDS(rbind(cwt_dat_out2, chi_dat2),
        here::here("data", "harvest", "cleaned_cyer_dat.rds"))


# as above but adjusts puget isbm for puget stocks (divides by two since 
# only northern fisheries will impact them)
cwt_dat_out3 <- cwt_dat_long2 %>% 
  filter(grepl("isbm", strata) | grepl("aabm_wcvi", strata),
         !grepl("isbm_nbc", strata)
  ) %>% 
  mutate(
    new_scaled_percent = ifelse(
      grepl("puget", strata) & indicator %in% c(
        "SAM_marked", "SSF_marked", "SPS_marked", "STL_marked", "SKY_marked",
        "SAM_unmarked", "SSF_unmarked", "SPS_unmarked", "STL_unmarked",
        "SKY_unmarked"
      ),
      0, #scaled_percent / 2,
      scaled_percent
    )
  ) %>% 
  group_by(year, indicator) %>% 
  summarize(
    focal_er = sum(new_scaled_percent)
  ) %>% 
  ungroup() %>% 
  mutate(
    clip = ifelse(grepl("unmarked", indicator), "N", "Y"),
    stock = str_split(indicator, "_")  %>%
      purrr::map(., head, n = 1) %>%
      unlist() 
  )
chi_dat3 <- chi_foo(cwt_dat_out3)
saveRDS(rbind(cwt_dat_out3, chi_dat3),
        here::here("data", "harvest", "cleaned_cyer_dat_adj.rds"))



## CLEAN PUGET SOUND CATCH/EFFORT DATA -----------------------------------------

# creel rec catch data from https://wdfw.wa.gov/fishing/reports/creel/puget-annual?sample_date=3&ramp=&catch_area=
paths <- list.files(path = here::here("data", "harvest", "puget_sound_catch"),
                    pattern = "\\.csv$", 
                    full.names = TRUE)
ps_dat <- lapply(paths, read.csv) %>% 
  bind_rows() %>% 
  janitor::clean_names() %>% 
  mutate(
    dttm = lubridate::mdy(sample_date),
    year = lubridate::year(dttm),
    month = lubridate::month(dttm),
    region = case_when(
      catch_area %in% c(
        "Area 4, Eastern portion",
        "Area 3, La Push"
      ) ~ "washington_coastal",
      catch_area %in% c(
        "Area 8-2, Ports Susan and Gardner",
        "Area 8-1, Deception Pass, Hope Island, and Skagit Bay",
        "Area 7, San Juan Islands",                             
        "Area 6, East Juan de Fuca Strait",
        "Bellingham Bay",
        "Area 5, Sekiu and Pillar Point",                       
        "Area 6-2, Eastern portion of Area 6",
        "Area 6-1, Western portion of Area 6",
        "Dungeness Bay"
      ) ~ "north_ps",
      TRUE ~ "south_ps"
    )
  ) %>% 
  filter(
    month > 4 & month < 11,
    !region == "washington_coastal"
  ) 

ps_dat %>% 
  group_by(year) %>%
  mutate(total_chinook = sum(chinook, na.rm = TRUE)) %>%
  group_by(region, year) %>%
  summarize(
    region_chinook = sum(chinook, na.rm = TRUE),
    total_chinook = first(total_chinook),  # same for entire year
    ppn_chinook = region_chinook / total_chinook,
    .groups = "drop"
  )
# approximately 50% of Puget Sound catch occurs south of arrays