## Prepare CTC harvest rate estimates using standard CTC mortality tables
# March 17, 2025

library(tidyverse)
library(readxl)

stock_key <- read.csv(
  here::here("data", "ctc_decoder", "ctc_stock_decoder.csv"),
  stringsAsFactors = FALSE
) 


## CLEAN CYER DATA -------------------------------------------------------------

cwt_dat_long2 <- readRDS(
  here::here(
    "data", "harvest", "cyer_dat_long.rds"
  )
)

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


# as above but removes puget isbm for puget stocks 
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