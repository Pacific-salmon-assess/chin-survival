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
) %>% 
  mutate(
    ctc_indicator = ifelse(ctc_indicator == "SRH/ELK", "SRH", ctc_indicator)
  )

# TODO: currently excludes terminal harvest; should this be included?
cyer_dat <- readRDS(here::here("data", "harvest", "cleaned_cyer_dat.rds")) %>% 
  rename(ctc_indicator = stock)


# ID tag codes to retain
kept_tags <- dat_tbl %>%
  # unnest and filter out stage 3
  unnest(cols = bio_dat) %>%
  filter(redeploy == "no") %>%
  pull(vemco_code)


## CLEAN AND EXPORT FOR MODEL FITTING ------------------------------------------

# TODO: decide whether to split upriver Columbia, lower Columbia and Puget Sound
# similarly to FR to fit hierarchical CJS
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
      agg_name == "Fraser 4.1" ~ "Fraser Sum. 4.1",
      acoustic_year %in% c("7719_2019", "7701_2019", "7692_2019", "7691_2019",
                           "7692_2019", "7690_2019") ~ "Fraser Sum. 4.1",
      acoustic_year %in% c("7696_2019", "5353_2022") ~ "Fraser Fall",
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
  left_join(., cyer_dat, by = c("year", "ctc_indicator")) 


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
            chin2 %>% 
              mutate(month = lubridate::month(date)) %>% 
              select(vemco_code = acoustic_year, month, year, acoustic_type, 
                     lat, lon, 
                     fl, lipid, year_day, adj_inj, ctc_indicator, isbm_cyer,
                     comment),
            by = "vemco_code") %>%
  mutate(
    redeploy = ifelse(acoustic_type %in% c("V13P", "V13"), "no", "yes"),
    terminal_p = case_when(
      stock_group %in% c("South Puget", "Up Col.") ~ 1,
      grepl("Fraser", stock_group) ~ 1,
      stock_group == "WCVI" & year %in% c("2021", "2022") ~ 1,
      TRUE ~ 0
    )
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


## PRELIM FIGURES --------------------------------------------------------------

ppn_foo <- function(group, response) {
  group_exp <- c("vemco_code", group)
  labs <- det_dat1 %>% 
    dplyr::select_at(group_exp) %>% 
    distinct() %>% 
    group_by_at(group) %>% 
    tally() %>% 
    mutate(year = as.factor(year))
  
  dat_out <- det_dat1 %>% 
    group_by_at(group) %>% 
    summarize(n = length(unique(vemco_code)),
              ppn = sum(.data[[response]] / n),
              se = sqrt((ppn * (1 - ppn)) / n),
              up = pmin(1, ppn + (1.96 * se)),
              lo = pmax(0, ppn - (1.96 * se)),
              .groups = "drop") %>% 
    mutate(year = as.factor(year))
  
  list("labs" = labs, "dat" = dat_out)
}


det_ppns <- ppn_foo(group = c("stock_group", "year"), response = "term_det") 
ggplot(data = det_ppns$dat, aes(y = ppn)) +
  geom_pointrange(aes(x = year, ymin = lo, ymax = up)) +
  facet_wrap(~stock_group) +
  ggsidekick::theme_sleek() +
  geom_text(data = det_ppns$lab, aes(x = year, y = -Inf, label  = n),
            position = position_dodge(width = 1),
            vjust = -0.5)

inj_ppns <- ppn_foo(group = c("injury", "year"), response = "term_det")
ggplot(data = inj_ppns$dat,  aes(x = year, y = ppn)) +
  geom_pointrange(aes(ymin = lo, ymax = up)) +
  facet_wrap(~injury) +
  ggsidekick::theme_sleek()+
  geom_text(data = inj_ppns$lab, 
            aes(x = year, y = -Inf, label  = n),
            position = position_dodge(width = 1),
            vjust = -0.5)

inj_ppns2 <- ppn_foo(group = c("comp_inj", "year"), response = "term_det")
ggplot(data = inj_ppns2$dat,  aes(x = year, y = ppn)) +
  geom_pointrange(aes(ymin = lo, ymax = up)) +
  facet_wrap(~comp_inj) +
  ggsidekick::theme_sleek()+
  geom_text(data = inj_ppns2$lab, 
            aes(x = year, y = -Inf, label  = n),
            position = position_dodge(width = 1),
            vjust = -0.5)


loc_ppns <- ppn_foo(group = c("hook_loc", "year"), response = "term_det")
ggplot(data = loc_ppns$dat,  aes(x = year, y = ppn)) +
  geom_pointrange(aes(ymin = lo, ymax = up)) +
  facet_wrap(~hook_loc) +
  ggsidekick::theme_sleek()+
  geom_text(data = loc_ppns$lab, 
            aes(x = year, y = -Inf, label  = n),
            position = position_dodge(width = 1),
            vjust = -0.5)

fin_ppns <- ppn_foo(group = c("fin_dam", "year"), response = "term_det")
ggplot(data = fin_ppns$dat,  aes(x = year, y = ppn)) +
  geom_pointrange(aes(ymin = lo, ymax = up)) +
  facet_wrap(~fin_dam) +
  ggsidekick::theme_sleek()+
  geom_text(data = fin_ppns$lab, 
            aes(x = year, y = -Inf, label  = n),
            position = position_dodge(width = 1),
            vjust = -0.5)

scale_ppns <- ppn_foo(group = c("scale_loss", "year"), response = "term_det")
ggplot(data = scale_ppns$dat,  aes(x = year, y = ppn)) +
  geom_pointrange(aes(ymin = lo, ymax = up)) +
  facet_wrap(~scale_loss) +
  ggsidekick::theme_sleek()+
  geom_text(data = scale_ppns$lab, 
            aes(x = year, y = -Inf, label  = n),
            position = position_dodge(width = 1),
            vjust = -0.5)

redep_ppns <- ppn_foo(group = c("redeploy", "year"), response = "det")
ggplot(data = redep_ppns$dat,  aes(x = year, y = ppn)) +
  geom_pointrange(aes(ymin = lo, ymax = up)) +
  facet_wrap(~redeploy) +
  ggsidekick::theme_sleek()+
  geom_text(data = redep_ppns$lab, 
            aes(x = year, y = -Inf, label  = n),
            position = position_dodge(width = 1),
            vjust = -0.5)


ggplot(det_dat1, aes(x = lipid_z, y = final_det)) +
  geom_point() +
  facet_wrap(~year) +
  ggsidekick::theme_sleek()
ggplot(det_dat1, aes(x = trough_time, y = det)) +
  geom_point() +
  facet_wrap(~year) +
  ggsidekick::theme_sleek()
ggplot(det_dat1, aes(x = mean_log_e, y = final_det)) +
  geom_point() +
  facet_wrap(~year) +
  ggsidekick::theme_sleek()
ggplot(det_dat1, aes(x = cyer_z, y = final_det)) +
  geom_point() +
  facet_wrap(~year) +
  ggsidekick::theme_sleek()

