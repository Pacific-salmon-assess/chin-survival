## Prepare CTC harvest rate estimates
# Aug 21, 2023
# Input data (details in readme)
# 1) CTC estimates for exploitation in marine areas (southern BC, Puget Sound; 
# coastal WA outstanding)
# 2) In-river harvest estimates for Columbia River summer and fall MUs
# TODO: WCVI AABM needs to be added for all years; reconcile Noel's data vs. 
# raw inputs


library(tidyverse)

stock_key <- read.csv(
  here::here("data", "ctc_decoder", "ctc_stock_decoder.csv"),
  stringsAsFactors = FALSE
) 

cyer_dat <- read.csv(
  here::here(
    "data", "harvest", "Interim_ERA_CatchYr_ExploitationRates_2023.csv"
    ),
  stringsAsFactors = FALSE
) %>% 
  janitor::clean_names() %>% 
  filter(
    mort_type == "TM"
  )


# CLEAN ------------------------------------------------------------------------

sus_stocks <- cyer_dat %>% 
  filter(is.na(escap)) %>% 
  pull(stock) %>% 
  unique()


# calculate 2023 mean exp rate for SUS fisheries as placeholder
mean_sus_exp_rates <- cyer_dat %>% 
  filter(catch_year < 2023) %>% 
  group_by(stock) %>% 
  summarize(
    mean_sus_isbm = mean(sus_isbm)
  )  

cyer_dat_long <- left_join(cyer_dat, mean_sus_exp_rates, by = "stock") %>% 
  mutate(
    sus_isbm = ifelse(catch_year == "2023", mean_sus_isbm, sus_isbm)
  ) %>% 
  select(stock, year = catch_year, seak_aabm:escap) %>% 
  pivot_longer(cols = c("seak_aabm", "can_aabm", "can_isbm", "sus_isbm", 
                        "term", "stray", "escap"),
               names_to = "strata", 
               values_to = "run_ppn") %>% 
  # rescale to correct for missing 2022 data and ensure all stock-year 
  # combinations sum to 100
  group_by(
    stock, year
  ) %>% 
  mutate(
    total_percent = sum(run_ppn)
  ) %>% 
  ungroup() %>% 
  mutate(
    scaled_percent = run_ppn / total_percent
  ) %>% 
  # drop missing SUS 2023 data that's empty
  filter(
    !(year == "2023" & stock %in% sus_stocks)
  )


# ALL 2023 data missing for SUS stocks; replace with 2019-2022 mean
mean_sus_stocks <- cyer_dat_long %>% 
  filter(stock %in% sus_stocks) %>% 
  group_by(stock, strata) %>% 
  summarize(mean_ppn = mean(run_ppn)) %>% 
  group_by(stock) %>% 
  mutate(
    total_percent = sum(mean_ppn),
    scaled_percent = mean_ppn / total_percent,
    year = 2023
  ) %>% 
  ungroup()


# Nicola and Chilliwack missing from data provided by Noel; use 2019-2022 data
# from published catch/escapement tables (cleaned in fraser-chinook-fsar)
cyer_pub <- readRDS(here::here("data", "harvest", "cyer_est_22.rds")) %>% 
  rename(stock = indicator)

# create dummy 2023 data based on mean percent run (currently these are the 
# 2022 SUS estimates as well)
cyer_2023 <- expand.grid(
  year = 2023,
  strata = unique(cyer_pub$strata),
  stock = c("CHI", "NIC")
) %>% 
  left_join(
    ., cyer_pub %>% select(stock, strata, mean_percent_run) %>% distinct(), 
    by = c("strata", "stock")
  ) %>% 
  filter(
    grepl("isbm", strata) | grepl("aabm_wcvi", strata),
    !grepl("nbc", strata),
    stock %in% c("CHI", "NIC")
  ) %>% 
  group_by(
    year, stock
  ) %>% 
  summarize(
    isbm_cyer = sum(mean_percent_run) / 100
  ) %>% 
  ungroup()


cyer_pub2 <- cyer_pub %>% 
  filter(
    grepl("isbm", strata) | grepl("aabm_wcvi", strata),
    !grepl("nbc", strata),
    stock %in% c("CHI", "NIC"),
    year > 2018
  ) %>% 
  group_by(
    year, stock
  ) %>% 
  summarize(
    isbm_cyer = sum(percent_run) / 100
  ) %>% 
  ungroup() %>% 
  rbind(., cyer_2023)


## combine to join in data_clean.R
cyer_dat_out<- rbind(
  cyer_dat_long %>% 
    select(stock, year, strata, scaled_percent),
  mean_sus_stocks %>% 
    select(stock, year, strata, scaled_percent)
) %>% 
  # exclude strata outside of sampling area
  filter(
    grepl("isbm", strata) | grepl("aabm_wcvi", strata)
  ) %>% 
  group_by(stock, year) %>% 
  summarize(
    isbm_cyer = sum(scaled_percent)
  ) %>% 
  rbind(., cyer_pub2) %>% 
  ungroup()

saveRDS(cyer_dat_out,
        here::here("data", "harvest", "cleaned_cyer_dat.rds"))


# OLD --------------------------------------------------------------------------

col_fall <- read.csv(
  here::here("data", "harvest", "Fall Chinook for Cam2.csv"),
  stringsAsFactors = FALSE
) %>% 
  janitor::clean_names() %>% 
  mutate(
    harvest = ifelse(type %in% c("Escapement", "Missing Fish"), "N", "Y")
  )

col_summer <- read.csv(
  here::here("data", "harvest", "2023 Spring JSR TABLES_020223_trim.csv"),
  stringsAsFactors = FALSE
) %>% 
  janitor::clean_names() 


# annual ctc harvest rate estimates 
marine_er <- ctc %>% 
  # select total mortality rather than landed catch
  filter(mort_type == "TM",
         !catch_year == "2022") %>% 
  pivot_longer(
    cols = aabm_wcvi_sport:terminal_canada_sport, 
    names_to = "fishery",
    values_to = "cyer"
  ) %>% 
  group_by(stock, catch_year) %>% 
  summarize(
    total_cyer = sum(cyer),
    .groups = "drop"
  ) %>% 
  left_join(., stock_key, by = "stock") 

ggplot(marine_er) +
  geom_boxplot(aes(x = as.factor(catch_year), y = total_cyer, fill = aggregate)) +
  facet_wrap(~aggregate)


# Columbia Fall harvest below Bonneville
col_fall_harvest <- col_fall %>% 
  filter(harvest == "Y",
         below_bon == "Y") %>% 
  group_by(return_year, mgmt_stock) %>% 
  summarize(total_catch_below_bon = sum(sum_of_n_fish),
            .groups = "drop") 
col_fall_return <- col_fall %>% 
  group_by(return_year, mgmt_stock) %>% 
  summarize(total_return = sum(sum_of_n_fish),
            .groups = "drop")
col_fall_er <- left_join(col_fall_harvest, col_fall_return, 
                        by = c("mgmt_stock", "return_year")) %>% 
  mutate(
    lower_river_er = total_catch_below_bon / total_return
  ) 

ggplot(col_fall_er) +
  geom_boxplot(aes(x = as.factor(return_year), y = lower_river_er))


# Columbia Summer harvest below Bonneville
col_summer_er <- col_summer %>% 
  mutate(
    lower_river_catch = sum(downstream_comm, downstream_sport, downstream_misc2,
                           downstream_treaty),
    lower_river_er = lower_river_catch / upriver_run) %>% 
  select(year, starts_with("lower_river")) 

ggplot(col_summer_er) +
  geom_boxplot(aes(x = as.factor(year), y = lower_river_er))



## Baranovs catch equation -----------------------------------------------------

## discrete parameters
# A #total mortality 
# S #total survival (S = 1 - A)
# u #fishing mortality
# v #natural mortality

## instantaneous parameters
# Z #total mortality (Z = -ln(1 - A))
# F #fishing mortality
# M #natural mortality

## conditional parameters (A = m + n - mn = u + v)
# m #fishing mortality (m = 1 - exp(F))
# n #natural mortality (n = 1 - exp(M))

## Type 2 fisheries u does not = m
# u = (FA)/Z
# v = (MA)/Z

## Baranovs equation
# average pop abundance during interval = Nbar = (NA)/Z
# total deaths = NA = ZNbar
# total natural deaths = Nv = MNbar
# catch = C = Nu = FNbar = (FNA)/Z

## available inputs
# A (survival estimated by tags)
# u exp rate estimated by:
# a) CWT recoveries/releases, b) Fraser run reconstruction, c) tag recoveries
# Nbar (escapement or CWTs in escapement)

# 1) use A and u to calculate F and M (for Z) 
A <- 0.3 # ~average for 19-21
Z <- -log(1 - A)
u <- 0.20 # recent average for fall 4_1
v <- A - u

#(Z / A) = (F / u) = (M / v)
FF <- (Z / A) * u
M <- Z - FF
v <- (M * A) / Z

# 2) use Baranovs equation to calculate abundance prior to fishery
Nbar <- 10000
NN <- ((Nbar * FF) * Z) / (FF * A)



