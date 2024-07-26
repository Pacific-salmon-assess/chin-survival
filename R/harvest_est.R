## Prepare CTC harvest rate estimates
# Aug 21, 2023
# Input data (details in readme)
# 1) CTC estimates for exploitation in marine areas (southern BC, Puget Sound; 
# coastal WA outstanding)
# 2) In-river harvest estimates for Columbia River summer and fall MUs


library(tidyverse)

stock_key <- read.csv(
  here::here("data", "harvest", "stock_key.csv"),
  stringsAsFactors = FALSE
) 

ctc <- read.csv(
  here::here("data", "harvest", "2023ERA_catchDistribution_CMZ_subset.csv"),
  stringsAsFactors = FALSE
) %>% 
  janitor::clean_names() 

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


# CLEAN ------------------------------------------------------------------------

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



