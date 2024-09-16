## Data Exploration and Summary Figures 
# Sep 13, 2024

library(tidyverse)


det_tbl <- readRDS(here::here("data", "det_history_tbl.RDS"))

chin <- readRDS(here::here("data", "cleanTagData_GSI.RDS")) %>% 
  filter(
    !charter == "msf",
    !is.na(acoustic)
  )


# export summary of stocks, stock groups and indicators to share
indicator_dat <- chin %>% 
  select(stock, cu_name, agg_name) %>% 
  distinct() %>% 
  arrange(agg_name)

write.csv(
  indicator_dat,
  here::here(
    "data", "ctc_stock_decoder.csv"
  ),
  row.names = FALSE
)
