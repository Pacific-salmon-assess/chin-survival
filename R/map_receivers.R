## Check location of terminal receivers for each stock group by plotting 
# receiver maps by year


library(tidyverse)
library(rmapshaper)
library(mapdata)
library(marmap)
library(sf)


w_can <- map_data("worldHires", region = c("usa", "canada")) %>%
  filter(long < -110) 
coast_plotting <- readRDS(here::here("data",
                                     "coast_major_river_sf_plotting.RDS"))



base_map <- ggplot() +
  geom_sf(data = coast_plotting, fill = "white", colour = "white") +
  labs(x = "", y = "") +
  ggsidekick::theme_sleek() +
  theme(strip.background = element_rect(colour="white", fill="white"),
        panel.background = element_rect(fill = "grey40"),
        legend.position = "top", 
        axis.text = element_blank()) +
  coord_sf(expand = FALSE)


## DEPLOYMENTS MAP -------------------------------------------------------------

chin <- readRDS(here::here("data", "cleanTagData_GSI.RDS"))
tag_dat <- readRDS(here::here("data", "surv_log_reg_data.rds"))

deploy_pts <- chin %>% 
  filter(acoustic_year %in% tag_dat$vemco_code)

deploy_map <- base_map + 
  coord_sf(xlim = c(-126.25, -124.5), ylim = c(48.4, 49.25)) + 
  geom_jitter(data = deploy_pts, 
              aes(x = lon, y = lat, fill = as.factor(year)),
              inherit.aes = FALSE, shape = 21,
              width = 0.025) +
  theme(
    legend.position = "none"
  )
# export for use in Rmd
saveRDS(deploy_map, here::here("figs", "deploy_map.rds"))


## RECEIVER LOCATIONS _---------------------------------------------------------

rec_all <- readRDS(here::here("data", 
                              "receivers_all.RDS"))$rec_all %>%
  dplyr::rename(latitude = station_latitude,
                longitude = station_longitude)
rec_all2 <- readRDS(here::here("data", 
                               "receivers_all.RDS"))$dup_rec_all %>%
  dplyr::rename(name = station_name,
                latitude = station_latitude,
                longitude = station_longitude) %>%
  mutate(
    project_name = fct_recode(project_name,
                              `NOAA Puget` = "PS",
                              Kintama = "Haro",
                              UBC = "temp_ubc",
                              `Feather River` = "feather",
                              SOBaD = "noaa_roegner")
  ) %>%
  droplevels() %>% 
  filter(
    !field_season == "2024"
  )


multi_year_map <- base_map + 
  coord_sf(xlim = c(-127.7, -122), ylim = c(46, 51)) + 
  geom_point(data = rec_all2, 
             aes(x = longitude, y = latitude),
             fill = "red",
             inherit.aes = FALSE, shape = 21) +
  facet_wrap(~field_season) +
  theme(
    legend.position = "none"
  )

png(here::here("figs", "multi-year-map.png"), 
    height = 5, width = 6.5, units = "in", res = 200)
multi_year_map
dev.off()

# export for use in Rmd
saveRDS(multi_year_map, here::here("figs", "receiver_map1.rds"))


# import PIT detections to include in receiver map
pit_dat <- readRDS(here::here("data", "cleanedPIT_detections.RDS")) %>% 
  mutate(region = "in_river",
         project_name = "pit",
         marine = "no") %>% 
  select(
    station_name = site, project_name, region, latitude, longitude, marine
  ) %>% 
  distinct()

rec <- rec_all  %>% 
  filter(
    !region == "alaska",
    !recover_year == "2024"
  ) %>% 
  select(station_name, project_name, region, latitude, longitude, marine) %>% 
  rbind(., pit_dat) %>% 
  mutate(
    project_name = ifelse(is.na(project_name), "Instream", project_name),
    # make sure regional definitions approximate those in 
    # chinTagging/prep_detection_histories
    region = case_when(
      grepl("BO", station_name) ~ "bonn",
      region %in% c("river", "columbia", "swwa_or") & 
        (longitude > -124.15 & longitude < -121.94 & 
           latitude < 46.3 & latitude > 45) ~ "lower_col",
      region == "river" ~ "in_river", #fish caught in terminal locations 
      marine == "no" ~ "in_river",
      region == "wcvi" & latitude < 48.405 & longitude > -125.5 ~ "nwwa",
      region == "fraser" ~ "in_river",
      region == "bark_snd" & latitude > 49.1 ~ "in_river",
      region == "bark_snd" & latitude < 49.1 ~ "wcvi",
      region == "disc_isl" & latitude > 49.25 ~ "nevi",
      TRUE ~ region
    )
  )

# use array key to plot array number by stock group 
array_tbl <- read.csv(here::here("data", 
                                 "surv_segment_key_2023.csv"))  %>% 
  # exclude stocks missing from CJS analysis
  filter(
    !stock_group %in% c("ECVI", "North Puget")
  ) %>%
  mutate(segment = array_num - 1,
         segment_name = str_replace(segment_name, " ", "\n")) %>% 
  distinct() %>%
  group_by(stock_group) %>% 
  group_nest()

# specify in-river receivers to drop 
cali <- rec %>% 
  filter(!region == "bonn",
         !(region == "in_river" & latitude > 44))
fraser <- rec %>% 
  filter(!region == "bonn",
         !(region == "in_river" & !(project_name %in% c("Instream", "Kintama"))))
col <- rec %>% 
  filter(!(region == "in_river" & 
             (latitude < 45 | latitude > 47)))
puget <- rec %>% 
  filter(!region == "bonn",
         !(region == "in_river" & 
             (latitude < 47 | latitude > 48.25)))
wa_or <- rec %>% 
  filter(!region == "bonn",
         !(region == "in_river"))
wcvi <- rec %>% 
  filter(!region == "bonn",
         !(region == "in_river" & 
             (latitude < 49.1 | longitude > -124)))
rec_dummy_list <- list(cali, fraser, col, puget, col, wa_or, wcvi)

# make plots
array_list <- purrr::pmap(
  list(array_tbl$stock_group, array_tbl$data, rec_dummy_list),
  function (x, y, z) {
    left_join(
      z, y, by = "region"
    ) %>% 
      mutate(segment_name = fct_reorder(segment_name, array_num),
             stock_group = x) 
  }
  ) 
cali_segs <- base_map +
  coord_sf(xlim = c(-126, -121), expand = FALSE) +
  geom_point(
    data = array_list[[1]],
    aes(x = longitude, y = latitude, fill = segment_name),
    shape = 21
  ) +
  scale_fill_viridis_d() +
  theme(#legend.position = "none",
        axis.ticks = element_blank()
        ) +
  facet_wrap(~stock_group)
seg_list <- purrr::map(
  array_list[-1],
  ~ base_map +
    coord_sf(xlim = c(-126, -119), ylim = c(45.5, 51)) +
    geom_point(
      data = .x,
      aes(x = longitude, y = latitude, fill = segment_name),
      shape = 21
    ) +
    scale_fill_viridis_d() +
    theme(#legend.position = "none",
          axis.ticks = element_blank()) +
    facet_wrap(~stock_group)
)

array_map2 <- cowplot::plot_grid(
  seg_list[[1]], seg_list[[2]], seg_list[[3]], seg_list[[4]], 
  seg_list[[5]], seg_list[[6]],  ncol = 3
  )
array_map <- cowplot::plot_grid(
  cali_segs, array_map2, rel_widths = c(0.16, 0.5)
  )
# export for use in Rmd
saveRDS(array_map, here::here("figs", "array_map.rds"))
