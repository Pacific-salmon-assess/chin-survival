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


## MAIN MAP --------------------------------------------------------------------

# deployment locations
rec_all <- readRDS(here::here("data", 
                              "receivers_all.RDS"))$rec_all
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


base_map <- ggplot() +
  geom_polygon(data = w_can, mapping = aes(x = long, y = lat, group = group), 
               color = "black", fill = "darkgrey") + 
  labs(x = "", y = "") +
  ggsidekick::theme_sleek() +
  theme(strip.background = element_rect(colour="white", fill="white"),
        legend.position = "top") +
  coord_fixed(xlim = c(-128.5, -122), ylim = c(48, 51), ratio = 1.3) 


multi_year_map <- base_map + 
  coord_fixed(xlim = c(-127.7, -122), ylim = c(46, 51), ratio = 1.3) + 
  geom_point(data = rec_all2, 
             aes(x = longitude, y = latitude),
             fill = "red",
             inherit.aes = FALSE, shape = 21) +
  facet_wrap(~field_season) +
  theme(
    legend.position = "none"
  )
# export for use in Rmd
saveRDS(multi_year_map, here::here("figs", "receiver_map1.rds"))

multi_year_map_region <- base_map + 
  coord_fixed(xlim = c(-127, -122), ylim = c(46, 50), ratio = 1.3) + 
  geom_point(data = rec_all2, 
             aes(x = longitude, y = latitude, fill = region),
             inherit.aes = FALSE, shape = 21) +
  scale_fill_discrete("Regions") +
  facet_wrap(~field_season)


## DEPLOYMENTS MAP -------------------------------------------------------------

chin <- readRDS(here::here("data", "cleanTagData_GSI.RDS"))
tag_dat <- readRDS(here::here("data", "cleaned_log_reg.rds"))

deploy_pts <- chin %>% 
  filter(acoustic_year %in% tag_dat$vemco_code)

deploy_map <- base_map + 
  coord_fixed(xlim = c(-126.25, -124), ylim = c(48, 49.25), ratio = 1.3) + 
  geom_point(data = deploy_pts, 
             aes(x = lon, y = lat, fill = as.factor(year)),
             inherit.aes = FALSE, shape = 21) +
  theme(
    legend.position = "none"
  )
# export for use in Rmd
saveRDS(deploy_map, here::here("figs", "deploy_map.rds"))




rec <- readRDS(
  here::here("data", "tagging_data", "for_glatos", "receivers_all.RDS")
)$rec_all %>% 
  mutate(
    region = case_when(
      region %in% c("columbia", "swwa_or") & 
        (station_longitude > -124.15 & station_longitude < -122.9 & 
           station_latitude < 46.3 & station_latitude > 45) ~ "col_est",
      marine == "no" ~ "in_river",
      region == "fraser" ~ "in_river",
      region == "bark_snd" & station_latitude > 49.1 ~ "in_river",
      region == "bark_snd" & station_latitude < 49.1 ~ "wcvi",
      region == "disc_isl" & station_latitude > 49.25 ~ "nevi",
      TRUE ~ region
    )
  ) %>% 
  filter(
    !region == "alaska"
  )

base_map <- readRDS(here::here("data", "generated_data", "base_map.RDS"))

base_map +
  coord_fixed(xlim = c(-126, -121), ylim = c(38, 51), ratio = 1.3) +
  geom_point(data = rec,
             aes(x = station_longitude, y = station_latitude, fill = region),
             shape = 21) +
  theme(legend.position = "right")


# use array key to plot array number by stock group 
array_tbl <- array_key %>% 
  group_by(agg) %>% 
  group_nest() 

plot_list <- purrr::map2(
  array_tbl$agg, array_tbl$data,
  function (x, y) {
    dd <- left_join(
      rec, y, by = "region"
    ) %>% 
      mutate(array_num = as.factor(array_num))
    base_map +
      coord_fixed(xlim = c(-126, -121), ylim = c(38, 51), ratio = 1.3) +
      geom_point(
        data = dd,
        aes(x = station_longitude, y = station_latitude, fill = array_num),
        shape = 21
      ) +
      scale_fill_viridis_d() +
      theme(legend.position = "right", axis.text = element_blank()) +
      labs(title = x)
  }
)

pdf(
  here::here("figs", "maps", "array_grouping_by_stock.pdf")
)
plot_list
dev.off()