## Check location of terminal receivers for each stock group by plotting 
# receiver maps by year


library(tidyverse)
library(rmapshaper)
library(mapdata)
library(marmap)
library(sf)
library(cowplot)
library(patchwork)
library(ggspatial)



coast_plotting <- readRDS(here::here("data",
                                     "coast_major_river_sf_plotting.RDS"))


base_map <- ggplot() +
  geom_sf(data = coast_plotting, color = "black", fill = "white") +
  labs(x = "", y = "") +
  theme_void() +
  theme(panel.background = element_rect(colour="black", fill="darkgrey"),
        legend.position = "top", 
        axis.text = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1)) +
  coord_sf(expand = FALSE)


## DEPLOYMENTS MAP -------------------------------------------------------------

# path to high res shapefiles for bc coast only
if (Sys.info()['sysname'] == "Windows") {
  coast_path <- "C:/Users/FRESHWATERC/OneDrive - DFO-MPO/General - Applied Salmon Ecology/Spatial Data/High_res_coastline"
} else {
  coast_path <- "/Users/cam/Google Drive/spatial/coastline_shapefiles/High_res_coastline"
}

bc_coast <- st_read(
  here::here(coast_path, "lpr_000b16a_e_dissolve_NAD83BCAlbers.shp")
) %>% 
  st_transform("+proj=longlat +datum=WGS84") %>% 
  st_make_valid() %>% 
  st_crop(
    xmin = -128, xmax = -121.5, ymin = 45.5, ymax = 51.5
  )


chin <- readRDS(here::here("data", "cleanTagData_GSI.RDS"))
tag_dat <- readRDS(here::here("data", "surv_log_reg_data.rds"))

deploy_pts <- chin %>% 
  filter(acoustic_year %in% tag_dat$vemco_code)

# primary map (now combined as inset with receiver maps)
deploy_map <- ggplot() +
  geom_sf(data = bc_coast, color = "black", fill = "darkgrey") +
  labs(x = "", y = "") +
  ggsidekick::theme_sleek() +
  coord_sf(expand = FALSE) + 
  coord_sf(xlim = c(-126.1, -125.15), ylim = c(48.5, 49.01)) +
  geom_jitter(data = deploy_pts, 
              aes(x = lon, y = lat),
              fill = "red", inherit.aes = FALSE, shape = 21, alpha = 0.6,
              width = 0.025) +
  theme(
    strip.background = element_rect(colour="white", fill="white"),
    legend.position = "none",
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
  ) 

# inset map
deploy_inset <- base_map +
  coord_sf(expand = FALSE, ylim = c(45, 51)) +
  geom_rect(aes(xmin = -126.2, xmax = -124.75, ymin = 48.25, ymax = 49.2),
            color = "red", fill = NA, linewidth = 1)

png(here::here("figs", "maps", "deploy_map.png"), height = 4.25, width = 5,
    units = "in", res = 250)
ggdraw() +
  draw_plot(deploy_map) +
  draw_plot(deploy_inset, x = 0.79, y = 0.16, width = 0.2, height = 0.2)
dev.off()


## RECEIVER LOCATIONS _---------------------------------------------------------

# NOTE: these locations only show VR and PIT arrays, not recaptures
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

rect_data <- data.frame(xmin = -126.05, xmax = -125.2, ymin = 48.45, ymax = 49)
multi_year_map <- ggplot() +
  geom_sf(data = coast_plotting, color = "black", fill = "white") +
  labs(x = "", y = "") +
  ggsidekick::theme_sleek() +
  geom_rect(data = rect_data,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            inherit.aes = FALSE,  
            color = "red", fill = NA) + 
  coord_sf(xlim = c(-127.7, -122), ylim = c(46, 51), expand = FALSE) +
  geom_point(data = rec_all2, 
             aes(x = longitude, y = latitude),
             fill = "blue", alpha = 0.5,
             inherit.aes = FALSE, shape = 21) +
  facet_wrap(~field_season) +
  scale_x_continuous(breaks = c(-127, -125, -123)) +
  theme(
    legend.position = "none",
    panel.background = element_rect(colour="black", fill="grey30"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
  )


world_sf <- rnaturalearth::ne_countries(scale = "medium", 
                                        returnclass = "sf") %>%
  filter(admin %in% c("United States of America", "Canada", "Mexico"))

# define your Lambert Conformal Conic projection
lambert_crs <- "+proj=lcc +lat_1=20 +lat_2=75 +lon_0=-135 +datum=WGS84 +units=m +no_defs"

# transform to Lambert
world_lambert <- st_transform(world_sf, crs = lambert_crs)
rect_ll <- st_as_sfc(
  st_bbox(c(xmin=-127.7, xmax=-122, ymin=46, ymax=51), crs = st_crs(4326))
)
rect_lambert <- st_transform(rect_ll, crs = lambert_crs)

# (optional) if you only want to plot the region around your map, crop first
bbox_ll <- st_bbox(c(xmin=-145, xmax=-120, ymin=20, ymax=75), crs = st_crs(4326))
world_clip <- st_crop(world_lambert, bbox_ll %>% st_transform(lambert_crs))

inset_world <- ggplot() +
  geom_sf(data = world_clip, fill = "grey30", color = NA) +
  geom_sf(data = rect_lambert, fill = NA, color = "blue", linewidth = 0.75) +
  theme_bw() +
  theme(
    legend.position = "none",
    panel.background = element_rect(colour = "black", fill = "white"),
    panel.border = element_rect(color = "black", fill = NA, size = 1)
  ) +
  scale_y_continuous(breaks = seq(30, 70, by = 10), expand = c(0, 0)) +
  scale_x_continuous(breaks = c(-140, -130, -120), expand = c(0, 0)) 
  

x0 <- 0.67;   y0 <- 0.25
w0 <- 0.025;   h0 <- 0.045
w1 <- 0.25;   h1 <- 0.62
x1 <- x0 - (w1 - w0) / 2
y1 <- y0 - (h1 - h0) / 2
 
png(here::here("figs", "maps", "multi-year-map.png"),
    height = 5, width = 6.5, units = "in", res = 200)
ggdraw() + 
  draw_plot(multi_year_map) +
  draw_plot(inset_world, 
            x = x1 + 0.1,
            y = y1,
            width = w1,
            height = h1) 
dev.off()


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
      # divide upstream and downstream fraser
      (region == "fraser" | (region == "river" & latitude > 49)) & 
        longitude > -122.5 ~ "fraser_2",
      (region == "fraser" | (region == "river" & latitude > 49)) & 
        (longitude < -122.5 & longitude > -123) ~ "fraser_1",
      (region %in% c("river", "columbia", "swwa_or") | marine == "no") & 
        (longitude > -124.15 & longitude < -121.94 & 
           latitude < 46.3 & latitude > 45) ~ "lower_col",
      region == "puget" & latitude < 47.3 ~ "in_river",
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
array_key <- read.csv(here::here("data", 
                                 "surv_segment_key_2023.csv")) %>% 
  # exclude stocks missing from CJS analysis
  # filter(
  #   !stock_group %in% c("ECVI", "North Puget", "WA_OR")
  # ) %>%
  mutate(segment = array_num - 1,
         segment_name = str_replace(segment_name, " ", "\n")) %>% 
  distinct() 
  
array_tbl <- array_key %>% 
  group_by(stock_group) %>% 
  group_nest()

# specify in-river receivers to drop 
cali <- rec %>% 
  filter(!region == "bonn",
         !(region == "in_river" & latitude > 44),
         !(region == "lower_col" & longitude > -123.7),
         !grepl("fraser", region))
fraser <- rec %>% 
  filter(!region == "bonn",
         !(region == "in_river" & !(project_name %in% c("Haro", "Instream", "Kintama"))),
         !(region == "lower_col" & longitude > -123.7))
col <- rec %>% 
  filter(!(region == "in_river" & 
             (latitude < 45 | latitude > 47)),
         !grepl("fraser", region))
puget <- rec %>% 
  filter(!region == "bonn",
         !(region == "in_river" & 
             (latitude < 47 | latitude > 48.25)),
         !(region == "in_river" & longitude > -12.45),
         !(region == "lower_col" & longitude > -123.7),
         !grepl("fraser", region))
wa_or <- rec %>%
  filter(!region == "bonn",
         !(region == "in_river"),
         !grepl("fraser", region))
wcvi <- rec %>% 
  filter(!region == "bonn",
         !(region == "in_river" & 
             (latitude < 49.1 | longitude > -124)),
         !(region == "lower_col" & longitude > -123.7),
         !grepl("fraser", region))
ecvi <- rec %>%
  filter(!region == "bonn",
         !(region == "in_river"),
         !grepl("fraser", region))
rec_dummy_list <- list(cali, ecvi, fraser, col, puget, puget, col, wa_or,
                       wcvi)

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
  theme(legend.position = "none",
        axis.ticks = element_blank()
        ) +
  facet_wrap(~stock_group)
seg_list <- purrr::map(
  array_list[-1],
  ~ base_map +
    coord_sf(xlim = c(-126, -121), ylim = c(45.5, 51)) +
    geom_point(
      data = .x,
      aes(x = longitude, y = latitude, fill = segment_name),
      shape = 21
    ) +
    scale_fill_viridis_d() +
    theme(legend.position = "none",
          axis.ticks = element_blank()) +
    facet_wrap(~stock_group)
)

array_map2 <- cowplot::plot_grid(
  seg_list[[3]], seg_list[[6]], 
  seg_list[[2]], seg_list[[5]], 
  ncol = 2
  )
array_map <- cowplot::plot_grid(
  cali_segs, array_map2, rel_widths = c(0.2, 0.3)
  )


png(here::here("figs", "maps", "array-map.png"),
    height = 5, width = 5, units = "in", res = 200)
array_map
dev.off()
