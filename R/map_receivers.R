## Map receiver locations for manuscript


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
  coord_sf(expand = FALSE) +
  ggsidekick::theme_sleek() +
  theme(panel.background = element_rect(colour="black", fill = "darkgrey"),
        legend.position = "top", 
        # axis.text = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1)) 
  


# path to high res shapefiles for bc coast only
if (Sys.info()['sysname'] == "Windows") {
  coast_path <- "C:/Users/FRESHWATERC/OneDrive - DFO-MPO/General - Applied Salmon Ecology/Spatial Data/coastline_shapefiles/High_res_coastline"
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
tag_dat <- readRDS(here::here("data", "surv_hts_data.rds"))


## DEPLOYMENTS MAP -------------------------------------------------------------


# location names
places <- tibble::tribble(
  ~name,               ~lon,     ~lat,
  "Juan de Fuca\nStrait", -124,  47.9,
  "Strait of\nGeorgia",   -123.5,  49.4,
  "Vancouver\nIsland",    -125.8,  49.8,
  "Puget Sound",         -122.6,  47.8,
  "Columbia\nRiver",      -122.2,  46,
  "Fraser\nRiver",        -122.2,  49.5
) %>% 
  sf::st_as_sf(coords = c("lon","lat"), crs = 4326)


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
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    axis.text = element_blank(),
    axis.ticks = element_blank()
  ) 

# inset map
deploy_inset <- base_map +
  coord_sf(expand = FALSE, ylim = c(45, 51), xlim = c(-128.5, -121.5)) +
  geom_rect(aes(xmin = -126.2, xmax = -124.75, ymin = 48.25, ymax = 49.2),
            color = "red", fill = NA, linewidth = 1) +
  geom_label(
    data = places,
    aes(geometry = geometry, label = name),
    stat = "sf_coordinates",
    size = 3,
    fill = "white",       # background color
    colour = "black"
  )

png(here::here("figs", "maps", "deploy_map.png"), height = 6.25, width = 5,
    units = "in", res = 250)
# cowplot::plot_grid(deploy_map, deploy_inset, nrow = 2, 
#                    rel_widths =  c(0.5, 0.7),
#                    rel_heights = c(0.4, 0.6))
# ggdraw() +
#   draw_plot(deploy_map) +
#   draw_plot(deploy_inset, x = 0.79, y = 0.16, width = 0.2, height = 0.2)
ggdraw() +
  draw_plot(deploy_inset) +
  draw_plot(deploy_map, x = 0.13, y = 0.03, width = 0.45, height = 0.45)
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



base_map2 <- ggplot() +
  geom_sf(data = coast_plotting, color = "black", fill = "white") +
  labs(x = "", y = "") +
  theme_void() +
  theme(panel.background = element_rect(colour="black", fill="darkgrey"),
        legend.position = "top", 
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1)) +
  coord_sf(expand = FALSE)


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
      (longitude > -124.05 & longitude < -121.94 & 
         latitude < 46.3 & latitude > 45) ~ "lower_col",
      (longitude > -124.15 & longitude < -123.65 & 
         latitude < 46.3 & latitude > 45) ~ "swwa_or",
      region == "puget" & latitude < 47.3 ~ "in_river",
      region == "river" ~ "in_river", #fish caught in terminal locations 
      marine == "no" ~ "in_river",
      region == "wcvi" & latitude < 48.405 & longitude > -125.5 ~ "nwwa",
      region == "fraser" ~ "in_river",
      region == "bark_snd" & latitude > 49.1 ~ "in_river",
      region == "bark_snd" & latitude < 49.1 ~ "wcvi",
      region == "disc_isl" & latitude > 49.25 ~ "nevi",
      TRUE ~ region
    ),
    receiver_type = ifelse(project_name == "pit", "PIT", "acoustic")
  )

# use array key to plot array number by stock group 
array_key <- read.csv(here::here("data", 
                                 "surv_segment_key_2023.csv")) %>% 
  mutate(segment = array_num - 1) %>% 
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
array_list1 <- purrr::pmap(
  list(array_tbl$stock_group, array_tbl$data, rec_dummy_list),
  function (x, y, z) {
    left_join(
      z, y, by = "region"
    ) %>% 
      mutate(segment_name = fct_reorder(segment_name, array_num),
             stock_group = x) 
  }
  ) 
array_list <- array_list1[c(1,3,4,6,7)]

shape_pal <- c(21, 23)
names(shape_pal) <- c("acoustic", "PIT")


cali_segs <- base_map2 +
  coord_sf(xlim = c(-126, -121), expand = FALSE) +
  geom_point(
    data = array_list[[1]],
    aes(x = longitude, y = latitude, fill = segment_name, shape = receiver_type),
  ) +
  scale_fill_viridis_d() +
  scale_shape_manual(values = shape_pal) +
  theme(axis.ticks = element_blank(),
        legend.position = "right",
        legend.key = element_rect(fill = "transparent", colour = NA),
        legend.title = element_blank(),
        legend.text = element_text(size = 6),
        legend.background = element_rect(fill = "white", colour = NA)
  ) +           
  guides(
    fill = guide_legend(
      keyheight = unit(0.25, "cm"),
      keywidth  = unit(0.25, "cm"),
      override.aes = list(shape = 21, colour = "black")
    ),
    shape = "none"
  ) +
  facet_wrap(~stock_group)
cali_legend <- get_legend(
  cali_segs
)
cali_inset <- ggdraw(
  cali_segs +
    theme(legend.position = "none")
) +
  draw_grob(cali_legend, x = 0.55, y = 0.79, width = 0.25, height = 0.25)

# vectors specifying legend location (2, 4, 1, 3)
x_vec <- c(0.55, #fraser
           0.67, #lcol
           0.61, # puget
           0.67) #ucol
y_vec <- c(0.74, 0.73, 0.74, 0.71)

seg_list <- purrr::pmap(
  list(array_list[-1], x_vec, y_vec),
  function (dat_in, x_vec, y_vec) {
    dat_in <- dat_in %>% 
      filter(
        !is.na(segment_name)
      ) %>% 
      droplevels()
    p_segs <- base_map2 +
      coord_sf(xlim = c(-126, -121), ylim = c(45.5, 51)) +
      geom_point(
        data = dat_in,
        aes(x = longitude, y = latitude, fill = segment_name, shape = receiver_type),
      ) +
      scale_fill_viridis_d() +
      scale_shape_manual(values = shape_pal) +
      theme(axis.ticks = element_blank(),
            legend.position = "right",
            legend.key = element_rect(fill = "transparent", colour = NA),
            legend.title = element_blank(),
            legend.text = element_text(size = 6),
            legend.background = element_rect(fill = "white", colour = NA)
      ) +           
      guides(
        fill = guide_legend(
          keyheight = unit(0.25, "cm"),
          keywidth  = unit(0.25, "cm"),
          override.aes = list(shape = 21, colour = "black")
        ),
        shape = "none"
      ) +
      facet_wrap(~stock_group)
    p_legend <- get_legend(
      p_segs
    )
    ggdraw(
      p_segs +
        theme(legend.position = "none")
    ) +
      draw_grob(p_legend, x = x_vec, y = y_vec, width = 0.2, height = 0.2)
  }
)

array_map2 <- cowplot::plot_grid(
  seg_list[[2]], seg_list[[4]], 
  seg_list[[1]], seg_list[[3]], 
  ncol = 2
  )
array_map <- cowplot::plot_grid(
  cali_inset, array_map2, rel_widths = c(0.2, 0.35)
  )


png(here::here("figs", "maps", "array-map.png"),
    height = 5.25, width = 5, units = "in", res = 200)
array_map
dev.off()
