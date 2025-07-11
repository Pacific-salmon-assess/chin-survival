### Compare Scaled Survival Rate Estimates for Chinook and Sockeye
## Used for WA-BC AFS Presentation
# Feb 24, 2025

library(tidyverse)

# scaled CJS posterior estimates from survival_cjs_hier.R
scaled_chin <- readRDS(
  here::here("data", "term_surv_rate.rds")
) %>% 
  mutate(
    study = "DFO",
    species = "chinook"
  )

scale_foo <- function (x, dist) {
  x^(1 / (dist / 100))
}

sox <- data.frame(
  median = c((6 / 252), (12 / 195)),
  se = c(
    sqrt(
      (6 / 252) * (1 - (6 / 252)) / 100
    ),
    sqrt(
      (12 / 195) * (1 - (12 / 195)) / 100
    )
  ),
  study = c("DFO", "UBC")
) %>% 
  mutate(
    stock_group = "Fraser", 
    low = median - se,
    up = median + se,
    species = "sockeye",
    data = "unscaled"
  ) 

scaled_sox <- sox %>% 
  mutate(
    median = scale_foo(median, 470),
    se = sqrt(median * (1 - median) / 100)
    ) %>% 
  mutate(
    low = median - se,
    up = median + se,
    data = "scaled"
  ) %>% 
  rbind(., sox)



surv_comp <- rbind(
  scaled_chin, scaled_sox %>% select(colnames(scaled_chin))
) %>% 
  mutate(
    xx = factor(
      paste(species, stock_group, study, sep = " "),
      labels = c("Cali", "Fraser", "Low Col.", "South Puget", "Up Col.",
                 "Fraser DFO", "Fraser UBC")
      )
  ) 

p1 <- ggplot(surv_comp %>% 
         filter(data == "unscaled",
                !study == "UBC")) +
  geom_pointrange(aes(x = xx, y = median, ymin = low, ymax = up, 
                      fill = species),
                  shape = 21) +
  labs(y = "Cumulative Survival Rate") +
  ggsidekick::theme_sleek() +
  theme(axis.title.x = element_blank(),
        legend.title = element_blank(),
        legend.position = "top") +
  lims(y = c(0, 1))

p2 <- ggplot(surv_comp %>% 
               filter(data == "scaled",
                      !study == "UBC")) +
  geom_pointrange(aes(x = xx, y = median, ymin = low, ymax = up, 
                      fill = species),
                  shape = 21) +
  labs(y = "Cumulative Survival\nRate per 100 km") +
  ggsidekick::theme_sleek() +
  theme(axis.title.x = element_blank(),
        legend.title = element_blank(),
        legend.position = "top") +
  lims(y = c(0, 1))

p3 <- ggplot(surv_comp %>% 
               filter(data == "unscaled")) +
  geom_pointrange(aes(x = xx, y = median, ymin = low, ymax = up, 
                      fill = species),
                  shape = 21) +
  labs(y = "Cumulative Survival Rate") +
  ggsidekick::theme_sleek() +
  theme(axis.title.x = element_blank(),
        legend.title = element_blank(),
        legend.position = "top") +
  lims(y = c(0, 1))

ggsave(
  filename = here::here("figs", "cjs", "sox_chin_comp_unscaled.png"),
  plot = p1,
  width = 6,
  height = 4,
  units = "in",
  device = "png"
)
ggsave(
  filename = here::here("figs", "cjs", "sox_chin_comp_scaled.png"),
  plot = p2,
  width = 6,
  height = 4,
  units = "in",
  device = "png"
)
ggsave(
  filename = here::here("figs", "cjs", "sox_chin_comp_unscaled_plusUBC.png"),
  plot = p3,
  width = 6,
  height = 4,
  units = "in",
  device = "png"
)
