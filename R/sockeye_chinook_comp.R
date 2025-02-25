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

scaled_sox <- data.frame(
  stock_group = "Fraser", 
  med = scale_foo(c((6 / 252), (12 / 195)), 470),
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
    lo = med - se,
    up = med + se,
    species = "sockeye"
  ) %>% 
  select(
    colnames(scaled_chin)
  )


surv_comp <- rbind(
  scaled_chin, scaled_sox
) %>% 
  mutate(
    xx = factor(
      paste(species, stock_group, study, sep = " "),
      labels = c("Cali", "Fraser", "Low Col.", "South Puget", "Up Col.",
                 "Fraser DFO", "Fraser UBC")
      )
  ) %>% 
  filter(!study == "UBC") %>% 
  ggplot(.) +
  geom_pointrange(aes(x = xx, y = med, ymin = lo, ymax = up, fill = species),
                  shape = 21) +
  labs(y = "Cumulative Survival\nRate per 100 km") +
  ggsidekick::theme_sleek() +
  theme(axis.title.x = element_blank(),
        legend.title = element_blank(),
        legend.position = "top") +
  lims(y = c(0, 1))

surv_comp2 <- rbind(
  scaled_chin, scaled_sox
) %>% 
  mutate(
    xx = factor(
      paste(species, stock_group, study, sep = " "),
      labels = c("Cali", "Fraser", "Low Col.", "South Puget", "Up Col.",
                 "Fraser DFO", "Fraser UBC")
    )
  ) %>% 
  ggplot(.) +
  geom_pointrange(aes(x = xx, y = med, ymin = lo, ymax = up, fill = species),
                  shape = 21) +
  labs(y = "Cumulative Survival\nRate per 100 km") +
  ggsidekick::theme_sleek() +
  theme(axis.title.x = element_blank(),
        legend.title = element_blank(),
        legend.position = "top") +
  lims(y = c(0, 1))

ggsave(
  filename = here::here("figs", "cjs", "sox_chin_comp.png"),
  plot = surv_comp,
  width = 6,
  height = 4,
  units = "in",
  device = "png"
)
ggsave(
  filename = here::here("figs", "cjs", "sox_chin_comp2.png"),
  plot = surv_comp2,
  width = 6,
  height = 4,
  units = "in",
  device = "png"
)
