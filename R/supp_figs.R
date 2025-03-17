### Survival Logistic Regression Supp Figures
## January 8, 025
## Supplementary figures showing input data to GLM

library(tidyverse)
library(ggplot2)


det_dat1 <- readRDS(here::here("data", "surv_log_reg_data.rds")) %>% 
  filter(
    stage_1 == "mature"
  ) %>% 
  mutate(
    stock_group = factor(
      stock_group, 
      levels = c(
        "Cali", "Low Col.", "Up Col.", "WA_OR", "WCVI", "ECVI", 
        "Fraser Spr. Yr.", "Fraser Sum. Yr.", "Fraser Sum. 4.1", "Fraser Fall", 
        "North Puget", "South Puget"
      ))
  )


## stock-specific traits
stock_pal <- c("#b30000", "#6baed6", "#08306b", "#fec44f", "#ccece6", "#238b45",
               "#bcbddc", "#807dba", "#54278f", "#3f007d", "#fde0dd", "#f768a1")
names(stock_pal) <- levels(det_dat1$stock_group)


png(here::here("figs", "supp", "stock_fl.png"), 
    height = 4, width = 5.5, units = "in", res = 200)
ggplot() +
  geom_boxplot(data = det_dat1,
               aes(x = stock_group, y = fl, fill = stock_group)) +
  scale_fill_manual(values = stock_pal) +
  labs(y = "Fork Length (cm)") +
  ggsidekick::theme_sleek() +
  theme(axis.title.x = element_blank(),
        legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_text(
    data = det_dat1 %>% 
      filter(!is.na(fl)) %>% 
      group_by(stock_group) %>% 
      summarize(
        n = length(unique(vemco_code))
      ),
    aes(x = stock_group, y = -Inf, label = n),
    vjust = -0.1
  )
dev.off()


png(here::here("figs", "supp", "stock_lipid.png"), 
    height = 4, width = 5.5, units = "in", res = 200)
ggplot() +
  geom_boxplot(data = det_dat1,
               aes(x = stock_group, y = lipid, fill = stock_group)) +
  scale_fill_manual(values = stock_pal) +
  labs(y = "Lipid Content (% wet weight)") +
  ggsidekick::theme_sleek() +
  theme(axis.title.x = element_blank(),
        legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_text(
    data = det_dat1 %>% 
      filter(!is.na(lipid)) %>% 
      group_by(stock_group) %>% 
      summarize(
        n = length(unique(vemco_code))
      ),
    aes(x = stock_group, y = -Inf, label = n),
    vjust = -0.1
  )
dev.off()


png(here::here("figs", "supp", "stock_date.png"), 
    height = 4, width = 5.5, units = "in", res = 200)
ggplot() +
  geom_boxplot(data = det_dat1,
               aes(x = stock_group, y = year_day, fill = stock_group)) +
  scale_fill_manual(values = stock_pal) +
  labs(y = "Tagging Date (year day)") +
  ggsidekick::theme_sleek() +
  theme(axis.title.x = element_blank(),
        legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_text(
    data = det_dat1 %>% 
      group_by(stock_group) %>% 
      summarize(
        n = length(unique(vemco_code))
      ),
    aes(x = stock_group, y = -Inf, label = n),
    vjust = -0.1
  )
dev.off()


png(here::here("figs", "supp", "stock_cyer.png"), 
    height = 4, width = 5.5, units = "in", res = 200)
ggplot() +
  geom_boxplot(data = det_dat1,
               aes(x = stock_group, y = focal_er, fill = stock_group)) +
  scale_fill_manual(values = stock_pal) +
  labs(y = "Exploitation Rate") +
  ggsidekick::theme_sleek() +
  theme(axis.title.x = element_blank(),
        legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_text(
    data = det_dat1 %>% 
      filter(!is.na(focal_er)) %>% 
      group_by(stock_group) %>% 
      summarize(
        n = length(unique(vemco_code))
      ),
    aes(x = stock_group, y = -Inf, label = n),
    vjust = -0.1
  )
dev.off()


## injury scores
png(here::here("figs", "supp", "injury_bar.png"), 
    height = 4, width = 4.5, units = "in", res = 200)
ggplot(det_dat1) +
  geom_bar(aes(x = adj_inj)) +
  labs(x = "Injury Score", y = "Tagged Individuals") +
  ggsidekick::theme_sleek()
dev.off()


