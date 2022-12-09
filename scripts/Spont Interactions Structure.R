library(tidyverse)
library(emmeans)
library(bbplot)
library(extrafont)
library(latex2exp)
library(patchwork)
library(arrow)
library(mgcv)

theme_set(
  bbc_style() +
    theme(
      text=element_text(family="Arial"),
      strip.text.x = element_text(size = 9, angle=0, hjust=0.5, margin = margin(t = 0, r = 0, b = 10, l = 0)),
      strip.text.y = element_text(size = 9, angle=0, hjust=0.5, margin = margin(t = 0, r = 0, b = 10, l = 0)),
      panel.spacing = unit(1, "lines"),
      axis.text.x = element_text(size=9, angle=0),
      axis.text.y = element_text(size=9),
      legend.text = element_text(size=9),
      legend.key.size = unit(0.5, 'cm'),
      legend.position = "right",
      axis.line=element_line(),
      axis.title.y = element_text(size=9, angle=90, margin = margin(t = 0, r = 5, b = 0, l = 5)),
      panel.border = element_blank(),
      strip.background = element_blank(),
      plot.title = element_text(size=9, face="plain", margin = margin(t = 0, r = 0, b = 5, l = 0))
    )
)


read_parquet("data/derived/graph/spont - graph.parquet") %>%
  filter(group != "DIS") %>%
  filter(shuffle == F) %>%
  filter(bin_width ==1) %>%
  mutate(
    group = factor(
      group, 
      levels=c("SAL", "CIT"), 
      labels=c("SAL",  "CIT"),
    )
  ) -> df_graph

swp_mod <- lm(swp ~ group, data=df_graph)
summary(swp_mod)

clust_mod <- lm(avg_clust ~ group, data=df_graph)
summary(clust_mod)

pl_mod <- lm(avg_path_len ~ group, data=df_graph)
summary(pl_mod)

dens_mod <- lm(density ~ group, data=df_graph)
summary(dens_mod)

pairs(swp_mod_ems)

