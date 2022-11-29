library(tidyverse)
library(patchwork)
library(bbplot)
library(latex2exp)
library(patchwork)
library(emmeans)
library(car)
library(arrow)

theme_set(
  bbc_style() +
    theme(
      text=element_text(family="Arial"),
      strip.text.x = element_text(size = 9, angle=0, hjust=0.5, margin = margin(t = 0, r = 0, b = 10, l = 0)),
      strip.text.y = element_text(size = 10, angle=0, hjust=0.5, margin = margin(t = 0, r = 0, b = 10, l = 0)),
      panel.spacing.x = unit(0.5, "lines"),
      panel.spacing.y = unit(2.5, "lines"),
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



######

data_dir <- "/Users/ruairiosullivan/repos/SSRI Interactions/data/derived"
fig_dir <- file.path("figs", "foot-shock")

neuron_types <- read_csv("data/derived/neuron_types.csv") %>%
  mutate(group = factor(group, levels=c("SAL", "CIT", "DIS"), ordered=T)) %>%
  filter(group %in% c("SAL", "CIT")) %>%
  mutate(neuron_type = factor(neuron_type, levels=c("SR", "SIR", "FF"), ordered = T))


fast_responders <- read_csv(
  file.path(data_dir, "fast_fs_foot_shock_unit_responders.csv")
) %>%
  left_join(neuron_types) %>%
  mutate(
    response=factor(
      fs_fast_response, 
      levels=c("inhibited", "non responder", "activated"),
      labels=c("Shock-\nInhibited", "Not\nShock-\nModulated", "Shock-\nActivated")
      )
  )

psth <- read_parquet(
  file.path(data_dir, "baseshock_counts_psth.parquet")
  ) %>%
  mutate(Block = factor(if_else(between(bin, left=0.05, right=0.2), "Pre", "Shock"), levels=c("Pre", "Shock"), ordered=T)) %>%
  mutate(group = factor(group, levels=c("SAL", "CIT", "DIS"), ordered=T)) %>%
  filter(group %in% c("SAL", "CIT")) %>%
  drop_na() 

mod <- lm(
  zcounts ~ window * neuron_type * group,
  data=psth
)


anova(mod)
emms <- emmeans(
  mod, ~ window | group | neuron_type
)

emms
prepost_by_type_by_group <- pairs(emms)
prepost_across_groups_by_type <- pairs(prepost_by_type_by_group, by="neuron_type")
prepost_across_groups_across_types <- pairs(prepost_across_groups_by_type, by=NULL)

prepost_by_type_by_group  # were neurons more active in post
prepost_across_groups_by_type  # within each neuron type, were effects different by drug
prepost_across_groups_across_types  # were the effects of drug larger among some neuron types

###### PLOTS

# Proportion 


fast_responders %>%
  drop_na() %>%
  group_by(group, neuron_type, response) %>%
  summarise(n = n()) %>%
  mutate(freq = ( n / sum(n) ) * 100) %>%
  mutate(neuron_type = factor(neuron_type, levels=c("SR", "SIR", "FF"), ordered = T)) %>%
  ggplot(
    aes(x=response, y=freq, fill=group)
  ) +
  geom_bar(
    stat="identity",  
    color="black", 
    width=0.5, 
    position=position_dodge(preserve = "single", width=0.59)
  ) +
  scale_fill_manual(values=c(SAL="grey", CIT="black")) +
  labs(y="") +
  facet_grid(rows=vars(neuron_type)) +
  guides(fill="none") +
  lims(y=c(0, 100)) + theme(axis.text.x = element_text(angle=22.5)) -> p_prop_neuron_types

p_prop_neuron_types

ggsave(file.path(fig_dir, "fast_fs_w1_prop_responders.png"), dpi=300, width=2.8, height=4)

########################## SECOND WINDOW


fast_responders_second_window <- read_csv(
  file.path(data_dir, "fast_fs_foot_shock_unit_responders_second_window.csv")
) %>%
  left_join(neuron_types) %>%
  mutate(
    response=factor(
      response_second_window, 
        levels=c("inhibited", "non responder", "activated"),
        labels=c("Shock-\nInhibited", "Not\nShock-\nModulated", "Shock-\nActivated")
        )
  )

psth <- read_parquet(
  file.path(data_dir, "baseshock_counts_psth_second_window.parquet")
) %>%
  left_join(neuron_types) %>%
  mutate(Block = factor(if_else(between(bin, left=0.5, right=0.8), "Pre", "Shock"), levels=c("Pre", "Shock"), ordered=T)) %>%
  filter(group %in% c("SAL", "CIT")) %>%
  drop_na() 

mod <- lm(
  zcounts ~ window * neuron_type * group,
  data=psth
)


anova(mod)
emms <- emmeans(
  mod, ~ window | group | neuron_type
)

emms
prepost_by_type_by_group <- pairs(emms)
prepost_across_groups_by_type <- pairs(prepost_by_type_by_group, by="neuron_type")
prepost_across_groups_across_types <- pairs(prepost_across_groups_by_type, by=NULL)

prepost_by_type_by_group  # were neurons more active in post
prepost_across_groups_by_type  # within each neuron type, were effects different by drug
prepost_across_groups_across_types  # were the effects of drug larger among some neuron types

###### PLOTS

# Proportion 


fast_responders_second_window %>%
  drop_na() %>%
  group_by(group, neuron_type, response) %>%
  summarise(n = n()) %>%
  mutate(freq = ( n / sum(n) ) * 100) %>%
  mutate(neuron_type = factor(neuron_type, levels=c("SR", "SIR", "FF"), ordered = T)) %>%
  ggplot(
    aes(x=response, y=freq, fill=group)
  ) +
  geom_bar(
    stat="identity",  
    color="black", 
    width=0.5, 
    position=position_dodge(preserve = "single", width=0.59)
  ) +
  scale_fill_manual(values=c(SAL="grey", CIT="black")) +
  labs(y="") +
  guides(fill="none") +
  facet_grid(rows=vars(neuron_type)) +
  lims(y=c(0, 100)) + theme(axis.text.x = element_text(angle=22.5)) -> p_prop_neuron_types

p_prop_neuron_types
ggsave(file.path(fig_dir, "fast_fs_w2_prop_responders.png"), dpi=300, width=2.8, height=4)

