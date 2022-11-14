library(tidyverse)
library(patchwork)
library(bbplot)
library(latex2exp)
library(lme4)
library(lmerTest)
library(patchwork)
library(emmeans)
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


###############


####### Load data

data_dir <- "/Users/ruairiosullivan/repos/SSRI Interactions/data/derived"

neuron_types <- read_csv("data/derived/neuron_types.csv") %>%
  mutate(group = factor(group, levels=c("SAL", "CIT", "DIS"), ordered=T)) %>%
  filter(group %in% c("SAL", "CIT")) %>%
  mutate(neuron_type = factor(neuron_type, levels=c("SR", "SIR", "FF"), ordered = T)) %>%
  mutate(width_basepost = width_basepost / 30)


read_csv(
  file.path(data_dir, "chal_responders.csv")
) %>%
  select(-neuron_type, -group) %>%
  left_join(select(neuron_types, group, neuron_type, neuron_id)) %>%
  mutate(
    response=factor(response, levels=c("non responder", "inhibited", "activated"))
  ) -> cit_responders


read_parquet(
  file.path(data_dir, "chal_binned.parquet")
  ) %>%
  filter(bin < 1800) %>%
  mutate(Block = factor(if_else(bin <= 0, "Pre", "Drug"), levels=c("Pre", "Drug"), ordered=T)) %>%
  select(-neuron_type, -Block) %>%
  left_join(select(neuron_types, group, neuron_type, neuron_id)) %>%
  filter(group %in% c("SAL", "CIT")) %>%
  drop_na() -> cit_counts


########### CIT MODELS

mod <- lm(
  zcounts ~ block * neuron_type * group,
  data=cit_counts
)

anova(mod)
emms <- emmeans(
  mod, ~ block | group | neuron_type
)

emms
prepost_by_type_by_group <- pairs(emms)
prepost_across_groups_by_type <- pairs(prepost_by_type_by_group, by="neuron_type")
prepost_across_groups_across_types <- pairs(prepost_across_groups_by_type, by=NULL)

prepost_by_type_by_group  # were neurons more active in post
prepost_across_groups_by_type  # within each neuron type, were effects different by drug
prepost_across_groups_across_types  # were the effects of drug larger among some neuron types




########## WAY #############

read_csv(
  file.path(data_dir, "way_responders.csv")
) %>%
  select(-neuron_type, -group) %>%
  left_join(select(neuron_types, group, neuron_type, neuron_id)) %>%
  mutate(
    response=factor(response, levels=c("non responder", "inhibited", "activated"))
  ) -> way_responders


read_parquet(
  file.path(data_dir, "way_binned.parquet")
) %>%
  filter(bin < 1800) %>%
  mutate(Block = factor(if_else(bin <= 0, "Pre", "Drug"), levels=c("Pre", "Drug"), ordered=T)) %>%
  select(-neuron_type, -Block) %>%
  left_join(select(neuron_types, group, neuron_type, neuron_id)) %>%
  filter(group %in% c("SAL", "CIT")) %>%
  drop_na() -> way_counts

############################ MOD

########### CIT MODELS

mod <- lm(
  zcounts ~ block * neuron_type * group,
  data=way_counts
)

anova(mod)
emms <- emmeans(
  mod, ~ block | group | neuron_type
)

emms
prepost_by_type_by_group <- pairs(emms)
prepost_across_groups_by_type <- pairs(prepost_by_type_by_group, by="neuron_type")
prepost_across_groups_across_types <- pairs(prepost_across_groups_by_type, by=NULL)

prepost_by_type_by_group  # were neurons more active in post
prepost_across_groups_by_type  # within each neuron type, were effects different by drug
prepost_across_groups_across_types  # were the effects of drug larger among some neuron types




###### Proportion plots

cit_responders %>%
  drop_na() %>%
  group_by(group, neuron_type, response) %>%
  summarise(n = n()) %>%
  mutate(freq = ( n / sum(n) ) * 100) %>%
  mutate(neuron_type = factor(neuron_type, levels=c("SR", "SIR", "FF"), ordered = T)) %>%
  filter(response == "inhibited") %>%
  mutate(response = factor(
    response, 
    levels=c("non responder", "inhibited", "activated"),
    labels=c("Not\nCitalopram-\nModulated", "Citalopram-\nInhibited", "Citalopram-\nActivated"),
    ordered = T
  )) %>%
  ggplot(
    aes(x=neuron_type, y=freq, fill=group)
  ) +
  geom_bar(
    stat="identity",  
    color="black", 
    width=0.5, 
    position=position_dodge(preserve = "single", width=0.8)
  ) +
  scale_fill_manual(values=c(SAL="grey", CIT="black")) +
  labs(y="") +
  lims(y=c(0, 100)) + theme(axis.text.x = element_text(angle=22.5))


way_responders %>%
  drop_na() %>%
  group_by(group, neuron_type, response) %>%
  summarise(n = n()) %>%
  mutate(freq = ( n / sum(n) ) * 100) %>%
  mutate(neuron_type = factor(neuron_type, levels=c("SR", "SIR", "FF"), ordered = T)) %>%
  filter(response == "activated") %>%
  mutate(response = factor(
    response, 
    levels=c("non responder", "inhibited", "activated"),
    labels=c("Not\nWay-\nModulated", "Way-\nInhibited", "Way-\nActivated"),
    ordered = T
  )) %>%
  ggplot(
    aes(x=neuron_type, y=freq, fill=group)
  ) +
  geom_bar(
    stat="identity",  
    color="black", 
    width=0.5, 
    position=position_dodge(preserve = "single", width=0.8)
  ) +
  scale_fill_manual(values=c(SAL="grey", CIT="black")) +
  labs(y="") +
  lims(y=c(0, 100)) + theme(axis.text.x = element_text(angle=22.5))
