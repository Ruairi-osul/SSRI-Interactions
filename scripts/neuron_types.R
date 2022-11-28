library(tidyverse)
library(patchwork)
library(bbplot)
library(latex2exp)
library(patchwork)
library(emmeans)
library(car)

theme_set(
  bbc_style() +
    theme(
      text=element_text(family="Arial"),
      strip.text.x = element_text(size = 9, angle=0, hjust=0.5, margin = margin(t = 0, r = 0, b = 10, l = 0)),
      strip.text.y = element_text(size = 10, angle=0, hjust=0.5, margin = margin(t = 0, r = 0, b = 10, l = 0)),
      panel.spacing.x = unit(0.2, "lines"),
      panel.spacing.y = unit(1.5, "lines"),
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

# dirs 
fig_dir <- file.path("figs", "spiketrain_props")

# Load Data

neuron_types <- read_csv("data/derived/neuron_types.csv") %>%
  mutate(group = factor(group, levels=c("SAL", "CIT", "DIS"), ordered=T)) %>%
  filter(group %in% c("SAL", "CIT")) %>%
  mutate(neuron_type = factor(neuron_type, levels=c("SR", "SIR", "FF"), ordered = T)) %>%
  mutate(width_basepost = width_basepost / 30)


# Preprocess Data

neuron_types_long <- neuron_types %>%
  pivot_longer(
    cols=c(mean_firing_rate, cv_isi_burst, width_basepost), 
    names_to="metric",
    values_to = "value"
    )


###########3 Neuron Type Stats

### Neuron Type Stats Preprocessing 
neuron_types_long_stats <- neuron_types_long %>%
  drop_na() %>%
  group_by(group, metric, neuron_type) %>%
  summarise(m=mean(value), se=sd(value)/n()) %>%
  ungroup()   

### Neuron Type Stats Plot
p_neuron_type_stats <- neuron_types_long_stats %>%
  ggplot(aes(x=group, y=m, fill=group, ymin=m - se, ymax=m + se)) + 
  geom_bar(
    stat="identity",  
    color="black", 
    width=0.5, 
    position=position_dodge(preserve = "single", width=0.2)
  ) +
  geom_errorbar(
    width=0.28, 
    color='#5c5c5c', 
    position=position_dodge(preserve = "single", width=0.8)
  ) +
  facet_grid(cols=vars(neuron_type), rows=vars(metric), scales="free") +
  scale_fill_manual(values=c(SAL="grey", CIT="black")) +
  labs(y="") +
  theme(
    axis.text.x = element_text(size=9, angle=45),
    strip.text.y = element_blank(),
    strip.text.x = element_text(size=12, margin = margin(b=18))
  )
p_neuron_type_stats

ggsave(file.path(fig_dir, "neuron_props_bars.png"), dpi=300, width=4, height=4)

#### NT Stats Stats

# CV
mod_cv_stats <- neuron_types_long %>%
  filter(metric == "cv_isi_burst") %>%
  drop_na() %>%
  lm(
    value ~ group * neuron_type, 
    data=.,
  )
aov_cv_stats <- anova(mod_cv_stats)
emms_cv_stats <- emmeans(mod_cv_stats, ~ group | neuron_type)
contrasts_cv_stats <- pairs(emms_cv_stats)
aov_cv_stats
contrasts_cv_stats

# MFR
mod_mfr_stats <- neuron_types_long %>%
  filter(metric == "mean_firing_rate") %>%
  drop_na() %>%
  lm(
    value ~ group * neuron_type, 
    data=.,
  )
aov_mfr_stats <- anova(mod_mfr_stats)
emms_mfr_stats <- emmeans(mod_mfr_stats, ~ group | neuron_type)
contrasts_mft_stats <- pairs(emms_mfr_stats)
aov_mfr_stats
emms_mfr_stats
contrasts_mft_stats

# waveform
mod_waveform_stats <- neuron_types_long %>%
  filter(metric == "width_basepost") %>%
  drop_na() %>%
  lm(
    value ~ group * neuron_type, 
    data=.,
  )
aov_waveform_stats <- anova(mod_waveform_stats)
emms_waveform_stats <- emmeans(mod_waveform_stats, ~ group | neuron_type)
contrasts_waveform_stats <- pairs(emms_waveform_stats)
aov_waveform_stats
emms_waveform_stats
contrasts_waveform_stats

########## Prop Neuron Types

## Prop Plot
p_prop_neuron_types <- neuron_types %>%
  drop_na() %>%
  group_by(group, neuron_type) %>%
  summarise(n = n()) %>%
  mutate(freq = (n / sum(n)) * 100 ) %>%
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
  guides(fill="none") +
  labs(y="")

p_prop_neuron_types
ggsave(file.path(fig_dir, "proportion_neuron_types.png"), dpi=300, width=2, height=1.5)

######### Volatility

### Preprocessing

nt_vol <- neuron_types %>%
  select(neuron_id, group, neuron_type) %>%
  left_join(read_csv("data/derived/spiketrain_stats_volitility.csv"))

nt_vol_long <- nt_vol %>%
  pivot_longer(
    cols=c(mean_firing_rate, cv_isi_burst), 
    names_to="metric",
    values_to = "value"
  )

nt_vol_long_stats <- nt_vol_long %>%
  drop_na() %>%
  group_by(group, metric, neuron_type) %>%
  summarise(m=mean(value), se=sd(value)/n()) %>%
  ungroup()   

### Volatility plot

p_vol <- nt_vol_long_stats %>%
  ggplot(aes(x=group, y=m, fill=group, ymin=m - se, ymax=m + se)) + 
  geom_bar(
    stat="identity",  
    color="black", 
    width=0.5, 
    position=position_dodge(preserve = "single", width=0.8)
  ) +
  geom_errorbar(
    width=0.28, 
    color='#5c5c5c', 
    position=position_dodge(preserve = "single", width=0.8)
  ) +
  facet_grid(cols=vars(neuron_type), rows=vars(metric), scales="free") +
  scale_fill_manual(values=c(SAL="grey", CIT="black")) +
  labs(y="") +
  theme(
    axis.text.x = element_text(size=9, angle=45),
    strip.text.y = element_blank(),
    strip.text.x = element_text(size=12, margin = margin(b=18))
  )

p_vol
### Volatility Stats

# Volatility Regularity
mod_cv_vol <- nt_vol_long %>%
  drop_na() %>%
  filter(metric == "cv_isi_burst") %>%
  lm(
    value ~ group * neuron_type, 
    data=.,
  )
aov_cv_vol <- anova(mod_cv_vol)
aov_cv_vol

# Volatility Rate
mod_mfr_vol <- nt_vol_long %>%
  drop_na() %>%
  filter(metric == "mean_firing_rate") %>%
  lm(
    value ~ group * neuron_type, 
    data=.,
  )
aov_mfr_vol <- anova(mod_mfr_vol)
aov_mfr_vol
