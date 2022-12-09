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


read_parquet("data/derived/corrs/spont - corr.parquet") %>%
  filter(group != "DIS") %>%
  mutate(
    nt_comb = factor(
      nt_comb, 
      levels=c("SR-SR", "SR-SIR", "SR-FF", "SIR-SIR", "SIR-FF", "FF-FF")
      ),
    group = factor(
      group, 
      levels=c("SAL", "CIT"), 
      labels=c("SAL",  "CIT"),
    ),
    distance = distance / 1000,
    is_neg = corr < 0,
    is_neg_sig = (corr < 0) & (sig == T),
    is_pos = corr > 0,
    is_pos_sig = (corr > 0) & (sig == T),
    mag = abs(corr),
    rec = if_else(corr < 0, 0, corr)
  ) %>%
  filter(bin_width == 1, shuffle==F) -> df_corr

prop_neg_mod <- glm(
  is_neg_sig ~ nt_comb + group + nt_comb:group + distance,
  data=filter(df_corr),
  family=binomial()
)
prop_neg_ems <- emmeans(
  prop_neg_mod, 
  specs = ~ group | nt_comb,
  type="response"
)
prop_neg_tab <- as_tibble(prop_neg_ems) %>%
  mutate(prop = str_c(round(prob * 100, 2),  round(SE * 100, 2), sep=" +-")) %>%
  select(group, nt_comb, prop) %>%
  pivot_wider(names_from=c(group), values_from=c(prop))
prop_neg_tab <- pairs(prop_neg_ems) %>% 
  as_tibble() %>%
  mutate(`Odds Ratio` = str_c(round(odds.ratio, 2),  round(SE, 2), sep=" +-")) %>%
  mutate(p = round(p.value, 4)) %>%
  select(nt_comb, `Odds Ratio`, p) %>%
  right_join(prop_neg_tab) %>%
  mutate(`Neuron Type Combination`=nt_comb) %>%
  select(`Neuron Type Combination`, SAL, CIT, `Odds Ratio`, p)
prop_neg_plot_all <- prop_neg_ems %>%
  as_tibble() %>%
  ggplot(aes(x=group, y=prob, fill=group, ymin=prob - SE, ymax=prob + SE)) + 
  geom_bar(stat="identity",  position=position_dodge(preserve = "single", width=0.4), color="black", width=0.5) +
  geom_errorbar(width=0.28, color='#5c5c5c', position=position_dodge(preserve = "single", width=0.8)) +
  scale_fill_manual(values=c(SAL="grey", CIT="black")) +
  facet_grid(cols=vars(nt_comb)) +
  guides(fill="none") +
  labs(y="Proportion of\nPairs With Negative\nCorrelations")
prop_neg_plot_sr <- prop_neg_ems %>%
  as_tibble() %>%
  filter(nt_comb == 'SR-SR') %>%
  ggplot(aes(x=group, y=prob, fill=group, ymin=prob - SE, ymax=prob + SE)) + 
  geom_bar(stat="identity",  position=position_dodge(preserve = "single", width=0.4), color="black", width=0.5) +
  geom_errorbar(width=0.28, color='#5c5c5c', position=position_dodge(preserve = "single", width=0.8)) +
  scale_fill_manual(values=c(SAL="grey", CIT="black")) +
  facet_grid(cols=vars(nt_comb)) +
  guides(fill="none") +
  labs(y="Proportion of\nPairs Negatively\nCorrelated") +
  lims(y=c(0, 0.3))


anova(prop_neg_mod, test="Chisq")
prop_neg_tab
prop_neg_plot_sr
# prop_neg_plot_all

# Prop pos

prop_pos_mod <- glm(
  is_pos ~ nt_comb + group + nt_comb:group + distance,
  data=df_corr,
  family=binomial()
)
prop_pos_ems <- emmeans(
  prop_pos_mod, 
  specs = ~ group | nt_comb,
  type="response"
)
prop_pos_tab <- as_tibble(prop_pos_ems) %>%
  mutate(prop = str_c(round(prob * 100, 2),  round(SE * 100, 2), sep=" +-")) %>%
  select(group, nt_comb, prop) %>%
  pivot_wider(names_from=c(group), values_from=c(prop))
prop_pos_tab <- pairs(prop_pos_ems) %>% 
  as_tibble() %>%
  mutate(`Odds Ratio` = str_c(round(odds.ratio, 2),  round(SE, 2), sep=" +-")) %>%
  mutate(p = round(p.value, 4)) %>%
  select(nt_comb, `Odds Ratio`, p) %>%
  right_join(prop_pos_tab) %>%
  mutate(`Neuron Type Combination`=nt_comb) %>%
  select(`Neuron Type Combination`, SAL, CIT, `Odds Ratio`, p)
prop_pos_plot_all <- prop_pos_ems %>%
  as_tibble() %>%
  ggplot(aes(x=group, y=prob, fill=group, ymin=prob - SE, ymax=prob + SE)) + 
  geom_bar(stat="identity",  position=position_dodge(preserve = "single", width=0.4), color="black", width=0.5) +
  geom_errorbar(width=0.28, color='#5c5c5c', position=position_dodge(preserve = "single", width=0.8)) +
  scale_fill_manual(values=c(SAL="grey", CIT="black")) +
  facet_grid(cols=vars(nt_comb)) +
  guides(fill="none") +
  labs(y="Proportion of\nPairs Negatively\nCorrelated")
prop_pos_plot_sr <- prop_pos_ems %>%
  as_tibble() %>%
  filter(nt_comb == 'SR-SR') %>%
  ggplot(aes(x=group, y=prob, fill=group, ymin=prob - SE, ymax=prob + SE)) + 
  geom_bar(stat="identity",  position=position_dodge(preserve = "single", width=0.4), color="black", width=0.5) +
  geom_errorbar(width=0.28, color='#5c5c5c', position=position_dodge(preserve = "single", width=0.8)) +
  scale_fill_manual(values=c(SAL="grey", CIT="black")) +
  facet_grid(cols=vars(nt_comb)) +
  guides(fill="none") +
  labs(y="Proportion of\nPairs Positively\nCorrelated") +
  lims(y=c(0, 0.65))

anova(prop_pos_mod, test = "Chisq")
prop_pos_tab
# prop_pos_plot_all
prop_pos_plot_sr

#### Mag neg

neg_mag_mod <- lm(
  mag ~ nt_comb + group + nt_comb:group + distance,
  data=filter(df_corr, is_neg_sig),
)
neg_mag_ems <- emmeans(
  neg_mag_mod,
  specs = ~ group | nt_comb,
)
neg_mag_tab <- as_tibble(neg_mag_ems) %>%
  mutate(`Spike Count Correlation` = str_c(round(emmean, 2),  round(SE, 2), sep=" +-")) %>%
  select(group, nt_comb, `Spike Count Correlation`) %>%
  pivot_wider(names_from=c(group), values_from=c(`Spike Count Correlation`))
neg_mag_tab <- pairs(neg_mag_ems) %>% as_tibble() %>%
  mutate(`Difference` = str_c(round(estimate, 2),  round(SE, 2), sep=" +-")) %>%
  mutate(p = round(p.value, 4)) %>%
  select(nt_comb, `Difference`, p) %>%
  right_join(neg_mag_tab) %>%
  mutate(`Neuron Type Combination`=nt_comb) %>%
  select(`Neuron Type Combination`, SAL, CIT, `Difference`, p)
neg_mag_plot_all <- neg_mag_ems %>%
  as_tibble() %>%
  ggplot(aes(x=group, y=emmean, fill=group, ymin=emmean - SE, ymax=emmean + SE)) +
  geom_bar(stat="identity",  position=position_dodge(preserve = "single", width=0.4), color="black", width=0.5) +
  geom_errorbar(width=0.28, color='#5c5c5c', position=position_dodge(preserve = "single", width=0.8)) +
  scale_fill_manual(values=c(SAL="grey", CIT="black")) +
  facet_grid(cols=vars(nt_comb)) +
  guides(fill="none") +
  labs(y="Magnitude of\nNegative Correlations") +
  theme(
    axis.text.x = element_text(size=10, angle=0),
  )
neg_mag_plot_sr <- neg_mag_ems %>%
  as_tibble() %>%
  filter(nt_comb == 'SR-SR') %>%
  ggplot(aes(x=group, y=emmean, fill=group, ymin=emmean - SE, ymax=emmean + SE)) +
  geom_bar(stat="identity",  position=position_dodge(preserve = "single", width=0.4), color="black", width=0.5) +
  geom_errorbar(width=0.28, color='#5c5c5c', position=position_dodge(preserve = "single", width=0.8)) +
  scale_fill_manual(values=c(SAL="grey", CIT="black")) +
  facet_grid(cols=vars(nt_comb)) +
  guides(fill="none") +
  labs(y="Magnitude of\nNegative Correlations") +
  theme(
    axis.text.x = element_text(size=10, angle=0),
  ) +
  lims(y=c(0, 0.35))

anova(neg_mag_mod, test = "Chisq")
pairs(neg_mag_ems)
neg_mag_tab
# neg_mag_plot_all
neg_mag_plot_sr

pos_mag_mod <- lm(
  rec ~ nt_comb + group + nt_comb:group + distance,
  data=filter(df_corr, is_pos_sig),
)
pos_mag_ems <- emmeans(
  pos_mag_mod,
  specs = ~ group | nt_comb,
)
pos_mag_tab <- as_tibble(pos_mag_ems) %>%
  mutate(`Spike Count Correlation` = str_c(round(emmean, 2),  round(SE, 2), sep=" +-")) %>%
  select(group, nt_comb, `Spike Count Correlation`) %>%
  pivot_wider(names_from=c(group), values_from=c(`Spike Count Correlation`))
pos_mag_tab <- pairs(pos_mag_ems) %>% as_tibble() %>%
  mutate(`Difference` = str_c(round(estimate, 2),  round(SE, 2), sep=" +-")) %>%
  mutate(p = round(p.value, 4)) %>%
  select(nt_comb, `Difference`, p) %>%
  right_join(pos_mag_tab) %>%
  mutate(`Neuron Type Combination`=nt_comb) %>%
  select(`Neuron Type Combination`, SAL, CIT, `Difference`, p)
pos_mag_plot_all <- pos_mag_ems %>%
  as_tibble() %>%
  ggplot(aes(x=group, y=emmean, fill=group, ymin=emmean - SE, ymax=emmean + SE)) +
  geom_bar(stat="identity",  position=position_dodge(preserve = "single", width=0.4), color="black", width=0.5) +
  geom_errorbar(width=0.28, color='#5c5c5c', position=position_dodge(preserve = "single", width=0.8)) +
  scale_fill_manual(values=c(SAL="grey", CIT="black")) +
  facet_grid(cols=vars(nt_comb)) +
  guides(fill="none") +
  labs(y="Spike Count\nCorrelation (R)") +
  theme(
    axis.text.x = element_text(size=10, angle=0),
  )
pos_mag_plot_sr <- pos_mag_ems %>%
  as_tibble() %>%
  filter(nt_comb == 'SR-SR') %>%
  ggplot(aes(x=group, y=emmean, fill=group, ymin=emmean - SE, ymax=emmean + SE)) +
  geom_bar(stat="identity",  position=position_dodge(preserve = "single", width=0.4), color="black", width=0.5) +
  geom_errorbar(width=0.28, color='#5c5c5c', position=position_dodge(preserve = "single", width=0.8)) +
  scale_fill_manual(values=c(SAL="grey", CIT="black")) +
  facet_grid(cols=vars(nt_comb)) +
  guides(fill="none") +
  labs(y="Spike Count\nCorrelation (R)") +
  theme(
    axis.text.x = element_text(size=10, angle=0),
  )


anova(pos_mag_mod)
pairs(pos_mag_ems)
pos_mag_tab
pos_mag_plot_all



###############

read_parquet("data/derived/corrs/spont - pcup.parquet") %>%
  filter(group != "DIS") %>%
  mutate(
      group = factor(
      group, 
      levels=c("SAL", "CIT"), 
      labels=c("SAL",  "CIT"),
    ),
  is_neg = cc < 0,
  is_neg_sig = (sig) & (cc < 0),
  is_pos_sig = (sig) & (cc > 0),
  mag = abs(cc),
  rec = if_else(cc < 0, 0, cc),
  neuron_type = factor(neuron_type, levels=c("SR", "SIR", "FF"))
  ) %>%
  filter(shuffle == F, bin_width == 1)-> df_pcup

pcup_mod <- lm(
  mag ~ group * neuron_type,
  data=filter(df_pcup)
)
anova(pcup_mod)
emms_pcup <- emmeans(
  pcup_mod, 
  specs= ~ group | neuron_type
)
pairs(emms_pcup)

pcup_plot_all <- emms_pcup %>%
  as_tibble() %>%
  ggplot(aes(x=group, y=emmean, fill=group, ymin=emmean - SE, ymax=emmean + SE)) + 
  geom_bar(stat="identity",  position=position_dodge(preserve = "single", width=0.4), color="black", width=0.5) +
  geom_errorbar(width=0.28, color='#5c5c5c', position=position_dodge(preserve = "single", width=0.8)) +
  scale_fill_manual(values=c(SAL="grey", CIT="black")) +
  facet_grid(cols=vars(neuron_type)) +
  labs(y="Population\nCoupling (R)") +
  theme(
    axis.text.x = element_text(size=10, angle=0),
  )

pcup_plot_sr <- emms_pcup %>%
  as_tibble() %>%
  filter(neuron_type == "SR") %>%
  ggplot(aes(x=group, y=emmean, fill=group, ymin=emmean - SE, ymax=emmean + SE)) + 
  geom_bar(stat="identity",  position=position_dodge(preserve = "single", width=0.4), color="black", width=0.5) +
  geom_errorbar(width=0.28, color='#5c5c5c', position=position_dodge(preserve = "single", width=0.8)) +
  scale_fill_manual(values=c(SAL="grey", CIT="black")) +
  facet_grid(cols=vars(neuron_type)) +
  labs(y="Population\nCoupling (R)") +
  theme(
    axis.text.x = element_text(size=10, angle=0),
  )
pcup_plot_sr

############# Partial correlation -> show reduction in SR-SR


read_parquet("data/derived/corrs/spont - pcorr.parquet") %>%
  filter(group != "DIS") %>%
  mutate(
    nt_comb = factor(
      nt_comb, 
      levels=c("SR-SR", "SR-SIR", "SR-FF", "SIR-SIR", "SIR-FF", "FF-FF")
    ),
    group = factor(
      group, 
      levels=c("SAL", "CIT"), 
      labels=c("SAL",  "CIT"),
    ),
    distance = distance / 1000,
    is_neg = pcorr < 0,
    is_pos = pcorr > 0,
    mag = abs(pcorr),
    rec = if_else(pcorr < 0, 0, pcorr)
  ) %>%
  filter(bin_width == 1, shuffle==F) -> df_pcorr

pcorr_mod <- lm(
  mag ~ nt_comb + group + nt_comb:group + distance,
  data=filter(df_pcorr),
)

anova(pcorr_mod)
emms_pcorr <- emmeans(
  pcorr_mod, 
  specs = ~ group | nt_comb,
  type="response"
)
pairs(emms_pcorr)

emms_pcorr %>%
  as_tibble() %>%
  ggplot(aes(x=group, y=emmean, fill=group, ymin=emmean - SE, ymax=emmean + SE)) +
  geom_bar(stat="identity",  position=position_dodge(preserve = "single", width=0.4), color="black", width=0.5) +
  geom_errorbar(width=0.28, color='#5c5c5c', position=position_dodge(preserve = "single", width=0.8)) +
  scale_fill_manual(values=c(SAL="grey", CIT="black")) +
  facet_grid(cols=vars(nt_comb)) +
  guides(fill="none") +
  labs(y="Partial Corr") +
  theme(
    axis.text.x = element_text(size=10, angle=0),
  )


prop_mod <- glm(
  sig ~ nt_comb + group + nt_comb:group + distance,
  data=filter(df_pcorr, bin_width==1, shuffle==F),
  family=binomial()
)

anova(prop_mod)
emms_prop <- emmeans(
  prop_mod, 
  specs = ~ group | nt_comb,
  type="response"
)
pairs(emms_prop)

pcorr_prop_all <- emms_prop %>%
  as_tibble() %>%
  ggplot(aes(x=group, y=prob, fill=group, ymin=prob - SE, ymax=prob + SE)) + 
  geom_bar(stat="identity",  position=position_dodge(preserve = "single", width=0.4), color="black", width=0.5) +
  geom_errorbar(width=0.28, color='#5c5c5c', position=position_dodge(preserve = "single", width=0.8)) +
  scale_fill_manual(values=c(SAL="grey", CIT="black")) +
  facet_grid(cols=vars(nt_comb)) +
  guides(fill="none") +
  labs(y="Proportion of\nPairs With Partial\nCorrelation")

pcorr_prop_sr <- emms_prop %>%
  as_tibble() %>%
  filter(nt_comb == 'SR-SR') %>%
  ggplot(aes(x=group, y=prob, fill=group, ymin=prob - SE, ymax=prob + SE)) + 
  geom_bar(stat="identity",  position=position_dodge(preserve = "single", width=0.4), color="black", width=0.5) +
  geom_errorbar(width=0.28, color='#5c5c5c', position=position_dodge(preserve = "single", width=0.8)) +
  scale_fill_manual(values=c(SAL="grey", CIT="black")) +
  facet_grid(cols=vars(nt_comb)) +
  guides(fill="none") +
  labs(y="Proportion of\nPairs With Partial\nCorrelation")


##### All together

prop_neg_tab
prop_pos_tab


(prop_neg_plot_sr + pcorr_prop_sr + pcup_plot_sr) / (pos_mag_plot_all)



  
  
  df_node
