library(tidyverse)
library(emmeans)
library(bbplot)
library(extrafont)
library(latex2exp)
library(patchwork)
library(arrow)
library(mgcv)

df_responders <- read_parquet("data/derived/graph/spont - responders.parquet") %>%
  filter(group %in% c("SAL", "CIT")) %>%
  mutate(
    response_fs_slow = factor(
      response_fs_slow, 
      levels=c("inhibited", "non responder", "activated"), 
      labels=c("Shock-Activated", "Not Shock-Modulated", "Shock-Inhibited"),
      ),
    fs_fast_response = factor(
      fs_fast_response, 
      levels=c("inhibited", "non responder", "activated"), 
      labels=c("Shock-Activated", "Not Shock-Modulated", "Shock-Inhibited"),
    ),
    response_second_window = factor(
      response_second_window, 
      levels=c("inhibited", "non responder", "activated"), 
      labels=c("Shock-Activated", "Not Shock-Modulated", "Shock-Inhibited"),
    ),
    group = factor(
      group, 
      levels=c("SAL", "CIT"), 
      labels=c("Shock-Activated",  "Shock-Inhibited"),
    )
    )

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


df_node <- read_parquet("data/derived/graph/spont - node.parquet") %>% 
  filter(group %in% c("SAL", "CIT")) %>%
  left_join(select(df_responders, -session_name, -group, -neuron_type), by="neuron_id") %>% 
  distinct()  %>%
  mutate(
    neuron_type = factor(neuron_type, levels=c("SR", "SIR", "FF")),
    group = factor(
      group, 
      levels=c("SAL", "CIT"), 
      labels=c("SAL",  "CIT"),
    )
    ) 
ylab_edge <- "Neuron Pair\nInteraction\nWeight"

#####################

read_parquet("data/derived/corrs/spont - corr.parquet") %>%
  filter(session_name != 'chronic_08') %>%
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
    is_neg = corr < -0.15,
    is_pos = corr > 0,
    mag = abs(corr),
    rec = if_else(corr < 0, 0, corr)
  ) -> df_corr


df_corr %>%
  filter(bin_width==1, shuffle==F) %>%
  filter(nt_comb == 'SR-SR') %>%
  ggplot(aes(x=corr)) +
  geom_histogram() +
  facet_grid(rows=vars(group)) +
  lims(x=c(-0.45, 0.45))


neg_mod <- glm(
  is_neg ~ nt_comb + group + nt_comb:group + distance,
  data=filter(df_corr, bin_width==1, shuffle==F),
  family=binomial()
)

anova(neg_mod)
emms_neg <- emmeans(
  neg_mod, 
  specs = ~ group | nt_comb,
  type="response"
)
pairs(emms_neg)

emms_neg %>%
  as_tibble() %>%
  ggplot(aes(x=group, y=prob, fill=group, ymin=prob - SE, ymax=prob + SE)) + 
  geom_bar(stat="identity",  position=position_dodge(preserve = "single", width=0.4), color="black", width=0.5) +
  geom_errorbar(width=0.28, color='#5c5c5c', position=position_dodge(preserve = "single", width=0.8)) +
  scale_fill_manual(values=c(SAL="grey", CIT="black")) +
  facet_grid(cols=vars(nt_comb)) +
  guides(fill="none") +
  labs(y="Proportion of\nPairs With Positive\nCorrelations")



anova(neg_mod)
emms_neg <- emmeans(
  neg_mod, 
  specs = ~ group | nt_comb,
  type="response"
)
pairs(emms_neg)

emms_neg %>%
  as_tibble() %>%
  ggplot(aes(x=group, y=prob, fill=group, ymin=prob - SE, ymax=prob + SE)) + 
  geom_bar(stat="identity",  position=position_dodge(preserve = "single", width=0.4), color="black", width=0.5) +
  geom_errorbar(width=0.28, color='#5c5c5c', position=position_dodge(preserve = "single", width=0.8)) +
  scale_fill_manual(values=c(SAL="grey", CIT="black")) +
  facet_grid(cols=vars(nt_comb)) +
  guides(fill="none") +
  labs(y="Proportion of\nPairs With Negative\nCorrelations")





mag_mod <- lm(
  mag ~ nt_comb + group + nt_comb:group + distance,
  data=filter(df_corr, bin_width==1, shuffle==F, corr > 0.05)
)
anova(mag_mod)
emms_mag <- emmeans(
  mag_mod,
  specs = ~ group | nt_comb
  )
emms_mag

em_tab <- as_tibble(emms_mag) %>%
  mutate(`Spike Count Correlation` = str_c(round(emmean, 2),  round(SE, 2), sep=" +-")) %>%
  select(group, nt_comb, `Spike Count Correlation`) %>%
  pivot_wider(names_from=c(group), values_from=c(`Spike Count Correlation`))

pairs(emms_mag) %>% as_tibble() %>%
  mutate(`Difference` = str_c(round(estimate, 2),  round(SE, 2), sep=" +-")) %>%
  mutate(p = round(p.value, 4)) %>%
  select(nt_comb, `Difference`, p) %>%
  right_join(em_tab) %>%
  mutate(`Neuron Type Combination`=nt_comb) %>%
  select(`Neuron Type Combination`, SAL, CIT, `Difference`, p)

emms_mag %>%
    as_tibble() %>%
    ggplot(aes(x=group, y=emmean, fill=group, ymin=emmean - SE, ymax=emmean + SE)) +
    geom_bar(stat="identity",  position=position_dodge(preserve = "single", width=0.4), color="black", width=0.5) +
    geom_errorbar(width=0.28, color='#5c5c5c', position=position_dodge(preserve = "single", width=0.8)) +
    scale_fill_manual(values=c(SAL="grey", CIT="black")) +
    facet_grid(cols=vars(nt_comb)) +
    guides(fill="none") +
    labs(y=ylab_edge) +
    theme(
      axis.text.x = element_text(size=10, angle=0),
    )

emms_mag

###############

read_parquet("data/derived/corrs/spont - pcup.parquet") %>%
  filter(group != "DIS") %>%
  # filter(session_name != 'chronic_08') %>%
  left_join(select(df_node, neuron_id, neuron_type)) %>%
  mutate(
      group = factor(
      group, 
      levels=c("SAL", "CIT"), 
      labels=c("SAL",  "CIT"),
    ),
  is_neg = cc < 0,
  mag = abs(cc),
  rec = if_else(cc < 0, 0, cc),
  is_sig = mag > 0.1
  ) -> df_pcup

pcup_mod <- lm(
  mag ~ group * neuron_type,
  data=filter(df_pcup, shuffle==F, bin_width==1)
)
anova(pcup_mod)
emms_pcup <- emmeans(
  pcup_mod, 
  specs= ~ group | neuron_type
)
pairs(emms_pcup)
emms_pcup %>%
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


#############

mod_ensemble <- glm(in_ensemble ~ neuron_type * group,
                    data= filter(df_node, bin_width == 1, shuffle==F),
                    family = binomial()
                    )
emms_ensemble <- emmeans(
  mod_ensemble,
  specs = ~ group | neuron_type,
  type="response"
)
pairs(emms_ensemble)

emms_ensemble %>%
  as_tibble() %>%
  ggplot(aes(x=group, y=prob, fill=group, ymin=prob - SE, ymax=prob + SE)) + 
  geom_bar(stat="identity",  position=position_dodge(preserve = "single", width=0.4), color="black", width=0.5) +
  geom_errorbar(width=0.28, color='#5c5c5c', position=position_dodge(preserve = "single", width=0.8)) +
  scale_fill_manual(values=c(SAL="grey", CIT="black")) +
  facet_grid(cols=vars(neuron_type)) +
  guides(fill="none") +
  labs(y="Proportion of\nPairs With Negative\nCorrelations")




#########

mod_test <- lm(clust ~ neuron_type * group,
                    data= filter(df_node, bin_width == 1, shuffle==F),
)
ems_test <- emmeans(
  mod_test,
  specs = ~ group | neuron_type
)
pairs(ems_test)

ems_test %>%
  as_tibble() %>%
  ggplot(aes(x=group, y=emmean, fill=group, ymin=emmean - SE, ymax=emmean + SE)) + 
  geom_bar(stat="identity",  position=position_dodge(preserve = "single", width=0.4), color="black", width=0.5) +
  geom_errorbar(width=0.28, color='#5c5c5c', position=position_dodge(preserve = "single", width=0.8)) +
  scale_fill_manual(values=c(SAL="grey", CIT="black")) +
  facet_grid(cols=vars(neuron_type)) +
  guides(fill="none") +
  labs(y="Proportion of\nPairs With Negative\nCorrelations")



###########

df_corr %>%
  filter(bin_width==1, shuffle==F) %>%
  group_by(session_name, nt_comb, group) %>%
  summarise(c = mean(corr), num=n()) %>%
  filter(nt_comb == 'SR-SR') %>%
  ggplot(aes(x=c, y=num)) +
  geom_point() +
  facet_grid(rows=vars(group))

df_corr %>%
  filter(bin_width==1, shuffle==F) %>%
  group_by(session_name, nt_comb, group) %>%
  summarise(c = mean(corr), num=n()) %>%
  filter(nt_comb == 'SR-SR') %>%
  filter(c < 0)
