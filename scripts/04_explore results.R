rm(list = ls())

# load libraries
library(here)
library(tidyverse)
library(brms)
library(marginaleffects)
library(gridExtra)
#library(corrplot)


# Load data and results ---------------------------------------------------

monitora <- readRDS(here("data", "monitora_2025.rds"))
data_wide <- readRDS(here("data", "data_wide.rds"))

pas <- read_rds(here("data", "pa_covars.rds")) %>%
  print()

species <- monitora$especies %>%
  print()

rm(monitora)
rm(data_wide)


#lambdapost <- readRDS(here("tss", "tss_lambdapost.rds"))
npost <- readRDS(here("tss", "tss_npost.rds"))
rpost <- readRDS(here("tss", "tss_rpost.rds"))
ppost <- readRDS(here("tss", "tss_ppost.rds"))


# explore detection p -----------------------------------------------------

ppost %>%
  group_by(pa,species,year) %>%
  summarise(p_mean = mean(p)) %>%
  ggplot(aes(x = p_mean)) +
  geom_histogram(bins = 20, 
                 color = "black", fill = "lightblue",
                 size = 0.2, alpha = 0.5) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  labs(x = "Detection prob.", y = "Freq") +
  theme(axis.text = element_text(size = 12)) +
  theme(
    axis.title.x = element_text(size = 14), # Adjust size for X-axis title
    axis.title.y = element_text(size = 14)  # Adjust size for Y-axis title
  )

ggsave(plot = last_plot(),
       scale = 4,
       units = "cm",
       width = 6,
       height = 4,
       dpi = 300,
       here("report", "images", "p_histogram.jpeg"))
# common recommendation: widths 6-8 cm (single column) or 16-18 cm (double column)
#maximum height of around 27 cm

# get detection summaries
p_summary <- ppost %>%
  group_by(pa,species) %>%
  summarise(p_mean = mean(p)) %>%
  ungroup() %>%
  summarise(pmean = round(mean(p_mean),2),
            pmin = round(min(p_mean),2),
            pmax = round(max(p_mean),2)) %>%
  print()

p_mean <- p_summary$pmean
p_min <- p_summary$pmin
p_max <- p_summary$pmax


# explore lambda ----------------------------------------

#lambdapost %>%
#  group_by(pa,species) %>%
#  summarise(mean_lambda = mean(lambda)) %>%
#  ungroup() %>%
#  ggplot(aes(x = mean_lambda)) +
#  geom_histogram(bins = 20,
#                 color = "black", fill = "lightblue",
#                 alpha = 0.6) +
#  theme_bw() +
#  theme(panel.grid.major = element_blank(),
#        panel.grid.minor = element_blank()) +
#  labs(x = "lambda médio", y = "Freq") +
#  xlim(0,100) # NB! tem um valor muito alto!


# explore n and pop trends ------------------------------------------------

npost %>%
  group_by(pa,species) %>%
  summarise(mean_n = mean(n)) %>%
  ungroup() %>%
  ggplot(aes(x = mean_n)) +
  geom_histogram(bins = 30,
                 color = "black", fill = "lightblue",
                 size = 0.2, alpha = 0.5) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  labs(x = "Average abundance", y = "Freq") +
  theme(axis.text = element_text(size = 12)) +
  theme(
    axis.title.x = element_text(size = 14), # Adjust size for X-axis title
    axis.title.y = element_text(size = 14)  # Adjust size for Y-axis title
  )

ggsave(plot = last_plot(),
       scale = 4,
       units = "cm",
       width = 6,
       height = 4,
       dpi = 300,
       here("report", "images", "n_histogram.jpeg"))
# common recommendation: widths 6-8 cm (single column) or 16-18 cm (double column)
#maximum height of around 27 cm

# check suspicious n (too large)
npost %>%
  group_by(pa,species) %>%
  summarise(mean_n = mean(n)) %>%
  ungroup() %>%
  arrange(desc(mean_n))
# suspicious: Cebus olivaceus from Monte Roraima and Guerlinguetus from Ipau-Anilzinho
# remove them in next iteration. I suspect they did not converge (check)

# check Dasyprocta iacki
npost %>%
  filter(pa == 123, species == "Dasyprocta iacki") %>%
  filter(sampled_year == "y") %>%
  group_by(species, year) %>%
  summarise(mean = mean(n, na.rm=TRUE),
            lower = quantile(n, prob = 0.025, na.rm=TRUE),
            upper = quantile(n, prob = 0.975, na.rm=TRUE)) %>%
  ungroup() %>%
  ggplot(aes(x = year, y = mean)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "lightblue", alpha = 0.5) +
  geom_line(color = "black", size = 0.2, alpha = 0.4) +
  geom_point(color = "black", size = 1, alpha = 0.6) +
  ylim(0,NA) +
  xlab("") +
  ylab("Relative abundance") +
  theme_bw() +
  theme(axis.text = element_text(size = 12)) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())


# vamos excluir temporariamente essa populacao da base de dados
# npost <- npost %>%
#   filter(! (pa == 189 & species == "Cebus olivaceus")) %>%
#   print()
# rpost <- rpost %>%
#   filter(! (pa == 189 & species == "Cebus olivaceus")) %>%
#   print()
# ppost <- ppost %>%
#   filter(! (pa == 189 & species == "Cebus olivaceus")) %>%
#   print()


# get relative abundance summaries
n_summary <- npost %>%
  group_by(pa,species) %>%
  summarise(n_mean = mean(n)) %>%
  ungroup() %>%
  summarise(nmean = round(mean(n_mean),2),
            nmin = round(min(n_mean),2),
            nmax = round(max(n_mean),2)) %>%
  print()

n_mean <- n_summary$nmean
n_min <- n_summary$nmin
n_max <- n_summary$nmax

# agrupar tendencias por uc ribbon
npost %>%
  filter(pa == 207) %>%
  filter(sampled_year == "y") %>%
  group_by(species, year) %>%
  summarise(mean = mean(n, na.rm=TRUE),
            lower = quantile(n, prob = 0.1, na.rm=TRUE),
            upper = quantile(n, prob = 0.9, na.rm=TRUE)) %>%
  ungroup() %>%
  ggplot(aes(x = year, y = mean)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "lightblue", alpha = 0.9) +
  geom_line(color = "black", size = 0.2, alpha = 0.4) +
  geom_point(color = "black", size = 1.5, alpha = 0.6) +
  ylim(0,NA) +
  #expand_limits(x = c(2014, 2024)) +
  scale_x_continuous(breaks = scales::breaks_pretty(n=5)) +
  xlab("") +
  ylab("Relative abundance") +
  #theme_classic() +
  theme_bw() +
  theme(axis.text = element_text(size = 12)) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  facet_wrap(~species, scales = "free", ncol = 3) +
  theme(strip.text = element_text(size = 14))

ggsave(plot = last_plot(),
       scale = 4,
       units = "cm",
       width = 8,
       height = 6,
       dpi = 300,
       here("results", "gurupi_trends_ribbon.jpeg"))
# common recommendation: widths 6-8 cm (single column) or 16-18 cm (double column)
#maximum height of around 27 cm


# agrupar por uc barras
npost %>%
  filter(sampled_year == "y") %>%
  filter(pa == 232) %>%
  #filter(species == "Dasyprocta fuliginosa") %>%
  group_by(species, year) %>%
  summarise(mean = mean(n, na.rm=TRUE),
            lower = quantile(n, prob = 0.025, na.rm=TRUE),
            upper = quantile(n, prob = 0.975, na.rm=TRUE)) %>%
  ungroup() %>%
  ggplot(aes(x = year, y = mean)) +
  geom_line(aes(x = year, y = mean), size = 0.2, alpha = 0.5) +
  geom_point(aes(x = year, y = mean),  shape = 21,
             color = "white", stroke = 3,
             fill = "black", size = 1.5) +
  geom_linerange(aes(x = year, y = mean,
                     ymin=lower, ymax=upper),
                 size = 0.2, alpha = 0.3) +
  ylim(0,NA) +
  #expand_limits(x = c(2014, 2024)) +
  scale_x_continuous(breaks = seq(2014, 2024, by = 2)) +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  labs(x = "", y = "Abund. estimada") +
  #theme_classic() +
  theme_bw() +
  theme(axis.text = element_text(size = 12)) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  facet_wrap(~species, scales = "free", ncol = 4) +
  theme(strip.text = element_text(size = 14))

ggsave(plot = last_plot(),
       scale = 4,
       units = "cm",
       width = 8,
       height = 6,
       dpi = 300,
       here("results", "cazumba_trends_linerange.jpeg"))
# common recommendation: widths 6-8 cm (single column) or 16-18 cm (double column)
#maximum height of around 27 cm


# explore beta1 -----------------------------------------------------------

# tabela r
tabela_r <- rpost %>%
  group_by(pa_name,species) %>%
  summarise(r_mean = round(mean(r),2),
            lower = round(quantile(r, prob = 0.025),2),
            upper = round(quantile(r, prob = 0.975),2),
            prob_increase = round(length(r[r>0])/n(),2),
            prob_decrease = round(length(r[r<0])/n(),2)) %>%
  ungroup() %>%
  print()

# add species covars
species <- species %>%
  rename(class = classe, order = ordem, species = especie, redlist = categoria) %>%
  select(class, order, species, redlist, mass, sociality, foraging_stratum, trophic_level) %>%
  print()

# add species taxonomy and traits
tabela_r <- tabela_r %>%
  left_join(species,
            by = "species") %>%
  print()

# add pa covariates
#pas <- read_rds(here("data", "pa_covars.rds")) %>%
#  print()

tabela_r <- tabela_r %>%
  left_join(pas,
            by = "pa_name") %>%
  print()

# add initial/average abundances
tabela_r <- tabela_r %>%
  left_join(npost %>%
              group_by(pa_name, species) %>%
              summarise(mean_n = mean(n)),
            by = c("pa_name", "species")) %>%
  print()

# add time series length
data_wide <- readRDS(here("data", "data_wide.rds"))
data_wide$indexes

tabela_r <- tabela_r %>%
  left_join(data_wide$indexes %>%
              select(pa_name, nyear),
            by = c("pa_name")) %>%
  print()

rm(data_wide)

# some duplications appeared while doing the left_joins above
tabela_r <- tabela_r %>% 
  distinct(pa_name, species, .keep_all = TRUE)

# save for later use
supplementary_mat_1 <- tabela_r %>%
  select(pa_name,
         order, species,
         redlist,
         r_mean, prob_decrease, prob_increase,
         nyear) %>%
  print()
saveRDS(supplementary_mat_1, here("results", "supplementary_mat_1.rds" ))

# quem tem > 95% de estar em declinio  
tabela_r %>%
  filter(prob_decrease > 0.95) %>%
  arrange(desc(prob_decrease)) %>%
  print(n=Inf)

# quais UCs com declinios  
tabela_r %>%
  filter(prob_decrease > 0.95) %>%
  group_by(pa_name) %>%
  count() %>%
  print(n=Inf)

# quais UCs com aumentos  
tabela_r %>%
  filter(prob_increase > 0.95) %>%
  group_by(pa_name) %>%
  count() %>%
  print(n=Inf)

# porcentagem pops de cada UC com declinios  
tabela_r %>%
  mutate(dec = case_when(prob_decrease > 0.95 ~ "decline",
                         .default = "no")) %>%
  group_by(pa_name) %>%
  summarise(n_pop = n(),
            n_dec = sum(dec == "decline"),
            prop_dec = round(sum(dec == "decline")/n(),1)) %>%
  arrange(desc(prop_dec)) %>%
  print(n=Inf)


# summarize r trends and save for quarto report
trends <- tabela_r %>%
  distinct(pa_name, species, .keep_all = TRUE) %>%
  summarise(n_declining = sum(prob_decrease > 0.95),
            percent_declining = round(100*sum(prob_decrease > 0.95)/n(),1),
            n_stable = sum(prob_decrease <= 0.95 & prob_decrease <= 0.95),
            percent_stable = round(100*sum(prob_decrease <= 0.95 & prob_increase <= 0.95)/n(),1),
            n_increasing = sum(prob_increase > 0.95),
            percent_increasing = round(100*sum(prob_increase > 0.95)/n(),1)) %>%
  print()
n_declining <- trends %>% pull(n_declining)
percent_declining <- trends %>% pull(percent_declining)
n_stable <- trends %>% pull(n_stable)
percent_stable <- trends %>% pull(percent_stable)
n_increasing <- trends %>% pull(n_increasing)
percent_increasing <- trends %>% pull(percent_increasing)
n_declining
percent_declining
n_stable
percent_stable
n_increasing
percent_increasing

mean_beta <- round(mean(tabela_r$r_mean),2)
min_beta <- round(min(tabela_r$r_mean),2)
max_beta <- round(max(tabela_r$r_mean),2)
mean_beta
min_beta
max_beta


result_summaries <- list(n_declining = n_declining,
                         percent_declining = percent_declining,
                         n_stable = n_stable,
                         percent_stable = percent_stable,
                         n_increasing = n_increasing,
                         percent_increasing = percent_increasing,
                         mean_beta = mean_beta,
                         min_beta = min_beta,
                         max_beta = max_beta,
                         p_mean = p_mean,
                         p_min = p_min,
                         p_max = p_max,
                         n_mean = n_mean,
                         n_min = n_min,
                         n_max = n_max)
result_summaries

# unify summaries in a single list and save
sample_summaries <- readRDS( here("results", "sample_summaries.rds"))
summaries <- c(sample_summaries, result_summaries)
saveRDS(summaries, here("results", "summaries.rds"))


# plot histogram of betas
tabela_r %>%
  mutate(trend = case_when(prob_decrease > 0.95 ~ "decreasing",
                           prob_increase > 0.95 ~ "increasing",
                           .default = "stable")) %>%
  mutate(trend = factor(trend, levels = c("increasing", "stable", "decreasing"))) %>%
  ggplot() +
  geom_histogram(aes(x = r_mean, fill = trend),
                 #bins = 20,
                 alpha = 0.8) +
  scale_fill_manual(values = c("increasing" = "blue", "stable" = "darkgrey", "decreasing" = "red")) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.1, alpha = 0.5) +
  theme_bw() +
  labs(x = expression(beta[1] ~ "(temporal trend)"),
       y = "Freq") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  #theme(axis.text = element_text(size = 12)) +
  theme(legend.position = c(0.2, 0.7)) +
  #theme(legend.position.inside = c(-0.25, 50)) +
  theme(legend.title = element_blank(),
        legend.text = element_text(size = 14)) +
  theme(
    axis.title.x = element_text(size = 14), # Adjust size for X-axis title
    axis.title.y = element_text(size = 14)  # Adjust size for Y-axis title
  )#+
  #labs(color = "Population trend")

ggsave(plot = last_plot(),
       scale = 4,
       units = "cm",
       width = 5,
       height = 5,
       dpi = 300,
       here("report", "images", "trend_histogram.jpeg"))
# common recommendation: widths 6-8 cm (single column) or 16-18 cm (double column)
#maximum height of around 27 cm

# plot ranking of growth rates
tabela_r %>%
  mutate(trend = case_when(prob_decrease > 0.95 ~ "dec",
                           prob_increase > 0.95 ~ "inc",
                           .default = "stable")) %>%
  #mutate(trend = case_when(upper < 0 ~ "dec",
  #                         lower > 0 ~ "inc",
  #                         .default = "stable")) %>%
  arrange(desc(r_mean)) %>%
  #filter(r_mean > -2.5) %>%
  rowid_to_column() %>%
  ggplot(aes(x = rowid, y = r_mean, color = trend)) +
  geom_hline(yintercept = 0, linewidth = 0.1, alpha = 0.3) +
  geom_vline(xintercept = nrow(tabela_r)/2, linewidth = 0.1, alpha = 0.3) +
  geom_linerange(aes(ymin = lower, ymax = upper),
                  linewidth = 0.4, alpha = 1) +
  scale_color_manual(values = c("inc" = "blue", "stable" = "darkgrey", "dec" = "red")) +
  labs(x = expression("Population (descending order of average" ~ beta[1] ~ ")"),
       #y = expression(beta[1] ~ "(temporal trend)"),
       y = expression(beta[1])) +
  theme_bw() +
  theme(legend.position="none") +
  #ylim(-0.8,0.8) +
  #scale_x_continuous(limits = c(1,NA))
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  theme(axis.text = element_text(size = 12)) +
  theme(
    axis.title.x = element_text(size = 14), # Adjust size for X-axis title
    axis.title.y = element_text(size = 14)  # Adjust size for Y-axis title
  )

ggsave(plot = last_plot(),
       scale = 4,
       units = "cm",
       width = 6,
       height = 4,
       dpi = 300,
       here("report", "images", "trend_rank_plot.jpeg"))
# common recommendation: widths 6-8 cm (single column) or 16-18 cm (double column)
#maximum height of around 27 cm


# declining pops
decliningSpecies <- tabela_r %>%
  filter(prob_decrease > 0.95) %>%
  distinct(pa_name,species) %>%
  mutate(decliningSpecies = paste(pa_name, species)) %>%
  pull(decliningSpecies)

# plot trends of significantly declining pops
npost %>%
  filter(paste(pa_name, species) %in% decliningSpecies) %>%
  filter(sampled_year == "y") %>%
  group_by(pa, species, year) %>%
  summarise(mean = mean(n, na.rm=TRUE),
            lower = quantile(n, prob = 0.025, na.rm=TRUE),
            upper = quantile(n, prob = 0.975, na.rm=TRUE)) %>%
  ungroup() %>%
  mutate(popID = paste(pa, species)) %>%
  ggplot(aes(x = year, y = mean, group = popID)) +
  geom_line(aes(x = year, y = mean), color = "darkblue", size = 0.2, alpha = 0.8) +
  #geom_smooth(aes(x = year, y = mean),
  #            color = "darkblue", se = FALSE, span = 0.8,
  #            size = 0.2, alpha = 0.5) +
  ylim(0,NA) +
  scale_x_continuous(breaks = seq(2014, 2024, by = 2)) +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  labs(x = "", y = "Abundance") +
  #theme_classic() +
  theme_bw() +
  theme(axis.text = element_text(size = 12)) +
  theme(axis.title = element_text(size = 16)) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

ggsave(plot = last_plot(),
       scale = 3,
       units = "cm",
       width = 8,
       height = 6,
       dpi = 300,
       here("report", "images", "declining_pops_trajectory.jpeg"))


# plotar somente pops con tendencia significativa
npost %>%
  filter(sampled_year == "y") %>%
  filter(paste(pa_name,species) %in% (tabela_r %>%
                                        #filter(order != "Primates") %>%
                                        filter(prob_decrease > 0.95 | prob_increase > 0.95) %>%
                                        mutate(decPops =paste(pa_name,species)) %>%
                                        pull(decPops))) %>%
  group_by(pa_name, species, year) %>%
  summarise(mean = mean(n, na.rm=TRUE),
            lower = quantile(n, prob = 0.025, na.rm=TRUE),
            upper = quantile(n, prob = 0.975, na.rm=TRUE)) %>%
  ungroup() %>%
  mutate(pa_name = case_when(#pa_name == "estacao_ecologica_da_terra_do_meio" ~ "Terra do Meio Ecological Station",
                             #pa_name == "estacao_ecologica_de_maraca" ~ "Maracá Ecological Station",
                             pa_name == "parque_nacional_da_amazonia" ~ "Amazônia National Park",
                             #pa_name == "parque_nacional_do_jau" ~ "Jaú National Park",
                             pa_name == "parque_nacional_do_juruena" ~ "Juruena National Park",
                             #pa_name == "parque_nacional_serra_da_cutia" ~ "Serra da Cutia National Park",
                             pa_name == "reserva_biologica_do_gurupi" ~ "Gurupi Biological Reserve",
                             pa_name == "reserva_extrativista_do_alto_tarauaca" ~ "Alto Tarauacá Extractive Reserve",
                             pa_name == "reserva_extrativista_do_cazumba_iracema" ~ "Cazumbá-Iracema Extractive Reserve")) %>%
  #mutate(pa_name = str_replace_all(pa_name, "_", " "),
  #       pa_name = str_to_title(pa_name)) %>%
  #arrange(pa_name, species) %>%
  #mutate(species = factor(species, levels=unique(species[order(pa_name)]), ordered=TRUE)) %>%
  ggplot(aes(x = year, y = mean)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "steelblue", alpha = 0.6) +
  geom_line(color = "black", size = 0.4, alpha = 0.4) +
  geom_point(color = "black", size = 2, alpha = 0.6) +
  ylim(0,NA) +
  scale_x_continuous(breaks = scales::breaks_pretty(n=5)) +
  xlab("") +
  ylab("Relative abundance") +
  theme_bw() +
  theme(axis.title.y = element_text(size = 24)) +
  theme(axis.text = element_text(size = 18)) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  facet_wrap(vars(pa_name,species),
             scales = "free_y",
             axes = "all",
             ncol = 2) +
  #theme(strip.background = element_blank()) +
  theme(panel.spacing = unit(2, "lines")) +
  theme(strip.text = element_text(size = 19),
        strip.background = element_rect(fill = "#FAFAFA"))

ggsave(plot = last_plot(),
       scale = 3.5,
       units = "cm",
       width = 12,
       height = 12,
       dpi = 300,
       here("report", "images", "significant_pop_trends_ribbon.jpeg"))



# explore potential factors affecting beta1 ----------------------------------------

# check forest loss (to be used in paper, no need for saving)
tabela_r %>% distinct(pa_name, loss) %>%
  summarise(min = min(loss),
            max = max(loss))


# plot relation between beta and log body mass
tabela_r %>%
  mutate(logmass = log(mass)) %>%
  ggplot(aes(x = logmass, y = r_mean, color = species)) +
  #geom_point(size = 1.5, alpha = 0.5) +
  geom_count(color = "grey30",) +
  #geom_smooth(method = "lm", se = FALSE,
  #            color = "black", size = 0.2, alpha = 0.7) +
  theme_bw() +
  theme(legend.position="none") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  theme(axis.text = element_text(size = 12)) +
  theme(axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14)) +
  labs(x = "Log Body mass", y = expression(beta))

ggsave(plot = last_plot(),
       scale = 4,
       units = "cm",
       width = 6,
       height = 4,
       dpi = 300,
       here("report", "images", "beta_vs_body_mass.jpeg"))


# plot relation between beta and average abundance
tabela_r %>%
  filter(mean_n < 50) %>% # keep outliers out
  ggplot(aes(x = mean_n, y = r_mean, color = species)) +
  geom_point(size = 2, color = "grey30", alpha = 0.5) +
  #geom_count(color = "grey30", alpha = 0.5) +
  geom_smooth(method = "loess",
              span = 1.5,
              se = FALSE,
              color = "black", size = 0.2, alpha = 0.7) +
  theme_bw() +
  theme(legend.position="none") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  theme(axis.text = element_text(size = 12)) +
  theme(axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14, angle = 0, vjust = 0.5)) +
  labs(x = "Average abundance", y = expression(beta[1]))

ggsave(plot = last_plot(),
       scale = 4,
       units = "cm",
       width = 6,
       height = 4,
       dpi = 300,
       here("report", "images", "beta_vs_average_abundance.jpeg"))


# plot relation between beta and monitoring duration
tabela_r %>%
  #filter(mean_n < 50) %>% # keep outliers out
  ggplot(aes(x = nyear, y = r_mean, color = species)) +
  geom_jitter(size = 2, width = 0.2, color = "grey30", alpha = 0.5) +
  theme_bw() +
  theme(legend.position="none") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  theme(axis.text = element_text(size = 12)) +
  theme(axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14, angle = 0, vjust = 0.5)) +
  labs(x = "Time series length (years)", y = expression(beta[1]))

ggsave(plot = last_plot(),
       scale = 4,
       units = "cm",
       width = 6,
       height = 4,
       dpi = 300,
       here("report", "images", "beta_vs_average_abundance.jpeg"))


# plot relation between beta and population inside PA
tabela_r %>%
  #ggplot(aes(x = pop_total, y = r_mean, color = pa)) +
  ggplot(aes(x = log(pop_total), y = r_mean, color = pa)) +
  #geom_point(size = 1, alpha = 0.5) +
  geom_count(color = "grey30") +
  #geom_smooth(method = "lm", se = FALSE,
  #            color = "black", size = 0.2, alpha = 1) +
  theme_bw() +
  theme(legend.position="none") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  theme(axis.text = element_text(size = 12)) +
  theme(axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14)) +
  labs(x = "Log population", y = expression(beta))


ggsave(plot = last_plot(),
       scale = 4,
       units = "cm",
       width = 6,
       height = 4,
       dpi = 300,
       here("report", "images", "beta_vs_population.jpeg"))



# plot relation between beta and forest loss
tabela_r %>%
  ggplot(aes(x = loss, y = r_mean, color = pa)) +
  #geom_point(size = 1, alpha = 0.5) +
  geom_count(color = "grey30") +
  #geom_smooth(method = "lm", se = FALSE,
  #            color = "black", size = 0.2, alpha = 1) +
  theme_bw() +
  theme(legend.position="none") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  theme(axis.text = element_text(size = 12)) +
  theme(axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14)) +
  labs(x = "Forest loss (%)", y = expression(beta))

ggsave(plot = last_plot(),
       scale = 4,
       units = "cm",
       width = 6,
       height = 4,
       dpi = 300,
       here("report", "images", "beta_vs_forest_loss.jpeg"))


# beta and conservation status boxplot
tabela_r %>%
  mutate(redlist = case_when(is.na(redlist) ~ "non-threatened",
                             .default = "threatened")) %>%
  ggplot(aes(x = redlist, y = r_mean)) +
  #geom_violin() +
  geom_boxplot(aes(fill = redlist)) +
  #geom_jitter(width = 0.02, size = 1, alpha = 0.4) +
  geom_hline(yintercept = mean(tabela_r$r_mean), linetype = "dashed", alpha = 0.5) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  labs(x = "Grau de ameaça", y = "Tendência") +
  guides(fill = guide_legend(title = "Threat category"))


# beta and conservation status histogram splitted
background_data <- tabela_r %>% select(-redlist)
tabela_r %>%
  mutate(redlist = case_when(is.na(redlist) ~ "non-threatened",
                             .default = "threatened")) %>%
  ggplot(aes(x = r_mean, fill = redlist)) +
  geom_histogram(data = background_data, aes(x = r_mean),
                 fill = "grey", alpha = 0.5) +
  geom_histogram(colour = "black", size = 0.2) +
  guides(fill = FALSE) +
  theme_bw() +
  facet_wrap(~ redlist) +
  theme(strip.text = element_text(size = 14)) +
  theme(axis.text = element_text(size = 12)) +
  theme(axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14)) +
  labs(x = expression(beta[t] ~ "(Temporal trend)"), y = "Count")

ggsave(plot = last_plot(),
       scale = 4,
       units = "cm",
       width = 6,
       height = 3,
       dpi = 300,
       here("report", "images", "beta_hist_redlist.jpeg"))

# beta and taxonomic order boxplot
tabela_r %>%
  ggplot(aes(x = order, y = r_mean)) +
  #geom_violin() +
  geom_boxplot(aes(fill = order)) +
  #geom_jitter(width = 0.02, size = 1, alpha = 0.4) +
  geom_hline(yintercept = mean(tabela_r$r_mean), linetype = "dashed", alpha = 0.5) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  labs(x = "Ordem", y = "Tendência") +
  guides(fill = guide_legend(title = "Order"))

# beta and taxonomic order histogram splitted
background_data <- tabela_r %>% select(-order)
tabela_r %>%
  ggplot(aes(x = r_mean, fill = order)) +
  geom_histogram(data = background_data, aes(x = r_mean),
                 fill = "grey", alpha = 0.5) +
  geom_histogram(colour = "black", size = 0.2) +
  guides(fill = FALSE) +
  theme_bw() +
  facet_wrap(~ order) +
  theme(strip.text = element_text(size = 14)) +
  theme(axis.text = element_text(size = 12)) +
  theme(axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14)) +
  labs(x = expression(beta[t] ~ "(Temporal trend)"), y = "Count")

ggsave(plot = last_plot(),
       scale = 4,
       units = "cm",
       width = 6,
       height = 6,
       dpi = 300,
       here("report", "images", "beta_hist_taxonomic_order.jpeg"))


# beta and PA type boxplot
tabela_r %>%
  ggplot(aes(x = pa_type, y = r_mean)) +
  #geom_violin() +
  geom_boxplot(aes(fill = pa_type)) +
  geom_hline(yintercept = mean(tabela_r$r_mean), linetype = "dashed", alpha = 0.5) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  labs(x = "Categoria de UC", y = "Tendência") +
  guides(fill = guide_legend(title = "PA type"))

# check confidence intervals
tabela_r %>%
  group_by(pa_type) %>%
  summarise(mean = mean(r_mean),
            lower = quantile(r_mean, prob = 0.025),
            upper = quantile(r_mean, prob = 0.975))

# beta and PA type histogram splitted
background_data <- tabela_r %>% select(-pa_type)
tabela_r %>%
  mutate(pa_type = case_when(pa_type == "strictly_protected" ~ "Strictly protected",
                             .default = "Sustainable use")) %>%
  ggplot(aes(x = r_mean, fill = pa_type)) +
  geom_histogram(data = background_data, aes(x = r_mean),
                 fill = "grey", alpha = 0.5) +
  geom_histogram(colour = "black", size = 0.2) +
  guides(fill = FALSE) +
  theme_bw() +
  facet_wrap(~ pa_type) +
  theme(strip.text = element_text(size = 14)) +
  theme(axis.text = element_text(size = 12)) +
  theme(axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14)) +
  labs(x = expression(beta[t] ~ "(Temporal trend)"), y = "Count")

ggsave(plot = last_plot(),
       scale = 4,
       units = "cm",
       width = 6,
       height = 3,
       dpi = 300,
       here("report", "images", "beta_hist_pa_type.jpeg"))


# summarise predictors ----------------------------------------------------

# summarise body mass, trophic_level, forest_loss and human_pop

# body mass
tabela_r %>%
  distinct(species, mass) %>%
  summarise(min_mass = min(mass),
            max_mass = max(mass),
            mean_mass = mean(mass))

tabela_r %>% distinct(species, mass) %>% arrange(mass)
tabela_r %>% distinct(species, mass) %>%
  ggplot() +
  #geom_histogram(aes(x = log(mass)), bins = 10) +
  geom_histogram(aes(x = mass), bins = 15)

# trophic level
tabela_r %>% 
  distinct(species, trophic_level) %>%
  group_by(trophic_level) %>% 
  count()

# forest loss
tabela_r %>%
  distinct(pa_name, loss) %>%
  mutate(loss = loss*100) %>%
  summarise(min_loss = min(loss),
            max_loss = max(loss),
            mean_loss = mean(loss))

# human population
tabela_r %>%
  distinct(pa_name, pop_total) %>%
  summarise(min_pop = min(pop_total),
            max_pop = max(pop_total),
            mean_pop = mean(pop_total))


# glmm beta1 --------------------------------------------------------------

# padronizar variaveis para GLMM
glmm_data <-  tabela_r %>%
  mutate(mass = as.numeric(scale(log(mass))),
         species = as.factor(dense_rank(species)),
         sociality = as.factor(sociality),
         vertical = as.factor(foraging_stratum),
         trophic = factor(trophic_level, levels = c("herbivore", "omnivore", "carnivore")),
         pa = as.factor(dense_rank(pa_name)),
         #pa_type = as.factor(dense_rank(pa_type)),
         pa_type = case_when(pa_type == "strictly_protected" ~ "strictly protected",
                        .default = "sustainable use"),
         pa_type = factor(pa_type, levels = c("strictly protected", "sustainable use")),
         sp = case_when(pa_type == "strictly_protected" ~ 1,
                        .default = 0),
         su = case_when(pa_type == "strictly_protected" ~ 0,
                        .default = 1),
         forest_loss = as.numeric(scale(loss)),
         #percent_loss_50km = as.numeric(scale(percent_loss_50km)),
         redlist = case_when(is.na(redlist) ~ "Not threatened",
                        .default = "Threatened"),
         pop_total = as.numeric(scale(log(pop_total))),
         #region = as.factor(region_ter_steege)
         ) %>%
  mutate(trophic = case_when(trophic %in% c("omnivore", "carnivore") ~ "secondary",
                             .default = "primary")) %>%
  mutate(trophic = factor(trophic, levels = c("primary", "secondary"))) %>%
  select(r_mean, pa, species, 
         order, mass, sociality, vertical, trophic, redlist,
         pa_type, forest_loss, pop_total) %>%
  print()

# checar correlacao entre preditores continuos
with(glmm_data, cor(forest_loss, pop_total))
with(tabela_r, cor(loss, pop_total))
glmm_data %>%
  distinct(pa, forest_loss, pop_total) %>%
  select(-pa) %>%
  cor()


# sintese massa, nivel trofico, PA loss, PA pop
tabela_r %>%
  distinct(species, mass) %>%
  summarize(mean_mass = mean(mass),
            min_mass = min(mass),
            max_mass = max(mass))

tabela_r %>%
  distinct(species, trophic_level) %>%
  group_by(trophic_level) %>%
  count()

tabela_r %>%
  distinct(pa_name, loss) %>%
  summarize(mean_loss = mean(loss),
            min_loss = min(loss),
            max_loss = max(loss))

tabela_r %>%
  distinct(pa_name, pop_total) %>%
  summarize(mean_pop = mean(pop_total),
            min_pop = min(pop_total),
            max_pop = max(pop_total))

  
# checar tutorial:
# https://ourcodingclub.github.io/tutorials/mixed-models/
#library(lme4)
#library(glmmTMB)
#mod1 <- glmmTMB(r_mean  ~ log(mass) + percent_loss + log(pop_total) + pa_type + (1|species) + (1|pa_name),
#                data = tabela_r)
#summary(mod1)
#plot(mod1)
#qqnorm(resid(mod1))
#qqline(resid(mod1)) 

# versao bayesiana
mod1 <- brm(r_mean  ~ #pa_type +
              forest_loss + pop_total +
              mass + trophic +#sociality + vertical +
              (1|species) + (1|pa),
            family = gaussian,
            chains = 3, iter = 5000, warmup = 2500, thin = 10,
            data = glmm_data)

summary(mod1, prob = c(0.95))
round(summary(mod1, prob = c(0.95))$fixed, 3)

# save results table for the paper
glmm_table <- round(summary(mod1, prob = c(0.95))$fixed, 3) %>%
  rownames_to_column(var = "Parameter") %>%
  rename(Mean = Estimate,
         "Lower CI" = "l-95% CI",
         "Upper CI" = "u-95% CI") %>%
  select(Parameter, Mean, "Lower CI", "Upper CI", Rhat) %>%
  tibble() %>%
  print()

saveRDS(glmm_table, here("results", "glmm_table.rds"))


#plot(mod1)
#ranef(mod1)
#coef(mod1)
conditional_effects(mod1)

# plot predicted responses for a given predictor
p1 <- plot_predictions(mod1, condition = "mass") +
  labs(x = "Body mass (standardized)", y = expression(beta[1])) +
  theme_bw() +
  theme(strip.text = element_text(size = 14)) +
  theme(axis.text = element_text(size = 14)) +
  theme(axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14, angle = 0, vjust = 0.5)) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

p2 <- plot_predictions(mod1, condition = "trophic") +
  labs(x = "Trophic level", y = expression(beta[1])) +
  theme_bw() +
  theme(strip.text = element_text(size = 14)) +
  theme(axis.text = element_text(size = 14)) +
  theme(axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14, angle = 0, vjust = 0.5)) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

p3 <- plot_predictions(mod1, condition = "forest_loss") +
  labs(x = "Forest loss (standardized)", y = expression(beta[1])) +
  theme_bw() +
  theme(strip.text = element_text(size = 14)) +
  theme(axis.text = element_text(size = 14)) +
  theme(axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14, angle = 0, vjust = 0.5)) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

p4 <- plot_predictions(mod1, condition = "pop_total") +
  labs(x = "Human population (standardized)", y = expression(beta[1])) +
  theme_bw() +
  theme(strip.text = element_text(size = 14)) +
  theme(axis.text = element_text(size = 12)) +
  theme(axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14, angle = 0, vjust = 0.5)) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

grid.arrange(p1, p2, p3, p4, ncol = 2) 

# save plot
jpeg(here("report", "images", "predictions_plot.jpeg"),
     units = "cm",
     res = 120,
     width = 24, height = 20)
grid.arrange(p1, p2, p3, p4, ncol = 2) 
dev.off() # Close the device and save the file


# miscellaneous stuff -----------------------------------------------------

# plotar somente pops em aumento
npost %>%
  filter(sampled_year == "y") %>%
  filter(paste(pa_name,species) %in% (tabela_r %>%
                                        #filter(order != "Primates") %>%
                                        filter(prob_increase > 0.95) %>%
                                        mutate(decPops =paste(pa_name,species)) %>%
                                        pull(decPops))) %>%
  group_by(pa_name, species, year) %>%
  summarise(mean = mean(n, na.rm=TRUE),
            lower = quantile(n, prob = 0.025, na.rm=TRUE),
            upper = quantile(n, prob = 0.975, na.rm=TRUE)) %>%
  ungroup() %>%
  ggplot(aes(x = year, y = mean)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "lightblue", alpha = 0.5) +
  geom_line(color = "black", size = 0.2, alpha = 0.4) +
  geom_point(color = "black", size = 1, alpha = 0.6) +
  ylim(0,NA) +
  #expand_limits(x = c(2014, 2024)) +
  #scale_x_continuous(breaks = scales::breaks_pretty(n=5)) +
  xlab("") +
  ylab("Relative abundance") +
  #theme_classic() +
  theme_bw() +
  theme(axis.text = element_text(size = 12)) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  facet_wrap(vars(species,pa_name),
             scales = "free", ncol = 2) +
  theme(strip.text = element_text(size = 14))

ggsave(plot = last_plot(),
       scale = 4,
       units = "cm",
       width = 5,
       height = 5,
       dpi = 300,
       here("results", "increasing_pop_trends_ribbon.jpeg"))


#----- explore average growth rate global and per PA

# tabela r
tabela_r_pa <- rpost %>%
  group_by(pa_name) %>%
  summarise(r_mean = round(mean(r),2),
            lower = round(quantile(r, prob = 0.025),2),
            upper = round(quantile(r, prob = 0.975),2),
            prob_increase = round(length(r[r>0])/n(),2),
            prob_decrease = round(length(r[r<0])/n(),2)) %>%
  ungroup() %>%
  print()

# plot ranking of growth rates
tabela_r_pa %>%
  mutate(trend = case_when(prob_decrease > 0.95 ~ "dec",
                           prob_increase > 0.95 ~ "inc",
                           .default = "stable")) %>%
  arrange(desc(r_mean)) %>%
  rowid_to_column() %>%
  ggplot(aes(x = rowid, y = r_mean, color = trend)) +
  geom_hline(yintercept = 0, linewidth = 0.1, alpha = 0.3) +
  #geom_vline(xintercept = nrow(tabela_r)/2, linewidth = 0.1, alpha = 0.3) +
  geom_point() +
  geom_linerange(aes(ymin = lower, ymax = upper),
                 linewidth = 0.4, alpha = 1) +
  scale_color_manual(values = c("inc" = "blue", "stable" = "darkgrey", "dec" = "red")) +
  labs(x = expression("Protected areas (descending order of average" ~ beta[1] ~ ")"),
       #y = expression(beta[1] ~ "(temporal trend)"),
       y = expression(beta[1])) +
  theme_bw() +
  theme(legend.position="none") +
  #ylim(-0.8,0.8) +
  #scale_x_continuous(limits = c(1,NA))
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  theme(axis.text = element_text(size = 12)) +
  theme(
    axis.title.x = element_text(size = 14), # Adjust size for X-axis title
    axis.title.y = element_text(size = 14)  # Adjust size for Y-axis title
  )

ggsave(plot = last_plot(),
       scale = 4,
       units = "cm",
       width = 6,
       height = 4,
       dpi = 300,
       here("report", "images", "trend_rank_plot.jpeg"))
# common recommendation: widths 6-8 cm (single column) or 16-18 cm (double column)
#maximum height of around 27 cm

# sample code to make figure comparing categories
ggplot(your_data, aes(x = Category, y = Value, fill = Period)) +
  geom_col(position = position_dodge(width = 0.9)) + # For bars, dodged by period
  geom_errorbar(aes(ymin = Lower_Error, ymax = Upper_Error),
                width = 0.2, position = position_dodge(width = 0.9)) + # For error bars
  facet_wrap(~ Period, scales = "free_x") + # Create separate panels for each period
  labs(y = "Increase in Lambda") + # Y-axis label
  theme_minimal() + # Or other theme for aesthetics
  theme(axis.title.x = element_blank()) # Remove x-axis title if desired  
  
