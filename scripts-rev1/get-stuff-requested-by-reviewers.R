# partindo do principio que objetos gerados em 01_fit_model estao no environment

#rm(list = ls())

#library(tidyverse)

res <- readRDS("results-rev1/res.rds")
#jags_data <- readRDS("data-rev1/jags_data.rds")



####################################################
###########===== REVIEWER 1 =======#################
####################################################


#----- assess goodness of fit

# posterior predictive check
# extract posterior samples
fit.obs <- res$sims.list$fit.obs
fit.rep <- res$sims.list$fit.rep
y.rep   <- res$sims.list$y.rep  # matrix of dim: (n_iterations x nobs)

# calculate Bayesian p-value
bayesian_p <- mean(fit.rep > fit.obs)
cat("Bayesian p-value:", bayesian_p, "\n")

# set up plotting layout
par(mfrow = c(1, 2), mar = c(5, 4, 4, 2))

# plot A: Chi-squared Discrepancy Plot (fit.obs vs. fit.rep)
plot(fit.obs, fit.rep, 
     xlab = "Discrepancy (observed data)", 
     ylab = "Discrepancy (replicated data)",
     main = paste("Freeman-Tukey test\nBayesian p-value =", round(bayesian_p, 3)),
     pch = 20, col = rgb(0.2, 0.2, 0.2, 0.4),
     asp = 1)
abline(0, 1, col = "red", lwd = 2, lty = 2) # 1:1 reference line

# plot B: Observed vs. Mean Replicated Counts
mean.y.rep <- colMeans(y.rep)

plot(jags_data$y, mean.y.rep, 
     xlab = "Observed counts (y)", 
     ylab = "Posterior Mean predicted counts",
     main = "Observed vs. Predicted counts",
     pch = 20, col = rgb(0, 0.3, 0.7, 0.5))
abline(0, 1, col = "red", lwd = 2, lty = 2) # 1:1 reference line

# reset plot margins
par(mfrow = c(1, 1))


#-----  specify the number of posterior samples

# NB! to do


# ----- display posterior mean along with 95 and 80% CI


# NB! por enquanto, selecionado para uma UC
# depois, fazer para spp selecionadas

# versao linerange
units_df |>
  mutate(n = apply(res$sims.list$n, 2, mean),
         year = year+2013) |>
  left_join(pop_attr |> select(pa, pop, species)) |>
  filter(pa == 20) |> # cazumba
  group_by(species, year) |>
  summarise(mean = mean(n, na.rm=TRUE),
            lower95 = quantile(n, prob = 0.1, na.rm=TRUE),
            upper95 = quantile(n, prob = 0.9, na.rm=TRUE),
            lower80 = quantile(n, prob = 0.2, na.rm=TRUE),
            upper80 = quantile(n, prob = 0.8, na.rm=TRUE)) |>
  ungroup() |>
  complete(year, nesting(species)) |>
  ggplot(aes(x = year, y = mean)) +
  geom_line(aes(x = year, y = mean), linewidth = 0.3, alpha = 0.8) +
  geom_point(aes(x = year, y = mean),  shape = 21,
             color = "white", stroke = 2,
             fill = "black", size = 4, alpha = 1) +
  geom_linerange(aes(x = year, y = mean,
                     ymin=lower95, ymax=upper95),
                 linewidth = 0.7, alpha = 1) +
  geom_linerange(aes(x = year, y = mean,
                     ymin=lower80, ymax=upper80),
                 linewidth = 3, alpha = 0.5) +
  #coord_cartesian(xlim = c(2014, 2025), expand = FALSE) +
  scale_x_continuous(breaks = scales::breaks_pretty(n=5)) +
  ylim(0,NA) +
  xlab("") +
  ylab("Abundância relativa") +
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
       "results-rev1/cazumba_trends_linerange.jpeg")


# versao ribbon
units_df |>
  mutate(n = apply(res$sims.list$n, 2, mean),
         year = year+2013) |>
  left_join(pop_attr |> select(pa, pop, species)) |>
  filter(pa == 20) |> # cazumba
  group_by(species, year) |>
  summarise(mean = mean(n, na.rm=TRUE),
            lower95 = quantile(n, prob = 0.1, na.rm=TRUE),
            upper95 = quantile(n, prob = 0.9, na.rm=TRUE),
            lower80 = quantile(n, prob = 0.2, na.rm=TRUE),
            upper80 = quantile(n, prob = 0.8, na.rm=TRUE)) |>
  ungroup() |>
  complete(year, nesting(species)) |>
  ggplot(aes(x = year, y = mean)) +
  geom_ribbon(aes(ymin = lower95, ymax = upper95), fill = "steelblue", alpha = 0.5) +
  geom_ribbon(aes(ymin = lower80, ymax = upper80), fill = "steelblue", alpha = 0.7) +
  geom_line(color = "black", linewidth = 0.2, alpha = 0.4) +
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
       "results-rev1/cazumba_trends_ribbon.jpeg")
# common recommendation: widths 6-8 cm (single column) or 16-18 cm (double column)
#maximum height of around 27 cm



#----- please clarify the calculation method for lpi

# escrever... para os codigos, ver 02_compute_lpi


#----- summary of the number of observations per species

data_long$ydata |>
  left_join(data_long$pop_indexes |>
              select(pa,pop,species),
            by = c("pa", "pop")) |>
  group_by(species) |>
  summarise(n = sum(y)) |>
  ungroup() |>
  arrange(desc(n)) |>
  mutate(species = fct_inorder(species)) |>
  #rowid_to_column() |>
  #print() |>
  ggplot(aes(x = species, y = log(n))) +
  geom_point() +
  #geom_line() +
  theme_minimal()


# the rhat and effective sample size must be in supplementary material

# set up plotting layout
par(mfrow = c(1, 2), mar = c(5, 4, 4, 2))


# parameter n
hist(res$Rhat$n,
     main = "R-hat, n",
     xlab = "R-hat", 
     ylab = "Freq")
abline(v = 1.1, lty = "dashed", lwd = 2, col = "red")

hist(res$n.eff$n,
     main = "Effective sample size, n",
     xlab = "n.eff", 
     ylab = "Freq")


# parameter p
hist(res$Rhat$p,
     main = "R-hat, p",
     xlab = "p", 
     ylab = "Freq")
abline(v = 1.1, lty = "dashed", lwd = 2, col = "red")

hist(res$n.eff$p,
     main = "Effective sample size, p",
     xlab = "n.eff", 
     ylab = "Freq")


# parameter r
hist(res$Rhat$r,
     main = "R-hat, r",
     xlab = "R-hat", 
     ylab = "Freq")
abline(v = 1.1, lty = "dashed", lwd = 2, col = "red")

hist(res$n.eff$r,
     main = "Effective sample size, r",
     xlab = "n.eff", 
     ylab = "Freq")


# parameter b_body_mass
hist(res$Rhat$b_body_mass,
     main = "R-hat, b_body_mass",
     xlab = "R-hat", 
     ylab = "Freq")
abline(v = 1.1, lty = "dashed", lwd = 2, col = "red")

hist(res$n.eff$b_body_mass,
     main = "Effective sample size, b_body_mass",
     xlab = "n.eff", 
     ylab = "Freq")


# parameter b_troph
hist(res$Rhat$b_troph,
     main = "R-hat, b_troph",
     xlab = "R-hat", 
     ylab = "Freq")
abline(v = 1.1, lty = "dashed", lwd = 2, col = "red")

hist(res$n.eff$b_threat,
     main = "Effective sample size, b_troph",
     xlab = "n.eff", 
     ylab = "Freq")


# parameter n
hist(res$Rhat$b_threat,
     main = "R-hat, b_threat",
     xlab = "R-hat", 
     ylab = "Freq")
abline(v = 1.1, lty = "dashed", lwd = 2, col = "red")

hist(res$n.eff$b_threat,
     main = "Effective sample size, b_threat",
     xlab = "n.eff", 
     ylab = "Freq")


# parameter b_forest_loss
hist(res$Rhat$b_forest_loss,
     main = "R-hat, b_forest_loss",
     xlab = "R-hat", 
     ylab = "Freq")
abline(v = 1.1, lty = "dashed", lwd = 2, col = "red")

hist(res$n.eff$b_forest_loss,
     main = "Effective sample size, b_forest_loss",
     xlab = "n.eff", 
     ylab = "Freq")


# parameter n
hist(res$Rhat$b_pop_size,
     main = "R-hat, b_pop_size",
     xlab = "R-hat", 
     ylab = "Freq")
abline(v = 1.1, lty = "dashed", lwd = 2, col = "red")

hist(res$n.eff$b_pop_size,
     main = "Effective sample size, b_pop_size",
     xlab = "n.eff", 
     ylab = "Freq")


# parameter b_mgmt
hist(res$Rhat$b_mgmt,
     main = "R-hat, b_mgmt",
     xlab = "R-hat", 
     ylab = "Freq")
abline(v = 1.1, lty = "dashed", lwd = 2, col = "red")

hist(res$n.eff$b_mgmt,
     main = "Effective sample size, b_mgmt",
     xlab = "n.eff", 
     ylab = "Freq")


# reset plot margins
par(mfrow = c(1, 1))


####################################################
###########===== REVIEWER 2 =======#################
####################################################

#----- clarify whether uncertainty was propagated into glmm

# NB! foi incluido no codigo jags, falta reescrever nos metodos


# consider showing the distributions of key continuous covariates

# set up plotting layout
par(mfrow = c(3, 2), mar = c(5, 4, 4, 2))

hist(pa_attr$loss*100,
     main = "PA percentage of forest loss",
     xlab = "Forest Loss (%)", 
     ylab = "Freq") 

hist(pa_attr$pop_total,
     main = "PA human population size",
     xlab = "Population", 
     ylab = "Freq") 

barplot(table(pa_attr$pa_type),
     main = "PA protected area type",
     ylab = "N") 

  
hist(pop_attr$mass,
     main = "Species body mass (g)",
     xlab = "Body mass (g)", 
     ylab = "Freq")

barplot(table(pop_attr$trophic_level),
        main = "Species trophic level",
        ylab = "N") 

barplot(table(pop_attr$categoria),
        main = "Species threat status",
        ylab = "N") 

# reset plot margins
par(mfrow = c(1, 1))



#----- consider complementing histogram in Fig. 2 with
# plot showing beta1 and 95 CIs

pop_attr |>
  mutate(rmean = apply(res$sims.list$r, 2, mean),
         l80 = apply(res$sims.list$r, 2, quantile, prob = 0.2),
         u80 = apply(res$sims.list$r, 2, quantile, prob = 0.8),
         l95 = apply(res$sims.list$r, 2, quantile, prob = 0.025),
         u95 = apply(res$sims.list$r, 2,  quantile, prob = 0.975)) |>
  arrange(desc(rmean)) %>%
  mutate(species = fct_inorder(species)) |>
  rowid_to_column() %>%
  ggplot(aes(x = rowid, y = rmean)) +
  geom_linerange(aes(ymin = l95, ymax = u95),
                 colour = "steelblue",
                 linewidth = 0.7, alpha = 0.8) +
  geom_linerange(aes(ymin = l80, ymax = u80),
                 colour = "steelblue",
                 linewidth = 1.5, alpha = 0.81) +
  geom_hline(yintercept = 0, linewidth = 0.3, alpha = 0.5) +
  #scale_color_manual(values = c("inc" = "blue", "stable" = "darkgrey", "dec" = "red")) +
  labs(x = "Population (rank in descending order of average r)",
       y = "Population trend (r)") +
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
       "results-rev1/trend_rank_plot.jpeg")
# common recommendation: widths 6-8 cm (single column) or 16-18 cm (double column)
#maximum height of around 27 cm



####################################################
###########===== REVIEWER 3 =======#################
####################################################


#----- I suggest presenting 80% CI in the plot showing pop-specific trends

# already done (see above)


#----- consider an approach that carries uncertainty to the second model

# already done (see above)


#----- how was uncertainty propagated into the lpi

# already done (see above). explain in the methods and in figure caption


#----- 