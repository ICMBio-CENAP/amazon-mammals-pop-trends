# Programa Monitora
# Analise de dados
# Modelo multi-regiao, multi-especies, N-mixture
# tss (time for space substitution)
# Elildo Carvalho Jr @ICMBio/CENAP 21 out 2026

# NB! o modelo esta salvo como arquivo txt na pasta 02_scripts
# caso seja preciso editar o modelo, fazer diretamente no arquivo txt

rm(list = ls())

#----- carregar pacotes
library(tidyverse)
library(jagsUI)
#library(abind)

# carregar dados previamente preparados
data_long <- readRDS("01_entrada/data_long.rds")


# 1. Definir combinações únicas de pa-pop
pa_pop_combos <- data_long$ydata |>
  distinct(pa, pop) |>
  mutate(pa_pop_id = row_number())

# 2. Definir unidades amostrais únicas
units_df <- data_long$ydata |>
  distinct(pa, pop, site, year) |>
  mutate(unit_id = row_number()) |>
  left_join(pa_pop_combos, by = c("pa", "pop")) |>
  mutate(group_id = as.integer(interaction(pa, pop, year, drop = TRUE)))

# 3. Definir combinações únicas de pa-pop-year
pa_pop_year_df <- units_df |>
  distinct(pa, pop, year, pa_pop_id) |>
  arrange(pa, pop, year) |>
  mutate(pa_pop_year_id = row_number())

# 4. Vincular pa_pop_year_id de volta em units_df
units_df <- units_df |>
  left_join(pa_pop_year_df |> select(pa, pop, year, pa_pop_year_id), by = c("pa", "pop", "year"))

# 5. Definir sites únicos e mapear a seu respectivo combo pa-pop
sites_df <- data_long$ydata |>
  distinct(pa, pop, site) |>
  mutate(site_id = row_number()) |>
  left_join(pa_pop_combos, by = c("pa", "pop"))

# 6. Mapear observações originais às unidades, sites e pa_pop_year
data_mapped <- data_long$ydata |>
  left_join(units_df |> select(pa, pop, site, year, unit_id, pa_pop_id, pa_pop_year_id), by = c("pa", "pop", "site", "year")) |>
  left_join(sites_df |> select(pa, pop, site, site_id), by = c("pa", "pop", "site"))

# 7. Obter ID numérico da PA
pa_lookup <- data_long$ydata |>
  distinct(pa) |>
  arrange(pa) |>
  mutate(pa_id = row_number())

data_mapped <- data_mapped |>
  left_join(pa_lookup, by = "pa")

# 8. Construir a lista completa para o JAGS
jags_data <- list(
  y = as.integer(data_mapped$y),
  nobs = nrow(data_mapped),
  
  # Dimensões
  npa = n_distinct(data_mapped$pa),
  n_pa_pop_combos = nrow(pa_pop_combos),
  n_total_sites = nrow(sites_df),
  n_combo_years = nrow(pa_pop_year_df), 
  n_units = nrow(units_df),
  n_groups = max(units_df$group_id),  
  
  # Mapear vetores para parâmetros e efeitos aleatórios
  pa_map = as.integer(pa_pop_combos |> left_join(pa_lookup, by = "pa") |> pull(pa_id)), 
  combo_year_map = as.integer(pa_pop_year_df$pa_pop_id), 
  combo_map = as.integer(sites_df$pa_pop_id),
  combo_unit = as.integer(units_df$pa_pop_id),
  site_unit = as.integer(units_df |> left_join(sites_df, by = c("pa", "pop", "site", "pa_pop_id")) |> pull(site_id)),
  year_unit = as.integer(units_df$year - min(units_df$year) + 1),
  group_id = as.integer(units_df$group_id),     
  
  # Mapear vetores no nível das observações
  combo_obs = as.integer(data_mapped$pa_pop_id),
  site_obs = as.integer(data_mapped$site_id),
  combo_year_obs = as.integer(data_mapped$pa_pop_year_id), 
  unit_obs = as.integer(data_mapped$unit_id),
  pa_obs = as.integer(data_mapped$pa_id)
)

str(jags_data)

#----- parametros a monitorar
params <- c(
  # essenciais
  'n',
  'p', 'r',
  #'n_sum',
  # goodness-of-fit
  'fit.obs', 'fit.rep', 'bayesian.p', 'y.rep'
  # nao essenciais (desligar se precisar de vlocidade)
  # 'lambda',
  # 'a0', 'b0',
  # 'mu.a0.global', 'mu.a0.pa',
  # 'mu.r.global', 'mu.r.pa',
  # 'mu.b0.global', 'mu.b0.pa'
)


#----- valores iniciais

# inicializar z com todos os valores == 1
z_init <- rep(1, nrow(sites_df))

# encontrar as maiores contagens
unit_max_y <- data_mapped |>
  group_by(unit_id) |>
  summarise(max_y = max(y, na.rm = TRUE)) |>
  ungroup() |>
  print()

# fazer join com  units_df para garantir que a ordem "matches" com 1:n_units,
# preencher observacoes faltantes com um valor inicial seguro (e.g., 1 or max_y + 1)
unit_inits <- units_df |>
  left_join(unit_max_y, by = "unit_id") |>
  mutate(init_n = pmax(max_y + 1, 1, na.rm = TRUE)) |> # N >= max observado
  print()

# definir a funcao para atualizar inits
inits <- function() {
  list(
    z = z_init,
    n = unit_inits$init_n,
    mu.a0.global = rnorm(1, 0, 0.5),
    mu.b0.global = rnorm(1, 0, 0.5)
  )
}

inits()


# #----- parametros mcmc
# ni <- 150000
# nt <- 800
# nb <- 50000
# nc <- 3

# para testes
ni <- 300
nt <- 1
nb <- 50
nc <- 3
 

# # para testes
# ni <- 500
# nt <- 1
# nb <- 100
# nc <- 3


# para testes
# ni <- 50000
# nt <- 60
# nb <- 20000
# nc <- 3


# check predicted size of posterior
((ni-nb)*3)/nt

# remover objetos para liberar espaco
rm(list = c("mod", "pa_pop_combos", "pa_year_df",
            "units_df", "sites_df", "data_mapped")) 

# chamar o jags
system.time(mod <- jagsUI::jags(data = jags_data,
                                inits = inits,
                                parameters.to.save = params,
                                #model.file = "02_scripts/model-olre-hier.txt",
                                #model.file = "02_scripts/model-olre-pa-level-random-effects.txt",
                                model.file = "02_scripts/model-20261118.txt",
                                n.chains = nc,
                                n.thin = nt,
                                n.iter = ni,
                                n.burnin = nb,
                                parallel = TRUE))

# salvar
#saveRDS(mod, "03_saida/res_olre.rds")

# elapsed/3600
123371.863/3600
# model fitting took 29.6 hours
# using fedora linux in asus vivobook intel core i5


#---- GOF -----

#mod <- readRDS("03_saida/res_olre.rds")

# posterior predictive check
# 1. Extract posterior samples
fit.obs <- mod$sims.list$fit.obs
fit.rep <- mod$sims.list$fit.rep
y.rep   <- mod$sims.list$y.rep  # matrix of dim: (n_iterations x nobs)

# 2. Calculate Bayesian p-value
bayesian_p <- mean(fit.rep > fit.obs)
cat("Bayesian p-value:", bayesian_p, "\n")

# Set up a 1-row, 2-column plotting layout
par(mfrow = c(1, 2), mar = c(5, 4, 4, 2))

# 3. Plot A: Chi-squared Discrepancy Plot (fit.obs vs. fit.rep)
plot(fit.obs, fit.rep, 
     xlab = "Discrepancy (observed data)", 
     ylab = "Discrepancy (replicated data)",
     main = paste("Freeman-Tukey test\nBayesian p-value =", round(bayesian_p, 3)),
     pch = 20, col = rgb(0.2, 0.2, 0.2, 0.4),
     asp = 1)
abline(0, 1, col = "red", lwd = 2, lty = 2) # 1:1 reference line

# 4. Plot B: Observed vs. Mean Replicated Counts
mean.y.rep <- colMeans(y.rep)

plot(jags_data$y, mean.y.rep, 
     xlab = "Observed counts (y)", 
     ylab = "Posterior Mean predicted counts",
     main = "Observed vs. Predicted counts",
     pch = 20, col = rgb(0, 0.3, 0.7, 0.5))
abline(0, 1, col = "red", lwd = 2, lty = 2) # 1:1 reference line

# reset plot margins
par(mfrow = c(1, 1))



hist(mod$Rhat$n)
abline(v = 1.1, lty = "dashed", lwd = 2, col = "red")

round(length(which(mod$Rhat$n > 1.1))/length(mod$Rhat$n), 2)

hist(mod$n.eff$n)
hist(mod$n.eff$p)

#traceplot(mod, parameter = "n")
#traceplot(mod, parameter = "n_sum")
#traceplot(mod, parameter = "p")

# check n
str(mod$sims.list$n)
hist(mod$sims.list$n)
#hist(mod$sims.list$n, xlim = c(0,100))


# check p
str(mod$sims.list$p)
hist(mod$sims.list$p)

# check r
str(mod$sims.list$r)
hist(mod$sims.list$r)


npost <- units_df |>
  mutate(n = apply(mod$sims.list$n, 2, mean)) |>
  group_by(pa,pop,year) |>
  summarise(mean_n = mean(n)) |>
  print()
hist(npost$mean_n)

units_df |>
  mutate(n = apply(mod$sims.list$n, 2, mean)) |>
  filter(pa == 27) |> # cazumba
  group_by(pop, year) |>
  summarise(mean = mean(n, na.rm=TRUE),
            lower = quantile(n, prob = 0.1, na.rm=TRUE),
            upper = quantile(n, prob = 0.9, na.rm=TRUE)) |>
  ungroup() |>
  complete(year, nesting(pop)) |>
  ggplot(aes(x = year, y = mean)) +
  geom_line(aes(x = year, y = mean), linewidth = 0.3, alpha = 0.8) +
  geom_point(aes(x = year, y = mean),  shape = 21,
             color = "white", stroke = 2,
             fill = "black", size = 2) +
  geom_linerange(aes(x = year, y = mean,
                     ymin=lower, ymax=upper),
                 linewidth = 0.2, alpha = 0.5) +
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
  facet_wrap(~pop, scales = "free", ncol = 3) +
  theme(strip.text = element_text(size = 14))


#----- alternative way of getting pop-year 
# 1. Extract posterior samples of N (dimensions: iterations x n_units)
n_samples <- mod$sims.list$n  # or extract from your MCMC object

# 2. Convert to a tidy format or compute column-wise across iterations
# Assuming columns of n_samples correspond 1-to-1 with row indices of units_df:
units_with_posterior <- units_df |>
  mutate(unit_index = row_number())

# 3. Compute the mean abundance per pa-pop-year for EACH MCMC iteration
# This preserves the full posterior distribution (including uncertainty for missing sites)
posterior_means <- map_dfr(1:nrow(n_samples), function(iter) {
  units_with_posterior |>
    mutate(n_val = n_samples[iter, unit_index]) |>
    group_by(pa, pop, year) |>
    summarise(mean_n = mean(n_val), .groups = "drop") |>
    mutate(iteration = iter)
})

# 4. Summarize into posterior mean and 95% credible intervals
pa_pop_year_summary <- posterior_means |>
  group_by(pa, pop, year) |>
  summarise(
    mean_abundance = mean(mean_n),
    lower_ci = quantile(mean_n, 0.025),
    upper_ci = quantile(mean_n, 0.975),
    .groups = "drop"
  )

print(pa_pop_year_summary)


pa_pop_year_summary |>
  filter(pa == 36) |> # cazumba
  ggplot(aes(x = year, y = mean_abundance)) +
  geom_line(aes(x = year, y = mean_abundance), linewidth = 0.3, alpha = 0.8) +
  geom_point(aes(x = year, y = mean_abundance),  shape = 21,
             color = "white", stroke = 2,
             fill = "black", size = 2) +
  geom_linerange(aes(x = year, y = mean_abundance,
                     ymin=lower_ci, ymax=upper_ci),
                 size = 0.2, alpha = 0.5) +
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
  facet_wrap(~pop, scales = "free", ncol = 3) +
  theme(strip.text = element_text(size = 14))
