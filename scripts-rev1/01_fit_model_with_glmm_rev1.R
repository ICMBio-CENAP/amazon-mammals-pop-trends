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
data_long <- readRDS("data-rev1/data_long.rds")

#str(data_long)

# 1. Definir combinações únicas de pa-pop e anexar atributos de PA e espécie
pa_attr <- data_long$pa_indexes |>
  select(pa, pa_type, pop_total, loss) |>
  distinct() |>
  print()

pop_attr <- data_long$pop_indexes |>
  select(pa, pop, pop_id, species, categoria, mass, trophic_level) |>
  distinct() |>
  print()

pa_pop_combos <- data_long$ydata |>
  distinct(pa, pop) |>
  mutate(pa_pop_id = row_number()) |>
  left_join(pa_attr, by = "pa") |>
  left_join(pop_attr, by = c("pa", "pop")) |>
  mutate(
    # Variáveis contínuas (padronizadas após transformação log quando apropriado)
    forest_loss_z = as.numeric(scale(loss)),
    pop_size_z    = as.numeric(scale(log1p(pop_total))),
    body_mass_z   = as.numeric(scale(log(mass))),
    
    # Variáveis categóricas e espécie (sem NAs para índices inteiros no JAGS)
    mgmt_f    = factor(replace_na(as.character(pa_type), "unknown")),
    troph_f   = factor(replace_na(as.character(trophic_level), "unknown")),
    threat_f  = factor(replace_na(as.character(categoria), "unlisted")),
    species_f = factor(species),
    
    mgmt_idx   = as.integer(mgmt_f),
    troph_idx  = as.integer(troph_f),
    threat_idx = as.integer(threat_f),
    sp_idx     = as.integer(species_f)
  ) |>
  print()

# Verificação de integridade das covariáveis de r[j]
stopifnot(
  nrow(pa_pop_combos) == n_distinct(data_long$ydata$pa, data_long$ydata$pop),
  !anyNA(pa_pop_combos[, c("forest_loss_z", "pop_size_z", "body_mass_z",
                           "mgmt_idx", "troph_idx", "threat_idx", "sp_idx")])
)

# 2. Definir unidades amostrais únicas
units_df <- data_long$ydata |>
  distinct(pa, pop, site, year) |>
  mutate(unit_id = row_number()) |>
  left_join(pa_pop_combos |> select(pa, pop, pa_pop_id), by = c("pa", "pop")) |>
  mutate(group_id = as.integer(interaction(pa, pop, year, drop = TRUE))) |>
  print()

# 3. Definir combinações únicas de pa-pop-year
pa_pop_year_df <- units_df |>
  distinct(pa, pop, year, pa_pop_id) |>
  arrange(pa, pop, year) |>
  mutate(pa_pop_year_id = row_number()) |>
  print()

# 4. Vincular pa_pop_year_id de volta em units_df
units_df <- units_df |>
  left_join(pa_pop_year_df |>
              select(pa, pop, year, pa_pop_year_id),
            by = c("pa", "pop", "year")) |>
  print()

# 5. Definir sites únicos e mapear a seu respectivo combo pa-pop
sites_df <- data_long$ydata |>
  distinct(pa, pop, site) |>
  mutate(site_id = row_number()) |>
  left_join(pa_pop_combos |> select(pa, pop, pa_pop_id), by = c("pa", "pop")) |>
  print()

# 6. Mapear observações originais às unidades, sites e pa_pop_year
data_mapped <- data_long$ydata |>
  left_join(units_df |>
              select(pa, pop, site, year, unit_id, pa_pop_id, pa_pop_year_id),
            by = c("pa", "pop", "site", "year")) |>
  left_join(sites_df |>
              select(pa, pop, site, site_id),
            by = c("pa", "pop", "site")) |>
  print()

# 7. Obter índice contínuo 1..npa da PA (evitando conflito com coluna pa_id de ydata)
pa_lookup <- data_long$ydata |>
  distinct(pa) |>
  arrange(pa) |>
  mutate(pa_idx = row_number()) |>
  print()

data_mapped <- data_mapped |>
  left_join(pa_lookup, by = "pa") |>
  print()

# 8. Construir a lista completa para o JAGS
jags_data <- list(
  y = as.integer(data_mapped$y),
  nobs = nrow(data_mapped),
  
  # Dimensões
  npa = nrow(pa_lookup),
  n_pa_pop_combos = nrow(pa_pop_combos),
  n_total_sites = nrow(sites_df),
  n_combo_years = nrow(pa_pop_year_df),
  n_units = nrow(units_df),
  n_groups = max(units_df$group_id),
  
  # Dimensões do sub-modelo de r[j]
  n_mgmt    = nlevels(pa_pop_combos$mgmt_f),
  n_troph   = nlevels(pa_pop_combos$troph_f),
  n_threat  = nlevels(pa_pop_combos$threat_f),
  n_species = nlevels(pa_pop_combos$species_f),
  
  # Mapear vetores para parâmetros e efeitos aleatórios
  pa_map = as.integer(pa_pop_combos |>
                        left_join(pa_lookup, by = "pa") |>
                        pull(pa_idx)),
  combo_year_map = as.integer(pa_pop_year_df$pa_pop_id),
  combo_map = as.integer(sites_df$pa_pop_id),
  combo_unit = as.integer(units_df$pa_pop_id),
  site_unit = as.integer(units_df |>
                           left_join(sites_df,
                                     by = c("pa", "pop", "site", "pa_pop_id")) |>
                           pull(site_id)),
  year_unit = as.integer(units_df$year - min(units_df$year) + 1),
  group_id = as.integer(units_df$group_id),
  
  # Covariáveis e mapeamento de espécie para r[j]
  forest_loss_j = pa_pop_combos$forest_loss_z,
  pop_size_j    = pa_pop_combos$pop_size_z,
  body_mass_j   = pa_pop_combos$body_mass_z,
  mgmt_j        = pa_pop_combos$mgmt_idx,
  troph_j       = pa_pop_combos$troph_idx,
  threat_j      = pa_pop_combos$threat_idx,
  sp_map        = pa_pop_combos$sp_idx,
  
  # Mapear vetores no nível das observações
  combo_obs = as.integer(data_mapped$pa_pop_id),
  site_obs = as.integer(data_mapped$site_id),
  combo_year_obs = as.integer(data_mapped$pa_pop_year_id),
  unit_obs = as.integer(data_mapped$unit_id),
  pa_obs = as.integer(data_mapped$pa_idx)
)

str(jags_data)

# save jags_data
saveRDS(jags_data, "data-rev1/jags_data.rds")


#----- parametros a monitorar
params <- c(
  # essenciais
  'n',
  'p', 'r',
  "b_forest_loss", "b_pop_size",  "b_body_mass",
  "b_mgmt",  "b_troph", "b_threat",
  #"eps_sp_r[sp_map[j]]", "eps_pa_r[pa_map[j]]",
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
# ni <- 250
# nt <- 1
# nb <- 10
# nc <- 3


# para testes
# ni <- 500
# nt <- 1
# nb <- 100
# nc <- 3


# para testes
ni <- 1500
nt <- 2
nb <- 500
nc <- 3


# check predicted size of posterior
((ni-nb)*3)/nt

# remover objetos para liberar espaco
rm(list = c("mod"#, 
            #"pa_pop_combos", "pa_year_df",
            #"units_df", "sites_df", "data_mapped"
            )) 

# model_file <- "scripts-rev1/mod-with-glmm.txt" # check exact path
# normalizePath(model_file, mustWork = TRUE)
# lines <- readLines(model_file)
# grep("mgmt", lines, value = TRUE)
# tail(lines, 15)

# chamar o jags
system.time(mod <- jagsUI::jags(data = jags_data,
                                inits = inits,
                                parameters.to.save = params,
                                model.file = "scripts-rev1/mod-with-glmm.txt",
                                n.chains = nc,
                                n.thin = nt,
                                n.iter = ni,
                                n.burnin = nb,
                                parallel = FALSE))

# salvar
saveRDS(mod, "results-rev1/res.rds")

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

#traceplot(mod, parameter = "n")
#traceplot(mod, parameter = "n_sum")
#traceplot(mod, parameter = "p")

# check n
str(mod$sims.list$n)
hist(mod$sims.list$n)
hist(mod$n.eff$n)
#hist(mod$sims.list$n, xlim = c(0,100))


# check p
str(mod$sims.list$p)
hist(mod$sims.list$p)
hist(mod$n.eff$p)

# check r
str(mod$sims.list$r)
hist(mod$sims.list$r)


# check r predictors
hist(mod$sims.list$b_body_mass, 
     main="Body mass effect", xlab = "Effect")
abline(v = 0, lty = "dashed", lwd = 2, col = "red")

hist(mod$sims.list$b_pop_size, 
     main="Human population size effect", xlab = "Effect")
abline(v = 0, lty = "dashed", lwd = 2, col = "red")

hist(mod$sims.list$b_forest_loss, 
     main="Forest loss effect", xlab = "Effect")
abline(v = 0, lty = "dashed", lwd = 2, col = "red")

str(mod$sims.list$b_mgmt)
hist(mod$sims.list$b_mgmt[,2], 
     main="Sustainable use", xlab = "Effect")
abline(v = 0, lty = "dashed", lwd = 2, col = "red")

boxplot(mod$sims.list$b_mgmt, 
        names = c("strict", "use"), 
        col = c("skyblue", "salmon"),
        main = "Group Comparison",
        ylab = "Value")

hist(mod$sims.list$b_troph[,2], 
     main="Secondary consumer effect", xlab = "Effect")
abline(v = 0, lty = "dashed", lwd = 2, col = "red")

boxplot(mod$sims.list$b_troph, 
        names = c("primary", "secondary"), 
        col = c("skyblue", "salmon"),
        main = "Group Comparison",
        ylab = "Value")

hist(mod$sims.list$b_threat[,2], 
     main="Threatened status effect", xlab = "Effect")
abline(v = 0, lty = "dashed", lwd = 2, col = "red")

boxplot(mod$sims.list$b_threat, 
        names = c("non-threatened", "threatened"), 
        col = c("skyblue", "salmon"),
        main = "Group Comparison",
        ylab = "Value")

npost <- units_df |>
  mutate(n = apply(mod$sims.list$n, 2, mean)) |>
  group_by(pa,pop,year) |>
  summarise(mean_n = mean(n)) |>
  print()
hist(npost$mean_n)

units_df |>
  mutate(n = apply(mod$sims.list$n, 2, mean),
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



units_df |>
  mutate(n = apply(mod$sims.list$n, 2, mean),
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


