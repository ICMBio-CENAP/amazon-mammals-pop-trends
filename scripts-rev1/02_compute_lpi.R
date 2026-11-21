
# compute Posterior LPI in R (Vectorized Across MCMC Draws)

calc_lpi_posterior <- function(n_draws, units_df, min_year = min(units_df$year),
                               offset = 1e-4, base_value = 100) {
  # n_draws: matrix [n_iter, n_units] from out$sims.list$n
  n_iter  <- nrow(n_draws)
  n_years <- max(units_df$year) - min_year + 1L
  n_combos <- max(units_df$pa_pop_id)
  
  yr_idx <- as.integer(units_df$year - min_year + 1L)
  
  # 1. Mean abundance across sites for each unit's (pa_pop_id, yr_idx)
  # Group unit indices by (pa_pop_id, yr_idx)
  unit_groups <- units_df |>
    mutate(u_row = row_number(), yr_idx = yr_idx) |>
    group_by(pa_pop_id, yr_idx) |>
    summarise(u_rows = list(u_row), .groups = "drop")
  
  # Array of mean abundance: [n_iter, n_combos, n_years]
  mean_abund <- array(NA_real_, dim = c(n_iter, n_combos, n_years))
  
  for (g in seq_len(nrow(unit_groups))) {
    j_i <- unit_groups$pa_pop_id[g]
    t_i <- unit_groups$yr_idx[g]
    cols <- unit_groups$u_rows[[g]]
    if (length(cols) == 1L) {
      mean_abund[, j_i, t_i] <- n_draws[, cols]
    } else {
      mean_abund[, j_i, t_i] <- rowMeans(n_draws[, cols, drop = FALSE])
    }
  }
  
  # 2. Consecutive-year log10 ratios per draw and combo: [n_iter, n_combos, n_years - 1]
  log10_d <- array(NA_real_, dim = c(n_iter, n_combos, n_years - 1L))
  for (t in seq_len(n_years - 1L)) {
    a_t0 <- mean_abund[, , t]
    a_t1 <- mean_abund[, , t + 1L]
    active <- !is.na(a_t0) & !is.na(a_t1)
    if (any(active)) {
      # Small offset avoids log10(0) if all sites in a combo-year have count 0
      num <- ifelse(active, a_t1 + offset, NA_real_)
      den <- ifelse(active, a_t0 + offset, NA_real_)
      log10_d[, , t] <- log10(num / den)
    }
  }
  
  # 3. Mean log10 change across active combos per transition & draw: [n_iter, n_years - 1]
  mean_d_t <- matrix(0, nrow = n_iter, ncol = n_years - 1L)
  for (t in seq_len(n_years - 1L)) {
    vals <- log10_d[, , t, drop = FALSE]
    # Mean across combos for each iteration ignoring NAs
    row_sums <- rowSums(vals, na.rm = TRUE)
    row_cnts <- rowSums(!is.na(vals))
    mean_d_t[, t] <- ifelse(row_cnts > 0L, row_sums / row_cnts, 0)
  }
  
  # 4. Accumulate LPI trajectory per draw: [n_iter, n_years]
  lpi_draws <- matrix(base_value, nrow = n_iter, ncol = n_years)
  for (t in 2:n_years) {
    lpi_draws[, t] <- lpi_draws[, t - 1L] * (10 ^ mean_d_t[, t - 1L])
  }
  
  # 5. Posterior summary by calendar year
  years_vec <- seq(min_year, max(units_df$year))
  summary_df <- tibble(
    year   = years_vec,
    median = apply(lpi_draws, 2, median),
    mean   = colMeans(lpi_draws),
    l95    = apply(lpi_draws, 2, quantile, probs = 0.025),
    u95    = apply(lpi_draws, 2, quantile, probs = 0.975)
  )
  
  list(lpi_draws = lpi_draws, summary = summary_df)
}

# Usage after fitting with jagsUI:
# lpi_res <- calc_lpi_posterior(mod$sims.list$n, units_df)
# print(lpi_res$summary)

lpi_res <- calc_lpi_posterior(mod$sims.list$n, units_df, base_value = 1)
print(lpi_res$summary)

# save
saveRDS(lpi_res, "results-rev1/lpi_res.rds")
# read
#lpi_res <- readRDS("results-rev1/lpi_res.rds")

# plot
lpi_res$summary |>
  ggplot(aes(x = year + 2013, y = mean)) +
  geom_ribbon(aes(ymin = l95, ymax = u95), alpha = 0.7, fill = "steelblue") +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey50") +
  geom_line(linewidth = 0.9) +
  labs(
    x = "Year",
    y = "Living Planet Index",
    title = "Posterior Living Planet Index with 95% Credible Interval"
  ) +
  theme_minimal()

