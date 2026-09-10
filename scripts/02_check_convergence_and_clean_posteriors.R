
# load libraries
library(here)
library(tidyverse)
library(jagsUI)
#library(monitoraSSM)

# load data and results
data_wide <- readRDS(here("data", "data_wide.rds"))
data_long <- readRDS(here("data", "data_long.rds"))
res <- readRDS(here("results", "res_tss_multiregion.rds"))


#----- check overall convergence

res$summary
hist(res$summary[,"Rhat"], main = "", xlab = "Rhat")
abline(v = 1.1, lty = "dashed", lwd = 2, col = "red")
# which parameters did not converged?
length(res$summary[which(res$summary[,"Rhat"] > 1.1),])
round(length(res$summary[which(res$summary[,"Rhat"] > 1.1),])/length(res$summary),2)
# traceplot check convergence
#traceplot(res, parameters = "n")
#traceplot(res, parameters = "p")
#traceplot(res, parameters = "a1")

# NB! the summary above includes "dummy nodes" (unobserved sites, pops ad years)
# these nodes are not linked to the likelihood so they do not affect parameters of interest
# however, they may still clutter convergence checks or summary tables
# so lets filter them out of the posteriors
# and after that we check parameter convergence again


#----- check convergence by parameter separately

#----- check n convergence

# first, let us create mask
# so that we retain only valid pops, sites and years in the posterior
npa <- data_long$npa
npop <- data_long$npop
nsite <- data_long$nsite
nyear <- data_long$nyear
npop_pa <- data_long$indexes$npops
nsite_pa <- data_long$indexes$nsite
start_year <- data_long$indexes$start_year-2013
end_year <- data_long$indexes$end_year-2013

valid_n_mask <- array(FALSE, dim = c(npa, npop, nsite, nyear))
for(pa in 1:npa) {
  valid_n_mask[pa, 1:npop_pa[pa], 1:nsite_pa[pa], start_year[pa]:end_year[pa]] <- TRUE
}

# check
str(valid_n_mask)
valid_n_mask[1,1,,] # 1st pop of 1st pa
valid_n_mask[2,1,,] # 1st pop of 2nd pa
valid_n_mask[22,1,,] # 1st pop of 22nd pa
valid_n_mask[34,1,,] # 1st pop of 34th pa

# use mask to "clean" the Rhat array
rhat_n <- res$Rhat$n
# set non-valid sites to NA using the inverted mask
rhat_n[!valid_n_mask] <- NA
rhat_n[1,1,,] # 1st pop of 1st pa
rhat_n[2,1,,] # 1st pop of 2nd pa
rhat_n[24,1,,] # 1st pop of 2nd pa
# check summary stats or max Rhat across valid sites only
max(rhat_n, na.rm = TRUE)
mean(rhat_n, na.rm = TRUE)
# filter for valid sites with potential convergence issues
which(rhat_n > 1.1, arr.ind = TRUE)

hist(rhat_n, breaks = 20, main = "", xlab = "Rhat")
abline(v = 1.1, lty = "dashed", lwd = 2, col = "red")
round(length(rhat_n[which(rhat_n > 1.1)])/length(rhat_n),2)

#traceplot(res, parameters = "n")

# store it back into the results list
res$Rhat$n <- rhat_n


# use mask to "clean" the n array
n_clean <- res$sims.list$n
# invert the mask so TRUE represents invalid/dummy sites
invalid_mask <- !valid_n_mask
# apply the 4D mask across the 5D array using a loop
for (i in 1:dim(n_clean)[1]) {
  # Extract the i-th MCMC draw and replace invalid sites with NA
  draw <- n_clean[i, , , , ]
  draw[invalid_mask] <- NA
  n_clean[i, , , , ] <- draw
}
n_clean
n_clean[1,,,,]
n_clean[1,1,,,]
n_clean[1,2,,,]
n_clean[1,22,1,,]
n_clean[1,34,1,,]

# store it back into the results list
res$sims.list$n <- n_clean 


#----- check p convergence

# first, let us create mask
# so that we retain only valid pops, sites and years in the posterior
valid_p_mask <- array(FALSE, dim = c(npa, npop, nyear))
for(pa in 1:npa) {
  valid_p_mask[pa, 1:npop_pa[pa], start_year[pa]:end_year[pa]] <- TRUE
}

# check
str(valid_p_mask)
valid_p_mask[1,,] # 1st pop of 1st pa
valid_p_mask[2,,] # 1st pop of 2nd pa
valid_p_mask[22,,] # 1st pop of 22nd pa
valid_p_mask[34,,] # 1st pop of 34th pa

# use mask to "clean" the Rhat array
rhat_p <- res$Rhat$p
# set non-valid sites to NA using the inverted mask
rhat_p[!valid_p_mask] <- NA
rhat_p[1,,] # 1st pop of 1st pa
rhat_p[2,,] # 1st pop of 2nd pa
rhat_p[24,,] # 1st pop of 2nd pa
# check summary stats or max Rhat across valid sites only
max(rhat_p, na.rm = TRUE)
mean(rhat_p, na.rm = TRUE)
# filter for valid sites with potential convergence issues
which(rhat_p > 1.1, arr.ind = TRUE)

hist(rhat_p, breaks = 20, main = "", xlab = "Rhat")
abline(v = 1.1, lty = "dashed", lwd = 2, col = "red")
round(length(rhat_p[which(rhat_p > 1.1)])/length(rhat_p),2)

#traceplot(res, parameters = "p")

# store it back into the results list
res$Rhat$p <- rhat_p


# use mask to "clean" the p array
p_clean <- res$sims.list$p
# invert the mask so TRUE represents invalid/dummy sites
invalid_mask <- !valid_p_mask
# apply the 4D mask across the 5D array using a loop
for (i in 1:dim(p_clean)[1]) {
  # Extract the i-th MCMC draw and replace invalid sites with NA
  draw <- p_clean[i, , , ]
  draw[invalid_mask] <- NA
  p_clean[i, , , ] <- draw
}
p_clean
p_clean[1,,,]
p_clean[1,1,,]
p_clean[1,2,,]
p_clean[1,22,,]
p_clean[1,34,,]

# store it back into the results list
res$sims.list$p <- p_clean 


#----- check beta_1 convergence (written as a1 here)

# first, let us create mask
# so that we retain only valid pops, sites and years in the posterior
valid_a1_mask <- array(FALSE, dim = c(npa, npop))
for(pa in 1:npa) {
  valid_a1_mask[pa, 1:npop_pa[pa]] <- TRUE
}

# check
str(valid_a1_mask)
valid_a1_mask[1,] # 1st pa
valid_a1_mask[2,] # 2nd pa
valid_a1_mask[22,] # 22nd pa
valid_a1_mask[34,] # 34th pa

# use mask to "clean" the Rhat array
rhat_a1 <- res$Rhat$a1
# set non-valid sites to NA using the inverted mask
rhat_a1[!valid_a1_mask] <- NA
rhat_a1[1,] # 1st pop of 1st pa
rhat_a1[2,] # 1st pop of 2nd pa
rhat_a1[24,] # 1st pop of 2nd pa
# check summary stats or max Rhat across valid sites only
max(rhat_a1, na.rm = TRUE)
mean(rhat_a1, na.rm = TRUE)
# filter for valid sites with potential convergence issues
which(rhat_a1 > 1.1, arr.ind = TRUE)

hist(rhat_a1, breaks = 20, main = "", xlab = "Rhat")
abline(v = 1.1, lty = "dashed", lwd = 2, col = "red")
round(length(rhat_a1[which(rhat_a1 > 1.1)])/length(rhat_a1),2)

#traceplot(res, parameters = "a1")

# store it back into the results list
res$Rhat$a1 <- rhat_a1


# use mask to "clean" the p array
a1_clean <- res$sims.list$a1
# invert the mask so TRUE represents invalid/dummy sites
invalid_mask <- !valid_a1_mask
# apply the 4D mask across the 5D array using a loop
for (i in 1:dim(a1_clean)[1]) {
  # Extract the i-th MCMC draw and replace invalid sites with NA
  draw <- a1_clean[i, , ]
  draw[invalid_mask] <- NA
  a1_clean[i, , ] <- draw
}
a1_clean
a1_clean[1,,]
a1_clean[1,2,]
a1_clean[1,22,]
a1_clean[1,34,]

# store it back into the results list
res$sims.list$a1 <- a1_clean 


# save the updated results
saveRDS(res, here("results", "res.rds"))


########################################
########################################
########################################
# below is the older version of convergence check... delete later!!!
# ...only after confirming the other scripts in the pipeline are working!
########################################
########################################
########################################


n_summary <- res$summary[grepl(pattern = "n", x = rownames(res$summary)), ]
hist(n_summary[,"Rhat"], breaks = 20, main = "", xlab = "Rhat")
abline(v = 1.1, lty = "dashed", lwd = 2, col = "red")
round(length(n_summary[which(n_summary[,"Rhat"] > 1.1),])/length(n_summary),2)

# check r (a1) convergence
r_summary <- res$summary[grepl(pattern = "a1", x = rownames(res$summary)), ]
hist(r_summary[,"Rhat"], breaks = 20, main = "", xlab = "Rhat")
abline(v = 1.1, lty = "dashed", lwd = 2, col = "red")

# check p convergence
p_summary <- res$summary[grepl(pattern = "p", x = rownames(res$summary)), ]
hist(p_summary[,"Rhat"], breaks = 20, main = "", xlab = "Rhat")
abline(v = 1.1, lty = "dashed", lwd = 2, col = "red")


#----- use Rhat to exclude pops that did not converged
# convergence was an issue only for the "n" parameter so lets focus on it
# filter summary rows for "n"
# create array with Rhat values
n_rhat_array <- with(as.data.frame.matrix(n_summary),
                     array(data = round(Rhat,2),
                           dim=dim(res$sims.list$n)[2:5],
                           dimnames = NULL))
str(n_rhat_array)
n_rhat_array[1,1,,]
# then we should use apply to find pops that did not converged well
# these should be removed later
# example
#apply(n_rhat_array, c(1,2), max, na.rm=TRUE) # max Rhat for n for a given population

# criar template para excluir pops e sites nonexistent do Rhat
n_template <- apply(data_wide$y, c(5,4,1,3), max, na.rm=TRUE)
max_values <- apply(n_template, c(1,2,3), max, na.rm=TRUE)
max_values[max_values == Inf | max_values == -Inf ] <- NA
max_values[!is.na(max_values)] <- 1
str(max_values)
max_values[2,3,]
for (i in 1:dim(n_template)[1]) {
  for (j in 1:dim(n_template)[2]) {
    for (k in 1:dim(n_template)[3]) {
      n_template[i,j,k,] <- (max_values[i,j,k])
    }
  }
}
n_template[1,3,,] # check
n_template[2,3,,] # check

# multiply n_rhat_array and n_template 
#so that nonexistent pops and sites become NAs
n_rhat_array <- n_rhat_array*n_template
n_rhat_array[1,3,,]
n_rhat_array[2,3,,]
hist(n_rhat_array[!is.na(n_rhat_array)], breaks = 20, main = "", xlab = "Rhat")
abline(v = 1.1, lty = "dashed", lwd = 2, col = "red")
round(length(which(n_rhat_array[!is.na(n_rhat_array)] > 1.1))/length(n_rhat_array[!is.na(n_rhat_array)]),2)



#----- check which populations did not converged for N parameter
# we should exclude them from further analysis
str(n_rhat_array) # dimension 1 corresponds to PA and 2 to population
rhat_meanN_check <- apply(n_rhat_array, c(1,2), mean, na.rm=TRUE) # check mean Rhat
rhat_meanN_check <- round(rhat_meanN_check,2)
head(rhat_meanN_check)
rhat_meanN_check[is.nan(rhat_meanN_check)] <- NA
head(rhat_meanN_check)
rhat_meanN_check

# how many individual parameters did not converge
nonConverg <- as.vector(which(rhat_meanN_check > 1.1))
allN <- as.vector(rhat_meanN_check)
allN <- allN[!is.na(allN)]
length(nonConverg)
length(nonConverg)/length(allN)
rhat_meanN_check[which(rhat_meanN_check > 1.1)]
hist(rhat_meanN_check[which(rhat_meanN_check > 1.1)])
min(rhat_meanN_check[which(rhat_meanN_check > 1.16)])
max(rhat_meanN_check[which(rhat_meanN_check > 1.16)])


# make a dataframe with mean Rhat by pop
data_dimnames<- dimnames(data_wide$y)
dimnames(rhat_meanN_check) <- list(data_dimnames[[5]],
                                   data_dimnames[[4]])
str(rhat_meanN_check)
rhat_meanN_check  <- as.data.frame.table(rhat_meanN_check)
head(rhat_meanN_check)
rhat_meanN_check <- tibble(rhat_meanN_check) %>%
  mutate(pa = as.numeric(as.character(Var1)),
         pop = as.numeric(as.character(Var2)),
         Rhat = as.numeric(as.character(Freq))
  ) %>%
  select(pa,pop,Rhat) %>%
  drop_na() %>%
  print()

# join with pop templates so that we know which pops to exclude
pop_N_convergence <- data_wide$pop_indexes %>%
  left_join(rhat_meanN_check, by = c("pa","pop")) %>%
  print()

# which pops did not converge?
pop_N_convergence %>%
  filter(Rhat > 1.1) %>%
  print(n=Inf)

