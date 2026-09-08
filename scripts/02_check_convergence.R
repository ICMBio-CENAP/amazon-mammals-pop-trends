
# load libraries
library(here)
library(tidyverse)
library(jagsUI)
#library(monitoraSSM)

# load data and results
data_wide <- readRDS(here("data", "data_wide.rds"))
data_long <- readRDS(here("data", "data_long.rds"))
#res <- readRDS(here("results", "res_tss.rds"))
res <- readRDS(here("results", "res_tss_multiregion.rds"))
#res <- readRDS(here("results", "tss_single_species.rds"))
#res <- readRDS(here("results", "tss_single_species_fixed_p.rds"))


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

#----- check convergence by parameter separately
# check n convergence
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

