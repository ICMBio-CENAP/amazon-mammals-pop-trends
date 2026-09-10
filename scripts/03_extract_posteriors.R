
# load libraries
library(here)
library(tidyverse)
#library(monitoraSSM)

# load data and results
data_wide <- readRDS(here("data", "data_wide.rds"))
data_long <- readRDS(here("data", "data_long.rds"))
#res <- readRDS(here("results", "res_tss_multiregion.rds"))
res <- readRDS(here("results", "res.rds"))

#----- create mask to filter the posterior
# sampled_year should tell if the pa was sampled on a given year
# regardless of which and how many sites were sampled on that year
# sampled_site should tell if a given site was ever sampled in a given pa
# regardless of in which years the site was sampled
# start_date and end_date should be used to crop the posterior
# so that only estimates from inside the monitoring period are kept
post_template <- tibble(pa = data_long$pa,
                        site = data_long$site) %>%
  distinct(pa,site) %>%
  mutate(sampled_site = "y") %>%
  full_join(tibble(pa = data_long$pa,
                   year = data_long$year) %>%
              distinct(pa,year) %>%
              mutate(sampled_year = "y"),
            by = "pa",
            relationship = "many-to-many") %>%
  full_join(data_long$indexes %>%
              select(pa,pa_name,start_year,end_year),
            by = "pa",
            relationship = "many-to-many") %>%
  select(pa,pa_name,site,sampled_site,year,sampled_year,start_year,end_year) %>%
  print()

# quick check
post_template %>%
  group_by(pa) %>%
  summarise(nsite = n_distinct(site),
            nyear = n_distinct(year)) %>%
  print()


#----- rename posteriors and remove objects from environment to save space
#lambda <- res$sims.list$lambda
n <- res$sims.list$n
p <- res$sims.list$p
r <- res$sims.list$a1
pop_indexes <- data_wide$pop_indexes
data_dimnames<- dimnames(data_wide$y)
n_samples <- res$mcmc.info$n.samples
rm(list=setdiff(ls(),
                c(#"lambda",
                  "n", "p", "r",
                  "post_template", "pop_indexes",
                  "n_samples", "data_dimnames")))


#----- get lambda posterior

# lambda posterior to vector format
#lambda <- res$sims.list$lambda
#str(lambda)
#dimnames(lambda) <- list(1:n_samples,
#                         data_dimnames[[5]],
#                         data_dimnames[[4]],
#                         data_dimnames[[1]],
#                         data_dimnames[[3]])
#str(lambda)
#lambdapost  <- as.data.frame.table(lambda)
#head(lambdapost)
#lambdapost <- tibble(lambdapost) %>%
#  mutate(iter = as.numeric(as.character(Var1)),
#         pa = as.numeric(as.character(Var2)),
#         site = as.numeric(as.character(Var3)),
#         pop = as.numeric(as.character(Var4)),
#         year = as.numeric(as.character(Var5)),
#         lambda = round(as.numeric(as.character(Freq)),2)
#  ) %>%
#  select(iter,pa,pop,site,year,lambda) %>%
#  print()

# keep only pa,and sites that are in actual data
#lambdapost <- lambdapost %>%
#  semi_join(post_template,
#            join_by(pa,site)) %>%
#  print()

# keep only pops that are in actual data
#lambdapost <- lambdapost %>%
#  semi_join(pop_indexes,
#            join_by(pa,pop)) %>%
#  print()

# calculate average lambda for iter, pa and pop
#lambdapost %>%
#  group_by(iter,pa,pop) %>%
#  summarise(lambda = mean(lambda)) %>%
#  ungroup() %>%
#  print()

# add pa_names and species names
#lambdapost <- post_template %>%
#  left_join(tibble(pop_indexes),
#            by = c("pa")) %>%
#  select(pa, pa_name,pop,species) %>%
#  drop_na(pa_name,species) %>%
#  left_join(lambdapost,
#            by = c("pa", "pop")) %>%
#  print()


# remove lambda to save space
#rm(lambda)

#----- get n posterior

# convert n to vector format
#n <- res$sims.list$n
dimnames(n) <- list(1:n_samples,
                    data_dimnames[[5]],
                    data_dimnames[[4]],
                    data_dimnames[[1]],
                    data_dimnames[[3]])
str(n)
npost  <- as.data.frame.table(n)
head(npost)
npost <- tibble(npost) %>%
  mutate(iter = as.numeric(as.character(Var1)),
         pa = as.numeric(as.character(Var2)),
         pop = as.numeric(as.character(Var3)),
         site = as.numeric(as.character(Var4)),
         year = as.numeric(as.character(Var5)),
         Nest = as.numeric(as.character(Freq))
  ) %>%
  select(iter,pa,pop,site,year,Nest) %>%
  print()

# keep only pa and pops that are in actual data
npost <- npost %>%
  semi_join(pop_indexes,
            join_by(pa,pop)) %>%
  print()

# keep only pa and sites that are in actual data
npost <- npost %>%
  semi_join(post_template,
            join_by(pa,site)) %>%
  print()

# remove n to save space
rm(n)
#rm(list = c("lambdapost", "p", "r")) # to save even more space

# add pa_names and species names
#npost <- post_template %>%
#  left_join(tibble(pop_indexes),
#            by = c("pa")) %>%
#  select(pa, pa_name,pop,species) %>%
#  drop_na(pa_name,species) %>%
#  left_join(npost,
#            by = c("pa", "pop")) %>%
#  print()

# add species names (and just in case remove again non-existent populations)
npost <- npost %>%
  left_join(pop_indexes,
            by = c("pa", "pop")) %>%
  filter(!is.na(species)) %>%
  print()

# add start and end year so we can crop (keep only years inside "core" monitoring)
npost <- npost %>%
  left_join(post_template %>%
              distinct(pa, pa_name, start_year, end_year),
            by = c("pa")) %>%
  print()

# crop years outside "core" monitoring window
npost <- npost %>%
  filter(year >= start_year,
         year <= end_year) %>%
  select(-c(start_year, end_year)) %>%
  print()


# add info on sampled/non-sampled years
npost <- npost %>%
  left_join(post_template %>%
              distinct(pa,year,sampled_year),
            by = c("pa","year")) %>%
  mutate(sampled_year = case_when(is.na(sampled_year) ~ "n",
                                  .default = sampled_year)) %>%
  print()

# reorder columns
npost <- npost %>%
  select(iter,Nest,pa,pa_name,species,site,year,sampled_year) %>%
  print()

# we only want annual estimates per species per pa
# lets simplify npost by summarizing N means across sites
npost <- npost %>%
  group_by(iter,pa,pa_name,species,year,sampled_year) %>%
  summarise(n = mean(Nest)) %>%
  ungroup() %>%
  print()


#----- get growth rate
#str(res$sims.list$r)
str(r)
dimnames(r) <- list(1:n_samples,
                    data_dimnames[[5]],
                    data_dimnames[[4]])
str(r)
rpost  <- as.data.frame.table(r)
head(rpost)
rpost <- tibble(rpost) %>%
  mutate(iter = as.numeric(as.character(Var1)),
         pa = as.numeric(as.character(Var2)),
         pop = as.numeric(as.character(Var3)),
         r = as.numeric(as.character(Freq))
  ) %>%
  select(iter,pa,pop,r) %>%
  print()

# keep only pops that are in actual data
rpost <- rpost %>%
  semi_join(pop_indexes,
            join_by(pa,pop)) %>%
  print()


rpost <- rpost %>%
  left_join(post_template %>%
              distinct(pa,pa_name),
            by = "pa") %>%
  left_join(pop_indexes,
            by = c("pa", "pop")) %>%
  select(iter, r, pa, pa_name, species) %>%
  drop_na(pa_name, species) %>%
  print()



#----- get detection p
#str(res$sims.list$p)
str(p)
dimnames(p) <- list(1:n_samples,
                    data_dimnames[[5]],
                    data_dimnames[[4]],
                    data_dimnames[[3]])
str(p)
ppost  <- as.data.frame.table(p)
head(ppost)
ppost <- tibble(ppost) %>%
  mutate(iter = as.numeric(as.character(Var1)),
         pa = as.numeric(as.character(Var2)),
         pop = as.numeric(as.character(Var3)),
         year = as.numeric(as.character(Var4)),
         p = as.numeric(as.character(Freq))
  ) %>%
  select(iter,pa,pop,year,p) %>%
  print()

# keep only pops that are in actual data
ppost <- ppost %>%
  semi_join(pop_indexes,
            join_by(pa,pop)) %>%
  print()

# add names
ppost <- ppost %>%
  left_join(post_template %>%
              distinct(pa, pa_name),
            by = "pa") %>%
  left_join(pop_indexes,
            by = c("pa", "pop")) %>%
  drop_na(pa_name, species) %>%
  select(iter,p,pa,pa_name,species,year) %>%
  print()

#save posteriors for later use
#posteriors <- list(npost = npost,
#                   lambdapost = lambdapost,
#                   rpost = rpost,
#                   ppost = ppost)
#saveRDS(posteriors, here("tss", "tss_posteriors.rds"))
#saveRDS(posteriors, here("tss", "tss_multiregion_posteriors.rds"))
saveRDS(npost, here("tss", "tss_npost.rds"))
saveRDS(rpost, here("tss", "tss_rpost.rds"))
saveRDS(ppost, here("tss", "tss_ppost.rds"))

