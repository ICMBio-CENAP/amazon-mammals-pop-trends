
# load libraries
library(here)
library(tidyverse)
#library(monitoraSSM)

# load data and results
data_wide <- readRDS(here("data", "data_wide.rds"))
data_long <- readRDS(here("data", "data_long.rds"))
#res <- readRDS(here("results", "res_tss_multiregion.rds"))
res <- readRDS(here("results", "res.rds"))

#----- rename posteriors and remove objects from environment to save space
n <- res$sims.list$n
p <- res$sims.list$p
r <- res$sims.list$a1
pa_indexes <- data_wide$indexes
pop_indexes <- data_wide$pop_indexes
data_dimnames<- dimnames(data_wide$y)
n_samples <- res$mcmc.info$n.samples
rm(list=setdiff(ls(),
                c("n", "p", "r",
                  "pa_indexes", "pop_indexes",
                  "n_samples", "data_dimnames")))


#----- get n posterior

# convert n to vector format
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
# (i.e., remove "dummy nodes")
npost <- npost %>%
  filter(! is.na(Nest)) %>%
  print()

# we only want annual estimates per species per pa
# lets simplify npost by summarizing N means across sites
npost <- npost %>%
  group_by(iter,pa,pop,year) %>%
  summarise(n = mean(Nest)) %>%
  ungroup() %>%
  print()

# remove n to save space
rm(n)

# add pa name, species names and start and end years
npost <- npost %>%
  left_join(pa_indexes[,c("pa", "pa_name", "start_year", "end_year")],
            by = c("pa")) %>%
  left_join(pop_indexes,
            by = c("pa", "pop")) %>%
  drop_na(pa, species) %>%
  print()

# add info on sampled/non-sampled years
# npost <- npost %>%
#   left_join(post_template %>%
#               distinct(pa,year,sampled_year),
#             by = c("pa","year")) %>%
#   mutate(sampled_year = case_when(is.na(sampled_year) ~ "n",
#                                   .default = sampled_year)) %>%
#   print()

# reorder columns
npost <- npost %>%
  #select(iter,n,pa,pa_name,species,year,sampled_year) %>%
  select(iter,n,pa,pa_name,species,year) %>%
  print()


#----- get per-individual detection
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

# keep only pa and pops that are in actual data
# (i.e., remove "dummy nodes")
rpost <- rpost %>%
  filter(! is.na(r)) %>%
  print()

rpost <- rpost %>%
  left_join(pa_indexes[,c("pa", "pa_name")],
            by = c("pa")) %>%
  left_join(pop_indexes,
            by = c("pa", "pop")) %>%
  select(iter, r, pa, pa_name, species) %>%
  drop_na(pa_name, species) %>%
  print()


#----- get per-species detection p
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

# keep only pa and pops that are in actual data
# (i.e., remove "dummy nodes")
ppost <- ppost %>%
  filter(! is.na(p)) %>%
  print()

# add names
ppost <- ppost %>%
  left_join(pa_indexes[,c("pa", "pa_name")],
            by = c("pa")) %>%
  left_join(pop_indexes,
            by = c("pa", "pop")) %>%
  select(iter, p, pa, pa_name, species, year) %>%
  drop_na(pa_name, species) %>%
  print()

#save posteriors for later use
saveRDS(npost, here("results", "tss_npost.rds"))
saveRDS(rpost, here("results", "tss_rpost.rds"))
saveRDS(ppost, here("results", "tss_ppost.rds"))

