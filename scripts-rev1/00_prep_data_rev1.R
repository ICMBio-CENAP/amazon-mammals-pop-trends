# Preparar dados para analises

# limpar ambiente de trabalho
rm(list = ls())

# carregar pacotes
library(tidyverse)
library(janitor)


# ler e ajustar dados
monitora <- readRDS("data-rev1/monitora_2025.rds") |>
  print()


pas <- monitora$ucs |>
  rename(pa = cnuc,
         pa_name = uc,
         category = categoria,
         biome = bioma) |>
  mutate(pa_name = make_clean_names(pa_name)) |>
  mutate(category = case_when(category == "protecao_integral" ~ "strictly_protected",
                              category == "uso_sustentavel" ~ "sustainable_use")) |>
  filter(biome == "amazonia") |> # NB! AMAZON ONLY
  select(-sigla) |>
  print()

samples <- monitora$amostragens |>
  rename(pa = cnuc,
         year = ano,
         date = data,
         site = ea,
         site_name = nome_ea,
         effort = esforco) |>
  group_by(pa, site, year) |>
  arrange(date) |>
  mutate(rep = dense_rank(date)) |>
  ungroup() |>
  relocate(rep, .before = year) |>
  filter(pa %in% pas$pa) |> # NB! AMAZON ONLY
  print()

records <- monitora$registros |>
  rename(pa = cnuc,
         site = ea,
         year = ano,
         date = data,
         class = classe,
         order = ordem,
         family = familia,
         genus = genero,
         species = especie,
         #id_level = nivel_taxon,
         grp_size = n_animais,
         dist_perpendicular = distancia,
         dist_trail = plaqueta,
         timestamp = hora_avistamento,
         pop = populacao) |>
  mutate(pop = str_to_lower(pop)) |>
  relocate(pop, .after = species) |>
  left_join(samples |>
              distinct(pa,site,date,rep),
            by = c("pa","site","date")) |>
  relocate(rep, .before = class) |>
  select(pa,site,year,rep,
         class,order,family,genus,species,#id_level,
         pop,grp_size) |>
  # remove NA species
  drop_na(species,pop) |>
  #filter(class == "Mammalia") |> # NB! MAMMALS ONLY
  print()

species <- monitora$especies |>
  rename(class = classe,
         order = ordem,
         family = familia,
         genus = genero,
         species = especie#,
         #red_list = categoria,
         #mass = massa,
         #generation_length_d = tempo_geracao_d,
         #sociality = socialidade,
         #trophic_level = nivel_trofico,
         #foraging_stratum = estrato_vertical,
         #activity_cycle = atividade
  ) |>
  filter(class == "Mammalia") |> # NB! MAMMALS ONLY
  print()

rm(monitora)


# keep only 2014-2024
samples <- samples |>
  filter(year < 2025)

records <- records |>
  filter(year < 2025)


# keep only PAs with >= 3 years of data
minimum_years <- samples |>
  group_by(pa) |>
  summarise(n_years = n_distinct(year)) |> 
  filter(n_years > 2) |>
  pull(pa) |>
  print()

pas <- pas |>
  filter(pa %in% minimum_years) |>
  print()

samples <- samples |>
  filter(pa %in% minimum_years) |>
  print()

records <- records |>
  filter(pa %in% minimum_years) |>
  print()



# checar se ha EAs amostradas num unico ano
samples |> arrange(pa, site, site_name) |>
  distinct(pa, site_name, year) |>
  group_by(pa, site_name) |>
  count() |>
  print(n=Inf)
# foram amostrados num unico ano:
# pa 68 ea_3
# pa 123 trilha_takuara
# pa 242 agua_branca
# pa 256 missilene_rodrigues
# pa 258 ea_3

# entao talvez ideal seja:
# primeiro filtrar EAs com >= 3 anos
# segundo fitrar UCs com >= 3 anos
# (o passo acima vai excluir UCs com > 3 anos, mas sem nenhuma EA com > 2 anos)


# keep only EAs with > 2 years of data
minimum_ea_years <- samples |>
  group_by(pa, site) |>
  summarise(n_years = n_distinct(year)) |> 
  filter(n_years > 2) |>
  mutate(paSite = paste(pa,site, sep = "_")) |>
  print()
# keep samples in minimum_ea_years
samples <- samples |>
  filter(paste(pa, site, sep = "_") %in% minimum_ea_years$paSite) |>
  print()
# keep only pas that remained in the updated samples object
pas <- pas |>
  filter(pa %in% samples$pa) |>
  print()
# keep only records from sites that remained in the updated minimum_ea_years
# and in the updated samples object
records <- records |>
  filter(paste(pa, site) %in% paste(samples$pa, samples$site)) |>
  print()

# checar novamente se ha EAs amostradas em < 3 anos
samples |> arrange(pa, site, site_name) |>
  distinct(pa, site_name, year) |>
  group_by(pa, site_name) |>
  count() |>
  print(n=Inf)


# keep only PAs with >= 2 sites
minimum_sites <- samples |>
  group_by(pa) |>
  summarise(n_sites = n_distinct(site)) |> 
  filter(n_sites > 1) |>
  pull(pa) |>
  print()

pas <- pas |>
  filter(pa %in% minimum_sites) |>
  print()

samples <- samples |>
  filter(pa %in% minimum_sites) |>
  print()

records <- records |>
  filter(pa %in% minimum_sites) |>
  print()



#----- define target species/populations
targets <- species |>
  # > 250 g
  filter(mass >= 250) |>
  # somente mamiferos
  filter(class == "Mammalia") |>
  # somente noturnos
  filter(activity_cycle != "nocturnal") |>
  # somente florestais, nativos
  filter(! genus %in% c("Ozotoceros", "Chrysocyon", "Lycalopex", "Cavia", "Galea")) |>
  # somente nativos
  filter(! genus %in% c("Sus", "Canis")) |>
  # somente terrestres
  filter(! genus %in% c("Lontra", "Pteronura", "Hydrochoerus")) |>
  arrange(class, order, family, species) |>
  distinct(class, order, family, species) |>
  print(n=Inf) |>
  pull(species)

head(targets)

# keep only target species in records
records <- records |>
  filter(species %in% targets) |>
  print()

# check selected species again
records |> group_by(class, order, family, species) |>
  summarise(n_pas = n_distinct(pa)) |>
  arrange(class, order, family, species) |>
  print(n=Inf)

# remove rare populations from raw_counts
# criteria: a population must be recorded in > 10% of samples in a given PA
# and in > 50% of EAs in a given PA
# species should be abundant enough to be reliably detected and counted multiple times
selected <- samples |>
  group_by(pa) |>
  summarise(n_days = n_distinct(site,year,rep),
            n_sites = n_distinct(site)) |>
  ungroup() |>
  full_join(records |>
              filter(class == "Mammalia") |>
              group_by(pa,species) |>
              summarise(n_days_sp = n_distinct(site,year,rep),
                        n_sites_sp = n_distinct(site))) |>
  mutate(prop_samples = n_days_sp/n_days,
         prop_sites = n_sites_sp/n_sites) |>
  select(pa,species,n_days,n_days_sp,prop_samples, prop_sites) |>
  mutate(pop = paste(pa, species)) |>
  arrange(desc(prop_samples)) |> #View()
  filter(prop_samples >= 1/10,
         #prop_sites >= 1/2
  ) |>
  print(n=Inf) |>
  pull(pop)


records <- records |>
  filter(paste(pa,species) %in% selected) |>
  print()

# check selected species again
records |>
  #filter(class == "Mammalia") |>
  filter(paste(pa,species) %in% selected) |>
  group_by(class, order, family, species) |>
  summarise(n_pas = n_distinct(pa)) |>
  arrange(class, order, family, species) |>
  print(n=Inf)

# let us check PAs with more than one squirrel
records |>
  filter(family == "Sciuridae") |>
  distinct(pa,species) |>
  group_by(pa) |>
  count() |>
  filter(n > 1) |>
  print(n=Inf)
# only Cazumba-Iracema and PN Jau have more than one squirrel

# let us check how many records of each species
records |>
  filter(pa %in% c(173,232), family == "Sciuridae") |>
  group_by(pa, species) |>
  count()

# since Cazumba have congeneric squirrels, let us group them as spadiceus (the more common)
# records <- records |>
#   mutate(species = case_when(pa == 232 & genus == "Hadrosciurus" ~ "Hadrosciurus spadiceus",
#          .default = species)) |>
#   print()
#   

# criteria: to be rare a population must be recorded in
# less than 2/3 of sites AND years
#too_rare <- records |>
#selected <- samples |>
#  group_by(pa) |>
#  summarise(n_sites_pa = n_distinct(site),
#            n_years_pa = n_distinct(year)) |>
#  full_join(records |>
#              group_by(pa,species) |>
#              summarise(n_sites_sp = n_distinct(site),
#                        n_years_sp = n_distinct(year))) |>
#  mutate(prop_sites = n_sites_sp/n_sites_pa,
#         prop_years = n_years_sp/n_years_pa) |>
#  select(pa, species,
#         prop_sites, prop_years) |>
#  arrange(prop_sites, prop_years) |>
#  filter(prop_sites >= 2/3 & prop_years >= 2/3) |>
#  arrange(species) |>
#  ungroup() |>
#  mutate(pop = paste(pa, species)) |>
#  print(n = Inf) |>
#  distinct(pop) |>
#  pull(pop)

#head(selected)

#records <- records |>
#  filter(paste(pa,species) %in% selected) |>
#filter(! pop %in% too_rare) |>
#  print()



# criar coluna pop com id_unico para cada populacao
records <- records |>
  arrange(pa, class, order, family, genus, species) |>
  group_by(pa) |>
  mutate(
    # calcular rank e formato como string de dois digitos
    pop_number = sprintf("%02d", dense_rank(paste(genus, species)))
  ) |>
  mutate(pop_id = paste(pa, pop_number, sep = "_")) |>
  select(-pop_number) |>
  #mutate(pop_number = as.numeric(pop_number)) |>
  ungroup() |>
  print()


#----- criar objeto y que sera input para modelo N-mixture
# primeiramente, criar raw_counts com contagens incluindo zeros

raw_counts <- samples |>
  select(pa,site,year,rep) |>
  full_join(records |> 
              mutate(n = 1),
            relationship = "many-to-many") |>
  select(pa,site,year,rep,pop_id,n) |>
  mutate(n = case_when(is.na(n) ~ 0,
                       .default = n)) |>
  # obter contagens por pop por rep
  group_by(pa, pop_id, year, site, rep, .drop=FALSE) |>
  summarise(count = sum(n, na.rm=TRUE)) |>
  ungroup() |>
  # agrupar por uc, e dentro do agrupamento
  # cruzar todas as possiveis combinacoes de pops com
  # `(ea, ano e rep)` ja existentes nos dados
  # e preencher n = 0 para obter contagens = zero
  group_by(pa) |>
  complete(pop_id,
           nesting(year, site, rep),
           fill = list(count = 0)) |>
  arrange(pa, pop_id, year, rep, site, count) |>
  ungroup() |>
  # o primeiro full_join nesta secao criou algumas pop_id == NA
  # elas prepresentam dias de amostragem sem nenhum registro das pops alvo
  # elas precisavam ser mantidas ate usarmos a funcao complete para
  # atribuir o zero para as opos alvo nos dias sem contagem
  # agora podemos remover:
  drop_na(pop_id) |>
  # finalmente, vamos adicionar coluna pop para representar id da populacao DENTRO da UC
  # ao contrario de pop_id, pop pode ter valores repetidos (reinicia contagem para cada uc)
  group_by(pa) |>
  mutate(pop = dense_rank(pop_id)) |>
  relocate(pop, .after = pop_id) |>
  ungroup() |>
  print()


# alguns sjustes adicionais
raw_counts <- raw_counts |>
  rename(y = count) |>
  mutate(pa_id = pa,
         pa = dense_rank(pa),
         year = year-2013) |>
  #select(cnuc,pop_id,ano,pa,pop,year,site,rep,y) |>
  select(pa_id,pop_id,year,pa,pop,year,site,rep,y) |>
  print()


#----- preparar os dados em long format para input no modelo JAGS

nobs <- nrow(raw_counts)
npop <- max(raw_counts$pop)
npa <- length(unique(raw_counts$pa))
nsite <- length(unique(raw_counts$site))
nyear <- length(unique(raw_counts$year))
nrep <- length(unique(raw_counts$rep))

# criar indices de uc (sera usado mais tarde para associar com a posterior)
pa_indexes <- raw_counts |>
  left_join(pas |> distinct(pa)) |>
  group_by(pa_id, pa) |>
  summarise(nsite = n_distinct(site),
            nyear = n_distinct(year),
            npops = n_distinct(pop_id),
            start = min(year),
            end = max(year),
            year_start = min(year) +2013,
            year_end = max(year) +2013) |>
  ungroup() |>
  tibble() |>
  print()

# add pa info
pa_info <- read_rds("data-rev1/pa_covars.rds") |>
  print()

pa_indexes <- pa_indexes |>
  left_join(pa_info |>
              rename(pa_id  = cnuc),
            by = "pa_id") |>
  relocate(pa_name, .after = pa_id) |>
  print()

# ajustes em pa_indexes
pa_indexes <- pa_indexes |>
  mutate(pa_type = factor(pa_type, levels = c("strictly_protected", "sustainable_use")),
         #forest_loss = as.numeric(scale(loss)),
         #pop_total = as.numeric(scale(log(pop_total)))
         ) |>
  print()


# criar indices de populacao/especie
pop_indexes <- raw_counts |>
  left_join(records |> distinct(genus, species, pop_id) |> drop_na(species)) |>
  distinct(pa_id,pa, pop_id, pop, genus, species) |>
  arrange(pa_id,pop_id) |>
  print()


# add species traits
pop_indexes <- pop_indexes |>
  left_join(species |>
              distinct(genus, species,
                       .keep_all = TRUE)) |> 
  select(pa_id, pa, pop_id, pop, 
         class, order, family, genus, species, 
         categoria, mass, generation_length_d, sociality,
         trophic_level, foraging_stratum, activity_cycle) |>
  print()
  

# ajustes nos traits
pop_indexes <- pop_indexes |>
  mutate(#mass = as.numeric(scale(log(mass))),
         sociality = as.factor(sociality),
         vertical = as.factor(foraging_stratum),
         trophic_level = case_when(trophic_level %in% c("omnivore", "carnivore") ~ "secondary",
                             .default = "primary"),
         trophic_level = factor(trophic_level, levels = c("primary", "secondary")),
         categoria = case_when(is.na(categoria) ~ "not_threatened",
                             .default = "threatened")) |>
  print()




# criar objeto daa_long para proxima etapa
data_long <- list(ydata = raw_counts,
                  nobs = nobs,
                  npop = npop,
                  npa = npa,
                  nsite = nsite,
                  nyesar = nyear,
                  nrep = nrep,
                  pa_indexes = pa_indexes,
                  pop_indexes = pop_indexes)

saveRDS(data_long, "data-rev1/data_long.rds")

