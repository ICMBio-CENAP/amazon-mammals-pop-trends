# Explorar os dados do programa Monitora
# obs: para verificar o significado das colunas, ver dicionario de dados em:
# https://docs.google.com/spreadsheets/d/1uMWgH0lF--7ZRHhXMReOhqdtvQ8pSDXpDRabuZBJNs0/edit?usp=sharing

# limpar ambiente de trabalho
rm(list = ls())

# carregar pacotes
library(tidyverse)


# ler e ajustar dados
monitora <- readRDS("01_entrada/monitora.rds") |>
  print()

# copiar objetos da lista monitora para o ambiente de trabalho
list2env(monitora, envir = .GlobalEnv) 

# checar
ucs
amostragens
registros
especies


#----- ucs

# quantas ucs, categorias, biomas
nrow(ucs)
ucs |>
  group_by(categoria) |>
  count()
ucs |>
  group_by(sigla) |>
  count()
ucs |>
  group_by(bioma) |>
  count()


#----- amostragens

# evolucao do programa ao longo do tempo
amostragens |>
  summarise(.by = c(cnuc, ano),
            n = n(),
            first_year = first(ano)) |>
  mutate(cnuc = fct_reorder(as.factor(cnuc), first_year, .fun = first)) |> 
  ggplot(aes(x = ano, y = cnuc, fill = n)) +
  geom_tile(size = 1, show.legend = FALSE) +
  scale_fill_gradient(low = "white", high = "darkblue") +
  #theme_classic() +
  theme_minimal() +
  labs(x = "", y = "UC")

# anos de amostragem por uc
amostragens |>
  summarise(.by = c(cnuc),
            n_anos = n_distinct(ano)) |>
  ggplot() +
  geom_histogram(aes(x = n_anos),
                 bins = 6,
                 colour = "black", fill = "steelblue",
                 linewidth = 0.2, alpha = 0.5) +
  theme_minimal() +
  labs(x = "Número de anos")

# distribuicao do esforco acumulado por uc
amostragens |>
  summarise(.by = cnuc,
            esforco_total = sum(esforco)) |>
  ggplot(aes(x = esforco_total)) +
  geom_histogram(colour = "black", fill = "steelblue",
                 linewidth = 0.2, alpha = 0.5) +
  theme_minimal() +
  labs(x = "Esforço acumulado (km)")


# distribuicao do esforco por amostragem individual
amostragens |>
  ggplot(aes(x = esforco)) +
  geom_histogram(colour = "black", fill = "steelblue",
                 linewidth = 0.2, alpha = 0.5) +
  theme_minimal() +
  labs(x = "Esforço (km)")


#----- registros

# quantas ordens familias, generos, especies na base de dados
registros |>
  group_by(classe) |>
  summarise(ordens = n_distinct(ordem),
            familias = n_distinct(familia),
            generos = n_distinct(genero),
            especies = n_distinct(genero, especie),
            populacoes = n_distinct(cnuc, genero, especie))

registros |>
  filter(classe == "Aves") |>
  group_by(classe, ordem, familia) |>
  summarise(n_ucs = n_distinct(cnuc),
            n_registros = n()) |>
  arrange(classe, ordem, familia) |>
  print(n=Inf)

registros |>
  filter(classe == "Mammalia") |>
  group_by(classe, ordem, familia) |>
  summarise(n_ucs = n_distinct(cnuc),
            n_registros = n()) |>
  arrange(classe, ordem, familia) |>
  print(n=Inf)

# quais as cinco especies de aves e mamiferos mais registradas
registros |>
  group_by(classe, genero, especie) |>
  count() |>
  group_by(classe) |>
  slice_max(order_by = n, n = 5, with_ties = FALSE) |>
  ungroup() |>
  print()

# idem, mas excluindo especie == NA
registros |>
  drop_na(especie) |>
  group_by(classe, genero, especie) |>
  count() |>
  group_by(classe) |>
  slice_max(order_by = n, n = 5, with_ties = FALSE) |>
  ungroup() |>
  print()
  

# numero total de registros por especie
registros |>
  group_by(genero, especie) |>
  count() |>
  ungroup() |>
  mutate(nome_cientifico = paste(genero, especie)) |> 
  mutate(nome_cientifico = fct_reorder(nome_cientifico, n, .desc = TRUE)) |> 
  ggplot(aes(x = nome_cientifico, y = n)) +
  geom_col() +
  theme_minimal() +
  scale_x_discrete(labels = NULL) +
  labs(x = "Espécie (ranking)", y = "Número de Registros")


# numero total de registros por populacao (especie + uc)
registros |>
  group_by(cnuc, genero, especie) |>
  count() |>
  ungroup() |>
  mutate(pop = paste(cnuc,genero, especie)) |> 
  mutate(pop = fct_reorder(pop, n, .desc = TRUE)) |> 
  ggplot(aes(x = pop, y = n)) +
  geom_col() +
  theme_minimal() +
  scale_x_discrete(labels = NULL) +
  labs(x = "População (ranking)", y = "Número de Registros")


# numero de ucs em que cada especie foi registrada
registros |>
  group_by(genero, especie) |>
  summarise(n_ucs = n_distinct(cnuc)) |>
  ungroup() |>
  mutate(nome_cientifico = paste(genero, especie)) |> 
  mutate(nome_cientifico = fct_reorder(nome_cientifico, n_ucs, .desc = TRUE)) |> 
  ggplot(aes(x = nome_cientifico, y = n_ucs)) +
  geom_col() +
  theme_minimal() +
  scale_x_discrete(labels = NULL) +
  labs(x = "Espécie (ranking)", y = "Número de UCs")


# percentagem de amostragens em que especie foi registrada
registros |>
  #drop_na(especie) |>
  group_by(cnuc, genero, especie) |>
  summarise(n_reps = n_distinct(ea, ano, data))|>
  ungroup() |>
  left_join(amostragens |>
              group_by(cnuc) |>
              count() |>
              ungroup(),
            by = "cnuc") |>
  mutate(percent_reps = n_reps/n) |>
  print() |>
  ggplot(aes(x = percent_reps)) +
  geom_histogram() +
  theme_minimal() +
  labs(x = "Percentagem das amostragens")


# checar serra do cipo
# obter cnuc
ucs |>
  filter(grepl('cipo', uc))
amostragens |>
  filter(cnuc == "0000.00.0148") |>
  print()
amostragens |>
  filter(cnuc == "0000.00.0148") |>
  summarise(n_ano = n_distinct(ano),
            n_eas = n_distinct(ea),
            n_amostras = n()) |>
  print()
registros |>
  filter(cnuc == "0000.00.0148") |>
  print()
