# Derive EBV's from IBA data ------------------------------------------------------------------
library(tidyverse)
library(lubridate)
library(vegan)
library(ade4)
library(adespatial)
library(data.table)
library(fundiversity)
library(cluster)
library(ape)
library(picante)

source("R/functions.R")
# read in data --------------------------------------------------------------------------------

# ------- Trait data
trait_M_se <- readRDS("data/IBA_data/all_traits_se.rds") 
trait_M_mg <- readRDS("data/IBA_data/all_traits_mg.rds") 

# ------ Community data
# ASV's in each cluster
asv_se <- fread("data/IBA_data/cleaned_noise_filtered_cluster_taxonomy_se.tsv") 
asv_mg <- fread("data/IBA_data/cleaned_noise_filtered_cluster_taxonomy_mg.tsv") 

# Reads per cluster 
cluster_counts_se <- fread("data/IBA_data/cleaned_noise_filtered_cluster_counts_se.tsv")
cluster_counts_mg <- fread("data/IBA_data/cleaned_noise_filtered_cluster_counts_mg.tsv")

# Spike in IDs 
spike_ins_mg  <- fread("data/IBA_data/biological_spikes_taxonomy_mg.tsv")
spike_ins_se  <- fread("data/IBA_data/biological_spikes_taxonomy_se.tsv")


# ------- Meta data
all_meta_se <- readRDS("data/tidydata/all_meta_se.rds")
all_meta_mg <- readRDS("data/tidydata/all_meta_mg.rds")

# functional diversity ------------------------------------------------------------------------

# ------------------------------
# ----- Family by trait matrices 
# ------------------------------

# Calculate pairwise distances per family then join to each species. 
# Swedish Traits
trait_M_se <- trait_M_se |> 
  mutate(feeding_niche = 
           case_when(str_detect(feeding_niche , "parasitoid") ~ "Parasitoid" , 
                     TRUE ~ feeding_niche)) |>
  distinct() |> 
  # Exclude families where feeding niche traits are resolved to subfamily rather than family
  filter(!family %in%  c("Staphylinidae","Braconidae","Ichneumonidae") )|> 
  column_to_rownames("family") |> 
  mutate(across(where(is.character) , as.factor))


# Malagasy
trait_M_mg<- trait_M_mg  |> 
  mutate(feeding_niche = 
           case_when(str_detect(feeding_niche , "parasitoid") ~ "Parasitoid" , 
                     TRUE ~ feeding_niche)) |>
  distinct() |> 
  # Exclude families where feeding niche traits are resolved to subfamily rather than family
  filter(!family %in%  c("Staphylinidae","Braconidae","Ichneumonidae") )|> 
  column_to_rownames("family") |> 
  mutate(across(where(is.character) , as.factor))


# ------------------------------
# ----- Project mixed trait data into continuous space
# ------------------------------

# Take the gower dissimilarity across all 3 functional traits
# Project into continuous space using PCOA

# Sweden
gower_se <- cluster::daisy(trait_M_se , metric = "gower") |> as.matrix()
pcoa_se <- ape::pcoa(gower_se , rn = attr(gower_se , which = "Labels")) # Use PCOA to summarise traits & take first 3 axes

# Madagascar
gower_mg <- cluster::daisy(trait_M_mg , metric = "gower") |> as.matrix()
pcoa_mg  <- ape::pcoa(gower_mg , rn = attr(gower_mg , which = "Labels")) # Use PCOA to summarise traits & take first 3 axes

#  ------ Checks
plot_pcoa(pcoa_obj=pcoa_se)
plot_pcoa(pcoa_obj=pcoa_mg)


# ------ Join to species data 
pc_axes_se <- data.frame(p1 = pcoa_se$vectors[,1], p2 = pcoa_se$vectors[,2], p3 = pcoa_se$vectors[,3]) |> 
  rownames_to_column("family")

pc_axes_mg <- data.frame(p1 = pcoa_mg$vectors[,1], p2 = pcoa_mg$vectors[,2], p3 = pcoa_mg$vectors[,3]) |> 
  rownames_to_column("family")

# ------------------------------
# ------- Matrices for functional diversity metric calculation 
# ------------------------------

# Trait by species 
# Swedish trait scores for each cluster
trait_scores_se <- asv_se |> 
  select(cluster , family = Family) |> 
  distinct() |> 
  _[pc_axes_se,on="family"] |> 
  _[,c("cluster","p1","p2","p3")] |> 
  drop_na() |> 
  arrange(cluster) |> 
  column_to_rownames("cluster") |> 
  as.matrix()

# Malagasy trait scores for each cluster
trait_scores_mg <- asv_mg |> 
  select(cluster , family = Family) |> 
  distinct() |> 
  _[pc_axes_mg,on="family"] |> 
  _[,c("cluster","p1","p2","p3")] |> 
  drop_na() |> 
  arrange(cluster) |> 
  column_to_rownames("cluster") |> 
  as.matrix()

# Sample by species
# Sweden 
scols <- grep("^P[0-9]+", colnames(cluster_counts_se), value = TRUE) # Sample IDs to sum read counts over
site_sample_M_se  <- cluster_counts_se |> 
  _[cluster %in% rownames(trait_scores_se)] |> 
  _[,(scols) := lapply(.SD, function(x) 1*(x > 0)),.SDcols=scols] |> 
  transpose(keep.names = "sampleID_NGI", make.names = "cluster") |> 
  column_to_rownames("sampleID_NGI") |> 
  as.matrix()

# Madagascar                  
scols <- grep("^P[0-9]+", colnames(cluster_counts_mg), value = TRUE) # Sample IDs to sum read counts over
site_sample_M_mg  <- cluster_counts_mg |> 
  _[cluster %in% rownames(trait_scores_mg)] |> 
  _[,(scols) := lapply(.SD, function(x) 1*(x > 0)),.SDcols=scols] |> 
  transpose(keep.names = "sampleID_NGI", make.names = "cluster") |> 
  column_to_rownames("sampleID_NGI") |> 
  as.matrix()

# ------------------------------
# ----- calculate functional diversity metrics
# ------------------------------

# Dispersion
options(future.globals.maxSize = 5.0 * 2e9) 
fun_disp_se <- fundiversity::fd_fdis(trait_scores_se , site_sample_M_se) |> select(sampleID_NGI = site , FDis)
fun_disp_mg <- fundiversity::fd_fdis(trait_scores_mg , site_sample_M_mg) |> select(sampleID_NGI = site , FDis)

# Eveness
fun_eve_se <- fundiversity::fd_feve(trait_scores_se, site_sample_M_se) |> select(sampleID_NGI = site , FEve)
fun_eve_mg <- fundiversity::fd_feve(trait_scores_mg, site_sample_M_mg) |> select(sampleID_NGI = site , FEve)


# ----  Bind with trap ID and weeks
fdisp_week_se <- all_meta_se |>
  select(sampleID_NGI, trapID, collecting_date) |> distinct() |> 
  mutate(week_year = week(collecting_date)) |> 
  right_join(fun_disp_se)|> drop_na(trapID)

fdisp_week_mg <- all_meta_mg |>
  select(sampleID_NGI, trapID, collecting_date) |> distinct() |> 
  mutate(week_year = week(collecting_date)) |> 
  right_join(fun_disp_mg)|> drop_na(trapID)


feve_week_se <- all_meta_se |>
  select(sampleID_NGI, trapID, collecting_date) |> distinct() |> 
  mutate(week_year = week(collecting_date)) |> 
  right_join(fun_eve_se) |> drop_na(trapID)

feve_week_mg <- all_meta_mg |>
  select(sampleID_NGI, trapID, collecting_date) |> distinct() |> 
  mutate(week_year = week(collecting_date)) |> 
  left_join(fun_eve_mg)|> drop_na(trapID)


# ----- Checks

ggplot(feve_week_se,aes(week_year,FEve))+
  geom_point()+
  geom_smooth()

ggplot(fdisp_week_se,aes(week_year,FDis))+
  geom_point()+
  geom_smooth()

ggplot(feve_week_mg,aes(week_year,FEve))+
  geom_point()+
  geom_smooth()

ggplot(fdisp_week_mg,aes(week_year,FDis))+
  geom_point()+
  geom_smooth()



# save ----------------------------------------------------------------------------------------


funList <- list(sweden     = list(fdisp     = fdisp_week_se , feve = feve_week_se) , 
                 madagascar = list(fdisp     = fdisp_week_mg , feve = feve_week_mg) ) 


saveRDS(funList , "data/tidydata/derived_ebvs/funList.rds")
