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

# ------ Community diversity EBVs
# Cluster occurrences
tax_counts_se <- readRDS("data/tidydata/tax_counts_se.rds")
tax_counts_mg <- readRDS("data/tidydata/tax_counts_mg.rds")


# ------ Genetic diversity EBV's
# ASV's in each cluster
asv_se <- fread("data/IBA_data/cleaned_noise_filtered_cluster_taxonomy_se.tsv") 
asv_mg <- fread("data/IBA_data/cleaned_noise_filtered_cluster_taxonomy_mg.tsv") 

# Reads per ASV
asv_counts_se <- fread("data/IBA_data/co1_asv_counts_se.tsv.gz")
asv_counts_mg <- fread("data/IBA_data/co1_asv_counts_mg.tsv.gz")


# Spike in IDs 
spike_ins_mg  <- fread("data/IBA_data/biological_spikes_taxonomy_mg.tsv")
spike_ins_se  <- fread("data/IBA_data/biological_spikes_taxonomy_se.tsv")


# EPA-NG annotations
tax_se <- fread("~/Downloads/asv_taxonomy_epang_se.tsv.gz") 
tax_mg <- fread("~/Downloads/asv_taxonomy_epang_mg.tsv.gz") 

# Chester tree
chester_tree <- read.tree("data/species_level_tree.nwk")

# ------- Meta data
all_meta_se <- readRDS("data/tidydata/all_meta_se.rds")
all_meta_mg <- readRDS("data/tidydata/all_meta_mg.rds")


# genetic diversity ---------------------------------------------------------------------------


# ----- Clean cluster and asv taxonomies
# Sweden
asv_se_clean <- asv_se |> 
  _[Phylum == "Arthropoda" ,] |> 
  _[!str_detect(Family , "unclassified*|_X*"),] |> 
  _[!str_detect(Order , "_X*"),] |> 
  dplyr::select(ASV_ID = ASV , cluster)

# Madagascar
asv_mg_clean <- asv_mg |> 
  _[Phylum == "Arthropoda" ,] |> 
  _[!str_detect(Family , "unclassified*|_X*"),] |> 
  _[!str_detect(Order , "_X*"),] |> 
  dplyr::select(ASV_ID = ASV , cluster)

# Cluster ids for spike ins
spike_cluster_se <- asv_se |> _[str_detect(Species , paste(spike_ins_se$Species,collapse="|")),] 
spike_cluster_mg <- asv_mg |> _[str_detect(Species , paste(spike_ins_mg$Species,collapse="|")),] 

# ------------------------------
# ------ Read number calibration
# ------------------------------

# ------ Get spike in ASV clusters
spikes_asv_se <- asv_se |> _[str_detect(Species , paste(spike_ins_se$Species,collapse="|")),] 
spikes_asv_mg <- asv_mg |> _[str_detect(Species , paste(spike_ins_se$Species,collapse="|")),] 

# ----- Filter ASV table to make calibration operation more efficient
# Sweden
asv_id    <- asv_counts_se$ASV_ID %in% asv_se_clean$ASV_ID # Only use ASV's found in the clean ASV data
sample_id <- all_meta_se$sampleID_NGI                      # Only use lysate samples that represent ecological communities 
asv_counts_sm_se <- asv_counts_se |>                       # Filter
  _[asv_id,colnames(asv_counts_se) %in% c("ASV_ID" ,sample_id) ,with=FALSE]


# Madagascar
asv_id    <- asv_counts_mg$ASV_ID %in% asv_mg_clean$ASV_ID # Only use ASV's found in the clean ASV data
sample_id <- all_meta_mg$sampleID_NGI                      # Only use lysate samples that represent ecological communities 
asv_counts_sm_mg <- asv_counts_mg |>                       # Filter
  _[asv_id,colnames(asv_counts_mg) %in% c("ASV_ID" ,sample_id) ,with=FALSE]

# ------ Get calibrated read numbers
# Sweden
scols_se <- grep("^P[0-9]+", colnames(asv_counts_sm_se), value = TRUE) # Sample IDs to sum read counts over
calib_counts_se <- calibrate_reads(asv_counts_sm_se , spikes_asv_se$ASV , sample_ids = scols_se)# Calibrate reads
setDT(calib_counts_se) ;set(calib_counts_se, j = "ASV_ID", value = asv_counts_sm_se$ASV_ID) # Reassign ASV_ID column

# Madagascar
scols_mg <- grep("^P[0-9]+", colnames(asv_counts_sm_mg), value = TRUE) # Sample IDs to sum read counts over
calib_counts_mg <- calibrate_reads(asv_counts_sm_mg , spikes_asv_mg$ASV , sample_ids = scols_mg) # Calibrate reads
setDT(calib_counts_mg) ; set(calib_counts_mg, j = "ASV_ID", value = asv_counts_sm_mg$ASV_ID) # Reassign ASV_ID column

# -----------------------------------------
# ------ Calculate shannon index as a measure of genetic diversity 
# -----------------------------------------

# Sweden
asv_counts_se_sm <- calib_counts_se |>   
  _[,RowSum := rowSums(.SD,na.rm = TRUE), .SDcols = scols_se] |>   # Calculate read counts per ASV
  _[,c("ASV_ID" , "RowSum")]                      # Keep ID and total read count

# Calculate shannon index per cluster
shn_ind_se <- asv_counts_se_sm |> 
  _[asv_se_clean , on = "ASV_ID", cluster := i.cluster] |>  # Join on cluster
  _[cluster %in% names(hiftx_se),] |> 
  _[,.(shn_ind = calc_shannon_index(RowSum)), by=cluster] |> 
  _[!str_detect(cluster , paste(spike_cluster_se$cluster,collapse="|")),] 

# Madagascar 
asv_counts_mg_sm <- calib_counts_mg |> 
  _[,RowSum := rowSums(.SD, na.rm = TRUE), .SDcols = scols_mg] |> 
  _[,c("ASV_ID" , "RowSum")]

# Calculate shannon index per cluster
shn_ind_mg <- asv_counts_mg_sm |> 
  _[asv_mg_clean , on = "ASV_ID", cluster := i.cluster] |> 
  _[cluster %in% names(hiftx_mg),] |> 
  _[,.(shn_ind = calc_shannon_index(RowSum)), by=cluster] |> 
  _[!str_detect(cluster , paste(spike_cluster_se$cluster,collapse="|")),] 



# ------ Cluster by sample lookups
tc_lu_se <- tax_counts_se |>  _[,c("cluster" , "sampleID_NGI")]
tc_lu_mg <- tax_counts_mg |>  _[,c("cluster" , "sampleID_NGI")]

# Join ASV's to cluster
sample_shn_ind_se <- shn_ind_se |> _[tc_lu_se , on = "cluster"]   
sample_shn_ind_mg <- shn_ind_mg |> _[tc_lu_mg , on = "cluster"]    

# ------ Average shannon index per sample 
mean_shn_se <- sample_shn_ind_se |>
  _[,.(mean_shn = mean(shn_ind , na.rm = TRUE)) , by = sampleID_NGI]

mean_shn_mg <- sample_shn_ind_mg |>
  _[,.(mean_shn = mean(shn_ind , na.rm = TRUE)) , by = sampleID_NGI]


# ------ Join with location data
g_shn_week_se <- all_meta_se |>
  dplyr::select(trapID , siteID,sampleID_NGI,collecting_date) |> distinct() |>
  full_join(mean_shn_se) |> drop_na(mean_shn) |> 
  mutate(week_year = week(collecting_date))

# ------ Join with location data
g_shn_week_mg <- all_meta_mg |>
  dplyr::select(trapID , siteID,sampleID_NGI,collecting_date) |> distinct() |>
  full_join(mean_shn_mg) |> drop_na(mean_shn) |> 
  mutate(week_year = week(collecting_date))


# ----- Checks
ggplot(g_shn_week_se , aes(week_year , mean_shn))+
  geom_point(alpha=.1)+
  geom_smooth()

ggplot(g_shn_week_mg , aes(week_year , mean_shn))+
  geom_point(alpha=.1)+
  geom_smooth()


# phylogenetic diversity ----------------------------------------------------------------------

# Data

# Clean & filter
sp_se <- tax_se |> 
  _[!Species %like% "unclassified"] |> 
  _[,c("Species")] |> 
  unique() |> 
  _[,Species := gsub(" " , "_" , Species)]

# Clean & filter
sp_mg <- tax_mg |> 
  _[!Species %like% "unclassified"] |> 
  _[,c("Species")] |> 
  unique() |> 
  _[,Species := gsub(" " , "_" , Species)]

# Tips to keep
tips_se <- chester_tree$tip.label[chester_tree$tip.label %in% sp_se$Species]
tips_mg <- chester_tree$tip.label[chester_tree$tip.label %in% sp_mg$Species]

# Pruned trees
tree_se <- keep.tip(chester_tree, tips_se)
tree_mg <- keep.tip(chester_tree, tips_mg)


##### ------ Phylogenetic diversity by sample

# lookups for cluster by species annotation
cluster_sp_se <- asv_se |>
  _[,c("Species" , "cluster")] |> 
  _[Species %in% str_replace(tips_se , "_" , " "),] |> unique()

cluster_sp_mg <- asv_mg |>
  _[,c("Species" , "cluster")] |> 
  _[Species %in% str_replace(tips_mg , "_" , " "),] |> unique()


# Site * species matrix for clusters with EPA-NG annotations
epang_M_se <- cluster_counts_se |>
  _[cluster_sp_se , on="cluster"] |> 
  _[,-"cluster"] |> 
  _[,Species := gsub(" " , "_" , Species)] |> 
  transpose(make.names="Species", keep.names = "sampleID_NGI") |> 
  column_to_rownames("sampleID_NGI")


epang_M_mg <- cluster_counts_mg |>
  _[cluster_sp_mg , on="cluster"] |> 
  _[,-"cluster"] |> 
  _[,Species := gsub(" " , "_" , Species)] |> 
  transpose(make.names="Species", keep.names = "sampleID_NGI")|> 
  column_to_rownames("sampleID_NGI")

pd_se <- pd(samp = epang_M_se, tree = tree_se , include.root=FALSE) |> rownames_to_column(var = "sampleID_NGI")
pd_mg <- pd(samp = epang_M_mg, tree = tree_se , include.root=FALSE) |> rownames_to_column(var = "sampleID_NGI")




pd_week_se <- all_meta_se |>
  dplyr::select(trapID , siteID,sampleID_NGI,collecting_date) |> distinct() |>
  full_join(pd_se) |> drop_na(PD) |> 
  mutate(week_year = week(collecting_date)) |> 
  drop_na(week_year)

pd_week_mg <- all_meta_mg |>
  dplyr::select(trapID , siteID,sampleID_NGI,collecting_date) |> distinct() |>
  full_join(pd_mg) |> drop_na(PD) |> 
  mutate(week_year = week(collecting_date)) |> 
  drop_na(week_year)



# ----- Checks
ggplot(pd_week_se , aes(week_year , PD))+
  geom_point(alpha=.1)+
  geom_smooth()

ggplot(pd_week_mg , aes(week_year , PD))+
  geom_point(alpha=.1)+
  geom_smooth()

# save ----------------------------------------------------------------------------------------


genList <- list(sweden     = list(g_shannon     = g_shn_week_se , pd = pd_week_se) , 
                madagascar = list(g_shannon     = g_shn_week_mg , pd = pd_week_mg) ) 


saveRDS(genList , "data/tidydata/derived_ebvs/genList.rds")

# bind all ebvs -------------------------------------------------------------------------------

commList <- readRDS("data/tidydata/derived_ebvs/commList.rds")
funList <- readRDS("data/tidydata/derived_ebvs/funList.rds")
genList <- readRDS("data/tidydata/derived_ebvs/genList.rds")


ebv_list <- list(sweden     = list(sr        = commList$sweden$sr, lcbd = commList$sweden$lcbd ,    # Community
                                   g_shannon = genList$sweden$g_shannon , pd = genList$sweden$pd   ,                    # Genetic
                                   fdisp     = funList$sweden$fdisp , feve = funList$sweden$feve) , # Traits
                 
                 madagascar = list(sr        = commList$madagascar$sr, lcbd = commList$madagascar$lcbd ,    # Community
                                   g_shannon = genList$madagascar$g_shannon , pd =genList$madagascar$pd ,                      # Genetic
                                   fdisp     = funList$madagascar$fdisp , feve = funList$madagascar$feve) ) # Traits

saveRDS(ebv_list , "data/tidydata/derived_ebvs/derived_ebvs.rds")
