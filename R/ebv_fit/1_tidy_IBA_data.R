# tidy IBA data -------------------------------------------------------------------------------
library(tidyverse)
library(data.table)

# import data ---------------------------------------------------------------------------------

# --------------- Spike in IDs 
spike_ins_mg  <- fread("data/IBA_data/biological_spikes_taxonomy_mg.tsv")
spike_ins_se  <- fread("data/IBA_data/biological_spikes_taxonomy_se.tsv")



#--------------- Cluster taxonomies
clusters_se <- fread("data/IBA_data/cleaned_noise_filtered_cluster_taxonomy_se.tsv") |>
                _[representative == 1 ,] |> 
                _[Phylum == "Arthropoda" ,] |> 
                _[!str_detect(Family , "unclassified*|_X*"),] |> 
                _[!str_detect(Order , "_X*"),] |> 
                _[!str_detect(Species , paste(spike_ins_se$Species,collapse="|")),] 


clusters_mg <- fread("data/IBA_data/cleaned_noise_filtered_cluster_taxonomy_mg.tsv") |>
                _[representative == 1 ,] |> 
                _[Phylum == "Arthropoda" ,] |> 
                _[!str_detect(Order , "_X*"),] |> 
                _[!str_detect(Species , paste(spike_ins_mg$Species,collapse="|")),] 


# ---------------- CLuster counts
counts_se <- fread("data/IBA_data/cleaned_noise_filtered_cluster_counts_se.tsv")
counts_mg <- fread("data/IBA_data/cleaned_noise_filtered_cluster_counts_mg.tsv")



# --------------- meta_data
se_seq_meta    <- fread("data/IBA_data/co1_sequencing_metadata_se.tsv") # Sequencing
mg_seq_meta    <- fread("data/IBA_data/co1_sequencing_metadata_mg.tsv")
se_sample_meta <- fread("data/IBA_data/samples_metadata_malaise_se.tsv")
mg_sample_meta <- fread("data/IBA_data/samples_metadata_malaise_mg.tsv")
se_site_meta   <- read_delim("data/IBA_data/sites_metadata_se.tsv" , locale=locale(encoding="UTF-16LE"),delim = "\t")
mg_site_meta   <- read_delim("data/IBA_data/sites_metadata_mg.tsv" , delim="\t")


# join metadata - sweden
all_meta_se <- full_join(se_seq_meta , se_sample_meta) |>
               full_join(se_site_meta) |> 
               filter(lab_sample_type == "sample" , str_detect(dataset,"lysate")) |> 
                dplyr::select(sampleID_NGI , trapID, placing_date, collecting_date, matches("lat|lon"),malaise_trap_type,siteID) |> 
                mutate(collecting_date = as.Date(collecting_date,format="%d/%m/%Y"),
                       placing_date    = as.Date(placing_date,format="%d/%m/%Y"))


# join metadata - madagascar
all_meta_mg <- full_join(mg_seq_meta , mg_sample_meta) |>
                full_join(mg_site_meta) |> 
                filter(lab_sample_type == "sample" , str_detect(dataset,"lysate")) |> 
                dplyr::select(sampleID_NGI , trapID, placing_date, collecting_date, matches("lat|lon"),malaise_trap_type,siteID) |> 
                mutate(collecting_date = as.Date(collecting_date,format="%d/%m/%Y"),
                       placing_date    = as.Date(placing_date,format="%d/%m/%Y"))



# tidy swedish data ---------------------------------------------------------------------------
counts_long_se <- counts_se |> 
                  _[cluster %in% clusters_se$cluster,] |> 
                  melt(id.vars = 1 , variable.name = "sampleID_NGI" , value.name = "read_count") |> 
                  _[read_count>0,] 



# Merge taxonomy with cluster IDS
tax_counts_se <- counts_long_se[clusters_se, on = "cluster"] |> 
                _[!str_detect(Order , "unclassified"),] |> 
                right_join(all_meta_se) |> 
                drop_na(cluster) |> 
                  _[collecting_date < '2020-01-01']
 

# tidy malagasy data --------------------------------------------------------------------------
counts_long_mg <- counts_mg |> 
                  _[cluster %in% clusters_mg$cluster,] |> 
                  melt(id.vars = 1 , variable.name = "sampleID_NGI" , value.name = "read_count") |> 
                  _[read_count>0,]


# Merge taxonomy with cluster IDS
tax_counts_mg <- counts_long_mg[clusters_mg, on = "cluster"] |> 
                 _[!str_detect(Order , "unclassified"),] |> 
                 right_join(all_meta_mg) |> 
                 drop_na(cluster) |> 
                 _[collecting_date < '2021-01-01']

# save data -----------------------------------------------------------------------------------

# site metadata + OTU observations
saveRDS(tax_counts_se , "data/tidydata/tax_counts_se.rds")
saveRDS(tax_counts_mg , "data/tidydata/tax_counts_mg.rds")

#meta data - useful for joins later
saveRDS(all_meta_se , "data/tidydata/all_meta_se.rds")
saveRDS(all_meta_mg , "data/tidydata/all_meta_mg.rds")
