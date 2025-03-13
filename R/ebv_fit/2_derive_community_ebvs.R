# Derive community EBV's from IBA data ------------------------------------------------------------------
library(tidyverse)
library(lubridate)
library(vegan)
library(ade4)
library(adespatial)
library(data.table)

source("R/functions.R")
# read in data --------------------------------------------------------------------------------

# ------ Community diversity EBVs
# Cluster occurrences
tax_counts_se <- readRDS("data/tidydata/tax_counts_se.rds")
tax_counts_mg <- readRDS("data/tidydata/tax_counts_mg.rds")

# ------- Meta data
all_meta_se <- readRDS("data/tidydata/all_meta_se.rds")
all_meta_mg <- readRDS("data/tidydata/all_meta_mg.rds")



# Species richness ----------------------------------------------------------------------------
#Get weekly species richness measures

# Sweden
se_richness <-  tax_counts_se |> 
  mutate(week_year  = lubridate::week(collecting_date)) |> 
  group_by(trapID, week_year,latitude_WGS84 , longitude_WGS84) |> 
  summarise(nOTU = n_distinct(cluster)) |>
  ungroup()

# Madagascar
mg_richness <-  tax_counts_mg |> 
  mutate(week_year = lubridate::week(collecting_date)) |> 
  group_by(trapID, week_year ,latitude_WGS84 , longitude_WGS84) |> 
  summarise(nOTU = n_distinct(cluster)) |>
  ungroup()

# ------------ Checks ---------- #
ggplot(se_richness , aes(week_year, nOTU))+
  geom_point(alpha=.1)+
  geom_smooth()

ggplot(mg_richness , aes(week_year, nOTU))+
  geom_point(alpha=.1)+
  geom_smooth()
# ------------ Checks ---------- #


# LCBD ----------------------------------------------------------------------------------------

#------- Calculate LCBD for each site at each time slice.



# Sweden
tax_counts_wide_se <- tax_counts_se |>
  dplyr::select(cluster, trapID, read_count, collecting_date) |>
  mutate(
    week_year = week(collecting_date),
    pres = 1*(read_count>0)) |>
  dplyr::select(-collecting_date,-read_count) |>
  pivot_wider(names_from = cluster,
              values_from = pres,
              values_fn = mean,
              values_fill=0) |> 
  drop_na(trapID)

# Remove low frequency taxa - these can add a lot of noise to beta diversity measures
hiftx_cs_se <- colSums(tax_counts_wide_se[,c(-1,-2)])
hiftx_se <- hiftx_cs_se[hiftx_cs_se>=5]

tax_counts_wide_se <- tax_counts_wide_se |> _[, c("trapID" , "week_year" , names(hiftx_se))]


# Madagascar
tax_counts_wide_mg <- tax_counts_mg |>
  dplyr::select(cluster,  trapID, read_count, collecting_date) |>
  mutate(
    week_year = week(collecting_date),
    pres = 1*(read_count>0)) |>
  dplyr::select(-collecting_date,-read_count) |>
  pivot_wider(names_from = cluster,
              values_from = pres,
              values_fn = mean,
              values_fill=0) |>
  drop_na(trapID) 

# Remove low frequency taxa - these can add a lot of noise to beta diversity measures
hiftx_cs_mg <- colSums(tax_counts_wide_mg[,c(-1,-2)])
hiftx_mg <- hiftx_cs_mg[hiftx_cs_mg>=5]

tax_counts_wide_mg <- tax_counts_wide_mg |> _[, c("trapID" , "week_year" , names(hiftx_mg))]



# ---------  monthly species * site matrices
# Sweden
comm_M_se <- tax_counts_wide_se |>
  split(tax_counts_wide_se$week_year)  |>
  purrr::map(~column_to_rownames(.x,"trapID")) 

# Madagascar
comm_M_mg <- tax_counts_wide_mg |>
  split(tax_counts_wide_mg$week_year)  |>
  purrr::map(~column_to_rownames(.x,"trapID")) 


# ---------  remove weeks that have fewer than forty observations
comm_M_se <- comm_M_se[-which(sapply(comm_M_se,nrow) < 40)]
comm_M_mg <- comm_M_mg[-which(sapply(comm_M_mg,nrow) < 40)]


# Calculate lcbd
lcbd_week_se <- purrr::map(comm_M_se , ~calc_lcbd(.x)) |> bind_rows(.id="week_year") 
lcbd_week_mg <- purrr::map(comm_M_mg , ~calc_lcbd(.x)) |> bind_rows(.id="week_year")

#join with trap ID - remove LCBD values <0 caused by samples with only 1-2  species detections. 
lcbd_week_mg <- all_meta_mg |>
  dplyr::select(trapID , siteID) |>
  distinct() |> full_join(lcbd_week_mg) |>
  filter(lcbd>=0) |>
  mutate(week_year = as.numeric(week_year))

lcbd_week_se <- all_meta_se |> 
  dplyr::select(trapID , siteID) |>
  distinct() |> full_join(lcbd_week_se) |>
  filter(lcbd>=0) |>
  mutate(week_year = as.numeric(week_year))



# ---- Checks
ggplot(lcbd_week_se , aes(as.numeric(week_year) , lcbd))+
  geom_point(alpha=.1)+
  geom_smooth()

ggplot(lcbd_week_mg , aes(as.numeric(week_year) , lcbd))+
  geom_point(alpha=.1)+
  geom_smooth()


# save ----------------------------------------------------------------------------------------

commList <- list(sweden     = list(sr        = se_richness, lcbd = lcbd_week_se) , 
                 madagascar = list(sr        = mg_richness, lcbd = lcbd_week_mg))

saveRDS(commList , "data/tidydata/derived_ebvs/commList.rds")
