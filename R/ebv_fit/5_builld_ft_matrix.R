# build data to feed to models ----------------------------------------------------------------
library(tidyverse)
library(lubridate)
library(raster)
library(GGally)
library(geosphere)

# prediction targets --------------------------------------------------------------------------

ebv_targets <- readRDS("data/tidydata/derived_ebvs/derived_ebvs.rds")


# meta_data -----------------------------------------------------------------------------------

all_meta_se <- readRDS("data/tidydata/all_meta_se.rds") |>
                dplyr::select(matches("long|lati|collecting_date|placing_date"),trap_ID = trapID) |>
                drop_na() |> distinct()  |> 
                mutate(sample_time = as.numeric(interval(min(placing_date) , max(collecting_date)))/(3600*24))
                
all_meta_mg <- readRDS("data/tidydata/all_meta_mg.rds") |> 
                dplyr::select(matches("long|lati|collecting_date|placing_date"),trap_ID = trapID) |>
                drop_na() |> distinct() |> 
                mutate(sample_time = as.numeric(interval(min(placing_date) , max(collecting_date)))/(3600*24))

# features ------------------------------------------------------------------------------------

# copernicus fractional cover  -------------------- #
fc_mg <- readRDS("data/tidydata/env_features/madagascar_copernicus_fractional_cover.rds")
fc_se <- readRDS("data/tidydata/env_features/sweden_copernicus_fractional_cover.rds")

# Extract 1km rasters
hab_cover_mg <- do.call(bind_cols, lapply(fc_mg$km_1, `[`, -1)) |> dplyr::select(-discrete_class)
hab_cover_se <- do.call(bind_cols, lapply(fc_se$km_1, `[`, -1)) |> dplyr::select(-discrete_class)
# copernicus fractional cover  -------------------- #


# climate -------------------- #
era5_se <- readRDS("data/tidydata/env_features/weekly_era5_se.rds") |> rename(week_year = week)
era5_mg <- readRDS("data/tidydata/env_features/weekly_era5_mg.rds") |> rename(week_year = week)
# climate -------------------- #


# Photoperiod & sampling time -------------------- #
photo_se <-  all_meta_se |> 
             mutate(day = yday(collecting_date), week_year=lubridate::week(collecting_date)) |> 
             mutate(photoperiod = daylength(latitude_WGS84 , yday(collecting_date))) |> 
             summarise(mean_photoperiod = mean(photoperiod) ,
                       sample_time = as.numeric(interval(min(placing_date) , max(collecting_date)))/(3600*24),
                       .by = c("week_year","trap_ID"))

photo_mg <-  all_meta_mg |> 
              mutate(day = yday(collecting_date), week_year=lubridate::week(collecting_date)) |> 
              mutate(photoperiod = daylength(latitude_WGS84 , yday(collecting_date))) |> 
              summarise(mean_photoperiod = mean(photoperiod) ,
                        sample_time = as.numeric(interval(min(placing_date) , max(collecting_date)))/(3600*24),
                        .by = c("week_year","trap_ID"))


# habitat characteristics ---------------------------------------------------------------------

# Decide which habitat features to include

# Swedish environmental data
all_sum_se <-hab_cover_se %>% 
  group_by(trap_ID) %>% 
  summarise_all(list(mean)) %>% ungroup() 

# Histograms
all_sum_se |>
  pivot_longer(-trap_ID) |> 
  ggplot(aes(value))+
  geom_histogram()+
  facet_wrap(~name)

# Malagasy environmental data
all_sum_mg <- hab_cover_mg %>% 
  group_by(trap_ID) %>% 
  summarise_all(list(mean)) %>% ungroup() 

# Histograms
all_sum_mg |>
  pivot_longer(-trap_ID) |> 
  ggplot(aes(value))+
  geom_histogram()+
  facet_wrap(~name)


# ---------------------------------------------------------------------------------------------

# Select useful habitat data --- Sweden 
hab_f_se <- hab_cover_se |> 
                dplyr::select(trap_ID,matches("forest|shrub|grass|crop|urban|water_cover|moss_cover")) |> 
                group_by(trap_ID) |> 
                summarise_all(list(mean)) |> ungroup()


# Select useful habitat data --- Madagascar
hab_f_mg <- hab_cover_mg |>
              dplyr::select(trap_ID,matches("forest|shrub|grass|crop")) |> 
              group_by(trap_ID) |> 
              summarise_all(list(mean)) |> ungroup()


# join features -------------------------------------------------------------------------------

features_se <- era5_se |> full_join(hab_f_se) |> left_join(photo_se)
features_mg <- era5_mg |> full_join(hab_f_mg) |> left_join(photo_mg)


# Join features to targets --------------------------------------------------------------------

# ---------------------------------------------------------------------------------------------
# ------- EBV's - Sweden
# ---------------------------------------------------------------------------------------------



# --- Community EBVS
# ---------------------------------------------------------------------------------------------


# Species richness
sr_M_se <- ebv_targets$sweden$sr |>
              dplyr::select(trap_ID = trapID , nOTU , week_year) |> 
              left_join(features_se) |> 
              arrange(trap_ID) |> 
              mutate(wk_sin = sin((2*pi*week_year)/52),
                     wk_cos = cos((2*pi*week_year)/52))

# LCBD
lcbd_M_se <- ebv_targets$sweden$lcbd |> 
              dplyr::select(trap_ID = trapID , lcbd , week_year) |> 
              left_join(features_se) |> 
              arrange(trap_ID) |> 
              mutate(wk_sin = sin((2*pi*week_year)/52),
                     wk_cos = cos((2*pi*week_year)/52))

# --- Trait diversity 
# ---------------------------------------------------------------------------------------------


# Dispersal
fdisp_M_se <- ebv_targets$sweden$fdisp|> 
              dplyr::select(trap_ID = trapID , FDis, week_year) |> 
              left_join(features_se) |> 
              arrange(trap_ID) |> 
              mutate(wk_sin = sin((2*pi*week_year)/52),
                     wk_cos = cos((2*pi*week_year)/52))


# Eveness
feve_M_se <- ebv_targets$sweden$feve|> 
              dplyr::select(trap_ID = trapID , FEve, week_year) |> 
              left_join(features_se) |> 
              arrange(trap_ID) |> 
              mutate(wk_sin = sin((2*pi*week_year)/52),
                     wk_cos = cos((2*pi*week_year)/52)) 

# --- Genetic EBVS
# ---------------------------------------------------------------------------------------------

# Shannon ASV index
shn_M_se <- ebv_targets$sweden$g_shannon |>
            dplyr::select(trap_ID = trapID , mean_shn, week_year) |> 
            left_join(features_se) |> 
            arrange(trap_ID) |> 
            mutate(wk_sin = sin((2*pi*week_year)/52),
                   wk_cos = cos((2*pi*week_year)/52))


# Phylogenetic diversity
pd_M_se <- ebv_targets$sweden$pd |>
  dplyr::select(trap_ID = trapID , PD, week_year) |> 
  left_join(features_se) |> 
  arrange(trap_ID) |> 
  mutate(wk_sin = sin((2*pi*week_year)/52),
         wk_cos = cos((2*pi*week_year)/52))




# ---------------------------------------------------------------------------------------------
# ------- EBV's - Madagascar
# ---------------------------------------------------------------------------------------------



# --- Community EBVS
# ---------------------------------------------------------------------------------------------

# Species richness
sr_M_mg <- ebv_targets$madagascar$sr |>
            dplyr::select(trap_ID = trapID , nOTU , week_year) |> 
            left_join(features_mg) |> 
            arrange(trap_ID) |> 
            mutate(wk_sin = sin((2*pi*week_year)/52),
                   wk_cos = cos((2*pi*week_year)/52))

# LCBD
lcbd_M_mg <- ebv_targets$madagascar$lcbd |> 
              dplyr::select(trap_ID = trapID , lcbd , week_year) |> 
              left_join(features_mg) |> 
              arrange(trap_ID) |> 
              mutate(wk_sin = sin((2*pi*week_year)/52),
                     wk_cos = cos((2*pi*week_year)/52))

# --- Trait diversity 
# ---------------------------------------------------------------------------------------------

# Dispersal
fdisp_M_mg <- ebv_targets$madagascar$fdisp|> 
                dplyr::select(trap_ID = trapID , FDis , week_year) |> 
                left_join(features_mg) |> 
                arrange(trap_ID) |> 
                mutate(wk_sin = sin((2*pi*week_year)/52),
                       wk_cos = cos((2*pi*week_year)/52))



# Eveness
feve_M_mg <- ebv_targets$madagascar$feve |> 
              dplyr::select(trap_ID = trapID , FEve , week_year) |> 
              left_join(features_mg) |> 
              arrange(trap_ID) |> 
              mutate(wk_sin = sin((2*pi*week_year)/52),
                     wk_cos = cos((2*pi*week_year)/52))

# --- Genetic EBVS
# ---------------------------------------------------------------------------------------------

# Shannon ASV index
shn_M_mg <- ebv_targets$madagascar$g_shannon |>
            dplyr::select(trap_ID = trapID , mean_shn , week_year) |> 
            left_join(features_mg) |> 
            arrange(trap_ID) |> 
            mutate(wk_sin = sin((2*pi*week_year)/52),
                   wk_cos = cos((2*pi*week_year)/52))

# phylogenetic diversity
pd_M_mg <- ebv_targets$madagascar$pd |>
  dplyr::select(trap_ID = trapID , PD, week_year) |> 
  left_join(features_mg) |> 
  arrange(trap_ID) |> 
  mutate(wk_sin = sin((2*pi*week_year)/52),
         wk_cos = cos((2*pi*week_year)/52))



# ---------------------------------------------------------------------------------------------
# -----  Full feature target lists
# ---------------------------------------------------------------------------------------------

# Sweden
full_M_se <- list(SR=sr_M_se , LCBD = lcbd_M_se , 
                  FD = fdisp_M_se , FE = feve_M_se , 
                  PD = pd_M_se, GSH = shn_M_se) |> 
             purrr::map(~as.data.frame(.x))


# Madagascar
full_M_mg <- list(SR=sr_M_mg , LCBD = lcbd_M_mg , 
                  FD = fdisp_M_mg , FE = feve_M_mg , 
                  PD = pd_M_mg,GSH = shn_M_mg) |> 
           purrr::map(~as.data.frame(.x))

# --- Quick checks 

# Check column number - all should be the same for each country
lapply(full_M_se, ncol) 
lapply(full_M_mg, ncol)

# Check order of columns - First column should be trap ID, second column should be prediction target
lapply(full_M_se, colnames)
lapply(full_M_mg, colnames)


# save --------------------------------------------------------------------------------------

saveRDS(list(sweden=full_M_se, madagascar=full_M_mg),"data/tidydata/full_M.rds")
saveRDS(list(sweden=full_M_se, madagascar=full_M_mg),"R/ebv_fit/7_dardel/full_M.rds")

