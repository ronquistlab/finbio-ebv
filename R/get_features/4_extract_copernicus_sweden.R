# get copernicus habitat data for sweden and sweden ---------------------------------------

library(tidyverse)
library(janitor)
library(lubridate)
library(raster)
library(rnaturalearth)
library(rnaturalearthdata)
library(sf)
library(terra)
library(data.table)

# ---------------------------------------------------------------------------------------------
# Documentation and data for copernicus land cover data found at:
#  https://doi.org/10.2909/c6377c6e-76cc-4d03-8330-628a03693042
# https://land.copernicus.eu/global/sites/cgls.vito.be/files/products/CGLOPS1_PUM_LC100m-V3_I3.4.pdf


IBA_locs <-  read_delim("data/IBA_data/sites_metadata_se.tsv" , 
                                           locale=locale(encoding="UTF-16LE"),delim = "\t") |>
  dplyr::select(x=longitude_WGS84 , y=latitude_WGS84 ,trap_id=trapID)

# boundaries for sweden
sweden       <- ne_countries(country = "sweden" , scale = "large")

buffers <- c("1km" = 1e3 , "2km" = 2e3 , "5km" = 5e3)

# get prediction layers  - Sweden  -------------------------------------------------------------

# read in cover data
ftype        <- raster("data/raw_envdata/copernicus/PROBAV_LC100_global_v3.0.1_2018-conso_Forest-Type-layer_EPSG-4326.tif")
sweden_ftype <- raster::crop(x = ftype , y = sweden) %>% mask( . , sweden) 

# Extract data in 1,2,5km radii around each trap location 
# Recode factors back to trap_IDs & code forest type  
# Pg. 29 table 5 of data manual - E = evergreen , D = deciduous , N = needle leaf , B = broad leaf
trap_key <- IBA_locs$trap_id ; names(trap_key) <- 1:length(IBA_locs$trap_id)
type_key <- c( "0"  = "unknown_forest" ,"1" = "ENF_forest" , "2" = "EBF_forest"  , "3" = "DNF_forest" , "4" = "DBF_forest"  , "5" = "mixed_forest") 

# loop over buffers
ftype_list <- list()
for(i in seq_along(buffers)){
  
  # extract data & recode
  lc_ex <- raster::extract(sweden_ftype ,  y = IBA_locs[,c("x","y")] , buffer = buffers[i] , df = TRUE) %>% as.data.frame()
  colnames(lc_ex) <- c("trap_ID" , "forest_type")
  ftype_list[[i]] <- lc_ex %>% mutate(trap_ID   = recode_factor(factor(trap_ID) , !!!trap_key) , 
                                      forest_type = recode_factor(factor(forest_type) , !!!type_key)) 
  
}

names(ftype_list) <- names(buffers)

# forest cover ------------------------------------------------------------

# Get forest cover for sweden
fcover        <- raster("data/raw_envdata/copernicus//PROBAV_LC100_global_v3.0.1_2018-conso_Tree-CoverFraction-layer_EPSG-4326.tif")
sweden_fcover <- raster::crop(x = fcover , y = sweden) %>% mask( . , sweden) 

plot(sweden_fcover)

# loop over buffers
fcover_list <- list()
for(i in seq_along(buffers)){
  
  # Extract data & recode
  lc_ex <- raster::extract(sweden_fcover ,  y = IBA_locs[,c("x","y")] , buffer = buffers[i] , df = TRUE) %>% as.data.frame()
  colnames(lc_ex) <- c("trap_ID" , "forest_cover")
  lc_ex <- lc_ex %>% mutate(trap_ID   = recode_factor(factor(trap_ID) , !!!trap_key)) 
  
  # Get fractional cover by type
  fcover_list[[i]] <- cbind(lc_ex , forest_type = ftype_list[[i]]$forest_type) %>% tibble() %>% 
    rownames_to_column() %>% 
    pivot_wider(names_from = "forest_type" , values_from = "forest_cover" , values_fill = 0) %>% 
    dplyr::select(-`NA`)
  
}

# name
names(fcover_list) <- names(buffers)


# crop cover ------------------------------------------------------------

# read in cover data
crop_cover        <- raster("data/raw_envdata/copernicus//PROBAV_LC100_global_v3.0.1_2018-conso_Crops-CoverFraction-layer_EPSG-4326.tif")
sweden_crop_cover <- raster::crop(x = crop_cover , y = sweden) %>% mask( . , sweden) 

# loop over buffers
crop_cover_list <- list()
for(i in seq_along(buffers)){
  
  # Extract data & recode
  lc_ex <- raster::extract(sweden_crop_cover ,  y = IBA_locs[,c("x","y")] , buffer = buffers[i] , df = TRUE) %>% as.data.frame()
  colnames(lc_ex) <- c("trap_ID" , "crop_cover")
  crop_cover_list[[i]] <- lc_ex %>% mutate(trap_ID   = recode_factor(factor(trap_ID) , !!!trap_key)) 
  
}

# names
names(crop_cover_list) <- names(buffers)



# grassland ---------------------------------------------------------------

# read in cover data
grass_cover        <- raster("data/raw_envdata/copernicus//PROBAV_LC100_global_v3.0.1_2018-conso_Grass-CoverFraction-layer_EPSG-4326.tif")
sweden_grass_cover <- raster::crop(x = grass_cover , y = sweden) %>% mask( . , sweden) 

plot(grass_cover)

# loop over buffers
grass_cover_list <- list()
for(i in seq_along(buffers)){
  
  # Extract data & recode
  lc_ex <- raster::extract(sweden_grass_cover ,  y = IBA_locs[,c("x","y")] , buffer = buffers[i] , df = TRUE) %>% as.data.frame()
  colnames(lc_ex) <- c("trap_ID" , "grass_cover")
  grass_cover_list[[i]] <- lc_ex %>% mutate(trap_ID   = recode_factor(factor(trap_ID) , !!!trap_key)) 
  
}

names(grass_cover_list) <- names(buffers)


# shrubs ------------------------------------------------------------------

# read in cover data
shrub_cover        <- raster("data/raw_envdata/copernicus//PROBAV_LC100_global_v3.0.1_2018-conso_Shrub-CoverFraction-layer_EPSG-4326.tif")
sweden_shrub_cover <- raster::crop(x = shrub_cover , y = sweden) %>% mask( . , sweden) 

plot(sweden_shrub_cover)

# loop over buffers
shrub_cover_list <- list()
for(i in seq_along(buffers)){
  
  # Extract data & recode
  lc_ex <- raster::extract(sweden_shrub_cover ,  y = IBA_locs[,c("x","y")] , buffer = buffers[i] , df = TRUE) %>% as.data.frame()
  colnames(lc_ex) <- c("trap_ID" , "shrub_cover")
  shrub_cover_list[[i]] <- lc_ex %>% mutate(trap_ID   = recode_factor(factor(trap_ID) , !!!trap_key)) 
  
}

names(shrub_cover_list) <- names(buffers)


# permanent water --------------------------------------------------------------------

# read in cover data
water_cover        <- raster("data/raw_envdata/copernicus//PROBAV_LC100_global_v3.0.1_2018-conso_SeasonalWater-CoverFraction-layer_EPSG-4326.tif")
sweden_water_cover <- raster::crop(x = water_cover , y = sweden) %>% mask( . , sweden) 

plot(sweden_water_cover)

# loop over buffers
water_cover_list <- list()
for(i in seq_along(buffers)){
  
  # Extract data & recode
  lc_ex <- raster::extract(sweden_water_cover ,  y = IBA_locs[,c("x","y")] , buffer = buffers[i] , df = TRUE) %>% as.data.frame()
  colnames(lc_ex) <- c("trap_ID" , "water_cover")
  water_cover_list[[i]] <- lc_ex %>% mutate(trap_ID   = recode_factor(factor(trap_ID) , !!!trap_key)) 
  
}

names(water_cover_list) <- names(buffers)

# seasonal water --------------------------------------------------------------------

# read in cover data
water_s_cover        <- raster("data/raw_envdata/copernicus//PROBAV_LC100_global_v3.0.1_2018-conso_SeasonalWater-CoverFraction-layer_EPSG-4326.tif")
sweden_water_s_cover <- raster::crop(x = water_s_cover , y = sweden) %>% mask( . , sweden) 

plot(sweden_water_s_cover)

# loop over buffers
water_s_cover_list <- list()
for(i in seq_along(buffers)){
  
  # Extract data & recode
  lc_ex <- raster::extract(sweden_water_s_cover ,  y = IBA_locs[,c("x","y")] , buffer = buffers[i] , df = TRUE) %>% as.data.frame()
  colnames(lc_ex) <- c("trap_ID" , "water_s_cover")
  water_s_cover_list[[i]] <- lc_ex %>% mutate(trap_ID   = recode_factor(factor(trap_ID) , !!!trap_key)) 
  
}

names(water_s_cover_list) <- names(buffers)

# Urban --------------------------------------------------------------------

# read in cover data
urban_cover        <- raster("data/raw_envdata/copernicus//PROBAV_LC100_global_v3.0.1_2018-conso_BuiltUp-CoverFraction-layer_EPSG-4326.tif")
sweden_urban_cover  <- raster::crop(x = urban_cover  , y = sweden) %>% mask( . , sweden) 

plot(sweden_urban_cover)

# loop over buffers
urban_cover_list <- list()
for(i in seq_along(buffers)){
  
  # Extract data & recode
  lc_ex <- raster::extract(sweden_urban_cover ,  y = IBA_locs[,c("x","y")] , buffer = buffers[i] , df = TRUE) %>% as.data.frame()
  colnames(lc_ex) <- c("trap_ID" , "urban_cover")
  urban_cover_list[[i]] <- lc_ex %>% mutate(trap_ID   = recode_factor(factor(trap_ID) , !!!trap_key)) 
  
}

names(urban_cover_list) <- names(buffers)

# snow --------------------------------------------------------------------

# read in cover data
snow_cover        <- raster("data/raw_envdata/copernicus//PROBAV_LC100_global_v3.0.1_2018-conso_Snow-CoverFraction-layer_EPSG-4326.tif")
sweden_snow_cover  <- raster::crop(x = snow_cover  , y = sweden) %>% mask( . , sweden) 

plot(sweden_snow_cover)

# loop over buffers
snow_cover_list <- list()
for(i in seq_along(buffers)){
  
  # Extract data & recode
  lc_ex <- raster::extract(sweden_snow_cover ,  y = IBA_locs[,c("x","y")] , buffer = buffers[i] , df = TRUE) %>% as.data.frame()
  colnames(lc_ex) <- c("trap_ID" , "snow_cover")
  snow_cover_list[[i]] <- lc_ex %>% mutate(trap_ID   = recode_factor(factor(trap_ID) , !!!trap_key)) 
  
}

names(snow_cover_list) <- names(buffers)

# moss --------------------------------------------------------------------

# read in cover data
moss_cover        <- raster("data/raw_envdata/copernicus//PROBAV_LC100_global_v3.0.1_2018-conso_MossLichen-CoverFraction-layer_EPSG-4326.tif")
sweden_moss_cover  <- raster::crop(x = moss_cover  , y = sweden) %>% mask( . , sweden) 

plot(sweden_moss_cover)

# loop over buffers
moss_cover_list <- list()
for(i in seq_along(buffers)){
  
  # Extract data & recode
  lc_ex <- raster::extract(sweden_moss_cover ,  y = IBA_locs[,c("x","y")] , buffer = buffers[i] , df = TRUE) %>% as.data.frame()
  colnames(lc_ex) <- c("trap_ID" , "moss_cover")
  moss_cover_list[[i]] <- lc_ex %>% mutate(trap_ID   = recode_factor(factor(trap_ID) , !!!trap_key)) 
  
}

names(moss_cover_list) <- names(buffers)


# discrete ----------------------------------------------------------------


# read in cover data
discrete_class        <- raster("data/raw_envdata/copernicus//PROBAV_LC100_global_v3.0.1_2018-conso_Discrete-Classification-map_EPSG-4326.tif")
sweden_discrete_class  <- raster::crop(x = discrete_class  , y = sweden) %>% mask( . , sweden) 

plot(sweden_discrete_class)

# loop over buffers
discrete_class_list <- list()
for(i in seq_along(buffers)){
  
  # Extract data & recode
  lc_ex <- raster::extract(sweden_discrete_class ,  y = IBA_locs[,c("x","y")] , buffer = buffers[i] , df = TRUE) %>% as.data.frame()
  colnames(lc_ex) <- c("trap_ID" , "discrete_class")
  discrete_class_list[[i]] <- lc_ex %>% mutate(trap_ID   = recode_factor(factor(trap_ID) , !!!trap_key)) 
  
}

names(discrete_class_list) <- names(buffers)

# check covariance --------------------------------------------------------

all_cover <- cbind(water = water_cover_list[[1]]$water_cover , 
                   grass = grass_cover_list[[1]]$grass_cover , 
                   shrub = shrub_cover_list[[1]]$shrub_cover , 
                   crop = crop_cover_list[[1]]$crop_cover , 
                   urban = urban_cover_list[[1]]$urban_cover,
                   snow = snow_cover_list[[1]]$snow_cover,
                   moss = moss_cover_list[[1]]$moss_cover,
                   fcover_list[[1]])


all_sum <- all_cover %>% 
  group_by(trap_ID) %>% 
  dplyr::select(-rowname) %>% 
  summarise_all(list(sum)) %>% ungroup() %>% 
  pivot_longer(matches(c("ENF" , "unknown" , "mixed" , "DNF" , "EBF" , "DBF" , "NA")) , values_to = "forest") %>% 
  dplyr::select(-name)


ggpairs(all_sum , columns = 2:9)



# save --------------------------------------------------------------------


saveList <- list(km_1 = list(fcover_list$`1km` , crop_cover_list$`1km` , shrub_cover_list$`1km`  , 
                             grass_cover_list$`1km` ,  water_cover_list$`1km` , water_s_cover_list$`1km` , 
                             urban_cover_list$`1km` , snow_cover_list$`1km` ,  moss_cover_list$`1km` , 
                             discrete_class_list$`1km`) , 
                 
                 km_2 = list(fcover_list$`2km` , crop_cover_list$`2km` , shrub_cover_list$`2km`  , 
                             grass_cover_list$`2km` ,  water_cover_list$`2km` , water_s_cover_list$`2km` , 
                             urban_cover_list$`2km` , snow_cover_list$`2km` ,  moss_cover_list$`2km` , 
                             discrete_class_list$`2km`), 
                 
                 km_5 = list(fcover_list$`5km` , crop_cover_list$`5km` , shrub_cover_list$`5km`  , 
                             grass_cover_list$`5km` ,  water_cover_list$`5km` , water_s_cover_list$`5km` , 
                             urban_cover_list$`5km` , snow_cover_list$`5km` ,  moss_cover_list$`5km` , 
                             discrete_class_list$`5km`))



# save
saveRDS(saveList  , "data/tidy_data/sweden_copernicus_fractional_cover.rds")




# mess  -------------------------------------------------------------------
