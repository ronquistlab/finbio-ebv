library(tidyverse)
library(lubridate)
library(xgboost)
library(patchwork)
library(ggh4x)
library(mgcv)
library(rnaturalearth)
library(rnaturalearthdata)
library(sf)
source("R/functions.R")

# prediction targets --------------------------------------------------------------------------

full_M <- readRDS("data/tidydata/full_M.rds")

# functions ---------------------------------------------------------------

# Plot ebv predictions - temporal
plot_preds_spatial <- function(model, ebv_df, country,
                               fill_var = NULL,    
                               week_year = 26,
                               sample_time = 7,
                               grid_res = 300,
                               quantiles = c(0, .95)) {
  
  
  trap_ID <- ebv_df$trap_ID[1]
  # Create a lat/lon grid over country's bounding box
  lon_seq <- seq(min(ebv_df$longitude_WGS84)-10, max(ebv_df$longitude_WGS84)+10, length.out = grid_res)
  lat_seq <- seq(min(ebv_df$latitude_WGS84)-10, max(ebv_df$latitude_WGS84)+10, length.out = grid_res)
  
  grid <- expand.grid(
    longitude_WGS84 = lon_seq,
    latitude_WGS84 = lat_seq
  )
  
  # Convert grid to sf and clip 
  grid_sf <- st_as_sf(grid, coords = c("longitude_WGS84", "latitude_WGS84"), crs = 4326)
  grid_in <- st_join(grid_sf, country, join = st_within) %>%  filter(!is.na(admin))  
  
  # Extract coords and drop geometry
  grid <- grid_in %>%
    mutate(
      longitude_WGS84 = st_coordinates(.)[, 1],
      latitude_WGS84 = st_coordinates(.)[, 2]
    ) %>%
    st_drop_geometry()
  
  # Add any required predictors
  grid$sample_time <- sample_time
  grid$week_year <- week_year
  grid$trap_ID <- trap_ID
  grid$ebv    <- fill_var
  # Predict
  grid$predicted <- predict(model, newdata = grid, type = "response" , exclude = "s(week_year)") 

  
  # Determine fill scale limits
  fill_limits <- quantile(ebv_df[fill_var], probs = quantiles, na.rm = TRUE)
  legend_breaks <- pretty(fill_limits, n = 3)
  
  
  # Get bounding box of the country (as numeric)
  bbox <- st_bbox(country)
  x_breaks <- pretty(c(bbox["xmin"], bbox["xmax"]), n = 4)
  
  
  # Plot with dynamic breaks
  gg <- ggplot() +
    geom_sf(data = country, fill = "grey90", color = "black") +
    geom_raster(data = grid, aes(x = longitude_WGS84, y = latitude_WGS84, fill = predicted), interpolate = TRUE) +
    scale_fill_viridis_c(name = "Prediction", option = "turbo", limits = fill_limits,breaks=legend_breaks) +
    labs(x = "Longitude", y = "Latitude") +
    scale_x_continuous(breaks = x_breaks) +
    theme_linedraw(base_size = 15) +
    theme(
      strip.text = element_text(size = 20), 
      legend.position = "bottom", 
      legend.title = element_blank(),
      legend.key.width = unit(1, "cm"),
      legend.text = element_text(size = 16),
      axis.text = element_text(size = 16),
      axis.title = element_text(size = 18)
    ) 
  
  return(gg)
}
# meta_data -----------------------------------------------------------------------------------

all_meta_se <- readRDS("data/tidydata/all_meta_se.rds") |>
  dplyr::select(matches("long|lati|collecting"),trap_ID=trapID) %>% 
  mutate(week_year = lubridate::week(collecting_date),trap_ID=factor(trap_ID))


all_meta_mg <- readRDS("data/tidydata/all_meta_mg.rds") |> 
  dplyr::select(matches("long|lati|collecting"),trap_ID=trapID) %>% 
  mutate(week_year = lubridate::week(collecting_date),trap_ID=factor(trap_ID))

# join ebvs to meta data --------------------------------------------------

ebvs_se <- lapply(full_M$sweden , function(x) select(x ,2, matches("trap|week_year|sample"))) 
ebvs_se <- lapply(ebvs_se , 
                  function(x){left_join(x,all_meta_se,relationship="many-to-many") %>%
                              mutate(trap_ID = factor(trap_ID)) %>% drop_na()})[-5]

# ----- Fit models
mod_sr_se   <- gam(nOTU ~ s(week_year,k=10,bs="cc") + s(longitude_WGS84, latitude_WGS84,k=20)+log(sample_time), data = ebvs_se$SR,  family = "nb" , method = "REML",select=TRUE)
mod_lcbd_se <- gam(lcbd ~ s(week_year,k=10,bs="cc") + s(longitude_WGS84, latitude_WGS84,k=20), data = ebvs_se$LCBD,  family = "betar" , method = "REML",select=TRUE)  
mod_FD_se   <- gam(FDis ~ s(week_year,k=10,bs="cc") + s(longitude_WGS84, latitude_WGS84,k=20) + log(sample_time), data = ebvs_se$FD,  family = "betar" , method = "REML",select=TRUE)  
mod_FE_se   <- gam(FEve ~ s(week_year,k=10,bs="cc") + s(longitude_WGS84, latitude_WGS84,k=20), data = ebvs_se$FE,  family = "betar" , method = "REML",select=TRUE)  
mod_GSH_se  <- gam(mean_shn ~ s(week_year,k=10,bs="cc") + s(longitude_WGS84, latitude_WGS84,k=20), data = ebvs_se$GSH,  family = "gaussian" , method = "REML",select=TRUE)  

# Stick in a list
modList_se <- list(mod_sr_se , mod_lcbd_se , mod_FD_se , mod_FE_se , mod_GSH_se)


# Madagascar --------------------------------------------------------------
ebvs_mg <- lapply(full_M$madagascar , function(x) select(x ,2, matches("trap|week_year|sample")))
ebvs_mg <- lapply(ebvs_mg , 
                  function(x){left_join(x, all_meta_mg, relationship = "many-to-many") %>%
                      mutate(trap_ID = factor(trap_ID)) %>% drop_na()})[-5]

# Fit models to each ebv ------------ 
families <- c("nb", "quasibinomial", "betar", "betar", "gaussian", "gaussian")

# ----- Get model formulas

# ----- Fit models
mod_sr_mg   <- gam(nOTU ~ s(week_year,k=10,bs="cc") + s(longitude_WGS84, latitude_WGS84,k=20)+log(sample_time), data = ebvs_mg$SR,  family = "nb", method = "REML", select = TRUE)
mod_lcbd_mg <- gam(lcbd ~ s(week_year,k=10,bs="cc") + s(longitude_WGS84, latitude_WGS84,k=20), data = ebvs_mg$LCBD,  family = "betar", method = "REML", select = TRUE)  
mod_FD_mg   <- gam(FDis ~ s(week_year,k=10,bs="cc") + s(longitude_WGS84, latitude_WGS84,k=20) + log(sample_time), data = ebvs_mg$FD,  family = "gaussian", method = "REML", select = TRUE)  
mod_FE_mg   <- gam(FEve ~ s(week_year,k=10,bs="cc") + s(longitude_WGS84, latitude_WGS84,k=20), data = ebvs_mg$FE,  family = "betar", method = "REML", select = TRUE)  
mod_GSH_mg  <- gam(mean_shn ~ s(week_year,k=10,bs="cc") + s(longitude_WGS84, latitude_WGS84,k=20), data = ebvs_mg$GSH,  family = "gaussian", method = "REML", select = TRUE)  

# Stick in a list
modList_mg <- list(mod_sr_mg, mod_lcbd_mg, mod_FD_mg, mod_FE_mg, mod_GSH_mg)


# ----------------------------------- 

# Make some predictions
# Temporal 
temp_preds_se <- lapply(1:length(modList_se) , function(x){
  
  newData <- expand.grid(
    week_year       = seq(1, 52, 1),
    sample_time     = 7,
    trap_ID         = levels(ebvs_se[[x]]$trap_ID),
    latitude_WGS84  = 60.0619,
    longitude_WGS84 = 15.87964
  )
  
  newData$preds  <- predict(modList_se[[x]], newdata = newData, type = 'response')
  newData$ebv <- names(ebvs_se[x])
  return(newData)  
    
})

# ----------------------------------- 
# Make some predictions
# Temporal 

temp_preds_mg <- lapply(1:length(modList_mg), function(x) {
  
  newData <- expand.grid(
    week_year       = seq(1, 52, 1),
    sample_time     = 7,
    trap_ID         = levels(ebvs_mg[[x]]$trap_ID),
    latitude_WGS84  = -18.33365,
    longitude_WGS84 = 47.34408
  )
  
  newData$preds  <- predict(modList_mg[[x]], newdata = newData, type = 'response')
  newData$ebv <- names(ebvs_mg[x])
  return(newData)  
})



# tidy and plot data -----------------------------------------------------

# Madagascar 
raw_data_mg <- ebvs_mg %>% 
  imap(~ .x %>%
         rename(value =1) %>%      # rename first column
         mutate(ebv= .y)) %>% 
  bind_rows()%>% 
  mutate(country = "Madagascar")

# Sweden
raw_data_se <- ebvs_se %>% 
  imap(~ .x %>%
         rename(value =1) %>%      # rename first column
         mutate(ebv= .y)) %>% 
  bind_rows() %>% 
  mutate(country = "Sweden")

# Get predictions in one df
temp_preds_all_mg <- temp_preds_mg %>% bind_rows() %>% mutate(country = "Madagascar")
temp_preds_all_se<- temp_preds_se %>% bind_rows() %>% mutate(country = "Sweden")

temp_preds_all <- bind_rows(temp_preds_all_mg , temp_preds_all_se)
raw_data_all <- bind_rows(raw_data_mg , raw_data_se)


temp_preds_sum <- temp_preds_all %>%
  summarise(preds_range = max(preds) - min(preds),.by = c("country","ebv"))

temp_preds_all$country <- factor(temp_preds_all$country , levels = c("Sweden" , "Madagascar"))
temp_preds_all$ebv <- factor(temp_preds_all$ebv , levels = c("SR" , "LCBD" , "FD","FE","GSH"))


raw_data_all$country <- factor(raw_data_all$country , levels = c("Sweden" , "Madagascar"))
raw_data_all$ebv <- factor(raw_data_all$ebv , levels = c("SR" , "LCBD" , "FD","FE","GSH"))

# Plot
p1 <- ggplot()+
  geom_point(data=raw_data_all , aes(week_year , value),alpha = 0.05,size=3)+
  geom_line(data=temp_preds_all , aes(x = week_year, y = preds), linewidth = 4, color = "blue") +
  theme_linedraw(base_size = 40)+
  theme()+
  facet_grid2(ebv~country,scales = "free_y",independent = "y")+
  labs(x="Week of the year" , y="EBV value")


library(patchwork)
tiff("plots/temporal_trends.tif" , width = 1200,height=1800,compression = "lzw")
p1
dev.off()

browseURL("plots/temporal_trends.tif")


# spatial trends ----------------------------------------------------------

sweden <- ne_countries(scale = "medium", country = "sweden", returnclass = "sf")
madagascar <- ne_countries(scale = "medium", country = "madagascar", returnclass = "sf")



# plot --------------------------------------------------------------------



# plot --------------------------------------------------------------------
ebv_names <- list("FDis" = "FD" , 
                  "FEve" = "FE" , 
                  "lcbd" = 'LCBD',
                  "nOTU" = "SR",
                  "mean_shn" = "GSH")


ebv_labeller <- function(variable,value){
  return(ebv_names[value])
}

p1_s_se <- plot_preds_spatial(mod_sr_se , ebvs_se$SR , country = sweden , fill_var = "nOTU") + facet_wrap(~ebv , labeller =  ebv_labeller)
p2_s_se <- plot_preds_spatial(mod_lcbd_se , ebvs_se$LCBD , country = sweden , fill_var = "lcbd")+ facet_wrap(~ebv , labeller =  ebv_labeller)
p3_s_se <- plot_preds_spatial(mod_FD_se , ebvs_se$FD , country = sweden , fill_var = "FDis")+ facet_wrap(~ebv , labeller =  ebv_labeller)
p4_s_se <- plot_preds_spatial(mod_FE_se , ebvs_se$FE , country = sweden , fill_var = "FEve")+ facet_wrap(~ebv , labeller =  ebv_labeller)
p5_s_se <- plot_preds_spatial(mod_GSH_se , ebvs_se$GSH , country = sweden , fill_var = "mean_shn")+ facet_wrap(~ebv , labeller =  ebv_labeller)

p1_s_mg <- plot_preds_spatial(mod_sr_mg , ebvs_mg$SR , country = madagascar , fill_var = "nOTU")+ facet_wrap(~ebv , labeller =  ebv_labeller)
p2_s_mg <- plot_preds_spatial(mod_lcbd_mg , ebvs_mg$LCBD , country = madagascar , fill_var = "lcbd")+ facet_wrap(~ebv , labeller =  ebv_labeller)
p3_s_mg <- plot_preds_spatial(mod_FD_mg , ebvs_mg$FD , country = madagascar , fill_var = "FDis")+ facet_wrap(~ebv , labeller =  ebv_labeller)
p4_s_mg <- plot_preds_spatial(mod_FE_mg , ebvs_mg$FE , country = madagascar , fill_var = "FEve")+ facet_wrap(~ebv , labeller =  ebv_labeller)
p5_s_mg <- plot_preds_spatial(mod_GSH_mg , ebvs_mg$GSH , country = madagascar , fill_var = "mean_shn")+ facet_wrap(~ebv , labeller =  ebv_labeller)

# plot --------------------------------------------------------------------
library(patchwork)


tiff("plots/spatial_trends.tif" , width = 2000,height=1000,compression = "lzw")
p1_s_se + p2_s_se + p3_s_se + p4_s_se + p5_s_se +
  p1_s_mg + p2_s_mg + p3_s_mg + p4_s_mg + p5_s_mg + 
  plot_layout(nrow = 2) 
dev.off()

browseURL("plots/spatial_trends.tif")

names(ebvs_se$SR)
