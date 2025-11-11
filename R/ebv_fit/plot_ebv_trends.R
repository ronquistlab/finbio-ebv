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
plot_preds_temp <- function(temp_preds, ebvs, response_var_name) {
  
  pList <- list()
  for (i in seq_along(ebvs)) {
    
    # Extract predictions and week-year
    df <- data.frame(
      preds = temp_preds[, i],
      week_year = seq(min(ebvs[[1]]$week_year), 52, 1),
      trap_ID = levels(ebvs[[i]]$trap_ID)
    ) %>%
      summarise(preds = mean(preds), .by = week_year) %>%
      mutate(ebv = response_var_name[i])
    
    # Extract the dataset containing the response
    response_data <- ebvs[[i]]
    
    
    # Create the plot
    pList[[i]] <- ggplot(response_data, aes(x = week_year, y = .data[[response_var_name[i]]])) +
      geom_point(alpha = 0.05,size=2) +
      geom_line(data = df, aes(x = week_year, y = preds), linewidth = 4, color = "blue") +
      theme_linedraw(base_size = 35)+
      theme(
        strip.text = element_text(size = 20))
    
  }
  return(pList)
  
}

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
  
  # Convert grid to sf and clip to Sweden
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
      legend.key.width = unit(1, "cm")
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

ebvs_se <- lapply(full_M$sweden , function(x) select(x ,2, matches("trap|week_year|sample"))) %>% 
ebvs_se <- lapply(ebvs_se , 
                  function(x){left_join(x,all_meta_se,relationship="many-to-many") %>%
                              mutate(trap_ID = factor(trap_ID)) %>% drop_na()})[-5]

# Fit models to each ebv ------------ 
families <- c("nb","quasibinomial","betar","betar","gaussian","gaussian")

# ----- Get model formulas
formulas <- lapply(1:length(ebvs_se), function(x) {
  
  ebv_df <- ebvs_se[[x]]
  response <- names(ebv_df)[1]              # Get name of target EBV
  ebv_df <- ebv_df  %>% mutate(trap_ID = factor(trap_ID))
  
  formula <- as.formula(paste(response, "~", 
                              paste('s(week_year,k=8,bs="cc") + 
                                             s(longitude_WGS84, latitude_WGS84, k = 20)+
                                             log(sample_time)')))
          })

# ----- Fit models
mod_sr_se   <- gam(nOTU ~ s(week_year,k=8,bs="cc") + s(longitude_WGS84, latitude_WGS84,k=20)+log(sample_time), data = ebvs_se$SR,  family = "nb" , method = "REML",select=TRUE)
mod_lcbd_se <- gam(lcbd ~ s(week_year,k=8,bs="cc") + s(longitude_WGS84, latitude_WGS84,k=20), data = ebvs_se$LCBD,  family = "betar" , method = "REML",select=TRUE)  
mod_FD_se   <- gam(FDis ~ s(week_year,k=8,bs="cc") + s(longitude_WGS84, latitude_WGS84,k=20) + log(sample_time), data = ebvs_se$FD,  family = "betar" , method = "REML",select=TRUE)  
mod_FE_se   <- gam(FEve ~ s(week_year,k=8,bs="cc") + s(longitude_WGS84, latitude_WGS84,k=20), data = ebvs_se$FE,  family = "betar" , method = "REML",select=TRUE)  
mod_GSH_se  <- gam(mean_shn ~ s(week_year,k=8,bs="cc") + s(longitude_WGS84, latitude_WGS84,k=20), data = ebvs_se$GSH,  family = "gaussian" , method = "REML",select=TRUE)  

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
mod_sr_mg   <- gam(nOTU ~ s(week_year,k=8,bs="cc") + s(longitude_WGS84, latitude_WGS84,k=20)+log(sample_time), data = ebvs_mg$SR,  family = "nb", method = "REML", select = TRUE)
mod_lcbd_mg <- gam(lcbd ~ s(week_year,k=8,bs="cc") + s(longitude_WGS84, latitude_WGS84,k=20), data = ebvs_mg$LCBD,  family = "betar", method = "REML", select = TRUE)  
mod_FD_mg   <- gam(FDis ~ s(week_year,k=8,bs="cc") + s(longitude_WGS84, latitude_WGS84,k=20) + log(sample_time), data = ebvs_mg$FD,  family = "gaussian", method = "REML", select = TRUE)  
mod_FE_mg   <- gam(FEve ~ s(week_year,k=8,bs="cc") + s(longitude_WGS84, latitude_WGS84,k=20), data = ebvs_mg$FE,  family = "betar", method = "REML", select = TRUE)  
mod_GSH_mg  <- gam(mean_shn ~ s(week_year,k=8,bs="cc") + s(longitude_WGS84, latitude_WGS84,k=20), data = ebvs_mg$GSH,  family = "gaussian", method = "REML", select = TRUE)  

# Stick in a list
modList_mg <- list(mod_sr_mg, mod_lcbd_mg, mod_FD_mg, mod_FE_mg, mod_GSH_mg)


# ----------------------------------- 

# Make some predictions
# Temporal 
temp_preds_se <- sapply(1:length(modList_se) , function(x){
  
  newData <- expand.grid(week_year      = seq(5,52,1) , 
                        sample_time     = 7,
                        trap_ID         = levels(ebvs_se[[x]]$trap_ID) , 
                        latitude_WGS84  = 60.0619 ,
                        longitude_WGS84 = 15.87964 )
  
  preds   <- predict(modList_se[[x]] , newdata = newData , type='response') 
    
})

# ----------------------------------- 
# Plot data 

ebv_names <- list(
  'FDis'="Fun Dispersal",
  'FEve'="Fun Eveness",
  'lcbd'="LCBD",
  'nOTU'="Species richness",
  'mean_shn' = "Genetic diversity"
)

ebv_labeller <- function(variable,value){
  return(ebv_names[value])
}

ebv_plots_se <- plot_preds_temp(temp_preds_se , ebvs_se , c("nOTU","lcbd","FDis","FEve","mean_shn"))

p1_se <- ebv_plots_se[[1]] + labs(x = "Week of the year" , y = "Species richness") + facet_wrap(~ebv , labeller =  ebv_labeller)
p2_se <-ebv_plots_se[[2]]  + labs(x = "Week of the year" , y = "LCBD")+ facet_wrap(~ebv , labeller =  ebv_labeller)
p3_se <-ebv_plots_se[[3]]  + labs(x = "Week of the year" , y = "Functional dispersal")+ facet_wrap(~ebv , labeller =  ebv_labeller)
p4_se <-ebv_plots_se[[4]]  + labs(x = "Week of the year" , y = "Functional eveness")+ facet_wrap(~ebv , labeller =  ebv_labeller)
p5_se <-ebv_plots_se[[5]]  + labs(x = "Week of the year" , y = "Genetic diversity")+ facet_wrap(~ebv , labeller =  ebv_labeller)

# ----------------------------------- 
# Make some predictions
# Temporal 
temp_preds_mg <- sapply(1:length(modList_mg), function(x) {
  
  newData <- expand.grid(
    week_year       = seq(1, 52, 1),
    sample_time     = 7,
    trap_ID         = levels(ebvs_mg[[x]]$trap_ID),
    latitude_WGS84  = -18.33365,
    longitude_WGS84 = 47.34408
  )
  
  preds <- predict(modList_mg[[x]], newdata = newData, type = 'response')
  
})
# ----------------------------------- 
# Plot data
ebv_plots_mg <- plot_preds_temp(temp_preds_mg, ebvs_mg, c("nOTU", "lcbd", "FDis", "FEve", "mean_shn"))

p1_mg <- ebv_plots_mg[[1]] + labs(x = "Week of the year", y = "Species richness")+ facet_wrap(~ebv , labeller =  ebv_labeller)
p2_mg <- ebv_plots_mg[[2]] + labs(x = "Week of the year", y = "LCBD")+ facet_wrap(~ebv , labeller =  ebv_labeller)
p3_mg <- ebv_plots_mg[[3]] + labs(x = "Week of the year", y = "Functional dispersal")+ facet_wrap(~ebv , labeller =  ebv_labeller)
p4_mg <- ebv_plots_mg[[4]] + labs(x = "Week of the year", y = "Functional eveness")+ facet_wrap(~ebv , labeller =  ebv_labeller)
p5_mg <- ebv_plots_mg[[5]] + labs(x = "Week of the year", y = "Genetic diversity")+ facet_wrap(~ebv , labeller =  ebv_labeller)

# combine plots -----------------------------------------------------------


library(patchwork)
tiff("plots/temporal_trends.tif" , width = 2000,height=1000,compression = "lzw")

(p1_se + p2_se + p3_se + p4_se + p5_se +
 p1_mg + p2_mg + p3_mg + p4_mg + p5_mg) +
  plot_layout(ncol = 5)

dev.off()

browseURL("plots/temporal_trends.tif")


# spatial trends ----------------------------------------------------------

sweden <- ne_countries(scale = "medium", country = "sweden", returnclass = "sf")
madagascar <- ne_countries(scale = "medium", country = "madagascar", returnclass = "sf")


# plot --------------------------------------------------------------------

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
  plot_layout(nrow = 2) + plot_annotation(tag_levels = 'A')
dev.off()

browseURL("plots/spatial_trends.tif")
