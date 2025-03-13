# format era5 ---------------------------------------------------------------------------------
library(tidyverse)
library(raster)
# data ----------------------------------------------------------------------------------------

flist <- list.files("data/tidydata/env_features/era5/" , full.names = TRUE)

IBA_locs_se  <- readRDS("data/tidydata/all_meta_se.rds") |>
                dplyr::select(longitude_WGS84 , latitude_WGS84,,trap_ID = trapID) |>
                drop_na() |> distinct() |> as.data.frame()

IBA_locs_mg <- readRDS("data/tidydata/all_meta_mg.rds") |>
              dplyr::select(longitude_WGS84 , latitude_WGS84,,trap_ID = trapID) |>
              drop_na() |> distinct()  |> as.data.frame()


#Lookups for site vs ID

sites_lk_se <- IBA_locs_se |> dplyr::select(trap_ID) |> mutate(ID = 1:n())
sites_lk_mg <- IBA_locs_mg |> dplyr::select(trap_ID) |> mutate(ID = 1:n())

# function ------------------------------------------------------------------------------------
# Function to extract daily data from an ERA5 NetCDF file for specific locations
extract_daily_era5 <- function(era5_file, locations) {
  
  # Extract the variable name from the file path by removing unnecessary parts
  varname <- str_remove_all(era5_file, c("data/tidydata/env_features/era5//20[0-9]+_era5_|.nc$|2m_")) 
  
  # Load the NetCDF file as a raster brick, which represents multiple raster layers (one per day)
  era5_daily_layers <- raster::brick(era5_file)
  
  # Extract raster values at specified locations and transform the data into a tidy format
  site_era5 <- raster::extract(era5_daily_layers, locations, df = TRUE) |> 
    # Convert from wide format (one column per day) to long format
    pivot_longer(-ID) |> 
    # Extract the numeric day from the layer names
    mutate(day = as.numeric(str_remove(name, "X"))) |> 
    # Keep only relevant columns: location ID, day, and extracted values
    dplyr::select(ID, day, value) |> 
    # Rename the "value" column with the variable name derived earlier
    rename(!!varname := value)
  
  # Return the processed data frame
  return(site_era5)
}

# Function to compute weekly mean values of an ERA5 variable
weekly_mean_era5 <- function(site_era5) {
  
  # Construct the name for the new variable representing the weekly mean
  new_var <- paste0("wm_", colnames(site_era5[3]))
  
  # Process the input data to compute weekly means
  wm_era5 <- site_era5 |> 
    # Shift days by 1 to align weeks correctly and calculate the week number
    mutate(day = day + 1, week = ceiling(day / 7)) |> 
    # Group data by location ID and week number
    group_by(ID, week) |> 
    # Calculate the mean of the third column (the ERA5 variable) for each group
    summarise(!!new_var := mean(cur_data_all()[[3]])) |> 
    # Remove grouping to return a standard data frame
    ungroup()
  
  # Return the processed data frame with weekly mean values
  return(wm_era5)
}

#Function to check variables 
check_vars <- function(era5_vars , y_var) {
  ggplot(era5_vars,aes(week,{{ y_var }},group=trap_ID))+
    geom_line()
}
# sweden --------------------------------------------------------------------------------------


wm_tmax_se <- extract_daily_era5(flist[[1]] , IBA_locs_se[,1:2]) |> 
              weekly_mean_era5() |> 
              left_join(sites_lk_se)


wm_tmin_se <- extract_daily_era5(flist[[2]] , IBA_locs_se[,1:2]) |> 
              weekly_mean_era5() |> 
              left_join(sites_lk_se)

wm_laih_se <- extract_daily_era5(flist[[3]] , IBA_locs_se[,1:2]) |> 
              weekly_mean_era5() |> 
              left_join(sites_lk_se)

wm_lail_se <- extract_daily_era5(flist[[4]] , IBA_locs_se[,1:2]) |> 
              weekly_mean_era5() |> 
              left_join(sites_lk_se)


wm_prec_se <- extract_daily_era5(flist[[5]] , IBA_locs_se[,1:2]) |> 
              weekly_mean_era5() |> 
              left_join(sites_lk_se)



# madagascar ----------------------------------------------------------------------------------


wm_tmax_mg <- extract_daily_era5(flist[[6]] , IBA_locs_mg[,1:2]) |> 
              weekly_mean_era5() |> 
              left_join(sites_lk_mg)


wm_tmin_mg <- extract_daily_era5(flist[[7]] , IBA_locs_mg[,1:2]) |> 
              weekly_mean_era5() |> 
              left_join(sites_lk_mg)

wm_laih_mg <- extract_daily_era5(flist[[8]] , IBA_locs_mg[,1:2]) |> 
              weekly_mean_era5() |> 
              left_join(sites_lk_mg)

wm_lail_mg <- extract_daily_era5(flist[[9]] , IBA_locs_mg[,1:2]) |> 
              weekly_mean_era5() |> 
              left_join(sites_lk_mg)


wm_prec_mg <- extract_daily_era5(flist[[10]] , IBA_locs_mg[,1:2]) |> 
              weekly_mean_era5() |> 
              left_join(sites_lk_mg)


# tidy data -----------------------------------------------------------------------------------

# Sweden
vars_seL     <- list(wm_tmax_se , wm_tmin_se, wm_laih_se , wm_lail_se , wm_prec_se)
era5_vars_se <- purrr::reduce(vars_seL, dplyr::left_join) |>
                dplyr::select(trap_ID, week, matches("temp|prec|lai")) |> 
                mutate(across(matches("temp") , ~.x-273.15))

#Madagascar
vars_mgL     <- list(wm_tmax_mg , wm_tmin_mg, wm_laih_mg , wm_lail_mg , wm_prec_mg)
era5_vars_mg <- purrr::reduce(vars_mgL, dplyr::left_join) |>
                dplyr::select(trap_ID, week, matches("temp|prec|lai")) |> 
                mutate(across(matches("temp") , ~.x-273.15))


# Check variables 
check_vars(era5_vars_se , wm_temp_max)
check_vars(era5_vars_se , wm_temp_min)
check_vars(era5_vars_se , wm_lai_low)
check_vars(era5_vars_se , wm_lai_high)
check_vars(era5_vars_se , wm_prec)

check_vars(era5_vars_mg , wm_temp_max)
check_vars(era5_vars_mg , wm_temp_min)
check_vars(era5_vars_mg , wm_lai_low)
check_vars(era5_vars_mg , wm_lai_high)
check_vars(era5_vars_mg , wm_prec)



# save ----------------------------------------------------------------------------------------

saveRDS(era5_vars_se,"data/tidydata/env_features/weekly_era5_se.rds")
saveRDS(era5_vars_mg,"data/tidydata/env_features/weekly_era5_mg.rds")
