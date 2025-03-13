# download ERA5 reanalysis climate data -------------------------------------------------------
library(ecmwfr)


# functions -----------------------------------------------------------------------------------

download_era5 <- function(dataset_short_name ,daily_statistic,  variable, year , month,day,target){
  
  request <- list(
    dataset_short_name = dataset_short_name ,
    product_type = "reanalysis",
    variable = variable,
    daily_statistic = daily_statistic,
    time_zone = "utc+00:00",
    frequency = "1_hourly",
    year = "2019",
    month = month,
    day = day,
    time = "00:00",
    data_format = "netcdf",
    download_format = "unarchived",
    target = target
  )
  
  file <- wf_request(
    request  = request,  # the request
    transfer = TRUE,     # download the file
    path     = "."       # store data in current working directory
  )
  
  
}

# requests ------------------------------------------------------------------------------------
months <- as.character(1:12)
days   <- as.character(1:365)

# temperature ---------------------------------------------------------------------------------


# 2m temp max 2019
download_era5(dataset_short_name = "derived-era5-single-levels-daily-statistics",
              variable = "2m_temperature",
              year = "2019",
              daily_statistic = "daily_maximum",
              month = months,
              day = days,
              target = "data/tidydata/env_features/era5/2019_era5_2m_temp_max.nc")

# 2m temp min 2019
download_era5(dataset_short_name = "derived-era5-single-levels-daily-statistics",
              variable = "2m_temperature",
              year = "2019",
              daily_statistic = "daily_minimum",
              month = months,
              day = days,
              target = "data/tidydata/env_features/era5/2019_era5_2m_temp_min.nc")

# 2m temp max 2020
download_era5(dataset_short_name = "derived-era5-single-levels-daily-statistics",
              variable = "2m_temperature",
              year = "2020",
              daily_statistic = "daily_maximum",
              month = months,
              day = days,
              target = "data/tidydata/env_features/era5/2020_era5_2m_temp_max.nc")

# 2m temp min 2020
download_era5(dataset_short_name = "derived-era5-single-levels-daily-statistics",
              variable = "2m_temperature",
              year = "2020",
              daily_statistic = "daily_minimum",
              month = months,
              day = days,
              target = "data/tidydata/env_features/era5/2020_era5_2m_temp_min.nc")



# precipitation -------------------------------------------------------------------------------


# precip 2019
download_era5(dataset_short_name = "derived-era5-single-levels-daily-statistics",
              variable = "total_precipitation",
              year = "2019",
              month = months,
              daily_statistic = "daily_mean",
              day = days,
              target = "data/tidydata/env_features/era5/2019_era5_prec.nc")

# precip 2020
download_era5(dataset_short_name = "derived-era5-single-levels-daily-statistics",
              variable = "total_precipitation",
              year = "2020",
              daily_statistic = "daily_mean",
              month = months,
              day = days,
              target = "data/tidydata/env_features/era5/2019_era5_prec.nc")

# Leaf area index ---------------------------------------------------------------------------------

download_era5(dataset_short_name = "derived-era5-single-levels-daily-statistics",
              variable = "leaf_area_index_low_vegetation",
              year = "2019",
              daily_statistic = "daily_mean",
              month = months,
              day = days,
              target = "data/tidydata/env_features/era5/2019_era5_lai_low.nc")


download_era5(dataset_short_name = "derived-era5-single-levels-daily-statistics",
              variable = "leaf_area_index_high_vegetation",
              year = "2019",
              daily_statistic = "daily_mean",
              month = months,
              day = days,
              target = "data/tidydata/env_features/era5/2019_era5_lai_high.nc")


download_era5(dataset_short_name = "derived-era5-single-levels-daily-statistics",
              variable = "leaf_area_index_low_vegetation",
              year = "2020",
              daily_statistic = "daily_mean",
              month = months,
              day = days,
              target = "data/tidydata/env_features/era5/2020_era5_lai_low.nc")


download_era5(dataset_short_name = "derived-era5-single-levels-daily-statistics",
              variable = "total_precipitation",
              year = "2020",
              daily_statistic = "daily_mean",
              month = months,
              day = days,
              target = "data/tidydata/env_features/era5/2020_era5_prec.nc")




