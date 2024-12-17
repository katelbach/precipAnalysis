################################################################################
################################################################################
#### Download station data of GeoSphere ########################################
################################################################################

# Loading libraries
library(tidyverse)
library(lubridate)
library(jsonlite)

# Loading functions
source("functions_extreme_precip.R")

# Downloading Metadata
md <- get_metadata(resolution = "hourly")

# Sub-selecting area around Vienna
x <- md[[1]] %>% filter(between(lon, 16.1,16.65) & between(lat,48.1,48.35)) %>% 
  filter(is.na(group_id)) %>% mutate(yr = year(valid_to)) %>% filter(yr == 2022)
x <- as_tibble(x)

# Downloading station data

path_out <- "../stations_data/"
for(i in 1:11)
{
  yr <- first(strsplit(x$valid_from[i], split = "-")[[1]])
  vf <- ifelse(as.integer(yr) < 1941, "1941-01-01T00:00:00", x$valid_from[i])
  download_ts(resolution = "hourly", parameter = "rr", valid_from = vf, 
              valid_to = "2024-11-30T00:00:00", id = x$id[i], path_out = path_out)
}



