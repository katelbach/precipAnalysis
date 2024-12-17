################################################################################
################################################################################
################ Analyzing station data ########################################
################################################################################

# Loading libraries
library(tidyverse)
library(lubridate)
library(sf)
library(zoo)
library(extRemes)
library(lmomco)
library(modifiedmk)
library(RColorBrewer)

# Load functions
source("functions_extreme_precip.R")

# Read data
path_out <- "../stations_data/"
x <- x %>% mutate(fname = paste0(path_out, id, ".csv"))
x <- x %>% mutate(df = pmap(.l = list(fname), .f = edit_timeseries))
x <- x %>% mutate(pmis = pmap_dbl(.l = list(df), .f = get_missing_values),
                  nyear = pmap_dbl(.l = list(df), .f = get_nyear_ts))
## Filter year 2024 for station analysis
x <- x %>% mutate(df = pmap(.l = list(df), .f = function(x) 
  x %>% mutate(year = year(date)) %>% filter(year != 2024)))

# Testing trend for Vienna Hohe-Warte
comp_trend(x$df[[3]], nyear = x$nyear[3]) 

# How many events over all stations larger than 50mm
x[c("name", "df")] %>% mutate(sel_events = pmap(.l = list(df), .f = function(df)
  df %>% filter(rr2 > 50))) %>% dplyr::select("name", "sel_events") %>% unnest(sel_events) 

### Get model for each station
x <- x %>% mutate(model_lm = pmap(.l = list(df, nyear), .f = function(df, nyear) get_model(x = df, nyear = nyear)),
                  model_mle = pmap(.l = list(df, nyear), .f = function(df, nyear) get_model(x = df, nyear = nyear, method = "MLE")))

# Return period for each model of Event August 2024
sapply(x$model_lm, get_p_event, q = 110.1)
sapply(x$model_mle, get_p_event, q = 110.1)


### Model for years 2004-2023 for each station
x <- x %>% mutate(am = pmap(.l = list(df), .f = function(x) 
  x %>% drop_na() %>% mutate(year = year(date)) %>% group_by(year) %>% summarize(am = max(rr2, na.rm = TRUE)) %>% 
    filter(between(year, 2004,2023)))) 

# Only use stations that have a full time series between 2004 and 2023
x2004 <- x %>% filter(nyear >= 20)
x2004 <- x2004 %>% mutate(m2004 = pmap(.l = list(am), .f = function(x) fevd(x = x$am, type = "GEV", method = "Lmoments")))
x2004 <- x2004 %>% mutate(rp2004 = pmap(.l = list(m2004), .f = function(model) ci(x = model, type = "return.level", return.period = c(20))))

### Create sf object with 20-year return period
xx <- x2004[c("name","lon", "lat", "rp2004")] %>% unnest(rp2004)
colnames(xx) <- c("name","lon", "lat","estimate")
xx <- xx %>% mutate(ci = rep(c("Lower CI", "Mean", "Upper CI"), times = 7))
xsf <- st_as_sf(xx, coords = c("lon", "lat"), crs = 4326)
xsf <- xsf %>% mutate(cut_estimate = cut(estimate, seq(20,120, by = 10), labels = seq(20,110, by = 10)))


#### List all INCA 2-hour precipitation files
fl <- list.files(path = "../data/INCA_Wien_2hour/", pattern = "*.nc$", full.names = TRUE)

# How many events are over 50mm, 75mm and 100mm in INCA data
ar <- do.call(rbind, lapply(fl, get_area))
ar %>% mutate(year = year(date)) %>% filter(year >= 2004) %>% 
  filter(area_100 > 0)

ar %>% mutate(year = year(date)) %>% filter(year >= 2004) %>% 
  filter(area_75 > 0)

ar %>% mutate(year = year(date)) %>% filter(year >= 2004) %>% 
  filter(area_50 > 0) %>% mutate(day = floor_date(date, "day")) %>% group_by(day) %>% 
  summarize(area50 = max(area_50))

### Annual Maxima time series for INCA data
am_inca <- brick(lapply(seq(2004,2023), get_max_year, flist = fl))

# 20 year return period for INCA data
inca20 <- calc(am_inca, calc_rp)

### Create stars object
stars_inca <- st_as_stars(inca20)
names(stars_inca) <- "estimate"
stars_inca <- stars_inca %>% mutate(cut_estimate = cut(estimate, seq(20,120, by = 10), 
                                                       labels = seq(20,110, by = 10)))
stars_inca <- st_set_dimensions(stars_inca, 3, values = c("Lower CI", "Mean", "Upper CI"), names = "ci")

### Shape file for Vienna
vienna <- read_sf("../data/shape_vienna.shp")

# Colour for plot
cols <- c("#F0F0F0", "#F6E8C3", "#DFC27D", brewer.pal(n = 7, "GnBu"))

# Plot of 20 year return period for INCA and stations with confidence intervals
plot_inca <- ggplot() +
  geom_stars(data = stars_inca["cut_estimate"], alpha = .5) +
  geom_sf(data = vienna, fill = NA, linewidth = 1) +
  geom_sf(data = xsf, aes(fill = cut_estimate), size = 3, shape = 21, col = "black", show.legend = FALSE) +
  scale_fill_manual(values = cols) +
  guides(fill = guide_legend(title = "", nrow = 1, label.position = "bottom")) +
  xlab("") +
  ylab("") +
  facet_wrap(~ci) +
  theme_bw() +
  theme(panel.grid = element_line(linewidth = .1), text = element_text(size = 16), 
        legend.position = "bottom", legend.spacing.x = unit(0, "pt"),  
        legend.spacing.y = unit(0, "pt"), legend.key.spacing = unit(1, "pt"), 
        legend.key.width = unit(36, "pt"), legend.box.just = "center", legend.box = "horizontal", 
        legend.text = element_text(hjust = -0.5, vjust = -0.1))