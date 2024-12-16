################################################################################
################################################################################
################################################################################
############### FUNCTIONS ######################################################
################################################################################
################################################################################

edit_timeseries <- function(fname)
{
  x <- read_csv(fname, show_col_types = FALSE)
  x <- x %>% dplyr::select(all_of(c("time", "rr")))
  # Get full date 
  start <- first(x$time)
  end <- last(x$time)
  df_date <- tibble(time = seq(start, end, by = "hour"))
  x <- left_join(df_date, x, by = "time")
  x <- x %>% mutate(rr2 = rollsumr(rr, fill = NA, k = 2), date = floor_date(time, "day"))
  x <- x %>% group_by(date) %>% summarize(rr = max(rr), rr2 = max(rr2)) %>% ungroup()
  x
}

get_missing_values <- function(x)
{
  x <- x %>% drop_na()
  fl <- length(seq(first(x$date), last(x$date), by = "day"))
  l <- dim(x)[1]
  round((1 - l/fl)*100, digits = 1)
}

get_nyear_ts <- function(x)
{
  x <- x %>% drop_na()
  nyear <- dim(x)[1]/365
  round(nyear, digits = 2)
}

get_threshold <- function(x, k = 4)
{
  x <- na.omit(x)
  l <- length(x)
  ind <- ceiling(l/365)*k
  th <- sort(x, decreasing = TRUE)[ind]
  th
}

get_model <- function(x, nyear, var_name = "rr2", method = "Lmoments")
{
  x <- x %>% drop_na()
  if(nyear < 25)
  {
    xx <- pull(x, !!as.name(var_name))
    th <- get_threshold(x = xx, k = 3)
    m <- fevd(x = na.omit(xx), threshold = th, type = "GP", method = method)
  }
  else
  {
    am <- x %>% mutate(year = year(date)) %>% group_by(year) %>% 
      summarize(am = max(!!as.name(var_name), na.rm = TRUE)) %>% pull(am)
    m <- fevd(x = am, type = "GEV", method = method)
  }
  m
}

get_p_event <- function(q, m) # Event and model is needed
{
  type <- m$type
  method <- m$method
  if(method == "Lmoments")
  {
    par <- findpars(m) 
  }
  if(method == "MLE")
  {
    par <- m$results$par
  }
  if(type == "GP")
  {
    p <- pevd(q = q, scale = par[1], shape = par[2], threshold = m$threshold, type = type)  
    rp <- 1 / (3 * (1 - p))
    ci <- ci(x = m, type = "par", return.period = rp, alpha = 0.5)
    p <- apply(ci, 2, function(par)
      pevd(q = q, scale = par[1], shape = par[2], threshold = m$threshold, type = type))
    rp <- 1/(1-p)
  }
  if(type == "GEV")
  {
    p <- pevd(q = q, loc = par[1], scale = par[2], shape = par[3], type = type) 
    rp <- 1/(1 - p)
    ci <- ci(x = m, type = "par", return.period = rp, alpha = 0.5)
    p <- apply(ci, 2, function(par) 
      pevd(q = q, loc = par[1], scale = par[2], shape = par[3], type = type))
    rp <- 1/(1-p)
  }
  rp
}

comp_trend <- function(x, nyear, var_name = "rr2")
{
  x <- x %>% drop_na()
  if(nyear < 25)
  {
    pot <- pull(x, !!as.name(var_name))
    th <- get_threshold(x = pot, k = 3)
    tt <- bbsmk(pot[pot >= th])
  }
  else
  {
    am <- x %>% mutate(year = year(date)) %>% group_by(year) %>% 
      summarize(am = max(!!as.name(var_name), na.rm = TRUE)) %>% pull(am)
    tt <- bbsmk(am) 
  }
  tt
}

get_metadata <- function(resolution = "daily") ### daily gets daily data, hourly 1h data, and min 10min data
{
  parsed_obj <- switch(resolution,
                       daily = httr::content(httr::GET("https://dataset.api.hub.zamg.ac.at/v1/station/historical/klima-v2-1d/metadata"), as="text") %>% 
                         fromJSON(),
                       hourly = httr::content(httr::GET("https://dataset.api.hub.zamg.ac.at/v1/station/historical/klima-v2-1h/metadata"), as="text") %>% 
                         fromJSON(),
                       min = httr::content(httr::GET("https://dataset.api.hub.zamg.ac.at/v1/station/historical/klima-v2-10min/metadata"), as="text") %>% 
                         fromJSON(),
                       tawes = httr::content(httr::GET("https://dataset.api.hub.zamg.ac.at/v1/station/historical/tawes-v1-10min/metadata"), as="text") %>% 
                         fromJSON())
  stations <- parsed_obj$stations
  stations <- stations %>% mutate(valid_from = gsub("\\+00:00", "", valid_from), valid_to = gsub("\\+00:00", "", valid_to), valid_to = gsub("2100", "2022", valid_to))
  stations <- stations %>% mutate(group_id = as.integer(group_id), id = as.integer(id))
  par <- parsed_obj$parameters[,1:4]
  out <- list(stations, par)
  names(out) <- c("metadata", "climatic_variables")
  out
}

download_ts <- function(resolution = "daily", parameter, valid_from, valid_to, id, path_out = NULL)
{
  url <- switch(resolution, 
                daily = paste0("https://dataset.api.hub.zamg.ac.at/v1/station/historical/klima-v2-1d?parameters=", parameter, "&start=",
                               valid_from, "&end=", valid_to,"&station_ids=",id,"&output_format=csv"),
                hourly = paste0("https://dataset.api.hub.zamg.ac.at/v1/station/historical/klima-v2-1h?parameters=", parameter, "&start=",
                                valid_from, "&end=", valid_to,"&station_ids=",id,"&output_format=csv"),
                min = paste0("https://dataset.api.hub.zamg.ac.at/v1/station/historical/klima-v2-10min?parameters=", parameter, "&start=",
                             valid_from, "&end=", valid_to,"&station_ids=",id,"&output_format=csv"),
                tawes = paste0("https://dataset.api.hub.zamg.ac.at/v1/station/historical/tawes-v1-10min?parameters=", parameter, "&start=",
                               valid_from, "&end=", valid_to,"&station_ids=",id,"&output_format=csv"))
  x <- read_csv(url, show_col_types = FALSE)
  x <- x %>% dplyr::select(!all_of(c("station")))
  if(is.null(path_out))
  {
    out <- x
  }
  else
  {
    fname_out <- paste0(path_out, id, ".csv")
    write_csv(x, fname_out)
    out <- "Saved"
  }
  out
}

get_date <- function(fname)
{
  b <- brick(fname)
  l <- dim(b)[3]
  day <- round(l/24, digits = 0)
  bn <- basename(fname)
  lbn <- last(strsplit(bn, split = "_")[[1]])
  dt <- gsub(".nc", "", lbn)
  yr <- substr(dt, 1, 4)
  m <- substr(dt, 5, 6)
  if(yr == "2023")
  {
    start <- ymd_hms(paste0(yr, "-", m, "-01 00:00:15"))
  }
  else
  {
    start <- ymd_hms(paste0(yr, "-", m, "-01 00:00:00")) 
  }
  end <- ymd_hms(paste0(yr, "-", m, "-", day," 24:00:00"))
  date <- seq(start, end, by = "hour")
  date
}

get_area <- function(fname)
{
  b <- brick(fname)
  arr <- as.array(b)
  res50 <- apply(arr, 3, FUN = function(x) sum(x >= 50))
  res75 <- apply(arr, 3, FUN = function(x) sum(x >= 75))
  res100 <- apply(arr, 3, FUN = function(x) sum(x >= 100))
  date <- get_date(fname)
  print(fname)
  df <- tibble(area_50 = res50, area_75 = res75, area_100 = res100, date = date)
  df
}
get_max_month <- function(fname, year)
{
  b <- brick(fname)
  if(year == 2003)
  {
    b <- calc(b, fun = function(x) ifelse(x > 75, 0, x)) 
  }
  out <- max(b, na.rm = TRUE)
  out
}

get_max_year <- function(flist, year)
{
  year_str <- paste0("_", year)
  fnames <- grep(year_str, flist, value = TRUE)
  l <- lapply(fnames, get_max_month, year = year)
  b <- brick(l)
  out <- max(b, na.rm = TRUE)
  out
}

calc_rp <- function(x, rp = 20)
{
  m <- fevd(x = x, type = "GEV", method = "Lmoments")
  cci <- ci(x = m, type = "return.level", return.period = rp)
  cci
}
