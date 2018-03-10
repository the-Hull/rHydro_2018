
# Load Packages -----------------------------------------------------------


library(rnrfa)
library(dplyr)
library(ggplot2)
library(magrittr)
library(purrr)
library(rgdal)
library(fs)
library(broom)
library(forcats)

source("./src/01_presentation.R")

# Pull meta data ----------------------------------------------------------



# Collect all station information
all_stations <- catalogue() # or rnfra::stationSummary 

# Initial exploration:
names(all_stations)

all_stations <- all_stations %>% 
      select(-`ma-station-id`, 
             -`maximum-gauging-stage-date-time`,
             -`maximum-gauging-flow-date-time`,
             -benchmark2,
             -categories) 


# Initial Checks ----------------------------------------------------------



# first look at data
str(all_stations)

# select columns where type change is necessary
char2num <- c("catchmentArea", "altitude", "maximum-gauging-flow", "lat", "lon")

all_stations[ ,char2num] %<>% sapply(as.numeric)
sapply(all_stations, class)

# plot 
exp_catchment_area.plot <- all_stations %>% 
      ggplot(aes(x = catchmentArea,
                 y = `maximum-gauging-flow`,
                 col = altitude)) +
      geom_point() +
      theme_bw() 
exp_catchment_area.plot



# Add country label -------------------------------------------------------

shp_files <- fs::dir_info(path = "./dat/GBR_adm",
                          glob = "*.shp") #from http://www.diva-gis.org/datadown

ogrListLayers("./dat/GBR_adm/GBR_adm0.shp")
gb_shp <- shp_files$path %>% 
      map(readOGR)


gb_countries <- gb_shp[[2]]

coordinates(all_stations) <- ~ lon + lat
proj4string(all_stations) <- proj4string(gb_countries)


all_stations<- all_stations@data %>% 
      mutate(country =  over(all_stations, gb_countries[ ,"NAME_1"])[[1]])

exp_catchment_area.plot <- all_stations %>% 
      ggplot(aes(x = catchmentArea,
                 y = `maximum-gauging-flow`,
                 col = country)) +
      geom_point(alpha = .3) +
      facet_wrap(~country, scales = "free_x") +
      theme_bw()
exp_catchment_area.plot




# First look at data retrieval --------------------------------------------

flow_data_subset <- gdf(all_stations$id[1:3], cl = NULL, metadata = TRUE) 

# Check units
purrr::modify_depth(flow_data_subset, 1, "meta") %>%
      map("units") %>%
      map_chr(as.character)

# repeat for full data set


system.time( flow_data <- gdf(all_stations$id, cl = NULL, metadata = TRUE) )
# save(flow_data, file = "./dat/flow_data_full.Rda")

# user  system elapsed 
# 2231.28   18.69 5930.53 

# Check units
units <- purrr::modify_depth(flow_data, 1, "meta") %>%
      map("units") %>%
      map(as.character)



# Get data ----------------------------------------------------------------


library(parallel)
# Use detectCores() to find out many cores are available on your machine
cl <- makeCluster(getOption("cl.cores", detectCores()-1))


# flow_data_cl <-  gdf(all_stations$id, cl = cl, metadata = T) 
# flow_data_cl <- flow_data_cl %>% 
#       set_names(all_stations$id)
# 
# system.time(flow_data_cl_nometa <-  gdf(all_stations$id, cl = cl, metadata = F))
# flow_data_cl_nometa <- flow_data_cl_nometa %>% 
#       set_names(all_stations$id)

# pull available rain data for catchments
rain_data <- cmr(id = all_stations$id ,metadata = F, cl = cl)
rain_data <- rain_data %>% 
      set_names(all_stations$id) %>% 
      discard(is.null) %>% 
      map(~tibble::tibble(date_time = zoo::index(.x) %>% as.POSIXct(),
                          year = lubridate::year(date_time),
                          month = lubridate::month(date_time),
                          precipitation_mm_month = zoo::coredata(.x)[ ,1])) %>% 
      map2(names(.), 
           ~mutate(.x, id = as.factor(.y))
      )



# prep flow data set
flow_data <- flow_data %>% 
      map("data") %>% 
      discard(is.null) %>%
      map(~tibble::tibble(date_time = zoo::index(.x) %>% as.POSIXct(),
                          year = lubridate::year(date_time),
                          month = lubridate::month(date_time),
                          year_day = lubridate::yday(date_time), 
                          flow_m_3s = zoo::coredata(.x)[ ,1])) %>% 
      map2(names(.), 
           ~mutate(.x, id = as.factor(.y))
      ) 


rain_data <- rain_data %>%  keep(.p = names(.) %in% names(flow_data))
flow_data <- flow_data %>% keep(.p = names(.) %in% names(rain_data))
# map2(keep(rain_data_test, names(.x) %in% names(.y)),
# left_join, by = c(id, year, month))

flow_rain_data <- flow_data %>% 
      map2(rain_data %>% map(select, -date_time),
           left_join, by = c("id", "year", "month")) 

flow_rain_data <- flow_rain_data %>% 
      map(~mutate(.x, water_year = lfstat::water_year(date_time,
                                                  assign = "end",
                                                  origin = "usgs")))
          
          

flow_rain_data_agg <- flow_rain_data %>% 
      map(function(x){
            x %>% 
                  mutate(n_days_in_month = lubridate::days_in_month(x$date_time)) %>% 
                  group_by(water_year, month) %>% 
                  summarize(date_time = mean(date_time, na.rm = T) %>% lubridate::floor_date(unit = "month"),
                            count_n = n(),
                            flow_m3_month = sum(flow_m_3s*60*60*24,
                                                na.rm = T),
                            precipitation_mm_month = mean(precipitation_mm_month,
                                                          na.rm = T),
                            n_days_in_month_missing = mean(n_days_in_month, na.rm = T) -
                                  count_n
                            
                  )
      })


fr <- flow_rain_data_agg %>% 
      bind_rows(.id = "id")

fr <- fr %>% 
      left_join(all_stations %>%
                      select(catchmentArea, river, altitude, country, hydrometricArea, id),
                by = "id")

qp.plot <- fr %>% 
      ggplot(aes(x = precipitation_mm_month,
                 y= flow_m3_month / (catchmentArea * 1e6),
                 color = as.factor(country))) +
      geom_point(alpha = .3) +
      facet_wrap(~month) 

qp.plot


qbox.plot <- fr %>% 
      ggplot(aes(x = month %>% as.factor(),
                 y= flow_m3_month / (catchmentArea * 1e6))) +
      geom_violin(fill = "grey", color = NA) +
      geom_boxplot(width = .1) +
      facet_grid(country~., scales = "free_y") +
      scale_y_log10() +
      theme_bw()
qbox.plot

          # flow_rain_data <- flow_rain_data %>%
          #       map(~select(.x, date_time, flow_m_3s, precipitaion_mm_month)) %>% 
          #       bind_rows(.id = "id")
          


# Load Catchment Descriptors ----------------------------------------------


read_cd3 <- function(data){
      readr::read_csv(data, skip = 8) %>% 
      rbind(., set_names(
            as.data.frame(
                  matrix(
                        colnames(.),
                        ncol = 2)),
            colnames(.))) %>% 
      set_names(nm = c("k", "v")) %>% 
      filter(!stringr::str_detect(.$k, "\\[.*")) %>% 
      tidyr::spread(key = "k", value = "v") %>% 
            mutate(id = readLines(data)[6] %>% as.numeric())
}


cd3_data <- fs::dir_info("./dat/CEH_FEH/", glob = "*.CD3",recursive = T)$path %>% 
      map_df(read_cd3) %>% 
      mutate_if(.predicate = function(x) any(!stringr::str_detect(x, ".*[a-zA-Z]")==T & !is.na(x)),
                as.numeric) %>% 
      select_if(.predicate = is.numeric)


uk_nested <- all_stations %>% 
      left_join(cd3_data %>% mutate(id = as.character(id)), by = "id") %>% 
      tidyr::nest(-country)


uk_nested <- uk_nested %>% 
      mutate(model = map(data, ~lm(`maximum-gauging-flow` ~ catchmentArea, data = .x)),
             predict = map2(model, data, predict))

# uk_unnest <- uk_nested %>% 
#       tidyr::unnest(data, predict)


models <- list(
      null_intercept = function(x){lm(`maximum-gauging-flow` ~ 1, data = x)},
      ca = function(x){lm(`maximum-gauging-flow` ~ catchmentArea, data = x)},
      ca_dplbar_bfihost = function(x){lm(`maximum-gauging-flow` ~ 
                                               catchmentArea + DPLBAR + BFIHOST, data = x)},
      ca_dplbar_bfihost_farl = function(x){lm(`maximum-gauging-flow` ~ 
                                                    catchmentArea + DPLBAR + BFIHOST + FARL, data = x)},
      ca_dpsbar_bfihost_farl_fpext = function(x){lm(`maximum-gauging-flow` ~ 
                                                          catchmentArea + DPSBAR + BFIHOST + FARL + FPEXT, data = x)}
      
      )

apply_model <- function(.model, ndf){
      
      ndf$model <- map(ndf$data, possibly(.model, NULL))
      
      return(ndf)
}

uk_nest <- models %>% 
      map_df(apply_model, uk_nested, .id = "id_model") %>% 
      # mutate(prediction = map2(model, data, predict),
      #        resid = map(model, resid)) %>% 
      select(id_model, country, model) %>% 
      mutate(param = map(model, tidy),
             assess = map(model, glance)) %>% 
      select(-model) 
      
uk_param <- uk_nest %>% select(-assess) %>% tidyr::unnest(param)
uk_assess <- uk_nest %>% select(-param) %>%  tidyr::unnest(assess)

uk_param <- uk_param %>% 
      mutate(estimate_adj = scale(estimate, center = T, scale = T),
             estimate_adj_factor = ifelse(p.value < 0.05,
                                          ifelse(estimate > 0, "positive", "negative"),
                                          "not sign.") %>% 
                   fct_relevel("positive", "negative"))


 param.plot <- uk_param %>% 
      ggplot(aes(y = estimate_adj,
                 x = term, fill = estimate_adj_factor)) +
      geom_bar(stat = "identity", position = "dodge",col = "gray20") +
      coord_flip() +
      facet_grid(country~id_model) +
      geom_hline(yintercept = 0, linetype = 2, col = "darkgrey") +
      theme_bw() +
      theme(panel.grid.major.x = element_blank(),
            panel.grid.minor.x = element_blank()) +
      scale_fill_manual(name = "Estimate", values = c("negative" = "darkred",
                                                      "positive" = "steelblue",
                                                      "not sign." = "grey95"))
param.plot

uk_assess <- uk_assess %>% group_by(country) %>% arrange(AIC, .by_group = T) %>% 
      rename(AIC_val = AIC) %>% 
      mutate(minAIC = ifelse(AIC_val == min(AIC_val, na.rm = T),"minimum", "> minimum"))


uk_aic_min <- uk_assess %>%  summarise(minAIC_val = min(AIC_val, na.rm = T)) %>% 
            right_join(uk_assess %>% select(country, id_model))

aic.plot <-  uk_assess %>% 
      ggplot(aes(x = id_model, y = AIC_val, group = 1)) +
      geom_line(linetype = 2) +
      geom_point(aes(fill = uk_assess$minAIC), size = 5, shape = 21, col = 'gray20') +
      stat_summary(geom = "text",
                   fun.y = function(x){min(x,na.rm=T)*0.97},
                   label = round(uk_assess$adj.r.squared,2)) +
      facet_grid(country~., scales = 'free_y') +
      theme_presi() +
      theme(panel.grid.major.y = element_line(color = "gray90"),
            panel.border = element_blank(),
            panel.spacing = unit(2, "lines"),
            axis.text.x = element_text(angle = 15, hjust = 1),
            plot.title = element_text(hjust = 0),
            plot.subtitle = element_text(colour = 'gray60')) +
      scale_y_continuous(expand = c(0.1,0)) +
      scale_fill_manual(name = "AIC comparison", 
                        values = c("> minimum" = "gray40",
                                   "minimum" = "steelblue")) +
      labs(x = "Model", y = "AIC") +
      ggtitle(label = "Model selection", subtitle = "via relative and absolute goodness of fit")
      # geom_bar(stat = "identity", position = 'dodge', col = 'gray20') +
      # geom_text(data = uk_aic_min,
      #           aes(x = id_model,
      #               y = minAIC_val,
      #               label = .97*round(uk_assess$adj.r.squared,2))) +
aic.plot

# test
