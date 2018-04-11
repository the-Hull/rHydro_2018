---
  title: "rHydro 2018"
author: "Alex Hurley <br> a.g.hurley@pgr.bham.ac.uk <br> https://aglhurley.rbind.io"
date: "March 24, 2018"
output: 
  html_document:
  theme: flatly
highlight: tango
toc: yes
toc_float: yes
code_folding: show
number_sections: true
df_print: paged
---
  
  
  ```{r setup, echo=FALSE, cache=FALSE, include=FALSE}
library(knitr)
# library(rmdformats)
library(magrittr)
# library(formatR)

## Global options
options(max.print="80")
opts_chunk$set(echo=FALSE,
               cache=TRUE,
               prompt=FALSE,
               tidy=FALSE,
               comment=NA,
               message=FALSE,
               warning=FALSE)
opts_knit$set(width=80, 
              fig.width = 10,
              fig.height = 10)
```



# Introduction

## Background

R has..

- a growing user base
- a plethora of packages for scientific purposes (stats, modelling, visualization). 
- a number of highly useful packages for hydrological work*flows* (e.g. `rnrfa`, `lfstat`, `airGR[teaching]`).
- increased its use as a general-purpose language (api's, file system management, `Rcpp`, `Rpy2`)
                                                   
                                                   <center>
                                                   
                                                   ![](../img/r_hydrology_packages.png)
                                                   </center>
                                                   
                                                   <br>
                                                   
                                                   ## Purpose
                                                   
                                                   This tutorial will showcase a set of R's most prominent and widely-applicable set of packages,  the `tidyverse`, as well as its programming philosophy, and how it can be applied in conjunction with hydrology work*flows*.
                                                   
                                                   <br>
                                                     
                                                     ## Aproach
                                                     
                                                     Using (peak flow and time series) data from the [National River Flow Archive](http://nrfa.ceh.ac.uk/) (NRFA), we will:
                                                     
                                                     - clean and tidy
                                                   - explore
                                                   - define and apply statistical models (answering questions, cf. below)
                                                   - streamline model selection
                                                   - communicate results
                                                   
                                                   
                                                   **Questions:**
                                                     
                                                     - Q1: What impact does catchment area have on peak flows, and does it differ between countries?
                                                     - Q2: Which other characteristics affect peak flows on catchment-level, and which set of characteristics make up  the most representative models?
                                                     
                                                     <br>
                                                     
                                                     # Why use the `tidyverse`?
                                                     
                                                     ## We are all.. 
                                                     
                                                     <center>
                                                     ![](../img/hydrology.png)
                                                   </center>
                                                     
                                                     <font color="#f2f2f2" size="+3">**.. humans**, **users**, **colleagues**, **learners**. </font>
                                                     
                                                     ---
                                                     
                                                     ## We all want to..
                                                     
                                                     <br>
                                                     
                                                     <center>
                                                     ![Data analysis pipeline (CITATION)](../img/tidy.png)
                                                   
                                                   </center>
                                                     <br>
                                                     
                                                     
                                                     ## Tidy philosphy
                                                     
                                                     <big> [The `tidyverse` philosophy mandates](http://tidyverse.tidyverse.org/articles/manifesto.html): </big>
                                                     
                                                     - Reuse existing data structures.
                                                   - Compose simple functions with the pipe.
                                                   - Embrace functional programming.
                                                   - Design for humans.
                                                   
                                                   <br>
                                                     
                                                     ## Concept: Piping
                                                     
                                                     <!-- $$ data~\rightarrow~f(data,~paramaters)~\rightarrow~output$$ -->
                                                     
                                                     > *data*  $\rightarrow$  **f**(*data*, parameters)  $\rightarrow$  outputs
                                                   
                                                   ```{r exp, echo = T, eval = F}
                                                   data %>% mean()
                                                   ```
                                                   <br>
                                                     `tidyverse` functions often replace `base` functions with pipe-compatible version.
                                                   
                                                   ```{r exp1, echo = T, eval = F}
                                                   names(data) <- c("A", "B")
                                                   data %>% set_names(c("A", "B"))
                                                   
                                                   ```
                                                   
                                                   <!-- That is, `rnorm(n = 10, mean = 5, sd = 3) %>% mean()` gives `r rnorm(n = 10, mean = 5, sd = 3) %>% mean() %>% round(2)` -->
                                                     
                                                     
                                                     
                                                     ---
                                                     
                                                     <big>**More tangible example with iris data set**: </big>
                                                     
                                                     ```{r iris_table, echo = F, rows.print=5}
                                                   library(dplyr)
                                                   
                                                   iris %>% 
                                                     select(Species, 
                                                            Sepal.Length, 
                                                            Sepal.Width, 
                                                            Petal.Length, 
                                                            Petal.Width) %>% # re-arrange data frame
                                                     head(20) 
                                                   
                                                   
                                                   # str(iris)
                                                   ```
                                                   
                                                   <br>
                                                     <big> **Application:** </big>
                                                     ```{r pipe_example, echo=T}
                                                   
                                                   
                                                   library(dplyr)
                                                   
                                                   iris %>% 
                                                     group_by(Species) %>% 
                                                     summarise_all(mean) %>% 
                                                     as.data.frame()
                                                   
                                                   
                                                   ```
                                                   
                                                   ---
                                                     
                                                     <br>
                                                     
                                                     ## Concept: Mapping
                                                     
                                                     
                                                     > *DATA* = (data1, data2, data3)  
                                                   > *DATA* $\rightarrow$  **map**(*DATA*, **f**, parameters) $\rightarrow$   **f**(*data[i]*, parameters) $\rightarrow$  outputs  
                                                   
                                                   ```{r ex2, eval=F, echo = T}
                                                   
                                                   # using purrr:map
                                                   DATA %>% map(summary)
                                                   
                                                   # for-loop
                                                   for(i in factor){
                                                     summary(DATA[i])
                                                   }
                                                   
                                                   # lapply
                                                   lapply(DATA, summary)
                                                   
                                                   ```
                                                   
                                                   <br>
                                                     
                                                     The `purrr::map` family of functions often offers additional, convenient functionalities.  
                                                   For example,
                                                   
                                                   
                                                   ```{r exp_lapply, echo=T, eval=F}
                                                   data <- lapply(list_of_files, read.csv)
                                                   data <- do.call(rbind, data)
                                                   ```
                                                   
                                                   can be replaced by 
                                                   
                                                   ```{r exp_df, echo=T, eval=F}
                                                   map_df(list_of_files, read.csv, .id = "origin")
                                                   ```
                                                   
                                                   
                                                   <big> **Application:** </big>
                                                     ```{r map_example, echo=T}
                                                   
                                                   
                                                   library(dplyr)
                                                   
                                                   iris %>% 
                                                     select(Species, Sepal.Length) %>% 
                                                     split(f = .$Species) %>% 
                                                     map(summary)
                                                   
                                                   
                                                   
                                                   ```
                                                   
                                                   <br>
                                                     
                                                     ---
                                                     
                                                     ## Concept: Nesting
                                                     
                                                     Data frames can store just about anything. Interestingly, on cell in a data frame can be made up of an entire list. 
                                                   
                                                   This makes for highly useful data management and analyses approaches:
                                                     
                                                     <big>Application: </big>
                                                     
                                                     ```{r exp_nest, rows.print = 3, echo = T}
                                                   library(tidyr)
                                                   
                                                   
                                                   # regular data frame
                                                   iris
                                                   
                                                   # nested data frame
                                                   nested_iris <- iris %>%
                                                     nest(-Species)
                                                   nested_iris
                                                   
                                                   # add a model
                                                   nested_iris <- 
                                                     nested_iris %>% 
                                                     mutate(model_1 = map(data, ~lm(Sepal.Length ~ Sepal.Width, data = .x)),
                                                            resid_1 = map(model_1, resid))
                                                   nested_iris
                                                   
                                                   
                                                   
                                                   ```
                                                   
                                                   <br>
                                                     
                                                     Think: hierarchichal data structures, meta data, multiple models, etc.  
                                                   The nesting can be reversed by via `unnest(nested_iris)`.
                                                   
                                                   ---
                                                     
                                                     <br>
                                                     
                                                     # Analyses: Data handling
                                                     
                                                     
                                                     
                                                     The next sections will:
                                                     
                                                     - load necessary packages
                                                   - download and read-in data using the dedicated UK NRFA package `rnfra`, and an archive of [catchment descriptors](http://nrfa.ceh.ac.uk/feh-catchment-descriptors), available [here](http://nrfa.ceh.ac.uk/winfap-feh-files).
                                                   - tidy data and prepare for modelling
                                                   
                                                   
                                                   ## Set-up
                                                   
                                                   Load necessary packages and custom plotting theme:
                                                     
                                                     ```{r pcks, echo = T, cache = F}
                                                   
                                                   # tidyverse packages
                                                   library(fs) # file management
                                                   library(dplyr) # data manipulation
                                                   library(purrr) # functional programming
                                                   library(magrittr) # data manipulation / tidy code
                                                   library(ggplot2) # data viz
                                                   library(broom) # tidy stats results 
                                                   library(forcats) # tidy factors
                                                   library(lubridate) # tidy dates
                                                   library(tibble)
                                                   library(tidyr)
                                                   
                                                   # Environmental and spatial data/analyses
                                                   library(rnrfa) # flow archive
                                                   library(rgdal) # spatial data sets
                                                   
                                                   # custom theme
                                                   source("./src/01_presentation.R")
                                                   
                                                   
                                                   
                                                   ```
                                                   
                                                   ---
                                                     
                                                     <br>
                                                     
                                                     
                                                     ## Peak Flows and Meta Data
                                                     
                                                     ---
                                                     
                                                     **Note:** *Data from the UK National River Flow Archive*.  
                                                   Please refer to http://nrfa.ceh.ac.uk/costs-terms-and-conditions prior to any re-use
                                                   
                                                   ---
                                                     
                                                     
                                                     ```{r station_info, echo = T, cache = T, message=FALSE, max.print = 20, rows.print = 4}
                                                   
                                                   
                                                   # Catchment and Peak Flow Data ----------------------------
                                                   
                                                   # Collect all station information using the rnrfa API
                                                   
                                                   all_stations <- catalogue() # rnfra::stationSummary available for meta data
                                                   
                                                   # remove some columns
                                                   all_stations <- all_stations %>% 
                                                     select(-`ma-station-id`, 
                                                            -`maximum-gauging-stage-date-time`,
                                                            -`maximum-gauging-flow-date-time`,
                                                            -benchmark2,
                                                            -categories) 
                                                   
                                                   
                                                   
                                                   # columns where class change is necessary
                                                   char2num <- c("catchmentArea", "altitude", "maximum-gauging-flow", "lat", "lon")
                                                   
                                                   all_stations[ ,char2num] %<>% map_dfc(as.numeric)
                                                   
                                                   
                                                   
                                                   
                                                   
                                                   # Next lines for adding country label based on lat/lon location
                                                   
                                                   # file info of shape files
                                                   shp_files <- fs::dir_info(path = "./dat/GBR_adm",
                                                                             glob = "*.shp") #from http://www.diva-gis.org/datadown
                                                   
                                                   
                                                   # read in UK admin. shape files
                                                   gb_shp <- shp_files$path %>% 
                                                     map(readOGR)
                                                   
                                                   # select country level admin. demarkation
                                                   gb_countries <- gb_shp[[2]]
                                                   
                                                   # adjust projection
                                                   coordinates(all_stations) <- ~ lon + lat
                                                   proj4string(all_stations) <- proj4string(gb_countries)
                                                   
                                                   
                                                   # add labels
                                                   all_stations<- all_stations@data %>% 
                                                     mutate(country =  over(all_stations, gb_countries[ ,"NAME_1"])[[1]])
                                                   
                                                   
                                                   # make table
                                                   all_stations %>% 
                                                     # select(id, country, river, location, catchmentArea, `maximum-gauging-flow`) %>%
                                                     as.data.frame()
                                                   
                                                   ```
                                                   
                                                   
                                                   
                                                   
                                                   
                                                   
                                                   
                                                   ```{r exp_plot, echo = T, fig.align='center', cache = F}
                                                   
                                                   # exploratory plot 
                                                   exp_catchment_area.plot <- all_stations %>% 
                                                     ggplot(aes(x = catchmentArea,
                                                                y = `maximum-gauging-flow`,
                                                                col = country,
                                                                text = paste("Name:",name))) +
                                                     geom_point(alpha = 0.5) +
                                                     scale_x_log10() +
                                                     scale_y_log10() +
                                                     # theme_presi() +
                                                     ggtitle("Relationship between peak flow and catchment area") +
                                                     labs(y = "log peak flow (cumecs)", x = "log Catchment Area (sq km)") +
                                                     facet_wrap(~country)
                                                   
                                                   
                                                   plotly::ggplotly(exp_catchment_area.plot)
                                                   
                                                   
                                                   ```
                                                   
                                                   ---
                                                     
                                                     <br>
                                                     
                                                     ## FEH Catchment Descriptors
                                                     
                                                     The [FEH catchment descriptor data](http://nrfa.ceh.ac.uk/winfap-feh-files) is stored in individual *\*.cd3* files per catchment (`id` as file name and unique identifier in file). They contain some preceding and trailing meta data in a non-rectangular format. We extract and skip lines as needed in a custom function to achieve a programmatic read-in.
                                                   
                                                   
                                                   ---
                                                     
                                                     **Note:** *Data from the UK National River Flow Archive*.  
                                                   Please refer to http://nrfa.ceh.ac.uk/costs-terms-and-conditions prior to any re-use
                                                   
                                                   ---
                                                     
                                                     ```{r cd3_function, echo=F}
                                                   
                                                   # custom function to read in files
                                                   read_cd3 <- function(data){
                                                     
                                                     readr::read_csv(data, skip = 8) %>%
                                                       rbind(., 
                                                             set_names(
                                                               as.data.frame(
                                                                 matrix(
                                                                   colnames(.),
                                                                   ncol = 2
                                                                 )
                                                               ),
                                                               colnames(.)
                                                             )
                                                       ) %>%
                                                       set_names(nm = c("k", "v")) %>%
                                                       filter(!stringr::str_detect(.$k, "\\[.*")) %>%
                                                       tidyr::spread(key = "k", value = "v") %>%
                                                       mutate(id = readLines(data)[6] %>% as.numeric())
                                                   }
                                                   
                                                   ```
                                                   
                                                   
                                                   
                                                   ```{r feh, echo=T, cache=T, eval = T, max.print = 20, rows.print = 5}
                                                   
                                                   
                                                   # uses custom function read_cd3, source code available in raw document.
                                                   
                                                   
                                                   cd3_data <- fs::dir_info("./dat/CEH_FEH/",
                                                                            glob = "*.CD3",
                                                                            recursive = T)$path %>%
                                                     map_dfr(read_cd3) %>%
                                                     mutate_if(.predicate = function(x){
                                                       any(!stringr::str_detect(x, ".*[a-zA-Z]")==T & !is.na(x))},
                                                       .funs = as.numeric) %>%
                                                     select_if(.predicate = is.numeric) %>%
                                                     mutate_each(funs(replace(., .<0, NA)))
                                                   
                                                   cd3_data %>% head(20) %>% as.data.frame()
                                                   
                                                   
                                                   
                                                   ```
                                                   
                                                   ---
                                                     
                                                     <br>
                                                     
                                                     ## Precipitation data
                                                     
                                                     ---
                                                     
                                                     **Note:** *Data from the UK National River Flow Archive*.  
                                                   Please refer to http://nrfa.ceh.ac.uk/costs-terms-and-conditions prior to any re-use
                                                   
                                                   ---
                                                     
                                                     ```{r precip, echo = T, cache = T, eval = T, rows.print = 5}
                                                   
                                                   library(parallel)
                                                   
                                                   cl <- makeCluster(getOption("cl.cores", detectCores()-1))
                                                   
                                                   # pull available rain data for catchments into list
                                                   # set names for each list according to id's
                                                   # drop empty elements
                                                   rain_data <- cmr(id = all_stations$id ,metadata = F, cl = cl) %>% 
                                                     set_names(all_stations$id) %>% 
                                                     discard(is.null) 
                                                   
                                                   # create tibble data structure in each element, using map_dfr...
                                                   # to combine list elements and add a column of id's
                                                   rain_data <- rain_data %>%
                                                     map_dfr(~tibble(date_time = zoo::index(.x) %>% as.POSIXct(),
                                                                     year = year(date_time),
                                                                     month = month(date_time),
                                                                     p_mm_month = zoo::coredata(.x)[ ,1]),
                                                             .id = "id")
                                                   
                                                   
                                                   rain_data %>% head(20) %>% as.data.frame()
                                                   
                                                   
                                                   ```
                                                   
                                                   ```{r precip_process, echo = T, rows.print = 5, cache = T}
                                                   
                                                   
                                                   # make tibble of mean of annual max. precip by
                                                   # 1) find annual maximum by grouping per yer
                                                   # 2) take mean of all years
                                                   # for all stations
                                                   
                                                   rain_data <- rain_data %>%
                                                     group_by(id, year) %>%
                                                     summarize(n_months = n(),
                                                               max_p_mm = ifelse(n_months == 12, 
                                                                                 max(p_mm_month, na.rm = T),
                                                                                 NA)) %>% 
                                                     group_by(id) %>%
                                                     summarize(n_year = n(),
                                                               mean_max_p_mm = mean(max_p_mm, na.rm = T)) %>% 
                                                     filter(n_year > 25)
                                                   
                                                   
                                                   rain_data %>% head(20) %>% as.data.frame()
                                                   
                                                   ```
                                                   
                                                   <br>
                                                     
                                                     A brief, visual check whether small `n_years` results in lower mean maximum precipitation: 
                                                     
                                                     ```{r precip_plot, echo = T, fig.align='center'}
                                                   
                                                   rain_data %>% 
                                                     ggplot(aes(x = n_year, y = mean_max_p_mm)) +
                                                     geom_point() +
                                                     theme_presi() +
                                                     labs(x = "Years with data (n)", y = "Av. annual max. P (mm)")
                                                   
                                                   ```
                                                   
                                                   **Conclusion:** data does not appear to have a bias toward small `n_years`
                                                   
                                                   ---
                                                     
                                                     <br>
                                                     
                                                     ## Finalized data set
                                                     
                                                     The final data set will contain 
                                                   
                                                   - peak flow 
                                                   - select meta data
                                                   - [FEH catchment descriptor data](http://nrfa.ceh.ac.uk/winfap-feh-files)
                                                   - precipitation (mean, annual maximum precipitation)
                                                   
                                                   ```{r join_data, echo=T,rows.print = 5, cache = F}
                                                   
                                                   
                                                   # some wrangling
                                                   all_stations <- all_stations %>% 
                                                     select(id, country, catchmentArea,  peak_flow_cumecs = `maximum-gauging-flow`)
                                                   
                                                   
                                                   
                                                   # join by id
                                                   all_stations <- all_stations %>% 
                                                     left_join(cd3_data %>% 
                                                                 mutate(id = as.character(id)),
                                                               by = 'id') %>% 
                                                     left_join(rain_data %>% 
                                                                 mutate(id = as.character(id)), by = 'id') 
                                                   
                                                   all_stations %>% head(20) %>% as.data.frame()
                                                   
                                                   # select variables for modelling (i.e. after checking for multi-colinearity)
                                                   all_stations <- all_stations %>% 
                                                     select(id, country, catchmentArea, peak_flow_cumecs,
                                                            DPSBAR, FARL, FPEXT, SAAR,SPRHOST,PROPWET)
                                                   
                                                   all_stations <- all_stations[complete.cases(all_stations), ]
                                                   
                                                   # number of observations per country
                                                   all_stations %>% dplyr::count(country) %>% as.data.frame()
                                                   
                                                   ```
                                                   
                                                   ---
                                                     
                                                     <br>
                                                     
                                                     # Analyses: Statistical Modelling
                                                     
                                                     
                                                     This sections outlines the definition and application of statistical models to answer our questions (recall: difference between countries; what controls peak flow). 
                                                   
                                                   In a first step we will apply a simple linear model, and then build upon `tidyverse`'s data-handling capabilities to expand our analyses. The chosen modelling framework is mainly for illustrative purposes and can easily be exchanged with GLM, LME, NLS, GAM, etc.
                                                   
                                                   ---
                                                   
                                                   <br>
                                                   
                                                   ## Q1: Does size matter?
                                                   
                                                   **Recall:** 
                                                   
                                                   Q1: What impact does catchment area have on peak flows, and does it differ between countries?
                                                   
                                                   ---
                                                   
                                                   ```{r resid_plot, fig.width = 8, fig.align='center', fig.height=5, echo = T}
                                                   
                                                   library(modelr)
                                                   
                                                   # transform peak flow and catchment area with log10 in new object
                                                   p_data <- all_stations %>% 
                                                   mutate(peak_flow_cumecs_log = log10(peak_flow_cumecs),
                                                   ca_log = log10(catchmentArea)) %>% 
                                                   filter(is.finite(peak_flow_cumecs_log),
                                                   is.finite(ca_log))
                                                   
                                                   # simple lm 
                                                   catchment_area_country.mod <- lm(peak_flow_cumecs_log ~ ca_log + country, data = p_data) 
                                                   catchment_area.mod <- lm(peak_flow_cumecs_log ~ ca_log, data = p_data)
                                                   
                                                   anova(catchment_area.mod, catchment_area_country.mod)
                                                   
                                                   
                                                   # check summary
                                                   catchment_area_country.mod %>% summary()
                                                   
                                                   
                                                   
                                                   # add residuals to new object
                                                   countries.resid <- p_data %>% 
                                                   add_residuals(catchment_area_country.mod, var = "resid_ca_c")
                                                   
                                                   # Quick peak at the residuals (don't forget QQ and Leverage Plots..)
countries.resid %>% 
  ggplot(aes(x = ca_log, y = resid_ca_c, col = country)) +
  geom_ref_line(h = 0, colour = "gray60", size = 1) +
  geom_point(size = 2, alpha = .6) +
  theme_presi() + 
  facet_wrap(~country)

```

**Conclusion:** We compared two models, and found that including country enhances our ability to represent the relationship between peak flows and catchment area: 
  
  - Peak flow increases approximately by an order of magnitude ($\beta_1$ = `r coef(catchment_area_country.mod)[2] %>% round(2)`) as we increase our catchment area by an order of magnitude (recall our  $log_{10}$ transformation)
- This relationship holds for all countries, but peak flows at a given catchment area tend to be largest for Wales ( > Scotland > Northern Irealand > England)

```{r tidy_mod_summary, echo = T}
catchment_area_country.mod %>% broom::tidy()

```


## Q2: Catchment characteristics

**Recall:**
  
  Q2: Which other characteristics affect peak flows on catchment-level, and which set of characteristics make up  the most representative models?
  
  ---
  
  Understanding peak flows is crucial from a management perspective (e.g. risk, resources). Additonal insight may be gleaned by applying "region" specific models individually. To this end, we apply a series of models to a nested data frame (by country).

Chosen catchment characteristics:
  
  - **DPSBAR**: Overall catchment steepness (mean Drainage Path Slope)
- **FARL**: Flood Attenuation by Reservoirs and Lakes (< 0.8 strong attenuation)
- **FPEXT**: Floodplain extent
- **PROPWET**: proportion of time catchment soils are wet
- **SAAR**: Average annual rainfall in the standard period (1961-1990) in millimetres
- **SPRHOST**: Standard percentage runoff (%), weighted by soil class across catchment

```{r alt_resid_plot, echo=F, eval=F}
countries.resid %>% 
  ggplot(aes(x = resid, fill = country)) +
  geom_ref_line(v = 0, colour = "black", size = 0.5) +
  geom_histogram() +
  theme_presi() + 
  facet_wrap(~country)


summary(catchment_area.mod)


```


```{r nested_models, echo = T, print.rows = 10, cache = F}



uk_nested <- p_data %>% 
  tidyr::nest(-country)


models <- list(
  null_intercept = function(x){lm(peak_flow_cumecs_log ~ 1, data = x)},
  ca = function(x){lm(peak_flow_cumecs_log ~ ca_log, data = x)},
  ca_meanMax_p = function(x){lm(peak_flow_cumecs_log ~ 
                                  ca_log + mean_max_p_mm, data = x)},
  ca_farl = function(x){lm(peak_flow_cumecs_log ~ 
                             ca_log + FARL, data = x)},
  ca_sprhost = function(x){lm(peak_flow_cumecs_log ~ 
                                ca_log + SPRHOST, data = x)},
  ca_fpext = function(x){lm(peak_flow_cumecs_log ~ 
                              ca_log + FPEXT, data = x)},
  ca_propwet = function(x){lm(peak_flow_cumecs_log ~ 
                                ca_log + PROPWET, data = x)},
  ca_farl_fpext = function(x){lm(peak_flow_cumecs_log ~ 
                                   ca_log + FARL + FPEXT, data = x)},
  ca_farl_propwet = function(x){lm(peak_flow_cumecs_log ~ 
                                     ca_log + FARL + PROPWET, data = x)},
  ca_farl_propwet_fpext = function(x){lm(peak_flow_cumecs_log ~ 
                                           ca_log + FARL + PROPWET + FPEXT, data = x)},
  
  ca_farl_fpext_sprhost_propwet = function(x){lm(peak_flow_cumecs_log ~ 
                                                   ca_log + FARL + FPEXT + SPRHOST + PROPWET, data = x)}
  
)

apply_model <- function(.model, ndf){
  
  ndf$model <- map(ndf$data, possibly(.model, NULL))
  
  return(ndf)
}



```



```{r apply_models, echo = T, rows.print = 7, cache = F}

uk_nest <- models %>% 
  map_df(apply_model, uk_nested, .id = "id_model") %>% 
  
  select(id_model, country, model) %>% 
  mutate(coefficients = map(model, tidy),
         performance = map(model, glance)) %>% 
  select(-model)

uk_nest %>% head(20)


```


```{r coef_plot_func, echo = F}
coef_plot <- function(){
  
  
  
  coef.plot <- uk_coefficients %>% 
    ggplot(aes(y = estimate,
               x = term, fill = estimate_factor)) +
    geom_bar(stat = "identity", position = "dodge",col = "gray20") +
    coord_flip() +
    facet_grid(country~id_model) +
    facet_grid(id_model~country, scales = "free_y") +
    
    geom_hline(yintercept = 0, linetype = 2, col = "darkgrey") +
    theme_bw() +
    theme(panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          legend.position = "top",
          axis.text.x = element_text(angle = 0)) +
    scale_fill_manual(name = "Estimate", values = c("negative" = "darkred",
                                                    "positive" = "steelblue",
                                                    "not sign." = "grey95"))
  return(coef.plot)
}



```



```{r model_assessment, echo=T, fig.height = 8, fig.align='center', fig.width = 10, cache = F}

# unnest data frame, dropping other nested column
uk_coefficients <- uk_nest %>%
  unnest(coefficients, .drop = T)

# add column with adjusted coefficients, useful for plotting, harder to interpret
# add columns with factors for plotting/classification
uk_coefficients <- uk_coefficients %>% 
  mutate(estimate_adj = scale(estimate, center = T, scale = T),
         estimate_adj_factor = ifelse(p.value < 0.05,
                                      ifelse(estimate > 0,
                                             "positive",
                                             "negative"),
                                      "not sign.") %>% 
           fct_relevel("positive", "negative"),
         estimate_factor = ifelse(p.value < 0.05,
                                  ifelse(estimate > 0,
                                         "positive",
                                         "negative"),
                                  "not sign.") %>% 
           fct_relevel("positive", "negative"),
         id_model = fct_relevel(as.factor(id_model), "null_intercept"))


# custom function based on gglplot2, code available in raw document
coef_plot()


```


```{r aic_plot_func, echo=F}

aic_plot <- function(){
  aic.plot <-  uk_performance %>% 
    ggplot(aes(x = id_model, y = dAIC, group = 1)) +
    geom_line(linetype = 2) +
    geom_point(aes(fill = uk_performance$minAIC), size = 5, shape = 21, col = 'gray20') +
    stat_summary(geom = "text",
                 fun.y = function(x){min(x,na.rm=T)*1.3  + 2},
                 label = round(uk_performance$adj.r.squared,3)) +
    facet_grid(country~., scales = 'free_y') +
    theme_presi() +
    theme(panel.grid.major.y = element_line(color = "gray90"),
          panel.border = element_blank(),
          panel.spacing = unit(2, "lines"),
          axis.text.x = element_text(angle = 45, hjust = 1),
          plot.title = element_text(hjust = 0),
          plot.subtitle = element_text(colour = 'gray60'),
          legend.position = "top") +
    scale_y_continuous(expand = c(0.1,0)) +
    scale_fill_manual(name = "AIC comparison", 
                      values = c("less support/discard" = "gray40",
                                 "keep" = "steelblue")) +
    labs(x = "Model", y = "delta AIC") +
    ggtitle(label = "Model selection", subtitle = "via relative and absolute goodness of fit")
  # geom_bar(stat = "identity", position = 'dodge', col = 'gray20') +
  # geom_text(data = uk_aic_min,
  #           aes(x = id_model,
  #               y = minAIC_val,
  #               label = .97*round(uk_assess$adj.r.squared,2))) +
  return(aic.plot)
}

```



```{r aic_plot, echo = T, fig.width=10, fig.height=10, cache = F}


uk_performance <- uk_nest %>%
  tidyr::unnest(performance, .drop = T)

uk_performance <- uk_performance %>%
  group_by(country) %>% 
  arrange(AIC, .by_group = T) %>% 
  rename(AIC_val = AIC) %>% 
  mutate(dAIC = AIC_val - min(AIC_val),
         minAIC = ifelse(dAIC < 10 ,"keep", "less support/discard")) %>% 
  filter(id_model != "null_intercept")


# custom function based on gglplot2, code available in raw document
aic_plot()

```



<!-- <script> -->
  <!--   var d = document.document.getElementsByTagName("img"); -->
  <!--   d.className += " img-responsive"; -->
  <!-- </script> -->