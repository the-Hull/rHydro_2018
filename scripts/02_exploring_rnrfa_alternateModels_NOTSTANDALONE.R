# uk_nested_ca <- uk_nested %>% 
#       mutate(data = map(data, 
#                         ~mutate(.x, 
#                                 max_flow_mm = (`maximum-gauging-flow` * 60 * 60 * 24) / catchmentArea * 1e6))
#             
#       )

all_stations <- all_stations %>% 
      mutate(max_flow_mm = (`maximum-gauging-flow` * 60 * 60 * 24) / (catchmentArea * 1e6))

load("./dat/rain_data_full.Rda")
names(rain_data) <- all_stations$id

rain_data_av <- rain_data %>% 
 
      set_names(all_stations$id) %>% 
      discard(is.null) %>% 
      map(~tibble::tibble(date_time = zoo::index(.x) %>% as.POSIXct(),
                          year = lubridate::year(date_time),
                          month = lubridate::month(date_time),
                          precipitation_mm_month = zoo::coredata(.x)[ ,1])) %>% 
      map2(names(.), 
           ~mutate(.x, id = as.factor(.y))
      ) %>% 
      map(filter, year >= "1990") %>% 
      keep(.p = ~length(unique(.x$year)) > 20) %>% 
      map_df(function(x){
           x %>%  group_by(id) %>% 
                  summarise(n_years = length(unique(year)),
                            precipitation_mm_mean = mean(precipitation_mm_month, na.rm = T))
      })

uk_nested <- all_stations %>% 
      left_join(cd3_data %>% mutate(id = as.character(id)), by = "id") %>% 
      left_join(rain_data_av, by = 'id') %>% 
      filter(is.finite(max_flow_mm)) %>% 
      tidyr::nest(-country)


models <- list(
      null_intercept = function(x){lm(max_flow_mm ~ 1, data = x)},
      dplbar_bfihost = function(x){lm(max_flow_mm ~
                                               DPLBAR + BFIHOST, data = x)},
      dplbar_bfihost_farl = function(x){lm(max_flow_mm ~
                                                    DPLBAR + BFIHOST + FARL, data = x)},
      dpsbar_bfihost_farl_fpext = function(x){lm(max_flow_mm ~
                                                          DPSBAR + BFIHOST + FARL, data = x)},
      saar_dpsbar_bfihost_farl = function(x){lm(max_flow_mm ~
                                                     SAAR + DPSBAR + BFIHOST + FARL, data = x)},
      pmean_dbspar_bfihost_farl = function(x){lm(max_flow_mm ~
      precipitation_mm_mean + DPSBAR + BFIHOST + FARL, data = x)}
      
      
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
                   fct_relevel("positive", "negative"),
             id_model = fct_relevel(as.factor(id_model), "null_intercept"))


param.plot <- uk_param %>% 
      ggplot(aes(y = estimate_adj,
                 x = term, fill = estimate_adj_factor)) +
      geom_bar(stat = "identity", position = "dodge",col = "gray20") +
      coord_flip() +
      facet_grid(country~id_model) +
      facet_grid(id_model~country, scales = "free_y") +
      
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
                   label = round(uk_assess$adj.r.squared,3)) +
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
