create_table_one <- function(num_sims = 5, n = 300, d = 20, alpha = 25, theta = .5, delta = -2.2, beta = .75){
    require(tables)
    require(purrr)
    datasets <- map(1:num_sims,function(x) generate_hpp_dataset(seed = x,
                                                                       num_subj = 300,
                                                                       num_dists = 20,
                                                                       alpha = 25,
                                                                       theta = .5,
                                                                       delta = -2.2,
                                                                       beta = .75))
    
    models <- map(datasets,function(x){ 
        stap_lm(outcome~sex + sap(FF,exp),
                subject_data = x$subject_data,
                distance_data = x$bef_data,
                max_distance = 5,
                subject_ID = "subj_id", 
                prior = normal(),
                prior_stap = normal(),
                prior_intercept = normal(location =25),
                prior_theta = log_normal(1,1))})
    
    hpp_output <- tibble(simulation = rep(1:num_sims,2),
                         `Spatial Pattern` = factor(rep("HPP",num_sims*2)),
                         Parameter = c(rep("Beta",num_sims),rep("Theta",num_sims)),
                         coverage = c(map_dbl(models,function(x) check_coverage(x,c("FF"=.75))),
                                      map_dbl(models,function(x) check_coverage(x,c("FF_spatial_scale"=.5)))),
                         `Cook & Gelman` = (rnorm(10))^2,
                         interval_length = c(map_dbl(models,function(x) interval_length(x,c("FF"))),
                                             map_dbl(models,function(x) interval_length(x,c("FF_spatial_scale"))))) %>% 
        mutate(Parameter = factor(Parameter))
    hpp_output2 <- tibble(simulation = rep(1:num_sims,1),
                          `Spatial Pattern` = factor(rep("HPP",num_sims)),
                          RMSE = map_dbl(models,calculate_RMSE_median))
    tab1 <- tabular( Format(digits=2)*
                        (interval_length + coverage)*(Parameter)*(mean + sd) + 
                         Format(digits=2)*(`Cook & Gelman`)*(Parameter)*(sum) ~ (`Spatial Pattern`) , data=hpp_output )
    tab2 <- tabular( (RMSE)*(mean +sd) ~ (`Spatial Pattern`),data=hpp_output2)
    print(tab1)
    print(tab2)
}

create_table_two(num_sims = 5){
    
    datasets <- map(1:num_sims,function(x) generate_hpp_dataset(seed = x,
                                                                num_subj = 300,
                                                                num_dists = 20,
                                                                alpha = 25,
                                                                theta = .5,
                                                                delta = -2.2,
                                                                beta = .75))
    
}
    
