#' create_table_one in AoAS paper
#'
#' @param num_sims number of simulations to run
#' @param n number of subjects to simulate
#' @param d number of distances to simulate
#' @param alpha intercept for generating outcome
#' @param theta true spatial scale under which datasets
#' @param delta simulated binary covariate regression effect
#' @param beta SAP effect
#' @param alpha_prior prior to be placed on intercept in model
#' @param beta_prior prior to be placed on SAP effect
#' @param theta_prior prior to be placed on spatial scale
#' @param delta_prior prior to be placed on simulated binary covariate effect
#' @param iter number of iterations for which to run the stap_glm or stapdnd_glmer sampler
#' @param warmup number of iterations to warmup the sampler
#' @param chains number of independent MCMC chains to draw
#' @param cores number of cores with which to run chains in parallel
#' @param file path to file to save tables to in .tex format
#' @return list of two table components table1_top -
#' which includes the coverage and diagnostic statistics broken down by parameter and
#' table1_bottom which has the RMSE for the entire model
#' @export
create_table_one <- function(num_sims = 5,
                             num_subj = 100L,
                             num_mesa_subj = 50L,
                             num_dists = 30,
                             alpha = 23,
                             theta = .5, delta = -2.2,
                             beta = .75, beta_bar = .85,
                             alpha_prior = rstap::normal(location = 25,autoscale =F),
                             beta_prior = rstap::normal(),
                             theta_prior = rstap::log_normal(location = 1 , scale = 1),
                             delta_prior = rstap::normal(),
                             iter = 2E3,
                             warmup = 1E3,
                             chains = 1,
                             cores = 1,
                             file = NULL){

    ## HPP
    datasets <- purrr::map(1:num_sims,function(x) generate_hpp_dataset(seed = x,
                                                                       num_subj = num_subj,
                                                                       num_dists = num_dists,
                                                                       alpha = alpha,
                                                                       theta = theta,
                                                                       delta = delta,
                                                                       beta = beta))


    models <- purrr::map(datasets,function(x){
        rstap::stap_glm(outcome~sex + sap(FF,exp),
                        subject_data = x$subject_data,
                        distance_data = x$bef_data,
                        max_distance = 10,
                        subject_ID = "subj_id",
                        prior = delta_prior,
                        prior_stap = beta_prior,
                        prior_intercept = alpha_prior,
                        prior_theta = theta_prior,
                        chains = chains,
                        cores = cores,
                        iter = iter,
                        warmup = warmup)})

    hpp_output <- tibble::tibble(simulation = rep(1:num_sims,2),
                                 `Spatial Pattern` = factor(rep("HPP",num_sims*2)),
                                 Parameter = c(rep("Beta",num_sims),rep("Theta",num_sims)),
                                 coverage = c(purrr::map_dbl(models,function(x) check_coverage(x,c("FF"=beta))),
                                              purrr::map_dbl(models,function(x) check_coverage(x,c("FF_spatial_scale"=theta)))),
                                 `Cook & Gelman` =  c(purrr::map_dbl(models,function(x) calculate_CG_stat(x,c("FF"=beta))),
                                                      purrr::map_dbl(models,function(x) calculate_CG_stat(x,c("FF_spatial_scale"=theta)))),
                                 interval_length = c(purrr::map_dbl(models,function(x) interval_length(x,c("FF"))),
                                                     purrr::map_dbl(models,function(x) interval_length(x,c("FF_spatial_scale"))))) %>%
        dplyr::mutate(Parameter = factor(Parameter))

    hpp_output2 <- tibble::tibble(simulation = rep(1:num_sims,1),
                                  `Spatial Pattern` = factor(rep("HPP",num_sims)),
                                  RMSE = purrr::map_dbl(models,calculate_RMSE_median))

    ## NHPP

    nhpp_datasets <- purrr::map(1:num_sims,function(x) generate_nhpp_dataset(seed = x,
                                                                             num_subj = num_subj,
                                                                             num_dists = num_dists,
                                                                             alpha = alpha,
                                                                             theta = theta,
                                                                             delta = delta,
                                                                             beta = beta))

    nhpp_models <- purrr::map(nhpp_datasets,function(x){
        rstap::stap_glm(outcome~sex + sap(FF,exp),
                        subject_data = x$subject_data,
                        distance_data = x$bef_data,
                        max_distance = 10,
                        subject_ID = "subj_id",
                        prior = rstap::normal(),
                        prior_stap = rstap::normal(),
                        prior_intercept = rstap::normal(location =25),
                        prior_theta = rstap::log_normal(1,1),
                        chains = chains,
                        cores = cores,
                        iter = iter,
                        warmup = warmup)})

    nhpp_output <- tibble::tibble(simulation = rep(1:num_sims,2),
                                  `Spatial Pattern` = factor(rep("NHPP",num_sims*2)),
                                  Parameter = c(rep("Beta",num_sims),rep("Theta",num_sims)),
                                  coverage =c(purrr::map_dbl(nhpp_models,function(x) check_coverage(x,c("FF"=beta))),
                                              purrr::map_dbl(nhpp_models,function(x) check_coverage(x,c("FF_spatial_scale"=theta)))),
                                  `Cook & Gelman` =  c(purrr::map_dbl(nhpp_models,function(x) calculate_CG_stat(x,c("FF"=beta))),
                                                       purrr::map_dbl(nhpp_models,function(x) calculate_CG_stat(x,c("FF_spatial_scale"=theta) ))),
                                  interval_length = c(purrr::map_dbl(nhpp_models,function(x) interval_length(x,c("FF"))),
                                                      purrr::map_dbl(nhpp_models,function(x) interval_length(x,c("FF_spatial_scale"))))) %>%
        dplyr::mutate(Parameter = factor(Parameter))


    nhpp_output2 <- tibble::tibble(simulation = rep(1:num_sims,1),
                                   `Spatial Pattern` = factor(rep("NHPP",num_sims)),
                                   RMSE = purrr::map_dbl(nhpp_models,calculate_RMSE_median))

    ## Matern

    matern_datasets <- purrr::map(1:num_sims,function(x) generate_matern_dataset(seed = x,
                                                                                 num_subj = num_subj,
                                                                                 num_dists = num_dists,
                                                                                 alpha = alpha,
                                                                                 theta = theta,
                                                                                 delta = delta,
                                                                                 beta = beta))

    matern_models <- purrr::map(matern_datasets,function(x){
        rstap::stap_glm(outcome~sex + sap(FF,exp),
                        subject_data = x$subject_data,
                        distance_data = x$bef_data,
                        max_distance = 10,
                        subject_ID = "subj_id",
                        prior = rstap::normal(),
                        prior_stap = rstap::normal(),
                        prior_intercept = rstap::normal(location =25),
                        prior_theta = rstap::log_normal(1,1),
                        chains = chains,
                        cores = cores,
                        iter = iter,
                        warmup = warmup)})

    matern_output <- tibble::tibble(simulation = rep(1:num_sims,2),
                                    `Spatial Pattern` = factor(rep("Matern",num_sims*2)),
                                    Parameter = c(rep("Beta",num_sims),rep("Theta",num_sims)),
                                    coverage =c(purrr::map_dbl(matern_models,function(x) check_coverage(x,c("FF"=beta))),
                                                purrr::map_dbl(matern_models,function(x) check_coverage(x,c("FF_spatial_scale"= theta)))),
                                    `Cook & Gelman` =  c(purrr::map_dbl(matern_models,function(x) calculate_CG_stat(x,c("FF"=beta) )),
                                                         purrr::map_dbl(matern_models,function(x) calculate_CG_stat(x,c("FF_spatial_scale"=theta)))),
                                    interval_length = c(purrr::map_dbl(matern_models,function(x) interval_length(x,c("FF"))),
                                                        purrr::map_dbl(matern_models,function(x) interval_length(x,c("FF_spatial_scale"))))) %>%
        dplyr::mutate(Parameter = factor(Parameter))


    matern_output2 <- tibble::tibble(simulation = 1:num_sims,
                                     `Spatial Pattern` = factor(rep("Matern",num_sims)),
                                     RMSE = purrr::map_dbl(matern_models,calculate_RMSE_median))

    ## MESA

    MESA_datasets <- purrr::map(1:num_sims,function(x) generate_mesa_dataset(seed = x,
                                                                             num_subj = num_mesa_subj,
                                                                             alpha = alpha,
                                                                             theta = theta,
                                                                             delta = delta,
                                                                             W = function(x) exp(-x)))

    MESA_models <- purrr::map(MESA_datasets,function(x){
        rstap::stapdnd_glmer(outcome~sex + sap_dnd_bar(FF,exp) + (visit_number|id),
                             subject_data = x$subject_data,
                             distance_data = x$bef_data,
                             max_distance = 10,
                             subject_ID = "id",
                             group_ID = "visit_number",
                             prior = rstap::normal(),
                             prior_stap = rstap::normal(),
                             prior_intercept = rstap::normal(location = 25),
                             prior_theta = rstap::log_normal(1,1),
                             chains = chains,
                             cores = cores,
                             iter = iter,
                             warmup = warmup)})


    MESA_output <- tibble::tibble(simulation = rep(1:num_sims,2),
                                  `Spatial Pattern` = factor(rep("MESA data",num_sims*2)),
                                  Parameter = c(rep("Beta",num_sims),rep("Theta",num_sims)),
                                  coverage =c(purrr::map_dbl(MESA_models,function(x) check_coverage(x,c("FF_dnd"=beta))),
                                              purrr::map_dbl(MESA_models,function(x) check_coverage(x,c("FF_spatial_scale"=theta)))),
                                  `Cook & Gelman` =  c(purrr::map_dbl(MESA_models,function(x) calculate_CG_stat(x,c("FF_dnd"=beta),TRUE)),
                                                       purrr::map_dbl(MESA_models,function(x) calculate_CG_stat(x,c("FF_spatial_scale"=theta),TRUE))),
                                  interval_length = c(purrr::map_dbl(MESA_models,function(x) interval_length(x,c("FF_dnd"))),
                                                      purrr::map_dbl(MESA_models,function(x) interval_length(x,c("FF_spatial_scale"))))) %>%
        dplyr::mutate(Parameter = factor(Parameter))


    MESA_output2 <- tibble::tibble(simulation = 1:num_sims,
                                   `Spatial Pattern` = factor(rep("MESA data",num_sims)),
                                   RMSE = purrr::map_dbl(MESA_models,calculate_RMSE_median))



    output <- rbind(hpp_output,nhpp_output,matern_output,MESA_output)

    output2 <- rbind(hpp_output2,nhpp_output2,matern_output2,MESA_output2)

    require(tables)
    tab1a <- tabular( Format(digits=2)*
                          (interval_length + coverage)*(Parameter)*(mean + sd) +
                          Format(digits=2)*(`Cook & Gelman`)*(Parameter)*(sum) ~ (`Spatial Pattern`) ,
                      data = output )

    tab1b <- tabular( Format(digits=2)*(RMSE)*(mean +sd) ~ (`Spatial Pattern`),
                      data=output2)



    if(!is.null(file)){
        Hmisc::latex(tab1a, file = file, booktabs = T )
        Hmisc::latex(tab1b, file = file, booktabs = T )
    }
    else{
        Hmisc::latex(tab1a,file= "~/Desktop/tab1a.tex", booktabs = T)
        Hmisc::latex(tab1b,file = "~/Desktop/tab1b.tex",booktabs = T)
    }
    print(tab1a)
    print(tab1b)

    return(list(top_table1 = tab1a,bottom_table1 = tab1b))
}
