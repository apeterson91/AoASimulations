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
#' @return list of two table components table1_top -
#' which includes the coverage and diagnostic statistics broken down by parameter and
#' table1_bottom which has the RMSE for the entire model
#' @export
create_table_one <- function(num_sims = 5,
                             n = 300,
                             num_MESA_subj = 100,
                             d = 20,
                             alpha = 25,
                             theta = .5, delta = -2.2,
                             beta = .75, beta_bar = .85,
                             sigma = 1,
                             alpha_prior = rstap::normal(location = 25,autoscale =F),
                             beta_prior = rstap::normal(),
                             theta_prior = rstap::log_normal(location = 1 , scale = 1),
                             delta_prior = rstap::normal(),
                             iter = 2E3,
                             warmup = 1E3,
                             cores = 1,
                             chains = 1,
                             seed = NULL,
                             file = NULL){
    if(is.null(seed))
        set.seed(24321)

    ## HPP
    datasets <- purrr::map(1:num_sims,function(x) generate_hpp_dataset(seed = x,
                                                                       num_subj = n,
                                                                       num_dists = d,
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
                                                                             num_subj_sim = 1,
                                                                             num_bef_sim = 1,
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
                                                                               num_subj_sim = 1,
                                                                               num_bef_sim = 1,
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
                                                                             num_subj = num_MESA_subj,
                                                                             alpha = alpha,
                                                                             beta = beta,
                                                                             delta = delta,
                                                                             beta_bar = beta_bar,
                                                                             theta = theta,
                                                                             sigma = sigma,
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
        Hmisc::latex(tab1a,file = file, booktabs = T )
        Hmisc::latex(tab1b,file = file, booktabs = T )
    }
    else{
        Hmisc::latex(tab1a,file= "~/Desktop/tab1a.tex", booktabs = T)
        Hmisc::latex(tab1b,file = "~/Desktop/tab1b.tex",booktabs = T)
    }
    print(tab1a)
    print(tab1b)

    return(list(top_table1 = tab1a,bottom_table1 = tab1b))
}

#' create_table_two
#'
#' @param num_sims
#'
#' @export
create_table_two <- function(num_sims = 5,
                             num_subj = 300,
                             num_dists = 20,
                             alpha = 25,
                             theta =  .5,
                             shape = 5,
                             delta = -2.2,
                             beta = .75,
                             cores = 2){

    # DLM

    ddatasets <- purrr::map(1:num_sims,function(x) generate_dlm_dataset(seed = x,))

    # Exponential
    edatasets <- purrr::map(1:num_sims,function(x) generate_hpp_dataset(seed = x,
                                                                num_subj = num_subj,
                                                                num_dists = num_dists,
                                                                alpha = alpha,
                                                                theta = theta,
                                                                delta = delta,
                                                                beta = beta,
                                                                W = function(x) exp(-x)))

    # Weibull

    wdatasets <- purrr::map(1:num_sims, function(x) generate_hpp_dataset(seed = x,
                                                                 num_subj = num_subj,
                                                                 num_dists = num_dists,
                                                                 alpha = alpha,
                                                                 theta = theta,
                                                                 shape = shape,
                                                                 delta = delta,
                                                                 beta = beta))
    # DLM - 1




    # DLM - 2

    # Model fitting

    # Fit Exponential under Exponential

    ee <- purrr::map(edatasets,function(x){
        rstap::stap_glm(outcome~sex + sap(FF,exp),
                subject_data = x$subject_data,
                distance_data = x$bef_data,
                max_distance = 10,
                subject_ID = "subj_id",
                prior = rstap::normal(),
                prior_stap = rstap::normal(),
                prior_intercept = rstap::normal(location =25),
                prior_theta = rstap::log_normal(1,1))})


    # Fit Exponential under Weibull

    ew <- purrr::map(edatasets,function(x){
        rstap::stap_glm(outcome~sex + sap(FF,wei),
                subject_data = x$subject_data,
                distance_data = x$bef_data,
                max_distance = 10,
                subject_ID = "subj_id",
                prior = rstap::normal(),
                prior_stap = rstap::normal(),
                prior_intercept = rstap::normal(location =25),
                prior_theta = rstap::log_normal(1,1))})

    # Fit Weibull under Exponential

    we <- purrr::map(wdatasets,function(x){
        rstap::stap_glm(outcome~sex + sap(FF,exp),
                subject_data = x$subject_data,
                distance_data = x$bef_data,
                max_distance = 10,
                subject_ID = "subj_id",
                prior = rstap::normal(),
                prior_stap = rstap::normal(),
                prior_intercept = rstap::normal(location =25),
                prior_theta = rstap::log_normal(1,1),
                cores = cores)})

    # Fit Weibull under Weibull

    ww <- purrr::map(wdatasets,function(x){
        rstap::stap_glm(outcome ~ sex + sap(FF,wei),
                subject_data = x$subject_data,
                distance_data = x$bef_data,
                max_distance = 10,
                subject_ID = "subj_id",
                prior = rstap::normal(),
                prior_stap = rstap::normal(),
                prior_intercept = rstap::normal(location = 25),
                prior_theta = rstap::log_normal(1,1))})

    term_exp  <- tibble::tibble(sim_id = 1:num_sims,
                                Simulated_Function = rep("Exponential",num_sims),
                                Modeled_Function = rep("Exponential",num_sims),
                                True_termination = uniroot(function(x) exp(-(x/theta)) - 0.05,interval = c(0,10))$root,
                                Stap_Termination = purrr::map_dbl(ee,function(a) rstap::stap_termination(a,max_value=100)[2])
    )

    term_expwei <- tibble::tibble(sim_id = 1:num_sims,
                                  Simulated_Function = rep("Exponential",num_sims),
                                  Modeled_Function = rep("Weibull",num_sims),
                                  True_termination = uniroot(function(x) exp(-(x/theta)) - 0.05,interval = c(0,10))$root,
                                  Stap_Termination = purrr::map_dbl(ee,function(a) rstap::stap_termination(a,max_value=100)[2])
    )

    term_weiexp <- tibble::tibble(sim_id = 1:num_sims,
                                  Simulated_Function = rep("Weibull",num_sims),
                                  Modeled_Function = rep("Exponential",num_sims),
                                  True_termination = uniroot(function(x){ exp(-(x/theta)^shape) - 0.05},interval = c(0,10))$root,
                                  Stap_Termination = purrr::map_dbl(ee,function(a) rstap::stap_termination(a,max_value=100)[2])
    )

    term_wei <- tibble::tibble(sim_id = 1:num_sims,
                                        Simulated_Function = rep("Weibull",num_sims),
                                        Modeled_Function = rep("Weibull",num_sims),
                                        True_termination = uniroot(function(x){ exp(-(x/theta)^shape) - 0.05},interval = c(0,10))$root,
                                        Stap_Termination = purrr::map_dbl(ww,function(a) rstap::stap_termination(a,max_value=100)[2])
    )
    out <- dplyr::bind_rows(term_exp,term_expwei,term_weiexp,term_wei) %>%
        dplyr::mutate(Termination_Difference=abs(True_termination - Stap_Termination)) %>%
        dplyr::group_by(Simulated_Function,Modeled_Function) %>%
        dplyr::summarise(mean_difference = mean(Termination_Difference)) %>%
        dplyr::ungroup() %>%
        tidyr::spread(Modeled_Function,mean_difference)

    as.tabular(out)
}


