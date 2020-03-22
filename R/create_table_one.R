#' create datasets for sample size for first simulation figure
#'
#' @param num_sims number of simulations to run
#' @param num_subj_processes number of subject processes to simulate
#' @param num_dists_processes number of distance processes to simulate
#' @param num_mesa_subj number of subjects to sample from MESA data for MESA analysis
#' @param alpha intercept for generating outcome
#' @param theta true spatial scale under which datasets
#' @param delta simulated binary covariate regression effect
#' @param beta SAP effect
#' @param alpha_prior prior to be placed on intercept in model, must be an rstap:: namespace object
#' @param beta_prior prior to be placed on SAP effect
#' @param theta_prior prior to be placed on spatial scale
#' @param delta_prior prior to be placed on simulated binary covariate effect
#' @param skip_MESA Boolean value that indicates whether to run MESA model or not
#' @param iter number of iterations for which to run the stap_glm or stapdnd_glmer sampler
#' @param warmup number of iterations to warmup the sampler
#' @param chains number of independent MCMC chains to draw
#' @param cores number of cores with which to run chains in parallel
#' @param file path to file to save tables to in .tex format
#' @return list of two table components table1_top -
#' which includes the coverage and diagnostic statistics broken down by parameter and
#' The remaining "raw" table component contains the pre-aggregation data-frames
#' @export
create_table_one <- function(num_sims = 5,
                             num_subj_processes = 11L,
                             num_dist_processses = 4L,
                             num_mesa_subj = 50L,
                             alpha = 23,
                             theta = .5, delta = -2.2,
                             beta = .75, beta_bar = .85,
                             sigma = 1,
                             alpha_prior = rstap::normal(location = 25, scale = 4, autoscale =F),
                             beta_prior = rstap::normal(location = 0, scale = 3, autoscale = F),
                             theta_prior = rstap::log_normal(location = 0, scale = 1),
                             delta_prior = rstap::normal(location = 0, scale = 3, autoscale = F),
                             skip_MESA = FALSE,
                             iter = 4E3,
                             warmup = 2E3,
                             chains = 1,
                             cores = 1,
                             seed = NULL,
                             file = NULL){


    if(is.null(seed)){
        seed <- 34234
        set.seed(seed)
    }

    ## NHPP
    ## num subj and num_dists are determined randomly
    nhpp_dataset <- generate_nhpp_table_one(seed = seed,
                                       alpha = alpha,
                                       theta = theta,
                                       delta = delta,
                                       beta = beta,
                                       num_subj_realizations = num_subj_processes,
                                       num_dist_realizations = num_dist_processses)

    nhpp_models <- purrr::map(1:num_sims,function(x){
        subject_data <- nhpp_dataset$subject_data
        subject_data$outcome <- subject_data$outcome + rnorm(n = nrow(subject_data),mean = 0,sd = sigma)
        rstap::stap_glm(outcome~sex + sap(FF,exp),
                        subject_data = subject_data,
                        distance_data = nhpp_dataset$bef_data,
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

    nhpp_output <- tibble::tibble(simulation = rep(1:num_sims,2),
                                  `Spatial Pattern` = factor(rep("NHPP",num_sims*2)),
                                  Parameter = c(rep("Beta",num_sims),rep("Theta",num_sims)),
                                  coverage =c(purrr::map_dbl(nhpp_models,function(x) check_coverage(x,c("FF"=beta),nonzero = F)),
                                              purrr::map_dbl(nhpp_models,function(x) check_coverage(x,c("FF_spatial_scale"=theta),nonzero = F))),
                                  Difference = c(purrr::map_dbl(nhpp_models,function(x) abs_diff(x,c("FF"=beta)) ),
                                                 purrr::map_dbl(nhpp_models,function(x) abs_diff(x,c("FF_spatial_scale"=theta)))),
                                  `Cook & Gelman` =  c(purrr::map_dbl(nhpp_models,function(x) calculate_CG_stat(x,c("FF"=beta))),
                                                       purrr::map_dbl(nhpp_models,function(x) calculate_CG_stat(x,c("FF_spatial_scale"=theta) ))),
                                  interval_length = c(purrr::map_dbl(nhpp_models,function(x) interval_length(x,c("FF"))),
                                                      purrr::map_dbl(nhpp_models,function(x) interval_length(x,c("FF_spatial_scale"))))) %>%
        dplyr::mutate(Parameter = factor(Parameter))


    num_bef <- nhpp_dataset$bef_data %>% dplyr::group_by(subj_id) %>% dplyr::count() %>% dplyr::pull(n)
    num_bef <- num_bef[1]
    ## HPP
    dataset <- generate_hpp_table_one(seed = seed,
                                     num_subj = nrow(nhpp_dataset$subject_data),
                                     num_dists = num_bef,
                                     alpha = alpha,
                                     theta = theta,
                                     delta = delta,
                                     beta = beta)


    models <- purrr::map(1:num_sims,function(x){
        subject_data <- dataset$subject_data
        subject_data$outcome <- subject_data$outcome + rnorm(n = nrow(subject_data),mean = 0,sd = sigma)
        rstap::stap_glm(outcome~sex + sap(FF,exp),
                        subject_data = subject_data,
                        distance_data = dataset$bef_data,
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
                                 coverage = c(purrr::map_dbl(models,function(x) check_coverage(x,c("FF"=beta),nonzero = F)),
                                              purrr::map_dbl(models,function(x) check_coverage(x,c("FF_spatial_scale"=theta),nonzero = F))),
								 Difference = c(purrr::map_dbl(models,function(x) abs_diff(x,c("FF"=beta))),
												purrr::map_dbl(models,function(x) abs_diff(x,c("FF_spatial_scale"=theta)))),
                                 `Cook & Gelman` =  c(purrr::map_dbl(models,function(x) calculate_CG_stat(x,c("FF"=beta))),
                                                      purrr::map_dbl(models,function(x) calculate_CG_stat(x,c("FF_spatial_scale"=theta)))),
                                 interval_length = c(purrr::map_dbl(models,function(x) interval_length(x,c("FF"))),
                                                     purrr::map_dbl(models,function(x) interval_length(x,c("FF_spatial_scale"))))) %>%
        dplyr::mutate(Parameter = factor(Parameter))


    if(!skip_MESA){
        ## MESA

        MESA_datasets <- purrr::map(1:num_sims,function(x) generate_mesa_dataset(seed = x,
                                                                                 num_subj = num_mesa_subj,
                                                                                 alpha = alpha,
                                                                                 theta = theta,
                                                                                 delta = delta,
                                                                                 beta = beta,
                                                                                 beta_bar = beta_bar,
                                                                                 K = function(x) exp(-x)))


        MESA_models <- purrr::map(MESA_datasets,function(x){
            rstap::stapdnd_glmer(outcome~sex + sap_dnd_bar(FF,exp) + (visit_number|id),
                                 subject_data = x$subject_data,
                                 distance_data = x$bef_data,
                                 max_distance = 10,
                                 subject_ID = "id",
                                 group_ID = "visit_number",
                                 prior = delta_prior,
                                 prior_stap = beta_prior,
                                 prior_intercept = alpha_prior,
                                 prior_theta = theta_prior,
                                 chains = chains,
                                 cores = cores,
                                 iter = iter,
                                 warmup = warmup)
        })
        MESA_output <- tibble::tibble(simulation = rep(1:num_sims,2),
                                  `Spatial Pattern` = factor(rep("MESA data",num_sims*2)),
                                  Parameter = c(rep("Beta",num_sims),rep("Theta",num_sims)),
                                  coverage =c(purrr::map_dbl(MESA_models,function(x) check_coverage(x,c("FF_dnd"=beta))),
                                              purrr::map_dbl(MESA_models,function(x) check_coverage(x,c("FF_spatial_scale"=theta)))),
							   	  Difference = c(purrr::map_dbl(MESA_models,function(x) abs_diff(x,c("FF_dnd"=beta))),
												purrr::map_dbl(MESA_models,function(x) abs_diff(x,c("FF_spatial_scale"=theta)))),
                                  `Cook & Gelman` =  c(purrr::map_dbl(MESA_models,function(x) calculate_CG_stat(x,c("FF_dnd"=beta),TRUE)),
                                                       purrr::map_dbl(MESA_models,function(x) calculate_CG_stat(x,c("FF_spatial_scale"=theta),TRUE))),
                                  interval_length = c(purrr::map_dbl(MESA_models,function(x) interval_length(x,c("FF_dnd"))),
                                                      purrr::map_dbl(MESA_models,function(x) interval_length(x,c("FF_spatial_scale"))))) %>%
            dplyr::mutate(Parameter = factor(Parameter))


        output <- rbind(hpp_output,nhpp_output,MESA_output)

    }
    else{
        output <- rbind(hpp_output,nhpp_output)
    }

    standard_error <- function(x) sd(x)/sqrt(length(x))
    require(tables)
    tab1a <- tabular( Format(digits=2)*
                          (interval_length + coverage + Difference)*(Parameter)*(mean + standard_error) +
                          Format(digits=2)*(`Cook & Gelman`)*(Parameter)*(sum) ~ (`Spatial Pattern`) ,
                      data = output )



    if(!is.null(file)){
		if(length(file)==0)
			Hmisc::latex(tab1a,file= "~/Desktop/tab1a.tex", booktabs = T)
		else
			Hmisc::latex(tab1a,file= paste0(file,"_tab1a.tex"), booktabs = T)
    }
    print(tab1a)

    return(list(top_table1 = tab1a,
                top_raw_output = output,
                num_subj = nrow(dataset$subject_data),
                num_bef = num_bef
                ))
}
