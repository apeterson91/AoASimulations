
#' @param num_sims number of simulations to run
#' @param num_subj number of subjects to sample from MESA data for MESA analysis
#' @param alpha intercept for generating outcome
#' @param theta true spatial scale under which datasets
#' @param delta simulated binary covariate regression effect
#' @param beta SAP effect
#' @param alpha_prior prior to be placed on intercept in model, must be an rstap:: namespace object
#' @param beta_prior prior to be placed on SAP effect
#' @param theta_prior prior to be placed on spatial scale
#' @param delta_prior prior to be placed on simulated binary covariate effect
#' @param iter number of iterations for which to run the stap_glm sampler
#' @param warmup number of iterations to warmup the sampler
#' @param chains number of independent MCMC chains to draw
#' @param cores number of cores with which to run chains in parallel
#' @param file path to file to save tables to in .tex format
#' @return list of two table components table1_top -
#' which includes the coverage and diagnostic statistics broken down by parameter and
#' The remaining "raw" table component contains the pre-aggregation data-frames
#' @export
run_simulation <- function(num_sims = 5,
                           data, 
                           num_subj = 100L,
                           pars = list(alpha = 23,
                                       shape_one = 1,
                                       scale_one = .5,
                                       shape_two = 1,
                                       scale_two = .3,
                                       delta = -2.2,
                                       beta = 1.2,
                                       sigma = 1),
                           modeled_K = c("exp","exp"),
                           priors = list(
                                         alpha = rstap::normal(location = 25,
                                                               scale = 4,
                                                               autoscale = FALSE),
                                          beta = rstap::normal(location = 0, 
                                                               scale = 3,
                                                               autoscale = FALSE),
                                            theta = rstap::log_normal(location = 0,
                                                                      scale = .75),
                                            delta = rstap::normal(location = 0,
                                                                  scale = 3,
                                                                  autoscale = FALSE)
                              ),
                              settings = list(iter = 5E3,
                                              warmup = 3E3,
                                              chains = 5,
                                              cores = 5),
                              seed = NULL)
{
  
    if(is.null(seed)){
        seed <- 34234
        set.seed(seed)
    }
  
  
  mdf <- generate_mesa_dataset(MESA = data,pars = pars, num_subj = num_subj)
  
  simdf <- purrr::map_dfr(1:num_sims,function(x) .fit_and_process_stap(sim_ix = x,
                                                                       data = mdf,
                                                                       modeled_K,
                                                                       pars,
                                                                       priors,
                                                                       settings
                                                                       ))
  
  return(simdf)
  
}

.fit_and_process_stap <- function(sim_ix, data, modeled_K, pars, priors,settings){
  
  
  ## generative function info 
  if(pars$shape_one == 1 && pars$shape_two == 1){
    spatial_exposure_f <- "exp"
    temporal_exposure_f <- "exp"
  }else if(pars$shape_one != 1 && pars$shape_two == 1){
    spatial_exposure_f <- "wei"
    temporal_exposure_f <- "exp"
  }
  else{
    spatial_exposure_f <- "wei"
    temporal_exposure_f <- "wei"
  } 
  
  ## modeled function
  if(modeled_K[1] == "exp" & modeled_K[2] == "exp"){
    f <- outcome ~ sex  + stap(FF,exp,cexp) + (visit_number|id)
    shape_one <- 1
    shape_two <- 1
    so_KLD <- NA
    st_KLD <- NA
    so_cg <- NA
    st_cg <- NA
    so_co <- NA
    st_co <- NA
    so_int <- NA
    st_int <- NA
    so_diff <- NA
    st_diff <- NA
  }else if(modeled_K[1] == "wei" && modeled_K[2] == "exp"){
    f <- outcome ~ sex  + stap(FF,wei,cexp) + (visit_number|id)
    shape_two <- 1
    st_KLD <- NA
    st_co <- NA
    st_int <- NA
    st_cg <- NA
    st_diff <- NA
  }
  else{
    f <- outcome ~ sex  + stap(FF,wei,cwei) + (visit_number|id)
  } 
  
  fit <- rstap::stap_glmer(f,
                             subject_data = data$subject_data %>% 
                               dplyr::mutate(outcome = outcome + 
                                               rnorm(dplyr::n(),
                                                     sd = pars$sigma)),
                             distance_data = data$bef_data %>% 
                               dplyr::select(id,visit_number,BEF,Distance),
                             time_data = data$bef_data %>% 
                               dplyr::select(id,visit_number,BEF,Time),
                             family = gaussian(),
                             group_ID = 'visit_number',
                             subject_ID = 'id',
                             prior = priors$delta,
                             prior_intercept = priors$alpha,
                             prior_stap = priors$beta,
                             prior_theta = priors$theta,
                             iter = settings$iter,
                             max_distance = 10,
                             max_time = max(data$bef_data$Time),
                             warmup = settings$warmup,
                             chains = settings$chains,
                             cores = settings$cores
                            )
  
  print(fit)
  post_pars_all <- as.matrix(fit)
  parnms <- colnames(post_pars_all)
  ix <- stringr::str_which(parnms,"FF")
  parnms <- parnms[ix]
  ps <- c(pars$beta,pars$scale_one,pars$scale_two)
  names(ps) <- parnms
  post_pars <- post_pars_all[,parnms]
  pss <- nrow(post_pars)
  prior_beta <- .get_prior_samples(priors$beta,pss)
  prior_theta <- .get_prior_samples(priors$theta,pss)
  KLD <- function(x,y) LaplacesDemon::KLD(x,y)$mean.sum.KLD
  
  ## handle possible shape parameters
  if(modeled_K[1] == "wei"){
    so_nm <- colnames(post_pars_all)[stringr::str_which(colnames(post_pars_all),
                                                        "spatial_shape")]
    so_truth <- pars$shape_one
    names(so_truth) <- so_nm
    shape_one <- post_pars_all[,so_nm,drop= TRUE]
    so_KLD <- KLD(prior_theta,shape_one)
    so_cg <- calculate_CG_stat(fit, so_truth, mer = TRUE)
    so_co <- check_coverage(fit,so_truth)
    so_int <- interval_length(fit,so_nm)
    so_diff <- abs_diff(fit,so_truth)
  }
  if(modeled_K[2] == "wei"){
    st_nm <- colnames(post_pars_all)[stringr::str_which(colnames(post_pars_all),
                                                        "temporal_shape")]
    st_truth <- pars$shape_two
    names(st_truth) <- st_nm
    shape_two <- post_pars_all[,st_nm,drop = TRUE]
    st_KLD <- KLD(prior_theta,shape_two)
    st_cg <- calculate_CG_stat(fit, st_truth, mer = TRUE)
    st_co <- check_coverage(fit,st_truth)
    st_int <- interval_length(fit,st_nm)
    st_diff <- abs_diff(fit,st_truth)
  }
  
  
  
  out <- dplyr::tibble(sim_ix = sim_ix,
                       Parameter = c("beta","theta_s",
                                     "theta_t","shape_s","shape_t"),
                       Truth = c(ps[1],ps[2],ps[3],
                                 pars$shape_one,pars$shape_two),
                       Est = c(apply(post_pars,2,median),
                               median(shape_one),median(shape_two)),
                       Spatial_f = spatial_exposure_f,
                       Temporal_f = temporal_exposure_f,
                       modeled_spatial_f = modeled_K[1],
                       modeled_temporal_f = modeled_K[2],
                       KLD = c(KLD(prior_beta,post_pars[,1]),
                               KLD(prior_beta,post_pars[,2]),
                               KLD(prior_theta,post_pars[,3]),
                               so_KLD,st_KLD
                               ),
                       CG = c(calculate_CG_stat(fit, ps[1], mer = TRUE),
                              calculate_CG_stat(fit, ps[2], mer = TRUE),
                              calculate_CG_stat(fit, ps[3], mer = TRUE),
                              so_cg,st_cg
                             ),
                       Coverage = c(check_coverage(fit, ps[1]),
                                    check_coverage(fit, ps[2]),
                                    check_coverage(fit, ps[3]),
                                    so_co,st_co
                                     ),
                       Interval = c(interval_length(fit, parnms[1]),
                                    interval_length(fit, parnms[2]),
                                    interval_length(fit, parnms[3]),
                                    so_int,st_int
                                    ),
                       abs_diff = c(abs_diff(fit, ps[1]),
                                    abs_diff(fit, ps[2]),
                                    abs_diff(fit, ps[3]),
                                    so_diff,st_diff
                                    ),
                       WAIC = rstap::waic(fit)
         )
  
  return(out)
}
  
.get_prior_samples <- function(prior,num_samples){
  dst <- prior$dist
  if(is.null(prior$scale))
    scl <- 1
  else
    scl <- prior$scale
  
  if(dst == "normal")
    return(rnorm(num_samples,mean = prior$location,sd = scl))
  else
    return(rlnorm(num_samples,meanlog = prior$location,sdlog = scl))
}
