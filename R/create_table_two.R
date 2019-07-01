#' create_table_two
#'
#' @param num_sims
#'
#' @export
create_table_two <- function(num_sims = 5,
                             num_subj = 100,
                             num_dists = 30,
                             alpha = 23,
                             theta =  .5,
                             shape = 5,
                             delta = -2.2,
                             beta = .75,
                             iter = 2E3,
                             warmup = 1E3,
                             chains = 1,
                             cores = 1,
                             file = NULL){

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

    ddatasets <- purrr::map(1:num_sims, function(x) generate_dlm_dataset(seed = x,
                                                                         alpha = alpha,
                                                                         beta = beta,
                                                                         delta = delta,
                                                                         W = function(x) (x<=theta)*1))

    # DLM - 2

    d2datasets <- purrr::map(1:num_sims, function(x) generate_dlm_dataset(seed = x,
                                                                         alpha = alpha,
                                                                         beta = beta,
                                                                         delta = delta,
                                                                         W = function(x) (x<=theta)*(1-x^2) ))

    # Fit DLM under DLM

    dlm_lists <- purrr::map(ddatasets,function(x){
        lag <- seq(from = .1,to = floor(max(x$bef_data$Distance)), by = .1)
        Conc <- suppressWarnings(x$bef_data %>%
            dplyr::mutate(bins = cut(Distance,breaks = c(0,lag),
                                     include.lowest = TRUE )) %>%
            dplyr::group_by(subj_id,bins) %>% dplyr::count() %>%
            dplyr::rename(count = n) %>%
            tidyr::spread(bins,count) %>%
            dplyr::ungroup() %>%
            dplyr::select(-subj_id))
        labs <- rep(0,ncol(Conc))
        names(labs) <- colnames(Conc)
        Conc <- as.matrix(tidyr::replace_na(Conc, replace = lapply(labs,identity)))
        out <- list(Conc = Conc[,1:(NCOL(Conc)-1)],
                    outcome = x$subject_data$outcome,
                    sex = x$subject_data$sex,
                    lag = lag)
    })

    cps <- numeric(num_sims) # changepoints
    for(i in 1:num_sims){
        assign("sex", dlm_lists[[i]]$sex, envir = globalenv())
        assign("outcome",dlm_lists[[i]]$outcome,envir = globalenv())
        assign("Conc", dlm_lists[[i]]$Conc, envir = globalenv())
        assign("lag",dlm_lists[[i]]$lag, envir = globalenv())
        fit <- dlmBE::dlm(outcome ~ sex + dlmBE::cr(lag,Conc))
        cps[i] <- lag[unlist(dlmBE::changePoint(fit))[1]+1]
    }


    dlm_lists <- purrr::map(d2datasets,function(x){
        lag <- seq(from = .1,to = floor(max(x$bef_data$Distance)), by = .1)
        Conc <- suppressWarnings(x$bef_data %>%
            dplyr::mutate(bins = cut(Distance,breaks = c(0,lag),
                                     include.lowest = TRUE )) %>%
            dplyr::group_by(subj_id,bins) %>% dplyr::count() %>%
            dplyr::rename(count = n) %>%
            tidyr::spread(bins,count) %>%
            dplyr::ungroup() %>%
            dplyr::select(-subj_id))
        labs <- rep(0,ncol(Conc))
        names(labs) <- colnames(Conc)
        Conc <- as.matrix(tidyr::replace_na(Conc, replace = lapply(labs,identity)))
        out <- list(Conc = Conc[,1:(NCOL(Conc)-1)],
                    outcome = x$subject_data$outcome,
                    sex = x$subject_data$sex,
                    lag = lag)
    })


    cp2s <- numeric(num_sims) # changepoints
    for(i in 1:num_sims){
        assign("sex", dlm_lists[[i]]$sex, envir = globalenv())
        assign("outcome",dlm_lists[[i]]$outcome,envir = globalenv())
        assign("Conc", dlm_lists[[i]]$Conc, envir = globalenv())
        assign("lag",dlm_lists[[i]]$lag, envir = globalenv())
        fit <- dlmBE::dlm(outcome ~ sex + dlmBE::cr(lag,Conc))
        cp2s[i] <- lag[unlist(dlmBE::changePoint(fit))[1]+1]
    }

    dlm_lists <- purrr::map(edatasets,function(x){
        lag <- seq(from = .1,to = floor(max(x$bef_data$Distance)), by = .1)
        Conc <- suppressWarnings(x$bef_data %>%
            dplyr::mutate(bins = cut(Distance,breaks = c(0,lag),
                                     include.lowest = TRUE )) %>%
            dplyr::group_by(subj_id,bins) %>% dplyr::count() %>%
            dplyr::rename(count = n) %>%
            tidyr::spread(bins,count) %>%
            dplyr::ungroup() %>%
            dplyr::select(-subj_id))
        labs <- rep(0,ncol(Conc))
        names(labs) <- colnames(Conc)
        Conc <- as.matrix(tidyr::replace_na(Conc, replace = lapply(labs,identity)))
        out <- list(Conc = Conc[,1:(NCOL(Conc)-1)],
                    outcome = x$subject_data$outcome,
                    sex = x$subject_data$sex,
                    lag = lag)
    })

    cpes <- numeric(num_sims) # changepoints
    for(i in 1:num_sims){
        assign("sex", dlm_lists[[i]]$sex, envir = globalenv())
        assign("outcome",dlm_lists[[i]]$outcome,envir = globalenv())
        assign("Conc", dlm_lists[[i]]$Conc, envir = globalenv())
        assign("lag",dlm_lists[[i]]$lag, envir = globalenv())
        fit <- dlmBE::dlm(outcome ~ sex + dlmBE::cr(lag,Conc))
        ci <- dlmBE::confint.dlMod(fit, coef=FALSE)
        non0 <- !(ci[, 1] <= 0.05 & ci[, ncol(ci)] >= 0.05)
        rslt <- lapply(dlmBE::lagIndex(fit),
               function(i) {
                   x <- non0[i]
                   which(x & !c(tail(x, -1), TRUE) & c(FALSE, head(x, -1)))
               })
        cpes[i] <- lag[unlist(rslt)[1]+1]
    }

    dlm_lists <- purrr::map(wdatasets,function(x){
        lag <- seq(from = .1,to = floor(max(x$bef_data$Distance)), by = .1)
        Conc <- suppressWarnings(x$bef_data %>%
            dplyr::mutate(bins = cut(Distance,breaks = c(0,lag),
                                     include.lowest = TRUE )) %>%
            dplyr::group_by(subj_id,bins) %>% dplyr::count() %>%
            dplyr::rename(count = n) %>%
            tidyr::spread(bins,count) %>%
            dplyr::ungroup() %>%
            dplyr::select(-subj_id))
        labs <- rep(0,ncol(Conc))
        names(labs) <- colnames(Conc)
        Conc <- as.matrix(tidyr::replace_na(Conc, replace = lapply(labs,identity)))
        out <- list(Conc = Conc[,1:(NCOL(Conc)-1)],
                    outcome = x$subject_data$outcome,
                    sex = x$subject_data$sex,
                    lag = lag)
    })

    cpws <- numeric(num_sims) # changepoints
    for(i in 1:num_sims){
        assign("sex", dlm_lists[[i]]$sex, envir = globalenv())
        assign("outcome",dlm_lists[[i]]$outcome,envir = globalenv())
        assign("Conc", dlm_lists[[i]]$Conc, envir = globalenv())
        assign("lag",dlm_lists[[i]]$lag, envir = globalenv())
        fit <- dlmBE::dlm(outcome ~ sex + dlmBE::cr(lag,Conc))
        ci <- dlmBE::confint.dlMod(fit, coef=FALSE)
        non0 <- !(ci[, 1] <= 0.05 & ci[, ncol(ci)] >= 0.05)
        rslt <- lapply(dlmBE::lagIndex(fit),
               function(i) {
                   x <- non0[i]
                   which(x & !c(tail(x, -1), TRUE) & c(FALSE, head(x, -1)))
               })
        cpws[i] <- lag[unlist(rslt)[1]+1]
    }

    # STAP Model fitting

    # Fit DLM under Exponential

    de <- purrr::map(ddatasets,function(x){
        rstap::stap_glm(outcome ~ sex + sap(FF,exp),
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

    # Fit DLM under Weibull

    dw <- purrr::map(ddatasets,function(x){
        rstap::stap_glm(outcome~sex + sap(FF,wei),
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

    # Fit DLM2 under Exponential

    d2e <- purrr::map(d2datasets,function(x){
        rstap::stap_glm(outcome~sex + sap(FF,exp),
                        subject_data = x$subject_data,
                        distance_data = x$bef_data,
                        max_distance = 10,
                        subject_ID = "subj_id",
                        prior = rstap::normal(),
                        prior_stap = rstap::normal(),
                        prior_intercept = rstap::normal(location = 25),
                        prior_theta = rstap::log_normal(1,1),
                        chains = chains,
                        cores = cores,
                        iter = iter,
                        warmup = warmup)})

    # Fit DLM2 under Weibull

    d2w <- purrr::map(d2datasets,function(x){
        rstap::stap_glm(outcome~sex + sap(FF,wei),
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
                prior_theta = rstap::log_normal(1,1),
                chains = chains,
                cores = cores,
                iter = iter,
                warmup = warmup)})


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
                prior_theta = rstap::log_normal(1,1),
                chains = chains,
                cores = cores,
                iter = iter,
                warmup = warmup)})

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
                cores = cores,
                chains = chains,
                cores = cores,
                iter = iter,
                warmup = warmup)})

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
                prior_theta = rstap::log_normal(1,1),
                chains = chains,
                cores = cores,
                iter = iter,
                warmup = warmup)})

    term_step_d <- tibble::tibble(sim_id = 1:num_sims,
                              Simulated_Function = "Step Function",
                              Modeled_Function = "DLM",
                              True_termination = theta,
                              Estimate_termination = cps,
                              )


    term_step_e <- tibble::tibble(sim_id = 1:num_sims,
                                Simulated_Function = rep("Step Function",num_sims),
                                Modeled_Function = rep("Exponential", num_sims),
                                True_termination = theta,
                                Estimate_termination = purrr::map_dbl(de,function(x) rstap::stap_termination(x,max_value=1E8)[2]))


    term_step_w <- tibble::tibble(sim_id = 1:num_sims,
                                Simulated_Function = rep("Step Function",num_sims),
                                Modeled_Function = rep("Weibull", num_sims),
                                True_termination = theta,
                                Estimate_termination = purrr::map_dbl(dw,function(x) rstap::stap_termination(x,max_value=1E8)[2]))


    term_q_d <- tibble::tibble(sim_id = 1:num_sims,
                              Simulated_Function = "Quadratic Step",
                              Modeled_Function = "DLM",
                              True_termination = theta,
                              Estimate_termination = cp2s)


    term_q_e <- tibble::tibble(sim_id = 1:num_sims,
                              Simulated_Function = rep("Quadratic Step",num_sims),
                              Modeled_Function = rep("Exponential", num_sims),
                              True_termination = theta,
                              Estimate_termination = purrr::map_dbl(d2e,function(x) rstap::stap_termination(x,max_value=1E8)[2]))


    term_q_w <- tibble::tibble(sim_id = 1:num_sims,
                                 Simulated_Function = rep("Quadratic Step",num_sims),
                                 Modeled_Function = rep("Weibull", num_sims),
                                 True_termination = theta,
                                 Estimate_termination = purrr::map_dbl(d2w,function(x) rstap::stap_termination(x,max_value=1E8)[2]))

    term_e_d <- tibble::tibble(sim_id = 1:num_sims,
                                 Simulated_Function = "Exponential",
                                 Modeled_Function = "DLM",
                                 True_termination = uniroot(function(x) exp(-(x/theta)) - 0.05,interval = c(0,10))$root,
                                 Estimate_termination = cpes)

    term_exp  <- tibble::tibble(sim_id = 1:num_sims,
                                Simulated_Function = rep("Exponential",num_sims),
                                Modeled_Function = rep("Exponential",num_sims),
                                True_termination = uniroot(function(x) exp(-(x/theta)) - 0.05,interval = c(0,10))$root,
                                Estimate_termination = purrr::map_dbl(ee,function(a) rstap::stap_termination(a,max_value=1E8)[2])
    )

    term_e_w <- tibble::tibble(sim_id = 1:num_sims,
                                  Simulated_Function = rep("Exponential",num_sims),
                                  Modeled_Function = rep("Weibull",num_sims),
                                  True_termination = uniroot(function(x) exp(-(x/theta)) - 0.05,interval = c(0,10))$root,
                                  Estimate_termination = purrr::map_dbl(ew,function(a) rstap::stap_termination(a,max_value=1E8)[2])
    )

    term_w_d <- tibble::tibble(sim_id = 1:num_sims,
                                  Simulated_Function = "Weibull",
                                  Modeled_Function = "DLM",
                                  True_termination = uniroot(function(x){ exp(-(x/theta)^shape) - 0.05},interval = c(0,10))$root,
                                  Estimate_termination = cpws)

    term_w_e <- tibble::tibble(sim_id = 1:num_sims,
                                  Simulated_Function = rep("Weibull",num_sims),
                                  Modeled_Function = rep("Exponential",num_sims),
                                  True_termination = uniroot(function(x){ exp(-(x/theta)^shape) - 0.05},interval = c(0,10))$root,
                                  Estimate_termination = purrr::map_dbl(we,function(a) rstap::stap_termination(a,max_value=1E8)[2])
    )

    term_wei <- tibble::tibble(sim_id = 1:num_sims,
                                        Simulated_Function = rep("Weibull",num_sims),
                                        Modeled_Function = rep("Weibull",num_sims),
                                        True_termination = uniroot(function(x){ exp(-(x/theta)^shape) - 0.05},interval = c(0,10))$root,
                                        Estimate_termination = purrr::map_dbl(ww,function(a) rstap::stap_termination(a,max_value=1E8)[2])
    )

    out <- dplyr::bind_rows(term_step_d,term_step_e,term_step_w,
                            term_q_d, term_q_e,term_q_w,
                            term_e_d,term_exp,term_e_w,
                            term_w_d,term_w_e,term_wei) %>%
        dplyr::mutate(Termination_Difference=abs(True_termination - Estimate_termination)) %>%
        dplyr::group_by(Simulated_Function,Modeled_Function) %>%
        dplyr::summarise(mean_difference = mean(Termination_Difference)) %>%
        dplyr::ungroup() %>%
        tidyr::spread(Modeled_Function,mean_difference)
}

