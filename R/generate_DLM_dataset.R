#' generates Distributed Lag Model or polynomial function effect
#'
#' @param seed an integer for set.seed()
#' @param num_subj number of subjects to simulate
#' @param num_dists number of distances/subject to simulate
#' @param alpha true intercept value
#' @param delta confounder effect
#' @param beta true bef effect
#' @param theta true spatial scale
#' @param W weight function
#'
#' @export
generate_dlm_dataset <- function(seed = NULL,
                                 num_subj = 100,
                                 num_dists = 30,
                                 bef_xlim = c(-1,1),
                                 bef_ylim=c(-1,1),
                                 subj_xlim = c(-1,1),
                                 subj_ylim = c(-1,1),
                                 alpha = 25,
                                 delta = -2.3,
                                 beta = .1,
                                 sigma = 1,
                                 K = function(x) {(x<=1.2)*1 }){
    if(!is.null(seed))
        set.seed(seed)
    else
        set.seed(3141)

    bef_locations <- tibble::tibble(x = runif(num_dists,
                                              min = bef_xlim[1],
                                              max = bef_xlim[2]),
                                    y = runif(num_dists,
                                              min = bef_ylim[1],
                                              max = bef_ylim[2]))
    subj_locations <- tibble::tibble(x = runif(num_subj,
                                               min = subj_xlim[1],
                                               max = subj_xlim[2]),
                                     y = runif(num_subj,
                                               min = subj_ylim[1],
                                               max = subj_ylim[2]))

    subj.xy <- as.matrix(subj_locations)
    feat.xy <- as.matrix(bef_locations)

    dists <- tibble::as_tibble(fields::rdist(subj.xy,feat.xy)) %>%
        dplyr::mutate(subj_id = 1:num_subj) %>%
        tidyr::gather(dplyr::contains("V"),key = "BEF",value="Distance") %>%
        dplyr::mutate(BEF = "FF")

    X <- dists %>% dplyr::group_by(subj_id) %>%
        dplyr::summarise(Exposure = sum(K(Distance))) %>%
        dplyr::ungroup() %>%
        dplyr::select(Exposure) %>%
        dplyr::pull()



    ## Simulate BMI
    ## -------------------------------------------------------------------


    female <- rbinom(num_subj, 1, 0.5)


    Z <- cbind(1, female)
    b <- c(alpha, delta)

    mu <- Z %*% b + X*beta
    y <- rnorm(n = num_subj, mean = mu, sd = sigma)

    subj_data <- tibble::tibble(subj_id = 1:nrow(subj.xy),
                                outcome = y,
                                sex = female)




    return(list(subject_data = subj_data,
                bef_data = dists,
                call = match.call(expand.dots = T)))
}

