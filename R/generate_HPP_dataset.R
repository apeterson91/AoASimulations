#' generates homogenous poisson process data
#'
#' @export
#' @param seed an integer for set.seed()
#' @param num_subj number of subjects to simulate
#' @param num_dists number of distances/subject to simulate
#' @param max_dist upper bound on distances to include
#' @param bef_xlim x coordinate bounds from which to simulate befs
#' @param bef_ylim y coordinate bounds from which to simulate bef positions
#' @param subj_xlim x coordinate bounds from which to simulate subject positions
#' @param subj_ylim y coordinate bounds from which to simulate subject positions
#' @param alpha true intercept value
#' @param delta confounder effect
#' @param beta true bef effect
#' @param theta true spatial scale
#' @param W weight function
generate_hpp_dataset <- function(seed,
                                 num_subj = 100L,
                                 num_dists = 30L,
                                 max_dist = 5L,
                                 bef_xlim = c(-1,1),
                                 bef_ylim= c(-1,1),
                                 subj_xlim = c(-1,1),
                                 subj_ylim = c(-1,1),
                                 alpha = 25,
                                 delta = -2.3,
                                 beta = .1,
                                 theta = 0.5,
                                 shape = NULL,
                                 sigma = 1,
                                 W = function(x,y) exp(-x)){

    if(!is.null(seed))
        set.seed(seed)
    if(length(num_dists)!=1 )
        stop("num_befs and num_dists should be vectors of the same length")
    if(!(length(subj_xlim)==2)&&length(subj_ylim)==2)
        stop("xlim and ylim should both be vectors of length 2")
    if(!(length(bef_xlim)==2)&&length(bef_ylim)==2)
        stop("xlim and ylim should both be vectors of length 2")

    bef_posdata <-  tibble::tibble(x = runif(n = num_dists,min = bef_xlim[1],
                                    max= bef_xlim[2]),
                        y = runif(n = num_dists, min = bef_ylim[1],
                                    max = bef_ylim[2]))

    subj_data <- tibble::tibble(x_coord = runif(n = num_subj,
                                      min = subj_xlim[1],
                                      max = subj_xlim[2]),
                        y_coord = runif(n = num_subj,
                                  min = subj_ylim[1],
                                  max = subj_ylim[2]),
                        subj_id = 1:num_subj)

    distances <- fields::rdist(as.matrix(subj_data)[,c("x_coord","y_coord")],
                               as.matrix(bef_posdata[,c("x","y")]))
    colnames(distances) <- paste0("BEF_",1:ncol(distances))
    distances <- dplyr::mutate(dplyr::as_tibble(distances),subj_id = 1:num_subj)
    distances <- tidyr::gather(distances,dplyr::contains("BEF"),key = "BEF",value="Distance")
    if(is.null(shape))
        X <- distances %>% dplyr::group_by(subj_id) %>%
        dplyr::summarise(Exposure = sum(W(Distance/theta))) %>%
        dplyr::ungroup() %>%
        dplyr::pull(Exposure)
    else
        X <- distances %>% dplyr::group_by(subj_id) %>%
        dplyr::summarise(Exposure = sum(exp(-(Distance/theta)^shape))) %>%
        dplyr::ungroup() %>%
        dplyr::pull(Exposure)

    sex <- rbinom(n = num_subj,size = 1,prob = .5)
    outcome <- alpha + sex*delta + X* beta + rnorm(n = num_subj,
                                                   mean = 0,
                                                   sd = sigma)
    subj_data <- dplyr::mutate(subj_data,outcome = outcome, sex = sex)

    return(list(subject_data = subj_data,
                bef_data = dplyr::mutate(distances,BEF="FF"),
                call = match.call(expand.dots = T)))
}



