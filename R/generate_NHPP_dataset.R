#' generates nonhomogenous poisson process data
#'
#' @param seed an integer for set.seed()
#' @param num_subj number of subjects to simulate
#' @param num_dists number of distances/subject to simulate
#' @param max_dist upper bound on distances to include
#' @param bef_xlim x coordinate bounds from which to simulate befs
#' @param bef_ylim y coordinate bounds from which to simulate bef positions
#' @param subj_xlim x coordinate bounds from which to simulate subject positions
#' @param subj_ylim y coordinate bounds from which to simulate subject positions
#' @param inten intensity function for nhpp as a function of spatial coords (x,y)
#' @param alpha true intercept value
#' @param delta confounder effect
#' @param beta true bef effect
#' @param theta true spatial scale
#' @param K exposure function
#'
#' @export
generate_nhpp_dataset <- function(seed = NULL,
                                 num_subj = 100L,
                                 num_dists = 30L,
                                 Lambda_xy = function(x,y){ 20 - x^2 -y^2 },
                                 max_dist = 5L,
                                 bef_xlim = c(-1,1),
                                 bef_ylim=c(-1,1),
                                 subj_xlim = c(-1,1),
                                 subj_ylim = c(-1,1),
                                 alpha = 25,
                                 delta = -2.3,
                                 beta = .1,
                                 theta = 0.5,
                                 sigma = 1,
                                 K = function(x) exp(-x)){

    if(!is.null(seed))
        set.seed(seed)
    if(!(length(subj_xlim)==2)&&length(subj_ylim)==2)
        stop("xlim and ylim should both be vectors of length 2")
    if(!(length(bef_xlim)==2)&&length(bef_ylim)==2)
        stop("xlim and ylim should both be vectors of length 2")

    bef_window <- spatstat::as.owin(list(xmin=bef_xlim[1],xmax=bef_xlim[2],
                           ymin=bef_ylim[1],ymax=bef_xlim[2]))

    subj_window <- spatstat::as.owin(list(xmin=subj_xlim[1],xmax=subj_xlim[2],
                                          ymin=subj_ylim[1],ymax=subj_xlim[2]))

    bef_points <- spatstat::rpoispp(lambda =Lambda_xy,
                                    win = bef_window,
                                    nsim = 100)

    subj_points <- spatstat::rpoispp(lambda = Lambda_xy,
                                     win = subj_window,
                                     nsim = 100)


    bef_posdata <-  tibble::tibble(x= unlist(purrr::map(bef_points,function(z) z[["x"]])),
                            y = unlist(purrr::map(bef_points,function(z) z[["y"]]))) %>%
        dplyr::sample_n(num_dists)


    subj_data <- tibble::tibble(x_coord = unlist(purrr::map(subj_points,function(z) z[["x"]])),
                                    y_coord = unlist(purrr::map(subj_points,function(z) z[["y"]]))) %>%
        dplyr::sample_n(num_subj) %>%
        dplyr::mutate(subj_id = 1:dplyr::n())


    distances <- fields::rdist(as.matrix(subj_data)[,c("x_coord","y_coord")],
                               as.matrix(bef_posdata[,c("x","y")]))
    colnames(distances) <- paste0("BEF_",1:ncol(distances))
    distances <- dplyr::mutate(dplyr::as_tibble(distances),subj_id = 1:dplyr::n())
    distances <- tidyr::gather(distances,dplyr::contains("BEF"),key = "BEF",value="Distance")
    X <- distances %>% dplyr::group_by(subj_id) %>%
        dplyr::summarise(Exposure = sum(K(Distance/theta))) %>%
        dplyr::ungroup() %>%
        dplyr::pull(Exposure)
    sex <- rbinom(n = nrow(subj_data),size = 1,prob = .5)
    eta <- alpha + sex*delta + X* beta
    outcome <- rnorm(n = nrow(subj_data),mean = eta, sd = sigma)
    subj_data <- dplyr::mutate(subj_data, outcome = outcome, sex = sex)

    return(list(subject_data = subj_data,
                bef_data = dplyr::mutate(distances,BEF="FF"),
                call = match.call(expand.dots = T)))
}
