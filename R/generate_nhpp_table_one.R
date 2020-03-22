#' generates nonhomogenous poisson process data without residual error for table one
#'
#' @param seed an integer for set.seed()
#' @param num_subj_realizations number of subjects position processes to simulate
#' @param num_dist_realizations number of bef processes to simulate
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
generate_nhpp_table_one <- function(
                                  seed = NULL,
                                  num_subj_realizations = 10L,
                                  num_dist_realizations = 4L,
                                  Lambda_xy = function(x,y){ (3 - x^2 -y^2) },
                                  bef_xlim = c(-1,1),
                                  bef_ylim=c(-1,1),
                                  subj_xlim = c(-1,1),
                                  subj_ylim = c(-1,1),
                                  alpha = 23,
                                  delta = -2.2,
                                  beta = .75,
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

    bef_points <- spatstat::rpoispp(lambda = Lambda_xy,
                                    win = bef_window,
                                    nsim = num_dist_realizations) ## Need to make sure enough data points are generated

    subj_points <- spatstat::rpoispp(lambda = Lambda_xy,
                                     win = subj_window,
                                     nsim = num_subj_realizations)


    bef_posdata <- purrr::map2_dfr(1:length(bef_points),bef_points,function(a,b){
        tibble::tibble(sim_id = a,
                       x = b$x,
                       y = b$y)})



    subj_data  <- purrr::map2_dfr(1:length(subj_points),subj_points,function(a,b){
        tibble::tibble(sim_id = a,
                       x_coord = b$x,
                       y_coord = b$y)})

    subj_data <- subj_data %>%
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
    sex <- rbinom(n = nrow(subj_data), size = 1, prob = .5)
    eta <- alpha + sex*delta + X* beta
    subj_data <- dplyr::mutate(subj_data, outcome = eta, sex = sex, X = X)

    return(list(subject_data = subj_data,
                bef_data = dplyr::mutate(distances,BEF="FF"),
                call = match.call(expand.dots = T)))
}
