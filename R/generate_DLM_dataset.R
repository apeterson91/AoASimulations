#' generates Distributed Lag Model or polynomial function effect
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
#' @param W weight function
generate_dlm_dataset <- function(seed,
                                 num_subj = 100L,
                                 num_dists = 50L,
                                 max_dist = 10L,
                                 bef_xlim = c(0L,1L),
                                 bef_ylim = c(0L,1L),
                                 subj_xlim = c(0L,1L),
                                 subj_ylim = c(0L,1L),
                                 alpha = 25L,
                                 delta = -2.3,
                                 beta = .1,
                                 kern = 20,
                                 theta = 40,
                                 residual_error = 1){

    n.loc <- num_dists
    loc <- 1:n.loc
    sigma.x <- residual_error
    mu.x <- log(0.2) - 0.5 * sigma.x^2
    phi <- 8/3


    xy <- expand.grid(x = loc, y = loc)  # grid of locations

    Sigma.L <- exp(-(0.99 / phi) * as.matrix(dist(xy))) + 0.1 * Matrix::Diagonal(nrow(xy))
    Sigma.L <- Matrix::t(as(chol(Sigma.L), "sparseMatrix"))

    Mu <- exp(c(mu.x + as.matrix(Sigma.L %*% rnorm(nrow(xy)))))
    ct <- rpois(nrow(xy), Mu)

    ## Simulate subject locations
    ## -------------------------------------------------------------------
    ## Subject locations sample uniformly over the space
    ##

    subj.xy <- xy[sample.int(nrow(xy), num_subj, replace = TRUE), ]

    feat.xy <- as.matrix(xy[rep(1:length(ct), ct), ])
    subj.xy <- as.matrix(subj.xy)
    rownames (feat.xy) <- rownames (subj.xy) <- NULL


    count.features <- function(xy, feature.xy, radii) {
        .dist <- function(x) sqrt(sum(x^2))
        dxy <- apply(sweep(feature.xy, 2, xy), 1, .dist)
        table(cut(dxy, radii, include.lowest = TRUE))
    }

    lag <- 1:num_dists
    Conc <- t(apply(subj.xy, 1, count.features,
                    feature.xy = feat.xy, radii = c(0, lag))
    )

    rescl <- function(x) x / max(x, na.rm = TRUE)

    age <- sample(18:80, num_subj, replace = TRUE,
           prob = (18:80 <= 60) + (18:80 > 60) * rescl(dnorm(18:80, 61, 9))
    )

    female <- rbinom(num_subj, 1, 0.5)
    X <- cbind(1,female)
    b <- c(alpha, delta)

    kern1 <- as.numeric(lag <= kern)
    theta1 <- kern1 / theta

    mu <- X %*% b + Conc %*% theta1
    y <- rnorm(num_subj, mu)

    subj_data <- tibble::tibble(id = 1:num_subj,
                                outcome = y,
                                sex = female)

    bef_data <- tibble::as_tibble(fields::rdist(subj.xy,feat.xy)) %>%
        dplyr::mutate(id = 1:num_subj) %>%
        tidyr::gather(dplyr::contains("V"),key="BEF",value="Distance") %>%
        dplyr::mutate(BEF = "FF")

    dlm_fit <- dlmBE::dlm(y ~ female + dlmBE::cr(lag,Conc))

    return(list(subject_data = subj_data,
                bef_data = bef_data,
                dlm_model = dlm_fit,
                call = match.call(expand.dots = T)))
}

