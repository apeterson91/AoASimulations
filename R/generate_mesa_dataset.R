#' generates stap dataset using MESA HFS-cohort distances
#'
#' @param seed an integer for set.seed()
#' @param num_subj number of subjects to simulate
#' @param max_dist upper bound on distances to include
#' @param alpha true intercept value
#' @param delta confounder effect
#' @param beta true bef effect
#' @param theta true spatial scale
#' @param W weight function
generate_mesa_dataset <- function(seed = NULL,
                                  num_subj = NULL,
                                  max_dist = 10,
                                  alpha = 23,
                                  delta = -2.2,
                                  beta = .65,
                                  beta_bar = .85,
                                  theta = .7,
                                  sigma = 1,
                                  W = function(x) exp(-x)){
    ### fake gender covariate

    if(!is.null(seed))
        set.seed(seed)
    else
        set.seed(3224314)

    if(is.null(num_subj)){
        num_subj = length(unique(MESA$id))
        MESA_df <- MESA
    }
    else{
        idno <- sample(MESA$id,size=num_subj,replace=F)
        MESA_df <- MESA %>% dplyr::filter(id %in% idno)
    }

    X <- MESA_df %>% dplyr::group_by(id,visit_number) %>%
        dplyr::summarise(Exposure = sum(W((Total_Kilometers / theta))) ) %>%
        dplyr::ungroup()

    Xbar <- X %>% dplyr::group_by(id) %>%
        dplyr::summarise(Mean_Exposure = mean(Exposure)) %>%
        dplyr::mutate(Sex = rbinom(n=dplyr::n(),size=1,prob=.5),
                         intercept = rnorm(n = dplyr::n(),sd = 1),
                         slope = rnorm(n=dplyr::n(),sd=1.2))

    X_diff <- X %>% dplyr::left_join(Xbar,by='id') %>%
        dplyr::mutate(X_diff = Exposure - Mean_Exposure,)

    eta <- alpha + delta * X_diff$Sex + X_diff$X_diff*beta + X_diff$Mean_Exposure*beta_bar +
        X_diff$intercept + X_diff$slope * X_diff$visit_number

    subj_data <- tibble::tibble(outcome = rnorm(n = nrow(X_diff),mean = eta,sd = sigma),
                              sex = X_diff$Sex,
                              id = X_diff$id,
                              visit_number = X_diff$visit_number
    )

    bef_data <- MESA_df %>%
        dplyr::mutate(visit_number = as.integer(visit_number),
               BEF = "FF")

    return(list(subject_data = subj_data,
                bef_data =  bef_data,
                call = match.call(expand.dots = T)))

}
