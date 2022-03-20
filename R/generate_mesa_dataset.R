#' generates stap dataset using MESA HFS-cohort distances
#'
#' @param seed an integer for set.seed()
#' @param MESA dataframe consisting of 4 columns: idno, visit_number, Distance, Time
#' @param num_subj number of subjects to simulate
#' @param prop_dist proportion of MESA distances to use -
#' set to 1 for all MESA distances, or smaller to speed model fitting time
#' @param max_dist upper bound on distances to include
#' @param alpha true intercept value
#' @param delta confounder effect
#' @param beta true bef effect
#' @param theta true spatial scale
#' @param K exposure function
#'
#' @export
generate_mesa_dataset <- function(seed = NULL,
                                  MESA = NULL,
                                  num_subj = 100,
                                  pars,
                                  K = function(d,t,
                                               shape_one = 1,
                                               scale_one = .5,
                                               shape_two = 1,
                                               scale_two = 1) {
                                      pweibull(d,
                                               shape = shape_one, 
                                               scale = scale_one,
                                               lower.tail = FALSE) *
                                      pweibull(t,
                                               shape = shape_two,
                                               scale = scale_two,
                                               lower.tail = TRUE)
                                      }
                                  )
{
                                  
    if(!is.null(seed))
        set.seed(seed)
    else
        set.seed(3224314)
    if(is.null(MESA) || !is.data.frame(MESA))
        stop("MESA data.frame required")
        

    idno <- MESA %>% 
        dplyr::distinct(id) %>% 
        dplyr::pull() %>%
        sample(., size = num_subj, replace = F)
    
    MESA_df <- MESA %>% 
        dplyr::filter(id %in% idno) %>%
        dplyr::group_by(id,visit_number) %>%
        dplyr::summarise(Exposure = sum(K(Distance,Time,
                                          shape_one = pars$shape_one,
                                          shape_two = pars$shape_two,
                                          scale_one = pars$scale_one,
                                          scale_two = pars$scale_two))) %>% 
        dplyr::left_join(MESA %>% 
                             dplyr::distinct(id) %>% 
                             dplyr::mutate(Sex = rbinom(n = dplyr::n(), size = 1, prob = 0.5),
                                           intercept = rnorm(n = dplyr::n(), sd = .25),
                                           slope = rnorm(n = dplyr::n(), sd = .35))
                         ) %>% 
        dplyr::ungroup()

    eta <- pars$alpha + 
       MESA_df$Sex * pars$delta + 
        MESA_df$Exposure * pars$beta + 
        MESA_df$intercept + 
        MESA_df$visit_number * MESA_df$slope 

    subj_data <- tibble::tibble(outcome = eta,
                                sex = MESA_df$Sex,
                                id = MESA_df$id,
                                intercept = MESA_df$intercept,
                                slope = MESA_df$slope,
                                visit_number = MESA_df$visit_number
    )

    bef_data <- MESA %>%
        dplyr::mutate(visit_number = as.integer(visit_number),
                      BEF = "FF")

    return(list(subject_data = subj_data,
                bef_data =  bef_data,
                call = match.call(expand.dots = T)))

}
