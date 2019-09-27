#' check coverage
#'
#' @param model stapreg object from rstap package
#' @param true_par a named vector with the true parameter values as elements
#' @param prob probability value in (0,1)
#' @param nonzero logical indicating zero should not be in the posterior interval
check_coverage <- function(model,true_par,prob = .90,nonzero = T){
    if(!all(names(true_par)%in%rownames(rstap::posterior_interval(model))))
        stop("Parameter"," ", names(true_par)," ","Not included in model, please try again")
    CI <- rstap::posterior_interval(model,par = names(true_par),prob = prob)
    if(nonzero)
        return( (CI[1]<=true_par && CI[2]>=true_par) && ( CI[2]<=0 || CI[1]>=0  ))
    else
        return((CI[1]<=true_par && CI[2]>=true_par))
}

#' calculate interval length
#'
#' @param model stapreg object from rstap package
#' @param par parameter names as vector of strings
#' @param prob scalar numeric in (0,1)
#'
interval_length <- function(model,par,prob = .95){
    ci <- rstap::posterior_interval(model,pars=par)
    return(median(ci[,2]-ci[,1]))
}

#' Calculate Absolute difference in true-median beta
#'
#' @param model stapreg object from rstap package
#' @param par parameter named vector parameter
#' with true parameter value as elements
abs_diff  <- function(model,par){
	nm <- names(par)
	if(!(nm %in% colnames(as.matrix(model))))
		stop("par must be a named vector estimated in the model")

	estimates <- apply(as.matrix(model)[,nm,drop = FALSE],2,median)

	return(abs(estimates - par))
}


#' calculate Median Root Mean Square Error
#'
#' @param model stapreg object from rstap package
#'
calculate_RMSE_median <- function(model){
    mean(sqrt((residuals(model))^2))
}

#' calculate C&G validation statistic
#'
#' @param model stapreg object from rstap package
#' @param par named vector
#'
calculate_CG_stat <- function(model,par,mer=FALSE){
    nm <- names(par)
    if(length(nm)>1)
        stop("par must be only 1 parameter")
    if(!mer){
        if(!(nm %in% names(coef(model))))
            stop("par must be a named vector estimated in the model")
        return(qnorm(mean(as.matrix(model)[,nm,drop=TRUE]>par), mean = 0, sd = 1)^2)
    }
    else{
        if(!(nm %in% names(coef(model)[[1]])))
            stop("par must be a coefficient name estimated in first group of model")
        return(qnorm(mean(as.matrix(model)[,nm,drop=TRUE]>par), mean=0, sd = 1)^2)
    }
}




