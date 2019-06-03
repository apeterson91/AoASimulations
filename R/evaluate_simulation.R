
library(rstap)

check_coverage <- function(model,true_par,prob = .90,nonzero = T){
    if(!all(names(true_par)%in%rownames(posterior_interval(model))))
        stop("Parameter"," ", names(true_par)," ","Not included in model, please try again")
    CI <- posterior_interval(model,par = names(true_par),prob = prob)
    if(nonzero)
        return( (CI[1]<=true_par && CI[2]>=true_par) && ( CI[2]<=0 || CI[1]>=0  ))
    else
        return((CI[1]<=true_par && CI[2]>=true_par))
}

interval_length <- function(model,par,prob = .95){
    ci <- posterior_interval(model,pars=par)
    return(median(ci[,2]-ci[,1]))
}


calculate_RMSE_median <- function(model){
    mean(sqrt((residuals(model))^2))
}





