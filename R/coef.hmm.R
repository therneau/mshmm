coef.hmm <- function(object, type=c("matrix", "raw"), ...) {
    type <- match.arg(type)
    if (type=="raw") object$coefficients
    else object$beta
}

