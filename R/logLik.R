##
## Get the log-likelihood
## ======================
## - `df` is the number of free parameters in the model
## - We include an `nstate` attribute
## - Do we want to come up with an `nobs` value?
logLik.icmsh <- function(object, ...){
    out <- unname(object$loglik["final"])
    attr(out, "n") <- object$n
    attr(out, "df") <- length(object$coef)
    attr(out, "nstate") <- object$nstate
    class(out) <- "logLik"
    out
}
