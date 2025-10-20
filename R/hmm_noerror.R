#
# An HMM response function that leads to a Markov model
#
hmm_noerror <- function(y, nstate, ...) {
    temp <- diag(nstate)
    if (any(y < 1 | y>nstate | floor(y) !=y))
        stop("y must be an integer between 1 and number of states")
    temp[,y, drop=FALSE]
}

# A maximization function that does no iteration   
hmmtest <- function(par, fn, ...) {
    fit <- fn(par)
    if (is.list(fit)) c(fit, list(coef=par))
    else list(loglik=fit, coef=par)
}
