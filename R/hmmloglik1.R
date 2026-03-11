# This next set of functions hmmloglik, hmmderiv, hmmboth will normally be
#  called via a maximization function; hmmloglik is also called by msh.fit  
#    hmmloglik: just the loglik, no derivatives
#    hmmgrad  : just the derivatives
#    hmmboth  : loglik and derivative
# These functions in turn call hmm1 (loglik) and hmm2 (loglik and deriv)
#   to do the actual work, they are called for each separate id, in parallel.
# Because of the call chain msh.fit -> hmmscore -> hmmboth -> hmm2 for 
#   instance, we need to pass everthing that hmm1/hmm2 need all the way down
#   the chain. I originally had these all functions' code within
#   the outer {} of hmm itself, which allows arguments to be found via lexical 
#   scoping, but code simply became too unweildy to be managable.
# "logfun" is a copy of hmm1 or hmm2, properly scoped, and passed down the chain
#  Those two get their first two arguments "who" (the id for which to compute
#  a loglik) and "beta" from here, all others are found in the msh.fit frame
# 
# A reminder: B is the matrix of parameters, of the same shape as cmap; it will
#  have zeros for covariates that are not used in a given linear predictor.
#  "param" is the vector of parameters from the maximizer, which does not
#   contain fixed coefs, "coefficients" in the result does have fixed coefs.
#  \beta is what we call the vector of parameters within the math documentation
grab <- function(x, what) 
    if (what %in% names(x)) x[[what]] else NULL

hmmloglik <- function(param, logfun) {
    B2 <- coef.to.B(param, cmap, B) #copy with updated parameters
    if (mc.cores > 1) {
        if (!fork)
            mcfit <- parLapply(hmm_cluster, unique(id), logfun, B=B2) 
        else mcfit <- mclapply(unique(id), logfun, B= B2,
                               mc.cores= mc.cores, mc.set.seed=FALSE)
    }
    else mcfit <- lapply(unique(id), logfun, B=B)
    
    if (any(sapply(mcfit, is.character))) {
        # failure
        words <- sapply(mcfit, function(x) ifelse(is.character(x), x, ""))
        if (any(words == "underflow")) {
            # if (debug > 1) browser()
            # assume a bad guess from a maximizer, return a bad hit
            return(NA)
        }
        words <- words[words!=""]
        stop(words[1])
    }

    tpar <- c(param)  # used for constraints
    if (!is.null(penmat)) loglik <- loglik - sum(tpar * (penmat %*% tpar))/2
    
    if (control$debug == -1) {
        # Hand back more stuff
        alpha <- sapply(mcfit, function(x) grab(x, "alpha"))
        offset <- sapply(mcfit, function(x) grab(x, "offset"))
        loglik <- sum(log(colSums(alpha)) + offset)    
        ecount <- sapply(mcfit, function(x) grab(x, "ecount"))
        rval <- list(alpha = alpha,
                     offset = offset,
                     loglik = loglik,
                     ecount=ecount)
        tpar <- param
        if (!is.null(penmat)) rval$penalty <- sum(param* (penmat %*% param))/2
        rval
    } else sum(unlist(mcfit)) # mcfit returns a single number
}

# This function is used by the score based iteration
hmmboth <- function(param, logfun){
    B2 <- coef.to.B(param, cmap, B) #copy with updated parameters
    if (mc.cores > 1) {
        if (!fork)
            mcfit <- parLapply(hmm_cluster, unique(id), logfun, B=B2)
        else mcfit <- mclapply(unique(id), logfun, B=B2,
                               mc.set.seed=FALSE, mc.cores=mc.cores)
    }
    else mcfit <- lapply(unique(id), logfun, B=B2)

    alpha <- sapply(mcfit, function(x) sum(grab(x, "alpha")))
    offset <- sapply(mcfit, function(x) grab(x, "offset"))
    if (any(sapply(mcfit, is.character))) {
        # failure
        words <- sapply(mcfit, function(x) ifelse(is.character(x), x, ""))
        if (any(words == "underflow")) {
            # assume a bad guess from a maximizer, return a bad hit
            # the routine won't use deriv or S in this case
            return(list(loglik=-Inf, deriv=NULL, S=NULL))
        }
        else {
            words <- words[words!=""]
            stop(words[1])
        }
    }
    else loglik <- sum(log(alpha) + offset)
    
    # total number of expm calls, total number that used the pade() function
    ecount <- rowSums(sapply(mcfit, function(x) grab(x, "ecount")))
    # hmm_count_of_calls was set to (0,0) before iteration 
    # This was a question of how often tied eigenvalues show up, across
    #  iterations
    hmm_count_of_calls <<- hmm_count_of_calls + ecount

    d.alpha <- sapply(mcfit, function(x) rowSums(grab(x, "deriv")))
    u <- d.alpha * rep(1/alpha, each=nparm)
    # this will be a matrix with nparm rows, one col per subject
    deriv <- rowSums(u)
    S <- tcrossprod(u)
    S2 <- tcrossprod(u - rowMeans(u))  # a better variance when deriv!=0 ?
    
    # S is an estimate of -Hessian, so add -(second deriv) of penalty
    if (!is.null(penmat)) {
        tpar <- param
        loglik <- loglik - sum(tpar *(penmat %*% tpar))/2
        deriv  <- deriv - c(param %*% penmat)
        S <- S + penmat
        S2 <- S2 + penmat
    }
    # if (!is.null(conmat)) {
    #     temp <-  c(conmat %*% c(1, param))  # vector of values
    #     loglik <- loglik + sum(log(pmax(temp,0)))
    #     first <- colSums(conmat[,-1]/ rep(temp, nparm))
    #     deriv <- deriv + first
    #     S <- S + outer(first, first)
    #     if (debug>1) {cat("in hmmboth\n"); browser()}
    # }
    
    list(loglik = loglik, deriv= deriv, S=S, S2=S2)
}

# This is used by optim
hmmgrad <- function(param, grfun) {
    B2 <- coef.to.B(param, cmap, B) #copy with updated parameters
    if (mc.cores > 1) {
        if (!fork)
            mcfit <- parLapply(hmm_cluster, unique(id), grfun, B=B2)
        else mcfit <- mclapply(unique(id), grfun, B=B2,
                               mc.set.seed=FALSE, mc.cores=mc.cores)
    }
    else mcfit <- lapply(unique(id), grfun, B=B2)
    
    if (any(sapply(mcfit, is.character))) {
        # failure
        words <- sapply(mcfit, function(x) ifelse(is.character(x), x, ""))
        if (any(words == "underflow")) {
            # assume a bad guess from a maximizer, return a bad hit
            loglik <- -2 * abs(initial.loglik)
            return(rep(0, length(param)))
        }
        else {
            words <- words[words!=""]
            stop(words[1])
        }
    }

    alpha <- sapply(mcfit, function(x) sum(grab(x, "alpha")))
    d.alpha <- sapply(mcfit, function(x) rowSums(grab(x, "deriv")))
    # this will be a matrix with npar rows, one col per subject
    u <- d.alpha * rep(1/alpha, each=nparm)

    if (is.null(penmat)) deriv <- rowSums(u)
    else deriv <- rowSums(u) - c(param %*% penmat)
    #if (!is.null(conmat)) {
    #    temp <- conmat %*% c(1, param)
    #    deriv <- deriv + colSums(conmat[,-1]/temp)
    #}
    deriv
}
