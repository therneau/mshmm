    
# This next set of functions hmmloglik, hmmderiv, hmmboth will normally be
#  called via a maximization function,  hmmloglik or hmmdebug are called
#  directly by hmm for the zeroth iteration.
#    hmmloglik: just the loglik, no derivatives
#    hmmgrad  : just the derivatives
#    hmmboth  : loglik and derivative, plus a bit more
#    hmmdb: for debugging, pass back everything
# These functions in turn call hmm1 (loglik) and hmm2 (loglik and deriv)
#   to do the actual work, they are called for each separate id, in parallel.
# Because of the call chain hmm -> hmmscore -> hmmboth -> hmm2 for 
#   instance, we need to pass everthing that hmm1/hmm2 need all the way down
#   the chain. I originally had these all functions' code within
#   the outer {} of hmm itself, which allows arguments to be found via lexical 
#   scoping, but code simply became too unweildy to be managable.
# "logfun" is a copy of hmm1 or hmm2, properly scoped, and passed down the chain
#  Those two get their first two arguments "who" (the id for which to compute
#  a loglik) and "beta" from here, all others are found in the hmm frame
# 
# A reminder: B is the matrix of parameters, of the same shape as cmap; it will
#  have zeros if some covariates that are not used for all transitions.
#  "param" is the label used by the maximizer, for the vector of parameters
#  "coefficients" is what it is labeled in the result, to match lm, glm, etc.
#  \beta is what we call the vector of parameters within the math documentation

hmmloglik <- function(param, B, cmap, id, mc.cores, fork, logfun) {
    B[cmap>0] <- param[c(cmap)]   #"param[cmap]" fails if cmap has 2 columns
    if (mc.cores > 1) {
        if (!fork)
            mcfit <- parLapply(hmm_cluster, unique(id), logfun, B=B) 
        else mcfit <- mclapply(unique(id), logfun, B= B,
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

    loglik <- sum(unlist(mcfit))
    tpar <- c(param)  # used for constraints
    if (!is.null(penmat)) loglik <- loglik - sum(tpar * (penmat %*% tpar))/2
    
    # This might be added back at a later date
    #if (!is.null(conmat)) {
    #    temp <-  conmat %*% tpar
    #    loglik <- loglik + sum(log(pmax(temp,0)))# -Inf if there are violations
    #}
    loglik
}

# hand back everything (debug).
# Called as the zero iteration rather than hmmloglik during the debugging phase
#  it returns a list with one element per id.
hmmdb <- function(param, B, cmap, id, mc.cores, fork, logfun) {
    B[cmap>0] <- param[c(cmap)]
    if (mc.cores > 1) {
        if (!fork)
            mcfit <- parLapply(hmm_cluster, unique(id), logfun, B= B) 
        else mcfit <- mclapply(unique(id), logfun, B= B,
                               mc.cores= mc.cores, mc.set.seed=FALSE)
    }
    else mcfit <- lapply(unique(id), logfun, B=B)

    alpha <- sapply(mcfit, function(x) grab(x, "alpha"))
    offset <- sapply(mcfit, function(x) grab(x, "offset"))
    if (any(sapply(mcfit, is.character))) {
        # failure
        words <- sapply(mcfit, function(x) ifelse(is.character(x), x, ""))
        if (any(words == "underflow")) {
            # assume a bad guess from a maximizer, return a bad hit
            loglik <- -Inf
            mcfit <- mcfit[which(words=="")]  # toss the bad results and go on
        }
        else {
            words <- words[words!=""]
            stop(words[1])
        }
    }
    else loglik <- sum(log(colSums(alpha)) + offset)
    
    ecount <- sapply(mcfit, function(x) grab(x, "ecount"))
    dd <- dim(mcfit[[1]]$deriv)
    rval <- list(alpha = alpha,
                 offset = offset,
                 loglik = loglik,
                 deriv = array(unlist(lapply(mcfit, function(x) grab(x, "deriv"))),
                               dim=c(dd, length(mcfit))),
                 ecount=ecount)
    tpar <- param
    if (!is.null(penmat)) rval$penalty <- sum(tpar* (penmat %*% tpar))/2
    #if (!is.null(conmat)) {
    #    temp <- tpar %*% conmat
    #    rval$constraint <- log(pmax(temp,0))
    #}
    rval
}

# This function is used by the score based iteration
hmmboth <- function(param, x, ytime, ystate, id, rindex, B, cmap, logfun) {
    B[cmap>0] <- param[c(cmap)]
    if (mc.cores > 1) {
        if (!fork)
            mcfit <- parLapply(hmm_cluster, unique(id), logfun, B=B)
        else mcfit <- mclapply(unique(id), logfun, B=B,
                               mc.set.seed=FALSE, mc.cores=mc.cores)
    }
    else mcfit <- lapply(unique(id), hmm2, B=B)

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
hmmgrad <- function(param, x, ytime, ystate, id, rindex, B, cmap, 
                    grfun) {
    B[cmap>0] <- param[c(cmap)]
    if (mc.cores > 1) {
        if (!fork)
            mcfit <- parLapply(hmm_cluster, unique(id), grfun, B=B)
        else mcfit <- mclapply(unique(id), grfun, B=B,
                               mc.set.seed=FALSE, mc.cores=mc.cores)
    }
    else mcfit <- lapply(unique(id), grfun, B=B)
    
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
