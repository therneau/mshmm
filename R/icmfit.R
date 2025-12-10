# Given all the data and setup, do the actual fit of an interval-censored
#  multi-state hazard model.  
# Compared to hmmfit, this does not have to worry about response functions
#  or initial estimates.  At this point, the same maximizers are used.

icmfit <- function(ytime, ystat, X, fid, otype, qmatrix, cmap, B, 
                  mfun, mpar, iter, constraint, penalty, mc.cores,      
                  control) {
    #Give this next variable a long name that won't be found in calling
    #  routines.  It is updated farther down the calling chain.
    hmm_count_of_calls <- c(0, 0)  #total calls to expm, number with tied eigens

    # set up for parallel
    if (mc.cores > 1 && !fork)
        hmm_cluster <- makeCluster(mc.cores) #start up parallel

    time1 <- proc.time()
    if (iter>0) {
        # compute and save the loglik for the initial parameters
        beta <- B.to.coef(B, cmap)
        fit0 <- hmmloglik()
        }


    time1 <- proc.time()

    if (missing(mfun)) {
        # This is a call with no iteration, use hmmloglik function
        initial.loglik <- numeric(0)
        fit <- do.call(mfun, mpar)
        initial.loglik <- penalty0 <- NULL
    }

    else {
        # Get an initial loglik
        # If the maximizer fails (NA or infinite) it returns the initial;
        #  set it to a dummy value to detect that
        initial.loglik <- numeric(0)
        initial.loglik <- hmmloglik(param)
        if (length(initial.loglik)==0)
            stop("unable to evalutate the likelihood at the intial parameters")

        # Compute the intial penalty too
        if (!is.null(penmat))
            penalty0 <- sum(param *(penmat %*% param))/2
        else penalty0 <- 0

        fit <- do.call(mfun, mpar)
    }
 
    
}    


hmmloglikc <- function(param, ...) {
    beta[cmap>0] <- param[c(cmap)]   # param[cmap] =bad if cmap has 2 columns
    if (mc.cores > 1) {
        if (!fork)
            mcfit <- parLapply(hmm_cluster, 1:nid, hmm1, beta=beta, ...)
        else mcfit <- mclapply(1:nid, hmm1, beta=beta, 
                               mc.set.seed=FALSE, mc.cores=mc.cores, ...)
    }
    else mcfit <- lapply(1:nid, hmm1, beta=beta, ...)
    
    if (any(sapply(mcfit, is.character))) {
        # failure
        words <- sapply(mcfit, function(x) ifelse(is.character(x), x, ""))
        if (any(words == "underflow")) {
            if (debug > 1) browser()
            # assume a bad guess from a maximizer, return a bad hit
            return(-2 * abs(initial.loglik))
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
    count <- sapply(mcfit, function(x) attr(x, "counts"))
    # this next line reaches back and changes the variable in parent of parent
    #  (the parent called the maximizer, which calls this)
    hmm_count_of_calls <<- hmm_count_of_calls + rowSums(count)
    loglik
}

# This function is used by the score based iteration
hmmboth <- function(param, ...) {
    beta[cmap>0] <- param[c(cmap)]
    if (mc.cores > 1) {
        if (!fork)
            mcfit <- parLapply(hmm_cluster, 1:nid, hmm2, beta=beta)
        else mcfit <- mclapply(1:nid, hmm2, beta=beta, 
                               mc.set.seed=FALSE, mc.cores=mc.cores)
    }
    else mcfit <- lapply(1:nid, hmm2, beta=beta)

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
    
    ecount <- rowSums(sapply(mcfit, function(x) grab(x, "ecount")))
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

hmmgrad <- function(param, ...) {
    beta[cmap>0] <- param[c(cmap)]
    if (mc.cores > 1) {
        if (!fork)
            mcfit <- parLapply(hmm_cluster, 1:nid, hmm2, beta=beta)
        else mcfit <- mclapply(1:nid, hmm2, beta=beta, 
                               mc.set.seed=FALSE, mc.cores=mc.cores)
    }
    else mcfit <- lapply(1:nid, hmm2, beta=beta)
    
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

    ecount <- rowSums(sapply(mcfit, function(x) grab(x, "ecount")))
    hmm_count_of_calls <<- hmm_count_of_calls + ecount
    
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
