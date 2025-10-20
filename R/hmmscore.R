# Automatically generated from hmmcode
hmmscore <- function(par, fn, gr, iter=30,
                     scale=2, shrink=1, debug=FALSE,
                     constraint, eps=1e-6, hessian =FALSE, ...) {
    if (scale <1.5) stop("invalid value for scale")
    if (iter < 1) stop("iteration count must be 1 or more")

    ncoef <- length(par)
    if (missing(constraint)) ncon <- 0
    else ncon <- nrow(constraint)

    # Some tracing information: the coefficients, loglik, and LM scale
    #  as we go along, and the active set if necessary.
    cmat <- matrix(0, iter+1, ncoef)
    cmat[1,] <- par
    logmat <- matrix(0, iter+1, 4,
               dimnames=list(0:iter, c("loglik", "LM scale", "lambda", "step")))
    lm <- 4  # starting point
    if (ncon>0) actmat <- matrix(0L, iter+1, ncon)
    
    # initial step
    fit <- gr(par)
    logmat[1,] <- c(fit$loglik, 0, 0,0)
    if (!is.list(fit) || is.null(fit$S))
        stop("hmmscore needs to call hmmboth")
    npar <- length(par)

    #  A simple solve() will blow up if the S +dmat is too close to
    #   singular; use a tempered svd that ignores the near zeros
    #  Solve using a generalized inverse -- the algorithm is essentially that
    #   of the ginv function in MASS
    pfun <- function(coef, deriv, S, dmat, shrink, tol=1e-9) {
        smat <- svd(S + dmat) 
        dpos <- (smat$d > max(smat$d[1]*tol, 0))
        dd <- ifelse(dpos, 1/smat$d, 0)
        # all the parentheses save a tiny bit of time 
        if (all(dpos)) x <- drop(smat$v %*% (dd*(t(smat$u) %*% deriv)))
        else if (!any(dpos)) stop("zero hessian in update") # impossible I think
        else x <-drop(smat$v[,dpos, drop=FALSE] %*%(dd[dpos] * 
                      (t(smat$u[,dpos, drop=FALSE]) %*% deriv)))
 
        coef + shrink * x
    }

    # The usual LM would start at diag(S *.001), we are more conservative
    lm <- 1/4
    log0 <- fit$loglik  # the initial loglik
    oldpar <- par
    halving <- 1
    for (i in 1:iter) {
        if (halving < 1/64) break # stuck: 6 iters in a row with no progress
        S <- fit$S
        dS <- diag(S)
 
        dmat <- diag((dS + max(dS)/1e3)*lm) # the smallest 

        if (ncon>0) {
            Ebeta <- as.vector(constraint[,-1] %*% par)
            boundary <- shrink*halving*(Ebeta + constraint[,1])
            temp <- solve.QP(S+dmat, -fit$deriv, -t(constraint[,-1]),
                             -boundary)
            newpar <- par- temp$solution*shrink*halving
            active <- sort(temp$iact)
            actmat[i+1, active] <- 1L  #which constraints are active
        }
        else {
            newpar <- pfun(par, fit$deriv, S, dmat, shrink * halving)
            active <- NULL  # no active constraints
            }
        
        cmat[i+1,] <- newpar
        newfit <- gr(newpar)
        #cat("score stop\n"); browser()
        
        if (newfit$loglik <= -2*abs(log0)) {
            # fn() returns -Inf when the computation fails due
            #  to overflow/underflow, a failure of the expm function.
            # In this case newfit$deriv will be null.  Step halving
            #  is our best solution.  This should be very rare.
            #  
            logmat[i+1,] <- c(newfit$loglik, lm, NA, -1) #bombed
            lm <- lm * scale
            halving <- halving / 2
        }
        else {
            # Use a quadratic approx to find a (possibly) improved guess
            # Equation 9.7.11 of Numerical Recipes in C, 1992
            gder <- sum(fit$deriv * (newpar -par))  # directional deriv
            lambda <- 0.5*gder/(fit$loglik + gder - newfit$loglik)
            logmat[i+1,] <- c(newfit$loglik, lm, lambda, halving)#standard step
  
            if (abs((newfit$loglik - fit$loglik)/fit$loglik) < eps) {
                trdata <- list(coef=cmat[1:(i+1),], loglik=logmat[1:(i+1),])
                if (newfit$loglik > fit$loglik) # very last step was better
                    rval <- c(newfit, list(coef=newpar, iter=i, converged=TRUE,
                                          trace=trdata))  #all done
                else {
                    logmat[i+1, 4] <- 0  # no step
                    rval <- c(fit, list(coef=par, iter=i, converged=TRUE,
                                       trace=trdata))
                }
                if (hessian) {  # compute the hessian matrix
                    # hessian contains the epsilon, usually
                    if (is.logical(hessian)) epsilon <- 1e-3
                    else  epsilon <- hessian
                    hmat <- matrix(0., ncoef, ncoef) 
                    for (i in 1:ncoef) {
                        delta <- rep(0., length(par))
                        delta[i] <- epsilon
                        temp1 <- gr(par + delta)
                        temp2 <- gr(par - delta)
                        hmat[i,] <- (temp2$deriv - temp1$deriv)/(2* epsilon)
                        if (debug>1) cat(" hessian ", i)
                    }
                    rval$hessian <- (hmat + t(hmat))/2  # symmetrize
                }
                return(rval)
            }

            if (debug > 1) browser()
            if (newfit$loglik < fit$loglik) {
                lm <- if (lm>scale) lm*scale else lm*scale^2 # increase penalty
                logmat[i+1, 4] <- 0  # no progress
                if (debug > 1) {cat("in hmmscore\n"); browser()}
                
                if (gder >0 & lambda>0 & lambda <1) { # backtrack
                    try2 <- (1-lambda)*par + lambda*newpar
                    fit2 <- gr(try2)  # we expect this to usually work
                    if (fit2$loglik > fit$loglik) { # shrinkage worked
                        logmat[i+1,4] <- fit2$loglik
                        fit <- fit2
                        par <- try2
                    } else halving <- halving /2  # and smaller steps
                }    
            }
            else { #successful step
                fit <- newfit
                par <- newpar
                lm <- lm/scale
                halving <- min(1, halving*1.5)  #recover slowly
            }
            if (debug) {
                cat("iter", i, "loglik", logmat[i+1,], "\n")
                if (ncon>0) cat("   active", active, '\n')
                }
        }
    }

    if (ncon>0) trdata <- list(coef=cmat, loglik=logmat, active=actmat)
    else trdata <- list(coef=cmat, loglik=logmat)

    rval <- c(fit, list(coef=par, iter= iter, converged = FALSE, trace=trdata))

    if (hessian) {  # compute the hessian matrix, same algorithm as optimHess
        # hessian contains the epsilon, usually
        if (is.logical(hessian)) epsilon <- 1e-3
        else  epsilon <- hessian
        hmat <- matrix(0., ncoef, ncoef) 
        for (i in 1:ncoef) {
            delta <- rep(0., length(par))
            delta[i] <- epsilon
            temp1 <- gr(par + delta)
            temp2 <- gr(par - delta)
            hmat[i,] <- (temp2$deriv - temp1$deriv)/(2* epsilon)
            if (debug>1) cat(" hessian ", i)
        }
        rval$hessian <- (hmat + t(hmat))/2  # symmetrize
    }
    rval
}        
