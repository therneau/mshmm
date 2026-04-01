predict.icmsh <- function(object, newdata, type=c("link", "rate", "pstate"),
                        se.fit=FALSE, times, absorb, ...) {
    if (!inherits(object, "icmsh")) stop("only valid for icmsh objects")
    type <- match.arg(type)
    Terms <- delete.response(terms(object))

    if (missing(newdata) || is.null(newdata)) 
        X <- model.matrix(object)
    else X <- model.matrix(object, data=newdata)

    if (!is.null(object$fit$accept)) 
        return(predict_hmcmc(object, X, type, se.fit, ...))

    eta <- X %*% coef(object, type="matrix")
    
    if (se.fit && type=="pstate") {
        warning("se.fit not available for pstate")
        se.fit <- FALSE
    }
    if (se.fit) {
        se <- 0* eta   # inherit the dimnames        
        ncoef <- length(object$coef)
        if (is.null(object$fit$hessian)) {
            # warning("std.err=TRUE and model has no hessian")
            vmat <- solve(object$fit$S)
        }
        else vmat <- solve(object$fit$hessian)
        for (j in 1:ncol(eta)) {
            Z <- matrix(0., nrow(X), length(object$coef))
            for (k in unique(object$cmap[,j])) {
                if (k>0) {
                    indx <- which(object$cmap[,j] ==k)
                    # length(indx) >1 is very odd -- two covariates forced to
                    #   have the same coef.  It will likely never happen, but...
                    if (length(indx) ==1) Z[,k] <- X[,indx]
                    else Z[,k] <- rowSums(X[,indx])
                }
            }
            se[,j] <- sqrt(diag(Z %*% vmat %*% t(Z)))
        }       
        if (type== "link") return(list(fit=eta, se.fit=se))
        else if (type=="rate") {
            # only return the rate parameters
            keep <- 1:sum(object$qmat > 0)            
            return(list(fit=exp(eta[,keep]), 
                        se.fit=exp(eta[keep])*se[,keep]))
        }
        else stop("std is not avaiable for type=", type) # future proof
    }

    if (type == "link") eta
    else if (type=="rate") {
        keep <- 1:sum(object$qmat > 0)        
        exp(eta[,keep])
    }
    else if (type=="pstate") {
        # For a probability in state curve the data is restricted to be a single
        #  set of sequential times, along with its set of time dependent
        #  covariates.
        if (missing(times)) stop("pstate requires a times argument")
        if (!is.numeric(times) || any(diff(times) <=0))
            stop("times must be numeric, in increasing order")
        ntime <- length(times)
        if (ntime-1 != nrow(X))
            stop("times vector must match dimension of new data")

        # The initial p(t) is assumed to hold at the smallest Y time
        nstate <- object$nstate
        pstate <- matrix(0., nrow= ntime, ncol= nstate)
        temp <- cumsum(object$nlp)
        p.param <- (temp[2]:temp[3])[-1]  # which eta columns for initial p?
        if (length(p.param ==0)) pstate[1,] <- object$pfun(nstate)
        else pstate[1,] <- object$pfun(nstate, eta[,p.param])
 
        rmat <- matrix(0.0, nstate, nstate)  #working matrix
        qmap <- which(object$qmatrix != 0)
        keep <- 1:length(qmap)   # linear predictors that map to rates
        dtime <- diff(times)
        for (i in 1:(ntime-1)) { 
            rmat[qmap] <- exp(eta[i,keep])
            diag(rmat) <- diag(rmat) - rowSums(rmat)
            pstate[i+1,] <- pstate[i,] %*% expm(dtime[i]* rmat)
        }
        dimnames(pstate) <- list(times, dimnames(object$qmatrix)[[1]])
        pstate
    }
}

predict_hmcmc <- function(object, X, type, se.fit, ...) {
    if (type=="pstate") stop("type = pstate not yet done for MCMC")

    # Get the entire array of predictions
    # It's easiest to do it for one column of coef(object, type="matrix")
    #  at at time
    beta <- coef(object, type="matrix")
    nsim <- nrow(object$fit$par)
    yhat <- array(0., dim=c(nrow(X), ncol(beta), nsim))
    
    for (j in 1:ncol(beta)) {
        map <- object$cmap[,j]
        btemp <- matrix(beta[,j], nrow=nrow(beta), ncol=nsim) 
        btemp[map>0, ] <- t(object$fit$par[,map])
        yhat[,j,] <- X %*% btemp
    }

    if (se.fit) se.sim <- drop(apply(yhat, 1:2, sd))
    meansim <- drop(apply(yhat, 1:2, mean))
    
    if (type=="rate") {
        keep <- 1:sum(object$qmat > 0)            
        if (se.fit) list(fit= exp(meansim[,keep]), 
                         se.fit= (se.sim * exp(meansim))[,keep])
        else meansim[,keep]
    }
    else if (type=='link') {
        if (se.fit) list(fit=meansim, se.fit= se.sim)
        else meansim
    }
    else stop("unknown type")
}
                         
             
      
