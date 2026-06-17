predict.cmsh <- function(object, newdata, type=c("link", "rate", "pstate"),
                        se.fit=FALSE, time, p0, absorb, ...) {
    if (!inherits(object, "cmsh")) stop("only valid for icmsh objects")
    type <- match.arg(type)
    Terms <- delete.response(terms(object))

    if (se.fit && type=="pstate") {
        warning("se.fit not available for pstate")
        se.fit <- FALSE
    }

    if (type %in% c("link", "rate")) {
        # compute eta
        if (missing(newdata) || is.null(newdata)) 
            X <- model.matrix(object)
        else X <- model.matrix(object, data=newdata)
        eta <- X %*% coef(object, matrix=TRUE)
    
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
                        # length(indx) >1 is very odd -- two covariates forced 
                        # to have the same coef.  It will likely never happen,
                        # but...
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
        } else { # no se
            if (type == "link") eta
            else if (type=="rate") {
                keep <- 1:sum(object$qmat > 0)        
                exp(eta[,keep])
            }
        }
        
    } else if (type=="pstate") {
        # For a probability in state curve the data is restricted to be a single
        #  set of sequential time, along with its set of time dependent
        #  covariates.
        nstate <- length(object$states)
        if (missing(time)) stop("type pstate requires a times argument")
        if (missing(p0)) stop("type pstate requires a p0 argument")
        if (length(p0) != length(object$states)) 
            stop("p0 must have an element for each state")
        if (any(p0<0 | p0 >1)) stop("p0 values must be between 0 and 1")
        p0 <- p0/sum(p0)
        if (!missing(absorb)) {
            if (!is.numeric(absorb) || any(absorb!= floor(absorb)) ||
                any(absorb<1 | absorb>nstate))
                stop("invalid value for absorb")
        } else absorb <- NULL

        if (is.vector(time)) dtime <-  data.frame(time= time)
        else dtime <- as.data.frame(time)
        # dtime is now a data frame with one or more columns
        if (any(sapply(dtime, function(x) !is.numeric(x) || any(diff(x) <=0))))
                stop("time must be numeric, in strictly increasing order")
        if (ncol(dtime) >1) { # all columns must be in lockstep
            for (i in 2:ncol(dtime))
                if (!isTRUE(all.equal(diff(dtime[,i]), diff(dtime[,1]))))
                    stop("all columns of 'time' must have idential increments")
        }
        
        # By default, model.matrix only grabs the state coefficients
        # We will get a new eta for each curve
        newdata <- as.data.frame(newdata)  # in case they gave a list
        ncurve <- nrow(newdata)
        ntime <- nrow(dtime)
        nstate <- length(object$states)
        pstate <- array(0, dim=c(ntime, ncurve, nstate))
        rmat <- matrix(0.0, nstate, nstate)  #working matrix
        qmap <- which(object$qmatrix != 0)
        cumhaz <- array(0, dim=c(ntime, ncurve, length(qmap)))
        delta <- diff(dtime[,1])
        beta <- coef(object, matrix=TRUE)

        for (k in 1:ncurve) {
            pstate[1,k,] <- p0
            dummy <- cbind(dtime, newdata[rep(k,ntime),,drop=FALSE]) #data frame
            eta <- model.matrix(object, data=dummy) %*% beta
            for (i in 1:(ntime-1)) { 
                rmat[qmap] <- exp(eta[i,])
                if (!is.null(absorb)) rmat[absorb,] <- 0
                diag(rmat) <- diag(rmat) - rowSums(rmat)
                pstate[i+1,k,] <- pstate[i,k,] %*% survexpm(rmat, delta[i])
                cumhaz[i+1,k,] <- cumhaz[i,k,] + rmat[qmap]* delta[i]
            }
        }

        dimnames(pstate) <- list(NULL, NULL,object$states)
        dimnames(cumhaz) <- list(NULL, NULL,
                               paste(row(rmat)[qmap], col(rmat)[qmap], sep=':'))

        if (ncurve==1) {
            pstate <- pstate[,1,]
            cumhaz <- cumhaz[,1,]
        }
        ret <- list(time= dtime[,1], pstate= pstate, cumhaz= cumhaz, 
                    tname=names(dtime), data=newdata, states=object$states)
        class(ret) <- "pstate.cmsh"
        ret
    }
}

    
