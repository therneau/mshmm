# Automatically generated from hmmcode
# This routine is used in the test suite, to call the internal derivative()
#  routine.
derivtest <- function(time, x, beta, qmat, cmap) {
    eta <- x %*% beta
    nstate <- ncol(qmat)
    ntransition <- sum(qmat > 0)
    if (ncol(eta) != ntransition) stop("ncol(beta) != ntransition")
    rindex <- which(qmat > 0)
    rmat <- matrix(0., nstate, nstate)
    rmat[qmat>0] <- exp(eta)
    diag(rmat) <- -rowSums(rmat)
    
    # The derivative oF Rmat wrt each coefficint
    nonzero <- (cmap >0)
    coeff <- unique(cmap[nonzero])
    dmat <- array(0., c(nstate, nstate, length(coeff)))
    rr <- row(cmap)[nonzero]
    cc <- col(cmap)[nonzero]
    kk  <- cmap[nonzero]
    for (i in seq_along(coeff)) {
        temp <- matrix(0, nstate, nstate)
        for (j in  which(kk == coeff[i])) {
            indx <- rindex[cc[j]]  #cc[j] = which linear predictor
            temp[indx] <- rmat[indx] * x[rr[j]]
        }
        
        diag(temp) <- -rowSums(temp)
        dmat[,,i] <- temp
    }
    derivative(rmat, time, dmat)
}

decomp <- function(rmat, time, eps=1e-6) {
    delta <- diff(sort(diag(rmat)))
    if (any(delta < eps) || any(rmat[row(rmat) > col(rmat)] >0)) 
        stop("invalid matrix")
    else .Call("cdecomp", rmat, time)
}
    
derivative <- function(rmat, time, dR, eps=1e-8) {
    ncoef <- dim(dR)[3]
    nstate <- nrow(rmat)
    dlist <- decomp(rmat, time)
    
    dmat <- array(0.0, dim=c(nstate, nstate, ncoef))
    vtemp <- outer(dlist$d, dlist$d,
                   function(a, b) {
                       ifelse(abs(a-b)< eps, time* exp(time* (a+b)/2),
                         (exp(a*time) - exp(b*time))/(a-b))})
    # any unique value of cmap appears on only one row of cmap,
    #  multiple times in that row if a coefficient is shared
    # two transitions can share a coef, but only for the same X variable
    for (i in 1:ncoef) {
        G <- dlist$Ainv %*% dR[,,i] %*% dlist$A
        V <- G*vtemp
        dmat[,,i] <- dlist$A %*% V %*% dlist$Ainv
    }
    dlist$dmat <- dmat
    dlist
    }
