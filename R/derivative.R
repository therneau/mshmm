# The user visible version of these is in survexpm.R
# This should never be called with a 'bad' rmat, that will have been screened
#  earlier.
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
        V <- G*vtemp   # elementwise multiplication
        dmat[,,i] <- dlist$A %*% V %*% dlist$Ainv
    }
    dlist$dmat <- dmat
    dlist
    }
