#
# exponential of a matrix, optionally with derivatives
# If derivatives are not needed and upper=FALSE just call the expm function
# 
# R:   a square transition matrix, row sums will be 0, each non-zero element 
#   R[i,j] is assumed to be exp(eta[i,j]), i.e., positive
# time: width of the time interval
# deriv: 0= no derive, 1=compute the derivatives wrt each theta, for theta =
#   each non-zero element in turn, 2 assume theta= exp(eta) and return 
#   derivatives wrt eta
# tol:  if the inverse condition number of the eigenmatrix is < tol, use the
#  pade method
#
survexpm <- function(R, time=1.0, deriv=FALSE, tol=1e-10, 
                     method=c("eigen", "pade")) {
    method <- match.arg(method)

    # R is a transition matrix, so the diagonal elements are 0 or negative
    # First a special case which should never arise, but ...
    if (length(R)==1) {  # a 1 by 1 matrix exp[t* exp(eta)]
        P <- matrix(exp(R*time),1,1) 
        if (deriv) 
            return(list(P=P, dmat= array(R*time*exp(R*time), dim=c(1,1,1))))
        else return(P) 
     }
   
    # The case where all values with R>0 lie in a single row or a single
    #  column have simple closed form solutions, but are so rare in this
    #  library as to not be worth the bother.  Very common in survival though.
    #if (!deriv) return(expm(R*time))  # just use expm funtion

    upper <- all(R[row(R) > col(R)] ==0) # upper triangular
    if (upper) {
        # if there is a tied eigenvalue, we have to use pade, but ignore zeros
        eigen <- diag(R)
        e0 <- eigen[abs(eigen) > .Machine$double.eps]
        #if (length(e0) >1 && any(diff(sort(e0)) < tol)) method <- "pade"
    }

    if (method=="eigen") {
        # get the eigen decomp
        storage.mode(R) <- "double"   # failsafe
        efit <- .Call("hmmeigen", R)
        # We currently do most of the work in R, in order to better understand
        #  the algorithm. The hmmeigen routine returns eigenvalues, left and
        #  right eigenvectors
        scale <- 1/colSums(Conj(efit$left)* efit$right) 
        iright <- scale* t(Conj(efit$left)) # inverse of right eigenvectors
        # If too close to singular the eigen algorithm is not accurate
        #  expm-eigen.c uses .Machine$double.eps, which I think is too forgiving
        #if (any(Mod(scale) > 1/tol)) method <- "pade"
        cond <- norm(iright)*norm(efit$right)  # more formally correct
        if (1/cond < tol) method <- "pade" 
    }

    if (deriv>0 ) { # the common case
        npos <- sum(R>0) # number of positive elments in R
        nstate <- nrow(R)
        # create the target array for the derivatives, for each eta that we
        #  need to consider. dR = derivative of the R matrix wrt an element,
        #  see 'matrix exponential, eigenvector formula' in the code vignette
        dR <- array(0, dim=c(nstate, nstate, npos))
        indx <- cbind(row(R)[R>0], col(R)[R>0])
        if (deriv==2) temp <- R[which(R>0)]
        else temp <- rep(1, npos)
        for (i in 1:npos) dR[indx[i,1], indx[i,], i] <- c(- temp[i], temp[i])
    }

    if (method == "eigen") {
        right <- efit$right
        P <- Re(right %*% diag(exp(time* efit$values)) %*% iright)
        if (deriv) {
            dmat <- array(0.0, dim=c(nstate, nstate, npos))
            vtemp <- outer(efit$values, efit$values, 
                           function(a, b) {
                               ifelse(abs(a-b)< tol, time* exp(time* (a+b)/2),
                               (exp(a*time) - exp(b*time))/(a-b))})
            for (i in 1:npos) {
                G <- iright %*% dR[,,i] %*% right
                V <- G*vtemp
                dmat[,,i] <- Re(right %*% V %*% iright)
            }
            list(P=P, deriv=dmat, method="eigen")
        } else P
    } else {
        # use the Pade method
        if (deriv) {
            fit <-pade(R*time, dR*time)
            fit$method <- "pade"
            fit
        }
        else pade(R*time)$P
    }    
}

#
# This is largely for testing
#
hmmeigen <- function(x) {
    if (!is.matrix(x) || !nrow(R)== ncol(x) || !is.numeric(x))
        stop("argument must be a square numeric matrix")
    storage.mode(x) <- "double"
    .Call("hmmeigen", x)
}
