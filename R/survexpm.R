survexpminit <- function(rmat) {
    # check the validity of the transition matrix, and determine if it
    #  is acyclic, i.e., can be reordered into an upper triangular matrix.
    if (!is.matrix(rmat) || nrow(rmat) != ncol(rmat) || any(diag(rmat) > 0) ||
        any(rmat[row(rmat) != col(rmat)] < 0))
        stop ("input is not a transition matrix")
    temp <- all.equal(rowSums(rmat), rep(0, ncol(rmat)), check.attributes=FALSE)
    if (!isTRUE(temp)) stop("input is not a transition matrix")
    nc <- ncol(rmat)
    lower <- row(rmat) > col(rmat)
    if (all(rmat[lower] ==0))  return(0)  # already in order
    
    # score each state by (number of states it follows) - (number it precedes)
    # DOI: 10.1007/978-3-662-48971-0_15 show that in general, determining if
    # a matrix can be permuted to upper triangular is hard and give an 
    # exponential time algorithm.  The crude algorithm below can get lucky if
    # the transition matrix is sparse, which many are, but it is best if the
    # user orders states in a way that makes it easy.
    ztemp <- 1*(rmat >0) # 0/1 matrix
    indx <- order(colSums(ztemp) - rowSums(ztemp))
    temp <- rmat[indx, indx]  # try that ordering
    browser()
    if (all(temp[lower]== 0)) indx  # it worked!
    else -1  # no ordering found: there is a loop in the states (or we failed)
}

survexpm <- function(R, time=1.0, setup, deriv=FALSE, all=FALSE, eps=1e-6) {
    # R is a transition matrix, so the diagonal elements are 0 or negative
    # First a special case which should never arise, but users...
    if (length(R)==1) {  # a 1 by 1 matrix
        P <- matrix(exp(R*time),1,1) 
        if (deriv) 
            return(list(P=P, dmat= array(time*exp(R*time), dim=c(1,1,1))))
        else return(P) 
     }
    
    # The case where all values with R>0 lie in a single row or a single
    #  column are important in the survival library, here they are so rare
    #  in the censored case to not be worth the bother
    npos <- sum(R>0) # number of positive elments in R
    nstate <- nrow(R)
    if (deriv>0) {
        # create the target array for the derivatives, for each eta that we
        #  need to consider. dR = derivative of the R matrix wrt an element,
        #  see 'matrix exponential, eigenvector formula' in the code vignette
        dR <- array(0, dim=c(nstate, nstate, npos))
        indx <- cbind(row(R)[R>0], col(R)[R>0])
        rpos <- which(R>0)
        if (deriv==1) {
            for (i in 1:npos) dR[indx[i,1], indx[i,], i] <- c(-1,1)*R[rpos[i]]
        } else {
            for (i in 1:npos) dR[indx[i,1], indx[i,], i] <- c(-1,1)
        }
    }

    if (!missing(setup) && ((setup[1]==0 && all(diff(sort(diag(R))) >eps))
        || (setup[1]>0 && all(diff(sort(diag(R[setup,setup]))) >eps)))) {
        # use the Kalbfleisch and Lawless decomposition
        #  start with a decomp of the upper triangular matrix
        if (setup[1]==0) dlist <- .Call("cdecomp", R, time)
        else dlist <- .Call("cdecomp", R[setup, setup], time)

        if (deriv) {
            dmat <- array(0.0, dim=c(nstate, nstate, npos))
            vtemp <- outer(dlist$d, dlist$d,
                           function(a, b) {
                               ifelse(abs(a-b)< eps, time* exp(time* (a+b)/2),
                               (exp(a*time) - exp(b*time))/(a-b))})
            for (i in 1:npos) {
                G <- dlist$Ainv %*% dR[,,i] %*% dlist$A
                V <- G*vtemp
                dmat[,,i] <- dlist$A %*% V %*% dlist$Ainv
            }
        }

        if (!all) dlist <- list(P =dlist$P)  # only keep this portion
        # undo the reordering, if needed
        if (setup[1] >0) {
            indx <- order(setup)
            dlist$P <- dlist$P[indx, indx]
            if (deriv) dlist$deriv <- dmat[indx,indx]
        } else {
            if (deriv) dlist$deriv <- dmat
        }
        if (length(dlist)==1) dlist$P else dlist
    } # end of K&L decomp
    else {
        # use the Pade method
        if (deriv) {
            temp <- pade(R*time, dR*time)
            if (!all) temp$nterm <- NULL
            temp
        }
        else pade(R*time)$P
    }    
}
