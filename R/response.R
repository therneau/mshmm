# Material for response functions

# Gaussian
igauss <- function(map, param= c("mean", "std")) {
    
hmmgauss <- function(y, nstate, eta, gradient=FALSE, statemap, weight=1) {
    hcheck(nstate, statemap)
    ngroup <- max(statemap)
    if (!is.matrix(eta)) stop("eta must be a matrix for hgauss")
    yprob <- matrix(0., nstate, length(y))
    if (gradient) ygrad <- array(0, dim=c(nstate, length(y), ncol(eta)))
    if (ncol(eta) == (ngroup +1)) {
        # common std for all, follows the means columns
        std <- exp(eta[,ngroup+1])
        for (i in which(statemap >0)) {
            j <- statemap[i]
            yprob[i,] <- weight* dnorm(y, eta[,j], sd=std, log=TRUE)
            if (gradient) {
                ygrad[i,,j] <- weight * (y-eta[,j])/std^2 # deriv wrt eta
                ygrad[i,,ngroup+1] <- weight * ((y-eta[,j])^2/std^2 -1)
            }
        }
    }
    else if (ncol(eta) == 2*ngroup) {
        # eta has the means for each, flllowed by std for each
        for (i in which(statemap>0)) {
            j <- statemap[i]
            std <- exp(eta[, ngroup+j])
            yprob[i,] <- weight * dnorm(y, eta[,j], sd=std, log=TRUE)
            if (gradient) {
                ygrad[i,,j] <- weight * (y-eta[,j])/std^2 # deriv wrt eta
                ygrad[i,,ngroup+j] <- weight * (((y-eta[,j])/std)^2 -1)
            }
        }
    }
    else stop("eta and statemap do not match dimensions")
    
    # convert back
    yprob <- exp(yprob)
    if (gradient) attr(yprob, "gradient") <- c(yprob) * ygrad
    yprob
}


# multivariate logit, first category is the reference
mlogit <- function(eta, gradient=FALSE) {
    m <- ncol(eta)
    denom <- 1 + rowSums(exp(eta))
    pi <- cbind(1, exp(eta)) /denom

    if (gradient) {
        if (nrow(eta) ==1) {  # only one subject
            pi2 <- drop(pi)
            #temp <- diag(pi) - outer(pi, pi) # next line is a touch faster
            temp <- diag(pi2) -  rep(pi2, m+1) * rep(pi2, each= m+1)
            dmat <- temp[,-1,drop=FALSE] 
            dim(dmat) <- c(1, dim(dmat))  # make it a "row" per subject
        }
        else {
            dmat <- array(0., dim=c(nrow(eta), m+1, m))
            dmat[,1,] <- -pi[,-1]/denom
            for (j in 1:m) { # for each column of eta
                dmat[,j+1,] <- -pi[,j+1]* pi[,-1]
                dmat[,j+1, j] <- pi[,j+1]*(1-pi[,j+1])
            }
        }
        attr(pi, "gradient") <- dmat
    }
    pi
}
