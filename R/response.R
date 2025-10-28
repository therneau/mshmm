# Material for response functions
# initial functions are called with the state map as the first argument, whose
#  length = number of states
#  unique peak for each unique value, no peak for NA values (usually death)
# marker arg = name of marker, used to create labels

# Gaussian
# return a list with 
#   args: will be passed to hmmgauss
#   pname: will be used to make labels for the linear predictors
#   npar: the number of parameters for this dist
#   param: the parameters dealt with by the initialize call
igauss <- function(map, marker, param) {
    umap <- unique(map[!is.na(map)])
    npeak <- length(umap)
    newmap <- ifelse(is.na(map), 0, match(map, umap))

    if (missing(param)) {
        # the usual case
        par <- c(paste0(marker, ".mean", umap, sep=''), 
                     paste0(marker, ".log(std)", umap))
        args <- list(map=newmap, npeak= npeak)
        list(args= args, pname= par, npar=2, param=1:2)
    }
    else {
        check <- match(param, c("mean", "std"), nomatch=0)
        if (any(check==0)) stop("unrecognized param argument: ", param)
        if (length(check)==2) {
            par <- c(paste0(marker, ".mean", umap, sep=''), 
                     paste0(marker, ".log(std)", umap))
            args <- list(map=newmap, npeak= npeak)
            list(args= args, pname=par, npar=2, param=1:2)
        } else if (check==1) {
            par <- paste0(marker, ".mean", umap, sep='')
            args <- list(map= newmap, npeak= npeak)
            list(ars=args, pname=par, npar=1, param=1)
        } else {
            par <- paste0(marker, ".log(std)", umap, sep='')
            args <- list(map= newmap, npeak= npeak)
            list(args= args, pname= par, npar=1, param=2)
        }
    }   
}                                     
    
hmmgauss <- function(y, eta, args, gradient=FALSE) {
    npeak <- args$npeak
    map   <- args$newmap  # for each state, newmap has 0 or the peak number
    nstate <- length(map)

    yprob <- matrix(0., nstate, length(y))
    if (gradient) ygrad <- array(0, dim=c(nstate, length(y), ncol(eta)))

    # eta will have 2*npeak columns
    mcount <- table(map[map!=0])
    for (i in 1:npeak) {
        j <- i+ npeak  # npeak means followed my npeak std
        yprob[map==i,] <- rep(dnorm(y, eta[,i], exp(eta[,j]), log=TRUE),
                              each= mcount[i])
        if (gradient) 
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
