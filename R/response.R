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
igaussian <- function(map, marker, param) {
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

# say that there were 5 states and 2 peaks, the first applies to states 1 and 2,
#  the second to 3 and 4, which state= death has no parameters for this 
#  biomarker (it's not measured on dead people). Then npeak =2, map will be
#  1,1,2,2,0, and eta will have 4 columns for mean, mean, log(std), log(std).
#  yprob will have 5 rows, one per state; 1 and 2 the same, 3 and 4 the same,
#  containing the probability for each observation in the columns.
#  Row 5 will be 0.
# The gradient will have dimensions of state, obs, eta, and returns the 
#  derivative of the density with respect to each eta. Again, in this case
#  margins (1:2,,) and (3:4,,) are duplicates, and ygrad[5,,] =0. The parent
#  hmm routine uses the chain rule (d f/ d eta)(d eta/ d beta) to get 
#  derivatives for the parameters beta.
# The code does computations on the log density scale until the last step.
#  (The parent routine needs density).
hmmgaussian <- function(y, eta, args, gradient=FALSE) {
    npeak <- args$npeak
    map   <- args$newmap  # for each state, newmap has 0 or the peak number
    nstate <- length(map)

    yprob <- matrix(0., nstate, length(y))
    if (gradient) ygrad <- array(0, dim=c(nstate, length(y), ncol(eta)))

    # eta will have 2*npeak columns
    mcount <- table(map[map!=0])
    for (i in 1:npeak) {
        j <- i+ npeak  # npeak means followed my npeak log(std values)
        std <- exp(eta[,j])
        yprob[map==i,] <- rep(dnorm(y, eta[,i], std, log=TRUE),
                              each= mcount[i])
        if (gradient) {
           gradient[map==i,,i] = (y- eta[,i])/std^2
           gradient[map==i,,j] <- (((y-eta[,i])/std)^2 -1)
        }
    }
    
    # convert from log to density
    yprob[map>0,] <- exp(yprob[map>0,])
    if (gradient) attr(yprob, "gradient") <- c(yprob) * ygrad
    yprob
}

# Logisitic, a little fatter tails
ilogistic <- igaussian
hmmlogisitc <-  function(y, eta, args, gradient=FALSE) {
    npeak <- args$npeak
    map   <- args$newmap  # for each state, newmap has 0 or the peak number
    nstate <- length(map)

    yprob <- matrix(0., nstate, length(y))
    if (gradient) ygrad <- array(0, dim=c(nstate, length(y), ncol(eta)))

    # eta will have 2*npeak columns
    mcount <- table(map[map!=0])
    for (i in 1:npeak) {
        j <- i+ npeak  # npeak means followed my npeak log(std values)
        std <- exp(eta[,j])* sqrt(3)/pi # the R logist has "scale" not "std"
        yprob[map==i,] <- rep(dlogis(y, eta[,i], std, log=TRUE),
                              each= mcount[i])
        if (gradient) {
           gradient[map==i,,i] = (y- eta[,i])/std^2
           gradient[map==i,,j] <- (((y-eta[,i])/std)^2 -1)
        }
    }
    
    # convert from log to density
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
