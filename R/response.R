# Automatically generated from hmmcode
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
hmmesetup <- function(emat, states, response, statemap=NULL, partial=NULL) {
    nstate <- nrow(emat)  # must have 1 row per state
    # If states is missing believe the row names of emat, if response the
    #  believe the column names.  That was the original version.
    # We now suggest they be included as an extra check.
    if (missing(states)) states <- rownames(emat)
    if (missing(response)) response <- colnames(emat)
    ny  <- length(response) # unique levels of y
    if (nstate != length(states) || any(rownames(emat) != states))
        stop("the pattern matrix rows do not match the states")
    if (any(is.na(match(colnames(emat), response))))
        stop("the pattern matrix columns do no match the response")

    if (missing(statemap)) statemap <- 1:nstate
    else {
        if (!is.integer(statemap)){
            statemap <- match(statemap, states)
            if (any(is.na(statemap)))
                stop("values in statemap do not match the states")
        }
    }    
    if (!identical(sort(unique(statemap[statemap!=0])), 1:nstate))
        stop("values in statemap do not match the states")

    if (any(is.na(emat)) || any(emat != floor(emat)) || any(emat < -1))
        stop("emat must be an array of integers")
    if (all(emat <1)) stop("no linear predictors are specified")
    temp <- table(row(emat), factor(emat, c(-1, 1:max(emat))),
                  useNA="no")
    if (any(temp[,1] != 1))
        stop("each row of emat must have exactly one reference state")
    if (any(colSums(temp)==0)) stop("there are unused linear predictors")
    lpmap <- lapply(split(emat, row(emat)), function(x) x[x>0])
    names(lpmap) <- states
    # lpmap is a list, one per row of emat, showing the linear predictors
    #  for each state
    temp <- sapply(lpmap, function(x) any(duplicated(x)))
    if (any(temp)) stop("cannot currently handle a repeated linear predictor in one row of the error matrix")
 
    # make the mapping array for mlogit, y value x linear predictor x state
    # The slice for each state shows the linear combination of eta values
    #  for that response
    ny <- length(response)
    np <- max(emat) + 1L  # the maximum number of non-zero phats for any state
    map <- array(0L, dim= c(np, ny, nstate))
    dimnames(map) <- list(paste0("p", 1:np), response, states)
    yindx <- match(colnames(emat), response) # someone might not match order
    for (i in 1:nstate) {
        for (j in 1:ncol(emat)) {
            if (emat[i,j] >0) map[emat[i,j]+1L, yindx[j], i] <- 1
            else if (emat[i,j] == -1) {
                map[1, yindx[j],i] <- 1  # only a reference
            }
        }
    }
    
    # add the partials, if any, to the map
    if (!missing(partial)) {
        if (!is.list(partial)) stop("the partial argument must be a list")
        pindex <- match(names(partial), response)
        if (any(is.na(pindex))) 
            stop("names of the elements of partial must match a response")
        if (any(names(partial) %in% colnames(emat)))
            stop("the name of partial element appears in the pattern matrix")
        if (any(is.na(match(unlist(partial), colnames(emat)))))
            stop("all elements of partial much refer to columns in the pattern matrix")
        # I can't think of any examples, but you perhaps could have a case 
        #  were "1 or 2" and "1 or 3" were both possible.  Give warning, since it
        #  is more likely to be a user's typo however.
        if (any(duplicated(unlist(partial))))
            warning("overlapping partial states")
        if (any(sapply(partial, length)) <1)
            stop("a partial argument must refer to at least 2 states")
        for (i in 1:length(pindex)) {
            j <- match(partial[[i]], response)
            # say that j=1:2, pindex[i] = 5
            # then map[5, ,state] = map[1,, state] + map[2,, state], across all
            #  the linear predictors
            map[, pindex[i], ] <- apply(map[,j,], c(1,3), sum)
        }
    }

    result <- list(lpmap=lpmap, emap= map, statemap= statemap) 
    class(result) <- "hmmesetup"
    result
}
hmmemat <- function(y, nstate, eta, gradient=FALSE, statemap, setup, weight) {
    if (!is.matrix(eta)) eta <- matrix(eta, nrow=length(y))
    if (!inherits(setup, "hmmesetup")) 
        stop("setup must be the result of hmmesetup")
    if (missing(statemap)) statemap <- setup$statemap
    if (length(statemap) != nstate) stop("wrong length for statemap")
    if (any(y<1 | y > dim(setup$emap)[1])) stop("y does not match setup")
    
    ny <- length(y)
    rmat <- matrix(0., nrow=nstate, ncol=ny)
    if (gradient) gmat <- array(0., dim=c(nstate, ny, ncol(eta)))
    # each row of rmat is a state. Do computaions
    for (i in 1:length(setup$lpmap)) {
        sindx <- which(statemap ==i)
        map <- setup$emap[,,i]
        if (length(setup$lpmap[[i]]) ==0) {
            # no errors for this state, e.g., death
            indx2 <- which(map[1,] == 1) # there should be only 1
            rmat[sindx, y %in% indx2] <- 1
            # gradient is zero for these cells
        }
        else {
            mtemp <- mlogit(eta[, setup$lpmap[[i]], drop=FALSE], gradient)
            map2 <- map[c(1, setup$lpmap[[i]] + 1L),]
            # Say that there are 5 unique y levels, but state 1 only can
            #  result in 3 of them = 2 linear predictors. Then mtemp will
            #  have 3 colums, reference cell first, and the gradient attribute
            #  2 columns (deriv is 0 for refrence cell). Hence the collapse
            #  of map to map2.
            # mtemp %*% map2 will produce a matrix with 1 col for each possible
            #  y value, we could then use cbind(1:ny, y) as a matrix subscript
            #  to obtain ny predictions for state i
            #    "yhat <- (mtemp %*% map2)[cbind(1:ny, y)]"
            # For each state, the final result is one yhat value per obs 
            #
            # but map2 is a 0/1 matrix with few 1s, can we avoid the matrix
            #  multiplication?  Yes, if map2[j,k] =1, then all obs with y=k will
            #  have mtemp[j,k] added to them.
            # Because the first index of an array varies fastest in R, if we need
            #  to add this increment to multiple states it suffices to 
            #  replicate each addition.
            #
            nr <- length(sindx) 
            i1 <- row(map2)[map2==1]
            i2 <- col(map2)[map2==1]
            i3 <- setup$lpmap[[i]]
            for (k in seq.int(along.with =i1)) {
                j <- which(y== i2[k])
                if (length(j) >0) { # if there are any matching y values
                    rmat[sindx, j] <- rmat[sindx,j] + rep(mtemp[j, i1[k]], each=nr)
                    if (gradient) 
                        gmat[sindx, j,i3 ] <-  gmat[sindx, j, i3] +
                            rep(attr(mtemp, "gradient")[j, i1[k], ], each=nr)
                }
            }
        }
    }

    if (!missing(weight)) {
        if (!is.numeric(weight) || length(weight) !=1 || weight <=0) 
            stop("invalid value for weight")
        rmat <- rmat^weight
        # if rmat is 0 then gmat is 0, avoid 0/0
        r2 <- c(ifelse(rmat==0, 1, rmat))
        if (gradient) gmat <- weight*r2^(weight-1)* gmat
    }
    if (gradient) attr(rmat, "gradient") <- gmat
    rmat
}
hmminit <- function(nstate, eta, gradient=FALSE) {
    if (is.vector(eta)) eta <- matrix(eta, nrow=1)
    extra <- nstate - (1 + ncol(eta))
    if (extra < 0) 
        stop("too many linear predictors for the number of states")
    temp <- mlogit(eta, gradient)

    if (extra ==0) temp[1,]
    else {
        if (gradient) {
            d2 <- matrix(0, nstate, ncol(eta))
            dmat <- attr(temp, "gradient")[1,,]
            d2[1:nrow(dmat),] <- dmat

            p <- c(as.vector(temp), rep(0., extra))
            attr(p, "gradient") <- d2
            p
        }
        else c(temp, rep(0., extra))
    }
}    
# a multinomial response function
hmulti <- function(y, nstate, eta, gradient=FALSE, statemap, weight=1) {
    if (!is.matrix(eta)) eta <- matrix(eta, nrow=length(y))
    temp <- mlogit(eta, gradient)
    if (nrow(temp) != length(y)) 
        stop("nrow(eta) != length(y)")
    stopifnot(is.matrix(statemap), nrow(statemap)== nstate,
              ncol(statemap)== (ncol(eta) +1))
    if (gradient) {
        gmat <- array(0., dim=c(nstate, length(y), ncol(eta)))
        gtemp <- attr(temp, "gradient")
    }
    rmat <- matrix(0., nrow= nstate, ncol =length(y))
        
    # This function is called once per subject: both y and nstate are short
    #  A loop is simpler than fancy indexing
    # (It looks like you could do it all at once, but gtemp[vector, vector, ]
    #  will grab too much.)
    for (i in 1:nstate) {
        indx1 <- match(y, statemap[i,], nomatch=0)
        for (obs in which(indx1>0)) {
            rmat[i, obs] <- temp[obs, indx1[obs]]
            if (gradient) gmat[i,obs,] <- gtemp[obs, indx1[obs],]
        }
    }

    if (weight!= 1) {
        rmat <- rmat^weight
        # if rmat is 0 then gmat is 0, avoid 0/0
        r2 <- c(ifelse(rmat==0, 1, rmat))
        if (gradient) gmat <- weight*r2^(weight-1)* gmat
    }
    if (gradient) attr(rmat, "gradient") <- gmat
    rmat
}
hmmncut <- function(y, nstate, eta, gradient=FALSE, cuts, statemap) {
    xd <- function(x, sd)  # x * dnorm(x)
        ifelse(is.finite(x), x * dnorm(x, 0, sd), 0)
        
    # This only accepts a single parameter
    if (is.matrix(eta) && ncol(eta) != 1) 
        stop("only a single linear predictor is allowed")
    estd <- exp(as.vector(eta))
    if (length(statemap) != nstate) stop("wrong length for statemap")
    if (!is.matrix(cuts) || ncol(cuts) !=2) 
        stop("cuts must be a 2 column matrix")
    if (any(statemap!= floor(statemap) | statemap <0))
        stop ("statemap must be integers >=0")
    if (any(statemap > nrow(cuts)))
        stop("statemap and cuts do not agree")
    yprob <- matrix(0., nstate, length(y))
    if (gradient) ygrad <- array(0., dim=c(nstate, length(y), 1))

    for (i in which(statemap>0)) {
        j <- statemap[i]
        yprob[i,] <- pnorm(cuts[j,2] -y, 0, estd) - pnorm(cuts[j,1] -y, 0, estd)
        if (gradient) 
            ygrad[i,,1] <- xd(cuts[j,1]-y, estd) - xd(cuts[j,2]-y, estd)
    }

    if (gradient) attr(yprob, "gradient") <- ygrad
    yprob
}
        
hmmlcut <- function(y, nstate, eta, gradient=FALSE, cuts, statemap) {
    # The canonical form of the distribution has variance pi^2/3
    #  We want eta=0 to be a variance of 1, so the variance parameter
    #  is rescaled.  
    
    xd <- function(x, sd)  # x * dnorm(x)
        ifelse(is.finite(x), x * dlogis(x, 0, sd), 0)
        
    # This only accepts a single parameter
    if (is.matrix(eta) && ncol(eta) != 1) 
        stop("only a single linear predictor is allowed")
    estd <- exp(as.vector(eta))* sqrt(3)/pi
    if (length(statemap) != nstate) stop("wrong length for statemap")
    if (!is.matrix(cuts) || ncol(cuts) !=2) 
        stop("cuts must be a 2 column matrix")
    if (any(statemap!= floor(statemap) | statemap <0))
        stop ("statemap must be integers >=0")
    if (any(statemap > nrow(cuts)))
        stop("statemap and cuts do not agree")
    yprob <- matrix(0., nstate, length(y))
    if (gradient) ygrad <- array(0., dim=c(nstate, length(y), 1))

    for (i in which(statemap >0)) {
        j <- statemap[i]
        yprob[i,] <- plogis(cuts[j,2] -y, 0, estd) - 
                     plogis(cuts[j,1] -y, 0, estd)
        if (gradient) 
            ygrad[i,,1] <- xd(cuts[j,1]-y, estd) - xd(cuts[j,2]-y, estd)
    }

    if (gradient) attr(yprob, "gradient") <- ygrad
    yprob
}
hcheck <- function(nstate, statemap) {
    if (length(statemap) != nstate) stop("wrong length for statemap")
    if (any(statemap!= floor(statemap) | statemap <0))
        stop ("statemap must be integers >=0")
    j <- sort(unique(statemap[statemap!=0]))
    if (any(j != seq(along=j)))
        stop("statemap values must be a group number, with no missing groups")
}

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

hmmbeta <- function(y, nstate, eta, gradient=FALSE, statemap, weight =1) {
    hcheck(nstate, statemap)
    ngroup <- max(statemap)
    if (!is.matrix(eta)) stop("eta must be a matrix for hbeta")
    yprob <- matrix(0., nstate, length(y))
    if (gradient) ygrad <- array(0, dim=c(nstate, length(y), ncol(eta)))
    if (ncol(eta) == 2*ngroup) {
        # exp(eta) contains shape1 and shape2
        for (i in which(statemap>0)) {
            j <- 2*statemap[i]-1  # columns of eta
            k <- j+1
            shape1 <- exp(eta[,j]); shape2 <- exp(eta[,k])
            yprob[i,] <- weight * dbeta(y, shape1, shape2, log=TRUE)
            if (gradient) {
                ygrad[i,,j] <- weight*(log(y) + digamma(shape1+shape2) - 
                                  digamma(shape1)) * shape1
                ygrad[i,,k] <- weight*(log(1-y) + digamma(shape1+shape2) - 
                                  digamma(shape2)) * shape2
            }
        }
    }
    else stop("eta and statemap do not match dimensions")

    # convert from log-density to density
    yprob <- exp(yprob)
    if (gradient) attr(yprob, "gradient") <- c(yprob) * ygrad
    yprob
}

hmmexp <- function(y, nstate, eta, gradient=FALSE, statemap, weight =1) {
    hcheck(nstate, statemap)
    ngroup <- max(statemap)
    if (!is.matrix(eta)) stop("eta must be a matrix for hbeta")
    yprob <- matrix(0., nstate, length(y))
    if (gradient) ygrad <- array(0, dim=c(nstate, length(y), ncol(eta)))
    if (ncol(eta) == ngroup) {
        # exp(eta) contains lambda
        for (i in which(statemap>0)) {
            j <- statemap[i]
            lambda <- exp(eta[,j])
            yprob[i,] <- weight * dexp(y, lambda, log=TRUE)
            if (gradient) {
                ygrad[i,,j] <- weight * (1/lambda -y)
             }
        }
    }
    else stop("eta and statemap do not match dimensions")
    
    # convert from log-density to density
    yprob <- exp(yprob)
    if (gradient) attr(yprob, "gradient") <- c(yprob) * ygrad
   yprob
}

hmmgamma <- function(y, nstate, eta, gradient=FALSE, statemap, weight=1) {
    hcheck(nstate, statemap)
    ngroup <- max(statemap)
    if (!is.matrix(eta)) stop("eta must be a matrix for hgamma")
    yprob <- matrix(0., nstate, length(y))
    if (gradient) ygrad <- array(0, dim=c(nstate, length(y), ncol(eta)))
    if (ncol(eta) == 2*ngroup) {
        # exp(eta) contains shape1 and scale
        for (i in which(statemap>0)) {
            j <- 2*statemap[i]-1  # columns of eta
            k <- j+1
            shape <- exp(eta[,j]); rate <- exp(-eta[,k])
            yprob[i,] <- weight* dgamma(y, shape=shape, rate=rate, log=TRUE)
            if (gradient) {
                f<- yprob[i,]
                ygrad[i,,j] <- weight*(log(rate) + log(y) - digamma(shape))* shape
                ygrad[i,,k] <- weight *(y*rate - shape)
            }
        }
    }
    else stop("eta and statemap do not match dimensions")
    
    # convert from log-density to density
    yprob <- exp(yprob)
    if (gradient) attr(yprob, "gradient") <- c(yprob) * ygrad
    yprob
}
