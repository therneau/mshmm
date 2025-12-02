# These two functions compute the loglik distribution for a single subject,
#  hmm1 = without derivatives, hmm2 = with derivative
# They are called via a chain of hmmfit -> maximizer -> (hmmloglik, hmmgrad, or
#  hmmboth) ->  (hmm1 or hmm2).  To avoid having to pass a stack of arguments
#  down the chain (X matrix, id, etc) we set their environment to hmmfit.
# If we included these at the bottom of the hmmfit.R file they would inherit
#  that automatically, but that makes the source file way too large.
# Note that this file name "hmmloglik.R" is after "hmmfit.R": they are run in
#  the right order by R CMD build.
# 

hmm1 <- function(who, beta) {
    rows <- which(id ==uid[who])  # the subjects of interest
    eta <- X[rows,] %*% beta
    # starting probability
    if (!is.null(iprob)) alpha <- iprob[who,]
    else if (is.null(p0fixed))
        alpha <- pfun(nstate, eta[1,b3], gradient=FALSE)
    else alpha <- p0fixed
    
    # Now the response functions for this set
    rneed <- (otype[rows]==3)
    rlist <- vector("list", ny)
    #    cat("in hmm1\n"); browser()
    for (k in 1:ny) {
        j <- b2map[[k]]  #columns of beta for this response
        indx <- rneed & !is.na(yobs[rows, k])
        yy <- yobs[rows[indx], k]
        if (length(yy) >0) {
            if (length(j)==0) 
                rlist[[k]] <-rfun[[k]](yy, nstate, gradient=FALSE)
            else rlist[[k]] <- rfun[[k]](yy, nstate, 
                eta[indx, j,drop=FALSE], gradient=FALSE)
        }
    }
    
    # Compute the collection of matrix exponentials for the subject
    # The upper routine sends back the array of results as a vector
    #  along with the number of times there were tied eigenvalues
    #  we'll send the ties back as an attribute
    # We don't need the last row for each subject.
    # The call to upper uses the Ward approx (nterm=0) rather than
    #  the Higham09.  The former seems to better match my pade routine.
    r2 <- length(rows)  # there should always be at least 2 rows per id
    if (length(rows) > 1) {  # but add a failsafe
        if (any(abs(eta[-r2,]) > .Machine$double.max.exp/2)) {
            # such a bad estimate that it may blow up the matrix exp
            if (debug >1) browser()
            return("underflow")
        }
        myexp <- .Call("upper", nstate, eta[-r2,,drop=FALSE], 
                       ytime[rows[-r2]], rindex, 1e-7, 0)
        ucount <- c(length(rows)-1, myexp$ties)
        Pmat <- array(myexp$P, dim=c(nstate, nstate, length(rows)-1))
        
        if (debug >2 & any(Pmat < -control$smallpos)) {
            cat ("stop1\n"); browser()}
        if (any(Pmat > (1+control$smallpos) | 
                Pmat < -control$smallpos)) return("underflow")
        Pmat <- pmax(Pmat, 0)  # we sometimes get tiny negative numbers
    } 
    
    # Now walk through the visits one by one
    offset <- 0   # watch out for underflow
    nc <- integer(ny)  # the number of non-censored & non-missing so far
    rmat <- matrix(0., nstate, nstate)

    for (jj in seq_along(rows)) {
        j <- rows[jj] # j is the index in the original data, jj in our subset               
        if (otype[j] ==1 ) { # interval censored outcome
            k <- ystat[j] # in this state, at this time
            alpha[-k] <- 0
        } else if (otype[j] == 2) {
            # exact event time (death)
            rmat[rindex] <- exp(eta[jj-1, b1]) #covariate just before this point
            k <- ystat[j]  # the exact state just entered
            alpha[k] <- sum(alpha*rmat[,k])
            alpha[-k] <- 0  # known to not be in another state
            alpha <- alpha * dtemp
            if (debug > 2) cat("A2: j=", j, "alpha=", alpha, "\n")
        }
        else if (otype[j]==3) {  # one or more markers observed
            temp <- rep(1, nstate)
            for (k in 1:ny) {
                if (!is.na(yobs[j,k])) {
                    nc[k] <- nc[k] +1
                    alpha <- alpha* rlist[[k]][,nc[k]]
                    temp <- temp * rlist[[k]][,nc[k]]
                }
            }
            if (!all(is.finite(alpha)) || sum(alpha) <=0) {
                if (debug > 1) browser()
                return("underflow") 
            }
            if (debug>2) cat("B: j=", j, "alpha=", alpha, "\n")
        }

        if (any(exactabsorb) & ystat[j]==0) {
            # can't be in one of the exact+absorbing states
            alpha[exactabsorb] <- 0
        }
        
        if (!is.null(entrytime) && entrytime[j] ==1) {
            # entry to the study
            temp <- sum(alpha*entry)
            if (temp <= 0) 
                return(paste("subject", uid[who],
                             "enters in an impossible state"))
            alpha <- alpha*entry /temp
            if (debug > 2) cat("A3: j=", j, "alpha=", alpha, "\n")
        }

        if (jj< r2) {  # if not the last obs
            # transition matrix
            alpha <- alpha %*% Pmat[,,jj]  # transition to next time point
            if (debug > 2) cat("C: j=", j, "alpha=", alpha, "\n")

            if (!all(is.finite(alpha)) || sum(alpha) <=0) {
                if (debug > 1) browser()
                return("underflow")
            }
            if (mean(alpha) < exp(-20)) { # beware underflow
                reset <- min(-20, log(mean(alpha)))
                if (debug > 2) cat(" offset=", offset,"reset=", reset, "\n")
                offset <- offset + reset
                alpha <- alpha * exp(-reset)
            }
        }
    }

    loglik <- offset + log(sum(alpha))
    attr(loglik, "counts") <- ucount
    loglik
}
environment(hmm1) <- environment(hmmfit)

# Helper function for the derivatives in hmm2
makeindex <- function(cmap, all=cmap) {
    nonzero <- (cmap > 0)
    parms <- sort(unique(all[all>0]))  # the parameter numbers for this group
    p <- ncol(cmap)    # number of linear predictors
    k <- match(cmap[nonzero], parms)  #parameter number
    nparm <- length(parms)
    rr <- row(cmap)[nonzero]  # which X to use
    cc <- col(cmap)[nonzero]  #which eta this is
    list(xindex= rr, tindex= cc + (k-1)*p, dim=c(p, nparm))
}

psetup <- function(rmat, rindex, nstate) {
    n.eta <- length(rindex)
    out <- array(0., c(nstate, nstate, n.eta))
    temp <- matrix(0., nstate, nstate)
    rr <- row(temp)[rindex]
    for (i in 1:n.eta) {
        temp2 <- temp
        exp.eta <- rmat[rindex[i]]  # elements of rmat are exp(eta)
        temp2[rindex[i]] <-  exp.eta
        temp2[rr[i], rr[i]] <- -exp.eta
        out[,,i] <- temp2
    }
    out
}

Rtrans <- vector("list", ny)  #one element per response function
if (bcount[2]) { #if there are response parameters
    tfun <- function(dmat, x, map) {
        tmat <- matrix(0., map$dim[1], map$dim[2])
        tmat[map$tindex] <- x[map$xindex]
        dmat %*% tmat
    }
    for (i in 1:ny) {
        if (length(b2map[[i]]) >0 && any(cmap[, b2map[[i]]] > 0)) {
            formals(tfun)[[3]] <- makeindex(cmap[,b2map[[i]], drop=FALSE],
                                            cmap[,b2])
            Rtrans[[i]] <- tfun
        }
    }
}
if (bcount[3]) { #if there are initial probability  parameters
    pitrans <- function(dmat, x, map) {
        tmat <- matrix(0., map$dim[1], map$dim[2])
        tmat[map$tindex] <- x[map$xindex]
        dmat %*% tmat
    }
    formals(pitrans)[[3]] <- makeindex(cmap[,b3, drop=FALSE])
}
cmap.b1 <- makeindex(cmap[,b1, drop=FALSE])

Ptrans <- function(alpha, dmat, x, map=cmap.b1) {
    tmat <- matrix(0., map$dim[1], map$dim[2])
    tmat[map$tindex] <- x[map$xindex]
    #treat dmat as though it were a matrix with first dim nstate*nstate
    dim(dmat) <- c(nstate*nstate, map$dim[1])
    dmat2 <- dmat %*% tmat  #transform
    t(rowsum(dmat2 * rep(alpha, nstate*map$dim[2]), rep(1:nstate, each=nstate),
             reorder=FALSE))
}    

if (!missing(exact)) {  
    dtemp <- col(qmatrix)[rindex]
    exactcol  <- which(dtemp== exact)  # list of linear predictors
    exacttrans <- function(R, x, map=cmap.b1) {
        rows <- which(qmatrix[,exact] > 0)  #non-zero elements of column d
        dmat <- matrix(0, nstate, map$dim[1])
        for (i in 1:length(rows)) 
            dmat[rows[i], exactcol[i]] <- R[rows[i], exact]

        tmat <- matrix(0., map$dim[1], map$dim[2])
        tmat[map$tindex] <- x[map$xindex]
        dmat %*% tmat
    }
}


hmm2 <- function(who, beta) {
    rows <- which(id ==uid[who])  # the subjects of interest
    eta <- X[rows,] %*% beta
    P.d  <- matrix(0., parmcount[1], nstate)
    R.d  <- matrix(0., nstate, parmcount[2])
    pi.d <- matrix(0., nstate, parmcount[3])

    # starting probability
    if (!is.null(iprob)) alpha <- iprob[who,]
    else if (is.null(p0fixed)) alpha <- pfun(nstate, eta[1,b3], gradient=TRUE)
    else alpha <- p0fixed
    if (bcount[3]) {
        pi.d <- pitrans(attr(alpha, 'gradient'), X[rows[1],])
        attr(alpha, 'gradient') <- NULL  # no longer needed
    }
    
    # Execute the response functions, over the uncensored obs
    rlist <- rgrad <- vector("list", ny)
    rneed <- (otype[rows]==3)  # observed markers
    for (k in 1:ny) {
        index <- rneed & !is.na(yobs[rows,k])
        j <- b2map[[k]]  #linear predictors for this response
        yy <- yobs[rows[index], k]
        if (length(yy) > 0) {
            temp <- rfun[[k]](yy, nstate, eta[index, j, drop=FALSE], 
                gradient= TRUE)
            rlist[[k]] <- temp
            rgrad[[k]] <- attr(temp, "gradient")
        }
    }

    # Walk through the observations one by one
    P.d  <- matrix(0., pcount[1], nstate)
    offset <- 0  # watch out for underflow
    ecount <- c(length(rows), 0)
    nc <- integer(ny)    #number non-censored so far
    r2 <- length(rows)
    rmat <- matrix(0., nstate, nstate)

    for (jj in seq_along(rows)) {
        j <- rows[jj]
        if (otype[j] ==1) { # interval censored
            k <- ystat[j]
            alpha[-k] <- 0
            P.d[,, -k] <- 0
        }
        else if (otype[j] == 2 & jj> 1) {
            # exact event time (death)
            k <- pstate[j]
            dtemp <- rmat[,k]  #rate at this point
            dtemp[k] <- 0      # this line should be redundant
            if (pcount[3]) pi.d <- pi.d * rep(dtemp, pcount[3])
            if (pcount[2]) R.d  <- R.d  * rep(dtemp, pcount[2])
            # Why the j-1 below?  A death density will depend on covariates
            #  measured prior to the death, not measured at the death
            # dtemp above already has this lag, since rmat is from prior iter
            if (pcount[1]) P.d  <- P.d *  rep(dtemp, each=pcount[1]) +
                               t(alpha * deathtrans(rmat, X[j-1,]))
            alpha <- alpha * dtemp
            if (debug > 2) {
                cat("\n death: alpha=", format(alpha), "\n")
                # if (pcount[3]) print(pi.d)
                # if (pcount[2]) print(R.d)
                # print(P.d)
            }
            if (debug > 2) cat("A2: j=", j, "alpha=", alpha, "\n")
        }
        else if (otype[j]==3) {  # marker(s) was observed
            for (k in 1:ny) {
                if (!is.na(yobs[j,k])) {
                    nc[k] <- nc[k] +1
                    temp <- rlist[[k]][,nc[k]]
                    if (pcount[3]) pi.d <- pi.d * temp 
                    if (pcount[1]) P.d  <- P.d * rep(temp, each=pcount[1])
                    if (pcount[2]) R.d  <- R.d * temp
                    if (!is.null(Rtrans[[k]])) { #if there are derivatives
                        dtemp <- Rtrans[[k]](rgrad[[k]][,nc[k],], X[j,])
                        R.d  <- R.d + alpha * dtemp
                    } 
                    if (debug>3) browser()
                    alpha <- alpha * temp
                    if (debug > 4) {
                        cat("\n response: alpha=", format(alpha), "\n")
                        #if (pcount[3]) print(pi.d)
                        #if (pcount[2]) print(R.d)
                        # print(P.d)
                    }
                }
            }
            if (!all(is.finite(alpha)) || sum(alpha) <=0) {
                if (debug>1) browser()
                return("underflow")
            }
            if (debug > 2) cat("B: j=", j, "alpha=", alpha, "\n")
        }

        if (any(exactabsorb) && ystat[j]==0) {
            alpha[exactabsorb] <- 0
            stop("need to fix derivatives")
        }

        if (!is.null(entrytime) && entrytime[j]==1) {
            # entry to the study
            temp <- sum(alpha*entry)
            if (temp <= 0) return(paste("subject", uid[who],
                                        "enters in an impossible state"))
            if (pcount[3]) pi.d <- (entry/temp)*( pi.d -
                                             alpha %*% (entry %*% pi.d)/temp )
            if (pcount[2]) R.d <- (entry/temp)* (R.d - 
                                                 alpha %*% (entry %*% R.d)/temp)
            # remember that P.d is (nparm, nstate)
            if (pcount[1]) P.d  <- (rep(entry, each=pcount[1])/temp) * 
                               (P.d - outer(c(P.d %*% entry), alpha) /temp)
            alpha <- alpha*entry /temp
            if (debug > 2) {
                cat("\n entry: alpha=", format(alpha), "\n")
                #if (pcount[3]) print(pi.d)
                # if (pcount[2]) print(R.d)
                # print(P.d)
            }
            if (debug > 2) cat("A3: j=", j, "alpha=", alpha, "\n")
        }

        
        if (jj < r2) { # not the last row
            # state matrix transformation P
            rmat[rindex] <- exp(eta[jj,b1])
            if (!all(is.finite(rmat))) {
                # a horrible beta can overflow
                if (debug > 1) 
                    save(rmat, beta, file=paste0("rfail", who, ".rda")) 
                return("underflow")
            }
            
            diag(rmat) <- diag(rmat) -rowSums(rmat)
            tder <- psetup(rmat, rindex, nstate)
            if (any(diff(sort(diag(rmat))) < 1e-6)) {
                ptemp <- pade(rmat *ytime[j], tder*ytime[j])
                ecount[2] <- ecount[2] +1
            }
            else ptemp <- derivative(rmat, ytime[j], tder)
            if (any(ptemp$P < -eps | ptemp$P >1)) {
                if (debug>1) 
                    save(ptemp, beta, file=paste0("pfail", who, ".rda"))
                return("underflow")
            }
            if (pcount[3]) pi.d <- t(ptemp$P) %*% pi.d 
            if (pcount[2]) R.d <-  t(ptemp$P) %*% R.d 
            if (pcount[1]) 
                P.d <-  P.d %*% ptemp$P + Ptrans(alpha, ptemp$dmat, X[j,])
            alpha <- drop(alpha %*% ptemp$P)   # ditch the dimensions
            if (debug > 4) {
                cat("\n j=", j, "jj=", jj, "alpha=", format(alpha), "\n")
                if (pcount[3]) print(pi.d)
                if (pcount[2]) print(R.d)
            }
            if (debug > 2) cat("C: j=", j, "alpha=", alpha, "\n")

            if (!all(is.finite(alpha)) || sum(alpha) <=0) {
                if (debug > 1) 
                    save(alpha, ptemp, beta, file=paste0("afail", who, "rda"))
                return("underflow")
            }
            if (mean(alpha) < exp(-20)) {
                offset <- offset -20
                alpha <- alpha * exp(20)
                pi.d <- pi.d * exp(20)
                R.d  <- R.d  * exp(20)
                P.d  <- P.d  * exp(20)
            }
        }
    }
    if (debug>3) browser()
    list(alpha=alpha, deriv= rbind(P.d, t(R.d), t(pi.d)), ecount=ecount,
         offset = offset)
}

environment(hmm2) <- environment(hmmfit)
