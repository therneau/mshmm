# These are simpler versions of hmm1 and hmm2 (found in hmmloglik.R), for the
#  case of simple interval censored: no hidden states, fixed initial probability
# We didn't really require these, e.g., call hmm1 with iprob present
#  and rfun = NULL.  This will be a tiny faster, but mostly it helped with
#  testing, development and just thinking things through to start simpler.
# This routine is called by mclapply(unique.id, ....) and computes the loglik
#  contribution for a single id.
# uid:    the id for which to compute the value
#  id:    vector of integer id values: 1,1,2,2,2,3, ... etc
# ytime:  time values
# ystate: 0= no observed state, 1:k the observed states
# eta   : matrix XB
# otype :  0= censored, 1= observed state, 2= exact & absorbing (death)
#   the id, ytime, ystate, eta, otype args are for all n rows of the data set
# iprob : starting vectors for each subject (m rows = number of subjects)
# beta  : coefficient vector
# absorb:  which states (if any) are absorbing, normally "death", integer vector
# qmat:   the matrix of valid transitions, 0= not, >0 valid
# istate:  starting state for each 
# 
# 
chmm1 <- function(uid, id, ytime, ystate, eta, otype, iprob) {
    rows <- which(id ==uid)  # the subjects of interest
    nstate <- ncol(iprob)
    # starting probability
    alpha <- iprob[uid,]
   
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
    
    # Now walk through time jj <- 2:number of rows
    # Pmat[,,1] is the transtion from time1 to time2
    offset <- 0   # watch out for underflow
    nc <- integer(ny)  # the number of non-censored so far

    for (jj in seq_along(rows)[-1]) {
        # transition matrix
        alpha <- alpha %*% Pmat[,,jj-1]  # transition to next time point
        if (debug > 2) cat("C: j=", jj-1, "alpha=", alpha, "\n")

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

        # observation
        j <- rows[jj] # j is the index in the original data, jj in our subset
        if (otype[j] ==1 ) { # interval censored outcome
            k <- ystat[j] # in this state, at this time
            alpha[-k] <- 0
        } else if (otype[j] == 2) {
            # exact event time (death)
            k <- ystat[jj] # the new state
            temp <- exp(eta[jj,])
            temp[k] <- 0  # temp = transition rates from other states
            alpha[k] <- sum(alpha*temp)
            alpha[-k]<- 0 # not in any other state
            if (debug > 2) cat("A2: j=", j, "alpha=", alpha, "\n")
        }
        else { # censored
        if (length(absorb) >0) {
            alpha[absorb] <- 0
        }

    loglik <- offset + log(sum(alpha))
    attr(loglik, "counts") <- ucount
    loglik
}

# Do the same, but with derivatives
# See hmmlogik.R for the Ptrans function
hmm2 <- function(uid, id, ytime, ystate, eta, beta, otype, iprob, qmatrix) {

    rows <- which(id ==uid)  # the subjects of interest
    nstate <- ncol(iprob)
    # starting probability
    alpha <- iprob[uid,]
    
    # working matrix for the derivatives
    P.d  <- matrix(0., ncol(eta), nstate)
    rmat <- matrix(, nstate, nstate)

    # Walk through the observations one by one
    offset <- 0  # watch out for underflow
    ucount <- c(length(rows)-1, 0)
    nc <- integer(ny)    #number non-censored so far
    r2 <- length(rows)
    rindex <- which(qmatrix >0)

    for (jj in seq_along(rows)[-1]) {
        # compute P and alpha
        rmat[qmatrix] <- exp(eta[jj-1,])
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
    }
    if (debug>3) browser()
    list(alpha=alpha, deriv= rbind(P.d, t(R.d), t(pi.d)), ecount=ecount,
         offset = offset)
}
