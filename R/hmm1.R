# These two functions compute the loglik distribution for a single subject,
#  hmm1 = without derivatives, hmm2 = with derivatives
# They are called via a chain of hmmfit -> maximizer -> (hmmloglik, hmmgrad, or
#  hmmboth) ->  (hmm1 or hmm2). The original hmm code defined all the routines
#  within the hmm function, which allows all the variables to be found by
#  inheritance (lexical scope), but the final .R file was just too unwieldy.
#
# A local copy of these functions is now made within hmmfit, which accomplishes
#  the same thing.  See 'scope' in the code vignette.
#
# Remember that B is a matrix of coefficients, of the same shape as cmap,
#  while the vector of coefficients is "param": the latter is what the maximizer
#  functions use.  
# In the math vignette, the vector of coefs is a Greek beta, so you will see
#  'beta' used a lot in the comments and description.
#
hmm1 <- function(who, B) {
    if (control$debug >0) cat("in hmm1\n")
    rows <- which(id == who)  # the subjects of interest
    eta <- X[rows,] %*% B
    # starting probability
    if (!is.null(iprob)) alpha <- iprob[who,]
    else if (is.null(p0fixed))
        alpha <- pfun(nstate, eta[1,e3], gradient=FALSE)
    else alpha <- p0fixed
    
    # Now the response functions for this set
    # rlist was set up by cmsh
    rneed <- (otype[rows]==3)
    if (any(rneed)) {
        rfun <- vector("list", nmarker) # results of calling the rfun() fcns
        for (k in 1:nmarker) {
            j <- rlist[[k]]$e2map  + nlp[1] #columns of B for this response
            indx <- rneed & !(is.na(ymarker[rows, k]))
            yy <- ymarker[rows[indx], k]
            if (length(yy) >0) {
                # all the response functions have gradient=FALSE as default
                if (length(j)==0) rfun[[k]] <- rlist[[k]]$rfun(yy) 
                else rfun[[k]] <- rlist[[k]]$rfun(yy, eta[indx, j,drop=FALSE]) 
             }
        }
    }

    # Compute the collection of matrix exponentials for the subject
    #  Perhaps eventually add a check for upper triangular and pass to survexpm
    r2 <- length(rows)  # there should always be at least 2 rows per id
    if (r2 >1) { # failsafe
        Pmat <- array(0, dim=c(nstate, nstate, r2-1))
        if (any(abs(eta[-r2,]) > .Machine$double.max.exp/2)) {
            # such a bad estimate that it may blow up the matrix exp
            if (control$debug >1) {cat("underflow "); browser()}
            return("underflow")
        }
        rmat <- matrix(0, nstate, nstate)
        for (i in (1:r2)[-r2]) {
            rmat[qmatrix>0] <- exp(eta[i, e1])
            diag(rmat) <- diag(rmat) - rowSums(rmat)
            Pmat[,,i] <- survexpm(rmat, dtime[rows[i]], deriv=FALSE) 
        }
        if (control$debug >2 & (any(Pmat < -control$smallpos)) ||
            any(Pmat > (1  + control$smallpos))) 
                {cat ("stop1\n"); browser()}
        Pmat <- pmax(Pmat, 0)  # we sometimes get tiny negative numbers
    } 

    # Now walk through the visits one by one
    offset <- 0   # watch out for underflow
    nc <- integer(nmarker)  # the number of otype=3 so far, per marker
    rmat <- matrix(0., nstate, nstate)

    for (jj in seq_along(rows)) {
        # j is the index in the original data, jj in our subset, e.g., eta, P
        j <- rows[jj] 
        if (otype[j] ==1 ) { # interval censored outcome
            k <- ystate[j] # in this state, at this time
            alpha[-k] <- 0
        } else if (otype[j] == 2) {
            # exact event time (death)
            k <- ystate[j]  # the exact state just entered
            # we need col k of rmat, using the prior covariates, otype is
            #  never 2 on the first row (jj=1)
            rmat[,k] <-0
            rmat[qmatrix>0] <- exp(eta[jj-1, e1]) 
            k <- ystate[j]  # the exact state just entered
            alpha[k] <- sum(alpha*rmat[,k])
            alpha[-k] <- 0  # known to not be in another state
            if (control$debug > 2) cat("A2: j=", j, "alpha=", alpha, "\n")
        }
        else if (otype[j]==3) {  # one or more markers observed
            for (k in 1:nmarker) {
                if (!is.na(ymarker[j,k])) {
                    nc[k] <- nc[k] +1
                    alpha <- alpha* rfun[[k]][,nc[k]]
                }
            }
            if (!all(is.finite(alpha)) || sum(alpha) <=0) {
                if (control$debug > 1) browser()
                return("underflow") 
            }
            if (control$debug>2) cat("B: j=", j, "alpha=", alpha, "\n")
        }
        else { # censored
            if (length(exactabsorb)>0) alpha[exactabsorb] <- 0
            # can't be in one of the exact+absorbing states
        }
        
        if (jj< r2) {  # if not the last obs
            # transition matrix
            alpha <- alpha %*% Pmat[,,jj]  # transition to next time point
            if (control$debug > 2) cat("C: j=", j, "alpha=", alpha, "\n")

            if (!all(is.finite(alpha)) || sum(alpha) <=0) {
                if (control$debug > 1) browser()
                return("underflow")
            }
            if (mean(alpha) < exp(-20)) { # beware underflow
                reset <- min(-20, log(mean(alpha)))
                if (control$debug > 2) cat(" offset=", offset,"reset=", reset, "\n")
                offset <- offset + reset
                alpha <- alpha * exp(-reset)
            }
        }
    }

    loglik <- offset + log(sum(alpha))
    loglik
}

hmm2 <- function(who,  B) {
    if (control$debug >1) cat("in hmm2\n")
    rows <- which(id == who)  # the subjects of interest
    eta <- X[rows,] %*% B
    P.d  <- matrix(0., nparm[1], nstate)
    R.d  <- matrix(0., nstate, nparm[2])   # might be 0 columns
    pi.d <- matrix(0., nstate, nparm[3])   # might be 0 columns

    # starting probability
    if (!is.null(iprob)) alpha <- iprob[who,]
    else if (is.null(p0fixed)) alpha <- pfun(nstate, eta[1,e3], gradient=TRUE)
    else alpha <- p0fixed
    if (nlp[3] >0) {
        pi.d <- pitrans(attr(alpha, 'gradient'), X[rows[1],])
        attr(alpha, 'gradient') <- NULL  # no longer needed
    }
    
    # Execute the response functions, over the censored obs
    rfun <- rgrad <- vector("list", nmarker)  # results from calling rfun()s
    # rfun will be a matrix with nstate rows and nk col, nk = number of
    #  rows which have marker k not missing
    # rgrad is array of (nstate, length(e2map), nk)= (state, lp, row)
    rneed <- (otype[rows]==3)  # observed markers
    # rlist was created in cmsh from marker2$response
    if (any(rneed)) {
        for (k in 1:nmarker) {
            index <- rneed & !is.na(ymarker[rows,k])
            j <- rlist[[k]]$e2map  + nlp[1] #linear predictors for this response
            yy <- ymarker[rows[index], k]
            if (length(yy) > 0) {
                if (length(j) ==0) rfun[[k]] <- rlist[[k]]$rfun(yy) # no lp
                else {
                    temp <- rlist[[k]]$rfun(yy, eta[index, j, drop=FALSE], 
                                            gradient= TRUE)
                    rfun[[k]] <- temp
                    rgrad[[k]] <- attr(temp, "gradient")
                }
            }
        }
    }

    # Walk through the observations one by one
    offset <- 0  # watch out for underflow
    nc <- integer(nmarker)    #number otype==3, so far, per marker
    r2 <- length(rows)
    rmat <- matrix(0., nstate, nstate)
    for (jj in seq_along(rows)) {
        j <- rows[jj]
        if (otype[j] ==1) { # interval censored
            if (control$debug >2) {cat("otype1 "); browser()}
            k <- ystate[j]
            alpha[-k] <- 0
            P.d[, -k] <- 0
        }
        else if (otype[j] == 2) {
            # exact event time (death)
            k <- ystate[j]        # transitioning to state k       
            klp <- (qmatrix[,k])  #linear preds for transitions to state k
            dtemp <- rmat[,k]  #rate at this point
            dtemp[k] <- 0      # we don't want -1*rowsum (r[k,k]) here
            if (nlp[3]) pi.d <- pi.d * rep(dtemp, nlp[3])
            if (nlp[2]) R.d  <- R.d  * rep(dtemp, nlp[2])
            if (control$debug >1) {cat("death "); browser()}
            # Why the j-1 below?  A death density will depend on covariates
            #  measured prior to the death, not measured at the death
            # dtemp above already has this lag, since rmat is from prior iter
            # The D matrix is zeros except for the death column,
            #  for derivatives see the discussion in the code vignette
            # 
            if (nlp[1] >0 ) {
                term1 <- P.d %*% dtemp # first term, for col k
                etemp <- eta.beta1(X[j-1,])[klp,] # eta to beta for these lp
                term2 <- alpha[klp>0] %*% (dtemp[klp>0] * etemp)
                P.d[,-k] <- 0; P.d[,k] <- c(term1) + c(term2)
            }
            alpha[k] <- sum(alpha * dtemp)
            alpha[-k] <- 0
            if (control$debug > 2) {
                cat("\n death: alpha=", format(alpha), "\n")
                # if (nlp[3]) print(pi.d)
                # if (nlp[2]) print(R.d)
                # print(P.d)
            }
            if (control$debug > 2) cat("A2: j=", j, "alpha=", alpha, "\n")
        }
        else if (otype[j]==3) {  # marker(s) were observed
            for (k in 1:nmarker) {
                if (!(is.na(ymarker[j,k]))) {
                    nc[k] <- nc[k] +1
                    temp <- rfun[[k]][,nc[k]]
                    if (nlp[3]) pi.d <- pi.d * temp 
                    if (nlp[1]) P.d  <- P.d * rep(temp, each=nlp[1])
                    if (nlp[2]) R.d  <- R.d * temp
                    if (!is.null(rgrad[[k]])) { #if there are derivatives
                        j <- rlist$e2map[[k]] # which linear predictors
                        # eta.b3(X[j,]) is d\eta/d\beta for all response lp
                        # rgrad[[k]] will nstate rows, length(j) cols
                        temp <- alpha %*% rgrad[[k]][,,nc[k]]
                        R.d <- R.d + temp %*% eta.b3(X[j,])[,rlist$e2map[[k]]]
                    } 
                    if (control$debug>3) browser()
                    alpha <- alpha * temp
                    if (control$debug > 4) {
                        cat("\n response: alpha=", format(alpha), "\n")
                        #if (nlp[3]) print(pi.d)
                        #if (nlp[2]) print(R.d)
                        # print(P.d)
                    }
                }
            }
            if (!all(is.finite(alpha)) || sum(alpha) <=0) {
                if (control$debug>1) browser()
                return("underflow")
            }
            if (control$debug > 2) cat("B: j=", j, "alpha=", alpha, "\n")
        } 
        else {  # censored
            if (length(exactabsorb)) {
                alpha[exactabsorb] <- 0
                P.d[,exactabsorb] <- 0
            }
        }

        if (jj < r2) { # not the last row
            # state matrix transformation P
            rmat[qmatrix>0] <-  exp(eta[jj, e1])
            if (!all(is.finite(rmat))) {
                # a horrible beta can overflow
                if (control$debug > 1) 
                    save(rmat, beta, file=paste0("rfail", who, ".rda")) 
                return("underflow")
            }
            
            diag(rmat) <- diag(rmat) -rowSums(rmat)
            ptemp <- survexpm(rmat, dtime[j], deriv=2)
            if (any(ptemp$P < -control$smallpos | ptemp$P >1)) {
                if (control$debug>1) 
                    save(ptemp, beta, file=paste0("pfail", who, ".rda"))
                return("underflow")
            }
            if (nlp[3]) pi.d <- t(ptemp$P) %*% pi.d 
            if (nlp[2]) R.d <-  t(ptemp$P) %*% R.d 

            if (nlp[1]) {
                # t1 = deriv of each of the nstate^2 elements of P (rows)
                #   wrt each parameter (col)
                t1 <- matrix(ptemp$deriv, nrow=nstate^2) %*% eta.beta1(X[j,])
                t2 <- matrix(alpha %*% matrix(t1, nrow=nstate), 
                             ncol=nstate, byrow=TRUE)
                P.d <-  P.d %*% ptemp$P + t2
            }
            alpha <- drop(alpha %*% ptemp$P)   # ditch the dimensions
            if (control$debug > 4) {
                cat("\n j=", j, "jj=", jj, "alpha=", format(alpha), "\n")
                if (nlp[3]) print(pi.d)
                if (nlp[2]) print(R.d)
            }
            if (control$debug > 2) cat("C: j=", j, "alpha=", alpha, "\n")

            if (!all(is.finite(alpha)) || sum(alpha) <=0) {
                if (control$debug > 1) 
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
    if (control$debug >3) browser()
    list(alpha=alpha, deriv= rbind(P.d, t(R.d), t(pi.d)),
         offset = offset)
}
