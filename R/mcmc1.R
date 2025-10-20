# A recasting of the routine from Jonathan.
#   sensible line widths for the code 
#   use the same argument names as other routines
#   retain his groups and tuning
#   I find it easier to have the burnin and main iteration as two 
#  separate chunks.  The downside to this is some replicated lines of code
#  and the danger of fixing only one of the copies, however the primary mcmc is
#  quite short.

mcmc1 <- function(par, fn, iter, burnin, prior.means= par, prior.std,
                  group, pvar, reset=100, trace=0, constraint, bfit=FALSE) {
    # par contains the initial coefficient values, which are also the
    #  prior means if prior.means is missing.
    npar <- length(par)
    if (length(prior.means) != npar) stop("wrong length for prior.means")

    #if (missing(pvar)) stop("pvar argument is required")
    if (length(prior.std) == 1) prior.std <- rep(prior.std, npar)
    if (length(prior.std) != npar) stop("wrong length for prior.std")
    prior.var <- diag(prior.std^2, npar, npar)

    if (missing(constraint)) ncon <- 0
    else {
        ncon <- nrow(constraint)
        if (ncol(constraint) != npar+1) stop("wrong dimension for constraints")
    }

    
    if (length(group) != npar) 
        stop("group must be the same length as the parameters")
    if (any(group!= floor(group) | group < 1))
        stop("group must be a vector of positive integers")
    ngroup <- max(group)
    temp <- factor(group, 1:ngroup)
    gcount <- table(temp)
    if (any(is.na(temp)) || any(gcount==0)) 
        stop("group must be integers from 1 to maxgroup")
    if (any(gcount==1))
        stop("each group must have at least 2 members")
    sigma <- vector("list", ngroup)  # pre-subset the prior variance
    if (missing(pvar)) {
        for (i in 1:ngroup) sigma[[i]] <- diag( sum(group==i) )
        }
    else {
        for (i in 1:ngroup) {
            sigma[[i]] <- pvar[group==i, group==i]
            diag(sigma[[i]]) <- 2*diag(sigma[[i]])  # more independent
        }
    }

    if (bfit) {
        bpar <- matrix(0., burnin, npar) # the burnin results
        blog <- double(burnin)
    }
    mpar <- matrix(0., iter, npar)  # the mcmc results
    tau  <- rep(.5, ngroup)  #used to tune the proposals during burnin
    # starting small is more conservative

    # initialize
    current  <-  par  # the current working estimate and posterior
    posterior <- fn(par) + dmvnorm(par, prior.means, prior.var, log=TRUE)
    # a constraint at 0 won't be exact, i.e., could have -1e-16 
    if (ncon >0) ceps <- 2* min(c(-.Machine$double.eps, 
                                  constraint %*% c(1, current)))

    # burnin iterations
    accept <- rep(0, ngroup)
    a2 <- matrix(0L, burnin, ngroup)  # details on who is kept
    a2[1,] <- 1L 
    if (bfit) {bpar[1,] <- par; blog[1] <- posterior}
    for (i in 2:burnin) {
        if (trace >0 && i%%trace==0) cat("burn-in ", i, '\n')
        for (j in 1:ngroup) {
            proposal <- current
            update <- rmvnorm(1, rep(0, gcount[j]), sigma[[j]] * tau[j])
            proposal[group==j] <- current[group==j] + update

            if (ncon >0) {
                # During burnin, skip over proposals that fail the 
                # constrataints.  (It helps get the variance matrices settled?)
                cfail <- any(constraint %*% c(1, proposal) < ceps) #failed?
                k <- 0
                while(cfail & k < 1000) {
                    update <- rmvnorm(1, rep(0, gcount[j]), sigma[[j]] * tau[j])
                    proposal[group==j] <- current[group==j] + update
                    cfail <- any(constraint %*% c(1, proposal) < ceps)
                   # if (k %%200 == 0) tau[j] <- tau[j]/2
                    k <- k +1
                }
            }   

            if (ncon>0 && any(constraint %*% c(1, proposal) < ceps))
                logpost <- -Inf
            else  logpost <- fn(proposal) + 
                dmvnorm(proposal, prior.means, prior.var, log=TRUE)
			
            logratio <- logpost - posterior
            if (log(runif(1,0,1)) < logratio) { #accept
                current <- proposal
                posterior <- logpost  
                accept[j]  <- accept[j] +1
                a2[i,j] <- 1L
            }
            if (bfit) {bpar[i,] <- current; blog[i] <- posterior}
        }
		
        # tuning --------------------------------------
        if (i%%reset ==0) {
            # update tau and reset the acceptance
            temp <- as.numeric(cut(accept/reset, c(-1, .1, .2, .3, .40, .50, 
                                        .7, .8, .9, 1)))
            # if accept < 10%, shrink tau by 10 fold, shrink a little if
            #  between .3 and .45, leave as is for .4 to .5, 
            # make it larger if too much acceptance.
            #  if too much acceptance
            if (trace <= reset) 
                cat("burn=", i, "accept=", round(accept/reset, 2), "\n")

            tau <- tau * c(.1, .2, .5, .8, 1, 1.25, 2, 3, 4)[temp]
            accept <- rep(0., ngroup)           
            # use the chain so far to reset the variance
            #  average with the prior one to avoid zeros (only an issue
            #  if the first iterations have very low acceptance due to
            #  a bad input variance)
            for (j in 1:ngroup) {
                k <- which(a2[,j] > 0L)
                if (length(k) > 10) {
                    cmat <- cov(bpar[k, group==j])
                    cmat <- ifelse(is.na(cmat), sigma[[j]], cmat)
                    sigma[[j]] <- (sigma[[j]] + cmat)/2
                }
                else if (length(k) ==0) {
                    # complete failure of the initial matrix
                    sigma[[j]] <- diag(diag(sigma[[j]])) * .01
                    if (trace > 0) cat("full reset i=", i, " j=", j, "\n")
                }
            }
        }
    }
   
    # The primary mcmc
    loglik <- double(iter)
    loglik[1] <- posterior   # current posterior
    mpar[1,] <- current
    accept <- rep(1, ngroup)  # starting proposal is "accepted"
    for (i in 2:iter) {
        if (trace>0 && i%%trace ==0) cat('iteration ', i,'\n')
 
        for (j in 1:ngroup) {
            proposal <- current  # working proposal
            update <- rmvnorm(1, rep(0, gcount[j]), sigma[[j]] * tau[j])
            proposal[group==j] <- proposal[group==j] + update
 
            if (ncon<0 && any(constraint %*% c(1,proposal) < ceps)) 
                logpost<- -Inf
            else logpost <- fn(proposal) + dmvnorm(proposal, prior.means,
                                                   prior.var, log=TRUE)
            logratio <- logpost - posterior
            if (log(runif(1,0,1)) < logratio) { #accept
                current <- proposal
                posterior <- logpost  
                accept[j]  <- accept[j] +1
            }
        }
        mpar[i,] <- current
        loglik[i] <- posterior
    }         
        
    rval <- list(loglik = loglik, par= mpar, accept=accept/iter, covlist=sigma,
                 tau=tau)
    if (bfit) {
        rval$bpar <- bpar
        rval$burn.ac <- a2
        rval$blog <- blog
    }
    rval  
}
