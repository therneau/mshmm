# A very simple MCMC fitter for the HMM code
mcmc0 <- function(par, fn, iter=0, burnin=1000, 
                 prior.std = 4, prop.std =.1, reset=300, trace=0) {

    npar <- length(par)
    nmean <- par  # par contains the prior means
    
    log0 <- fn(par) # initial log likelihood
    loglik <- double(iter)
    accept <- matrix(0L, burnin, npar)
    newpar <- matrix(0L, iter, npar)
    
    prior.std <- rep(prior.std, length=npar)
    prop.std <-  rep(prop.std, length=npar)
    if (trace>0) pfile <- file("mcmctrace", "w")

    # Do the burnin, which will be "burnin" full loops of npar adjustments each
    tpar <- par
    best <- log0 + sum(dnorm(tpar, nmean, prior.std, log=TRUE))
    if (trace >0) cat("burn in\n", file=pfile)
    for (i in 1:burnin) {
        update <- rnorm(npar, tpar, prop.std)  # the update for each term
        for (j in 1:npar) {
            temp <- tpar
            temp[j] <- update[j]
            newlk <- fn(temp) + sum(dnorm(temp, nmean, prior.std, log=TRUE))
            if (newlk > best || exp(newlk - best) > runif(1)) {
                accept[i,j] <- 1
                best <- newlk
                tpar[j] <- update[j]
            }
        }
        if (trace>0 && i%%trace==0)
            cat("i=", i, " loglik=", format(best), "\n", file=pfile)
        if (reset >0 && i%%reset==0) {
            # tune up the proposal std to get a better acceptance rate
            #   model: acceptance rate approx = exp(-k*sd)
            #   solve for k* newsd = 0.9, which gives 40% acceptance
            # If the current sd is just right we will tend to oscillate:
            #   a phat < .4 decreases the sd too much, etc.  So damp it
            #   a little.  This also stops explosion if phat=0 or 1
            #   
            phat <- (.4 + 2* colMeans(accept[seq(i-reset, i),]))/3
            k <- -log(phat)/prop.std
            prop.std <- .92/k 
            if (trace>0) {
                temp <- format(rbind(phat, k, prop.std), digits=2)
                cat(" phat    ", temp[1,], "\n", 
                    "k       ", temp[2,], "\n",
                    "prop.std", temp[3,], "\n", file=pfile)
            }
        }
    }       

    # Now do the iterations. 
    # It is no longer valid to update the proposal distribution
    accept <- matrix(0L, iter, npar)
    for (i in 1:iter) {
        update <- rnorm(npar, tpar, prop.std)  # the update for each term
        utemp  <- runif(npar, 0, 1)
        for (j in 1:npar) {
            temp <- tpar
            temp[j] <- update[j]
            newlk <- fn(temp) + sum(dnorm(temp, nmean, prior.std, log=TRUE))
            if (newlk > best || exp(newlk - best) > utemp[j]) {
                accept[i,j] <- 1
                best <- newlk
                tpar[j] <- update[j]
            } 
        }
        loglik[i] <- best
        newpar[i,] <- tpar
        if (trace>0 && (i%%trace)==0)
            cat("i=", i, " loglik=", format(best), "\n", file=pfile)
    }       
    
    if (trace>0) close(pfile)
    list(loglik = loglik, accept=accept, par=newpar)
}   
    
