# The hmm function has done all of the formula and parameter parsing, and
#  calls this function to do the work. 
# id   : subject id of 1,1,1,.. 2,2,2 etc. Data sorted by time within subject
# ytime: vector of times, only used to create the length of time between rows
# ystate: integer, 0= censor, 1,2,.. match the states in qmatrix
# X     : matrix of covariates for the linear predictors
# iprob : initial probability for each id, matrix with one column per state
# B      : matrix of coefficients, one column per linear predictor
# cmap  : integer matrix that maps covariates, parameters, linear predictors
#           one row per covariate, one column per linear predictor
# nlp   : length 3 vector: number of linear predictors for the transitions, the
#          markers, and the initial state.  The second and third can be 0
# ymarker: data frame containing the marker variables
# rlist  : the emission distributions, one per marker
# qmatrix: nstate by nstate matrix, 0= not a valid transition, >0 = valid
# mc.cores: number of cores for mclapply
# control:  see icmsh.control()
# mfun   : the function to use for maximization
# mfattr : list describing characteristics of mfun
# mpar   : optional list of parameters for mfun
# iter   : maximum number of iterations
# iexact : which states are exact, if any, i.e., the time of entry is known
# conmat : constraint matrix (can be null)
# penmat : penalty matrix (can be null)
msh.fit <- function(id, ytime, ystate, X, iprob, B,
                      cmap, nlp, ymarker, rlist, qmatrix,
                      mc.cores, control, mfun, mfattr, mpar, iter,
                      iexact, conmat, penmat) {

    qmatrix[qmatrix>0] <- 1:nlp[1]  # which nlp for each transition
    nmarker <- length(ymarker)
    nstate <- nrow(qmatrix)
    dtime  <- diff(ytime)  #the time interval to the next visit
    dimnames(X) <- NULL    # make X %*% B a bit faster

    # a 0 row in qmap = an absorbing state (you never leave)
    absorb <- (rowSums(qmatrix>0) ==0)
    # an exact state that is absorbing: any obs is known to not be in an
    #  'exactabsorb' state at a censoring time (there is no way to get there)
    exactabsorb <- iexact[absorb[iexact]]

    # the otype variable: 1= interval censored = is in a known state at
    #  this time point; 2= exact = known state and we know exactly when it
    #  was entered (e.g. death), 3= one or more markers, 0 = none of these
    # The likelihood for an exact outcome depends on covariate values before
    #   that time point, we have to treat 'exact' on obs 1 of a subject as
    #   an ordinary "we know their state at this time", i.e. a 1
    # What to do with an obs where markers are measured, but the state is
    #  known, i.e., ystate>0?  I think it is data dependent, so we warn the
    #  user. The current hmm1/hmm2 code will ignore the markers.
    temp1 <- (ystate %in% iexact) & duplicated(id)
    if (nmarker >0 ) {
        temp3 <- rowSums(sapply(ymarker, function(x) !is.na(x))) >0
         # this, it turns out, complains too much
#        if (any(temp3 & ystate>0))
#            warning("marker variables present for an obs with known state")
    } else temp3 <- 0
    otype <- ifelse(temp1, 2L,
                    ifelse(ystate>0, 1L, 3L*temp3))

    # nlp has the number of linear predictors for the rates, markers, and init
    #  state, each linear predictor is a column of cmap.
    # Much of the code is built around a central idea, which is that the 
    #  underlying computations of transition probability, markers, and intial 
    #  state are all in terms of linear predictors, e.g., the matrix exponential
    #  routines are called with a set of linear predictors eta, and they return 
    #  derivatives wrt eta. 
    # These become derivatives wrt the parmameters via matrices created by
    #  the eta.beta1, eta.beta2 and eta.beta3 functions; for the transition 
    #  matrix linear predictors, response function lp, and initial state lp
    #  respectively.  (These have disjoint portions of the parameter vector).
    # Because hmm1/hmm2 can be called lots of times (once per subject per
    #  iteration) we set up some indices to make them faster/simpler.
    # nparm is the number of iterated parameters for each of the three
    # e1, e2, e3 are the columns of eta for each
    nparm <- rep(0L,3)
    if (nlp[1] > 0) { #should always be true
        ctemp <- cmap[, 1:nlp[1], drop=FALSE]
        eta.beta1 <- derivfun(ctemp)   # chain rule for eta to param
        nparm[1] <- length(unique(ctemp[ctemp>0]))
        e1 <- 1:nlp[1]  # the columns of eta for transition matrix
    }
    if (nlp[2] >0) {
        e2 <- nlp[1] + 1:nlp[2] # cols of cmap for the markers
        ctemp <- cmap[, e2, drop=FALSE]
        eta.beta2 <- derivfun(ctemp)
        nparm[2] <- length(unique(ctemp[ctemp>0]))
    }
    if (nlp[3] >0) { # cols of cmap and eta for initial probability
        e3 <- nlp[1] +nlp[2] + 1:nlp[3]
        ctemp <- cmap[, e3, drop=FALSE]
        eta.beta3 <- derivfun(ctemp)
        nparm[3] <- length(unique(ctemp[ctemp>0]))
    }
 
    # Set up parallel, Windows can't fork, others can
    if (mc.cores >1 && control$makecluster) {
        fork <- FALSE
        hmm_cluster <- makeCluster(mc.cores) #start up parallel
        }
    else fork <- TRUE

    # Set up copies of the hmm1 and hmm2 functions to have the scope
    #  of this function. See "scope" in the code vignette for details.
    # We do the same with the hmmloglik, hmmgrad, and hmmboth routines.
    # Arguments to msh.fit that haven't yet been touched, and are needed 
    #  downstream are 'forced' so that they have been copied to the environment.
    # (A small bit of slight-of-hand is that "msh.fit" is later in the
    #   alphabet, so the "hmm" routines compile first and are available).
    # I don't really need the 'x' suffix below, but it helps me remember what
    #   I am doing.  
    force(iprob); 
    hmm1x <- hmm1; hmm2x <- hmm2
    hmmloglikx <- hmmloglik
    hmmgradx <- hmmgrad; hmmbothx <- hmmboth
    environment(hmm1x) <- environment()
    environment(hmm2x) <- environment()
    environment(hmmloglikx) <- environment()
    environment(hmmgradx)   <- environment()
    environment(hmmbothx)   <- environment()

    param <- B.to.coef(B, cmap)
    if (control$detail) {
        # special case-- return all the detail about a fit: per subject
        # alpha and derivatives
        # no iteration
        dfit <- hmmbothx(param, logfun=hmm2x, detail=TRUE)
        return(c(list(param=param), dfit))
    }
    # get the initial loglik and penalty
    initial.loglik <- hmmloglikx(param, logfun=hmm1x)

    if (length(initial.loglik) ==0) 
        stop("unable to evaluate at the intial parameters")

    if (!is.null(penmat))
        penalty0 <- sum(param *(penmat %*% param))/2
    else penalty0 <- 0

    if (iter ==0) {
        # no need to iterate, create return
        if (is.list(initial.loglik)) {
            # this occurs with certain debug options, which return everthing
            rval <- list(param=param, loglik=initial.loglik$loglik, 
                         penalty=penalty0,
                         fit= initial.loglik)
        } else rval <- list(param=param, loglik= initial.loglik,
                            penalty=penalty0)
        return(rval)
    }
       
    clist <- list() # build the do.call argument list
    clist[[mfattr$par]] <- param
    if (mfattr$deriv) {
        # most common call, iteration using gradients
        if (is.null(mfattr$gfun)) {
            # the maximizer expects gradients as an attribute
            clist[[mfattr$fn]] <- hmmbothx
        } else {
            clist[[mfattr$fn]] <- hmmloglikx
            clist[[mfattr$gfun]] <- hmmgradx
        }
        clist[["logfun"]] <- hmm2x
    } else {
        clist[[mfattr$fn]] <- hmmloglikx
        clist[["logfun"]] <- hmm1x
    }

    # figure out where to put the iter argument into the mfun call
    if (is.character(mfattr$iter)) clist[[mfattr$iter]] <- iter  #simple arg
    else { # add it to an mpar argument (optim uses control$maxit)
        tname <- names(mfattr$iter) # what's it called
        temp <- mpar[tname] 
        if (is.null(temp)) {
            temp <- list() # if null, add it as a list
            temp[[mfattr$iter]] <- iter
            clist[[tname]] <- temp
        } else if (is.null(temp[[mfattr$iter]])) {
            temp[[mfattr$iter]] <- iter
            mpar$tname <- temp
        }
    }
    fit <- do.call(mfun, clist)

    # find the fitted coefs in the output, and the loglik
    # optim uses par and value, hmmscore coef and loglik
    nfit <- names(fit)
    indx <- pmatch(c("coef", "par", "log", "value"), nfit, nomatch=0)

    param <- if (indx[1] >0) fit[[indx[1]]] 
             else if (indx[2]>0) fit[[indx[2]]] 
             else {
                 zz <- seq_along(nfit)[-indx]
                 fit[[zz[1]]]
             }
    loglik <-  if (indx[3] >0) fit[[indx[3]]] 
               else if (indx[4]>0) fit[[indx[4]]] else NULL

    # Compute the penalties
    if (!is.null(penmat)) {  #matrix ones later
        penalty <- sum(param * (penmat %*% param))/2
        pderiv  <-  c(param %*% penmat)
    }
    else penalty <- 0
    rval <- list(param=param, 
                 loglik= c(initial= initial.loglik, final= loglik),
                 penalty=penalty, iter=fit$iter, fit=fit)
    rval
}
