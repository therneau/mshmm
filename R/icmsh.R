# Fit an interval censored mult-state hazard model
icmsh <- function(formula, data, subset, weights,
                id, qmatrix, icoef, intercept,
                scale=TRUE, penalty, constraint,
                mfun=hmmscore, mpar=list(gr="hmmboth"), 
                statedata, death="death", exact=death,
                mc.cores= getOption("mc.cores", 1L), control,
                debug=0, fork= (.Platform$OS.type=="unix")) {
    Call <- match.call()
    time0 <- proc.time()

    # create a call to model.frame() that contains the formula (required)
    #  and any other of the relevant optional arguments
    #  but don't evaluate it just yet
    indx <- match(c("formula", "data", "subset", "weights", "id"),
                  names(Call), nomatch=0)
    if (indx[1] ==0) stop("a formula argument is required")
    if (indx[5] ==0) stop("an id argument is required")
    tform <- Call[c(1,indx)]  # only keep the arguments we wanted
    tform$na.action <- quote(stats::na.pass)  # NA done by hand, later
    tform[[1L]] <- quote(stats::model.frame)  # change the function called

    # The transitions matrix qmatrix should be square, non-negative, with
    #  the states as dimnames
    if (missing(qmatrix)) stop("the qmatrix argument is required")
    if (!is.matrix(qmatrix) || (nrow(qmatrix) != ncol(qmatrix))) 
        stop("qmatrix must be a square matrix")
    nstate <- nrow(qmatrix)
    temp <- dimnames(qmatrix)
    if (length(temp[[1]]) ==0) 
      stop("the dimnames of qmatrix must be the state names")
    else{
        statenames <- temp[[1]]
        if (length(temp[[2]]) > 0 && any(temp[[2]] != temp[[1]]))
            stop("row and column names for qmatrix must be identical")
    }
    if (any(diag(qmatrix) != 0))
        stop("the diagonal of qmatrix should be 0")
    if (any(qmatrix < 0)) stop("qmatrix elements must be >=0")
    qmap <- which(qmatrix != 0)   
    ntransitions <- length(qmap)

    if (missing(death)) {
        # don't complain if the default value of "death" is not in the list,
        #  death might not be a state for this model
        if (is.na(match(death, statenames))) death <- NULL
    } else {
        ideath <- match(death, statenames)
        if (any(is.na(ideath)))
            stop("death argument contains a state not in qmatrix")
    }
    if (missing(exact)) {
        # don't complain if the default value of "death" is not in the list,
        #  death might not be a state for this model
        if (is.na(match(exact, statenames))) exact <- NULL
    } else {
        iexact <- match(exact, statenames)
        if (any(is.na(iexact)))
            stop("exact argument contains a state not in qmatrix")
    }
    iexact <- unique(c(ideath, iexact))  # death is always an exact state

    # Is there state data?
    if (!missing(statedata)) { # check that it is okay
        if (!inherits(statedata, "data.frame"))
            stop("statedata must be a data frame")
        if (is.null(statedata$state)) 
            stop("statedata data frame must contain a 'state' variable")
        indx <- match(statedata$state, statenames)
        if (any(is.na(indx))) stop("statedata does not contain all the states")
        statedata <- statedata[indx,]  # same row order as the states
        # Statedata might have rows for states that are not in the data set,
        #  for instance if the hmm call had used a subset argument.  Any of
        #  those are eliminated by the above line.
    } else statedata <- data.frame(state=statenames)
     
    # if the formula is a list, do the first level of processing on it,
    #  which is to pick off the list of variable names
    if (is.list(formula)) {
        if (length(formula)==1 && is.formula(formula[[1]])) {
            # a list with only one formula
            multiform <- FALSE
            parse1 <- NULL
            dformula <- formula[[1]]
        } else {  
            multiform <- TRUE
            dformula <- formula[[1]]   # the default formula for transitions   
            parse1 <- parsecovar1(formula[-1])
        }
    } else {
        multiform <- FALSE   # formula is not a list of expressions
        parse1 <- NULL
        dformula <- formula
    }

    # create the master formula, used for model.frame
    # the term.labels + reformulate + environment trio is used in [.terms;
    #  if it's good enough for base R it's good enough for me
    if (!is.null(parse1)) {
        # right hand side of the default formula
        tlab <- attr(delete.response(terms(dform)), "term.labels")
        # rhs of everything else
        tlab <- c(tlab, unlist(lapply(parse1$rhs, function(x){
                        attr(terms.formula(x), "term.labels")})))
        newform <- reformulate(unique(tlab), dformula[[2]])
        environment(newform) <- environment(dformula)
        formula <- newform
    }

    # Evaluate the expanded formula to create the model frame
    tform$formula <- formula
    mf <- eval(tform, parent.frame())
    Terms <- terms(mf)
    if (nrow(mf) ==0) stop("data has 0 rows")
    termnames <- attr(Terms, 'term.labels')

        
    # Grab the response and validate it
    Y <- model.response(mf)
    if (!inherits(Y, "Surv") || attr(Y, "type" != "mright"))     
        stop("response must be of the form Surv(time, state)")
    else {
        ystate <- attr(Y, "states")
        if (length(ystate)!= length(statenames) || any(ystate != statenames))
            stop("response must have the same states as qmatrix, in the same order")
        ytime <- Y[,1]
        ystat <- Y[,2] # 0= censored, 
    }   
    
    # check that the data is sorted by time within subject, all rows for a
    # subject need to be contiguous. 
    Y <- model.response(mf)
    id <- model.extract(mf, "id")
    id2 <- match(id, unique(id))
    if (is.matrix(Y)) index <- order(id, Y[,1])
    else index <- order(id, Y)
    if (any(diff(index) != 1)) stop("data not sorted by time within each id")

    # Get the first pass of the X matrix, and from that information create
    # cmap and tmap
    X <- model.matrix(Terms, mf)
    xassign <- attr(X, "assign")
    # the next two are included in the output object, and used by later
    #  model.frame and model.matrix calls
    xlevels <- .getXlevels(delete.response(Terms), mf)
    contrasts <- attr(X, "contrasts")

    parse2 <- parsecovar2(parse1, statedata, dformula, Terms, qmatrix,
                          statenames, colnames(X), xassign)
    cmap <- parse2$cmap # coefficients for the transitions
    tmap <- parse2$tmap # terms for the transitions

    # Now deal with missings, which we couldn't do before
    # Y and id can't be missing,
    # missing rate variables are filled in with lvcf
    toss <- is.na(Y) | is.na(id)
    if (any(toss)) 
        stop("missing response or id value, which would lead to an invalid timeline for that subject")

    for (i in 1:ncol(X)) {
        if (any(is.na(X[,i])))
            X[,i] <- lvcf(id, X[,i])
    }
    if (any(is.na(X))) # this will only occur on the first row of an id
        stop("initial observation for a subject has missing values")

    weights <- model.weights(mf)
    if (length(weights) >0) stop("weights are not yet supported")

    # Initialize the coefficients.  We do this before scaling X, since the
    #  user's view of coefficients is always on the original scale.
    # First icoef, then any overrides from options
    #  We allow initial values from a prior model that has fewer terms
    # 
    B <- 0*cmap
    if (!missing(icoef)) {
        if (!missing(intercept)) stop("only one of intercept or icoef allowed")
        if (is.matrix(icoef)) {
            # allow for partial matching, so that a smaller model can feed a 
            #  larger
            rmatch <- match(row.names(icoef), row.names(imat))
            cmatch <- match(col.names(icoef), col.names(imat))
            if (any(is.na(rmatch))) 
                stop("icoef has covarates not in the model")
            if (any(is.na(cmatch)))
                stop("icoef has linear predictors not in the model")
            B[rmatch, cmatch] <- icoef
        } else if (is.numeric(icoef)) {
            if (length(icoef) != nparam) stop("wrong length for icoef")
            B <- coef.to.B(icoef, cmap)
        } else stop("icoef must be a numeric vector or matrix")
    }   

    # Import any init() or fixed() from the options
    # not yet done

    # Standardize the X matrix.  In the extrememly rare case that there
    # is a linear predictor that does not involve the intercept, e.g. a user had
    # factor(group)-1, we can't do so.
    if (Xassign[1]!=0 || any(cmap[1,rcol] ==0)) {
        if (!missing(scale) && scale)
            warning("not possible to scale the data")
        scale <- FALSE
    }
    if (scale) {
        rvar <- 2:ncol(X)  # don't rescale the intercept
        Xscale <- c(0, colMeans(X[,rvar]))
        Xscale <- c(1, apply(X[,rvar], 2, sd))
        for (i in rvar) X[,i] <- (X[,i]- Xmean[i])/Xscale[i]
        # we have XB = (X T^{-1}) (T B) where T is a transformation matrix
        # see rescaling in the code vignette for more detail
        btrans <- diag(Xscale)   # T matrix, transforms B
        btrans[1, rvar] <- Xmean[rvar]
        xtrans <- diag(1/Xscale)      # T inverse, transforms X
        xtrans[1, rvar] <- -(Xmean/Xscale)[rvar]
        B <- btrans %*% B  # the coefs were in terms of unscaled X
    }
    
    if (!missing(intercept)) {
        # I expect this to be a common option for a new fit
        #  If intercepts start at a sensible value, the iteration usually
        # succeeds.
        if (!is.null(names(intercept))) {
            # match by name
            indx <- match(names(intercept), colnames(cmap), nomatch=0)
            if (any(indx==0)) 
                stop("intercept name not found: ",
                     paste(names(intercept)[indx==0], collapse=' '))
            B[1, indx] <- intercept
        } else {
            # assume the right order
            if (length(intercept) == ncol(cmap)) 
                B[1,] <- intercept
            else stop("wrong length for intercept")
        }
    }

    # preprocess constraint and penalty
    if (!missing(constraint)) {
        constraint <- hmmconstraint(constraint, Terms)
        constraint <- constraint %*% xtran # the are written for untransformed X
    } else constraint <- NULL
    if (!missing(penalty)) {
        penalty <- hmmconstraint(penalty, Terms)
        penalty <- penalty %*% xtran # the are written for untransformed X
        penmat <- crossprod(penalty)
    } else penmat <- NULL
    
 
    fit <- icmfit(Y, X, id, cmap, qmat, beta, coef, death, exact, control)

    uid <- unique(id)
    nid <- length(uid)
    rindex <- which(qmatrix >0)

    #Give this next variable a long name that won't be found in calling
    #  routines.  It is updated farther down the calling chain.
    hmm_count_of_calls <- c(0, 0)  #total calls to expm, number with tied eigens

    if (mc.cores > 1 && !fork)
        hmm_cluster <- makeCluster(mc.cores) #start up parallel

    # get a component from a list, but don't fail if it isn't there
    grab <- function(x, what) 
        if (what %in% names(x)) x[[what]] else NULL

    hmmloglik <- function(param, ...) {
        beta[cmap>0] <- param[c(cmap)]   # param[cmap] =bad if cmap has 2 columns
        if (mc.cores > 1) {
            if (!fork)
                mcfit <- parLapply(hmm_cluster, 1:nid, hmm1, beta=beta)
            else mcfit <- mclapply(1:nid, hmm1, beta=beta, 
                                   mc.set.seed=FALSE, mc.cores=mc.cores)
        }
        else mcfit <- lapply(1:nid, hmm1, beta=beta)
        
        if (any(sapply(mcfit, is.character))) {
            # failure
            words <- sapply(mcfit, function(x) ifelse(is.character(x), x, ""))
            if (any(words == "underflow")) {
                if (debug > 1) browser()
                # assume a bad guess from a maximizer, return a bad hit
                return(-2 * abs(initial.loglik))
            }
            words <- words[words!=""]
            stop(words[1])
        }

        loglik <- sum(unlist(mcfit))
        tpar <- c(param)  # used for constraints
        if (!is.null(penmat)) loglik <- loglik - sum(tpar * (penmat %*% tpar))/2
       
        # This might be added back at a later date
        #if (!is.null(conmat)) {
        #    temp <-  conmat %*% tpar
        #    loglik <- loglik + sum(log(pmax(temp,0)))# -Inf if there are violations
        #}
        count <- sapply(mcfit, function(x) attr(x, "counts"))
        # this next line reaches back and changes the variable in parent of parent
        #  (the parent called the maximizer, which calls this)
        hmm_count_of_calls <<- hmm_count_of_calls + rowSums(count)
        loglik
    }

    # hand back everything (debug)
    hmmdb <- function(param, ...) {
        beta[cmap>0] <- param[c(cmap)]
        if (mc.cores > 1) {
            if (!fork)
                mcfit <- parLapply(hmm_cluster, 1:nid, hmm2, beta=beta)
            else mcfit <- mclapply(1:nid, hmm2, beta=beta, 
                                   mc.set.seed=FALSE, mc.cores=mc.cores)
        }
        else mcfit <- lapply(1:nid, hmm2, beta=beta)

        alpha <- sapply(mcfit, function(x) grab(x, "alpha"))
        offset <- sapply(mcfit, function(x) grab(x, "offset"))
        if (any(sapply(mcfit, is.character))) {
            # failure
            words <- sapply(mcfit, function(x) ifelse(is.character(x), x, ""))
            if (any(words == "underflow")) {
                # assume a bad guess from a maximizer, return a bad hit
                loglik <- -Inf
                mcfit <- mcfit[which(words=="")]  # toss the bad results and go on
            }
            else {
                words <- words[words!=""]
                stop(words[1])
            }
        }
        else loglik <- sum(log(colSums(alpha)) + offset)
        
        ecount <- sapply(mcfit, function(x) grab(x, "ecount"))
        hmm_count_of_calls <<- hmm_count_of_calls + rowSums(ecount)
        dd <- dim(mcfit[[1]]$deriv)
        rval <- list(alpha = alpha,
             offset = offset,
             loglik = loglik,
             deriv = array(unlist(lapply(mcfit, function(x) grab(x, "deriv"))),
                           dim=c(dd, length(mcfit))),
             ecount=ecount)
        tpar <- param
        if (!is.null(penmat)) rval$penalty <- sum(tpar* (penmat %*% tpar))/2
        #if (!is.null(conmat)) {
        #    temp <- tpar %*% conmat
        #    rval$constraint <- log(pmax(temp,0))
        #}
        rval
    }
    # This function is used by the score based iteration
    hmmboth <- function(param, ...) {
        beta[cmap>0] <- param[c(cmap)]
        if (mc.cores > 1) {
            if (!fork)
                mcfit <- parLapply(hmm_cluster, 1:nid, hmm2, beta=beta)
            else mcfit <- mclapply(1:nid, hmm2, beta=beta, 
                                   mc.set.seed=FALSE, mc.cores=mc.cores)
        }
        else mcfit <- lapply(1:nid, hmm2, beta=beta)

        alpha <- sapply(mcfit, function(x) sum(grab(x, "alpha")))
        offset <- sapply(mcfit, function(x) grab(x, "offset"))
        if (any(sapply(mcfit, is.character))) {
            # failure
            words <- sapply(mcfit, function(x) ifelse(is.character(x), x, ""))
            if (any(words == "underflow")) {
                # assume a bad guess from a maximizer, return a bad hit
                # the routine won't use deriv or S in this case
                return(list(loglik=-Inf, deriv=NULL, S=NULL))
            }
            else {
                words <- words[words!=""]
                stop(words[1])
                }
        }
        else loglik <- sum(log(alpha) + offset)
        
        ecount <- rowSums(sapply(mcfit, function(x) grab(x, "ecount")))
        hmm_count_of_calls <<- hmm_count_of_calls + ecount
        
        d.alpha <- sapply(mcfit, function(x) rowSums(grab(x, "deriv")))
        u <- d.alpha * rep(1/alpha, each=nparm)
        # this will be a matrix with nparm rows, one col per subject
        deriv <- rowSums(u)
        S <- tcrossprod(u)
        S2 <- tcrossprod(u - rowMeans(u))  # a better variance when deriv!=0 ?
        
        # S is an estimate of -Hessian, so add -(second deriv) of penalty
        if (!is.null(penmat)) {
            tpar <- param
            loglik <- loglik - sum(tpar *(penmat %*% tpar))/2
            deriv  <- deriv - c(param %*% penmat)
            S <- S + penmat
            S2 <- S2 + penmat
            }
       # if (!is.null(conmat)) {
       #     temp <-  c(conmat %*% c(1, param))  # vector of values
       #     loglik <- loglik + sum(log(pmax(temp,0)))
       #     first <- colSums(conmat[,-1]/ rep(temp, nparm))
       #     deriv <- deriv + first
       #     S <- S + outer(first, first)
       #     if (debug>1) {cat("in hmmboth\n"); browser()}
       # }
            
        list(loglik = loglik, deriv= deriv, S=S, S2=S2)
    }

    hmmgrad <- function(param, ...) {
        beta[cmap>0] <- param[c(cmap)]
        if (mc.cores > 1) {
            if (!fork)
                mcfit <- parLapply(hmm_cluster, 1:nid, hmm2, beta=beta)
            else mcfit <- mclapply(1:nid, hmm2, beta=beta, 
                                   mc.set.seed=FALSE, mc.cores=mc.cores)
        }
        else mcfit <- lapply(1:nid, hmm2, beta=beta)
        
        if (any(sapply(mcfit, is.character))) {
            # failure
            words <- sapply(mcfit, function(x) ifelse(is.character(x), x, ""))
            if (any(words == "underflow")) {
                # assume a bad guess from a maximizer, return a bad hit
                loglik <- -2 * abs(initial.loglik)
                return(rep(0, length(param)))
            }
            else {
                words <- words[words!=""]
                stop(words[1])
                }
        }

        ecount <- rowSums(sapply(mcfit, function(x) grab(x, "ecount")))
        hmm_count_of_calls <<- hmm_count_of_calls + ecount
        
        alpha <- sapply(mcfit, function(x) sum(grab(x, "alpha")))
        d.alpha <- sapply(mcfit, function(x) rowSums(grab(x, "deriv")))
        # this will be a matrix with npar rows, one col per subject
        u <- d.alpha * rep(1/alpha, each=nparm)

        if (is.null(penmat)) deriv <- rowSums(u)
        else deriv <- rowSums(u) - c(param %*% penmat)
        #if (!is.null(conmat)) {
        #    temp <- conmat %*% c(1, param)
        #    deriv <- deriv + colSums(conmat[,-1]/temp)
        #}
        deriv
    }
    # fill in the parameters
    mpar$par <- param   # the initial parameters
    if (is.null(mpar$fn)) mpar$fn  <- hmmloglik
    else mpar$fn <- get(mpar$fn)
    if (!is.null(mpar$gr) && is.character(mpar$gr)) mpar$gr <- get(mpar$gr)
    if (!is.null(conmat) && is.null(mpar[["constraint"]]))
        mpar$constraint <- conmat

    time1 <- proc.time()

    if (missing(mfun)) {
        # This is a call with no iteration, using the default hmmloglik function
        initial.loglik <- numeric(0)
        fit <- do.call(mfun, mpar)
        initial.loglik <- penalty0 <- NULL
    }

    else {
        # Get an initial loglik
        # If the maximizer fails (NA or infinite) it returns the initial;
        #  set it to a dummy value to detect that
        initial.loglik <- numeric(0)
        initial.loglik <- hmmloglik(param)
        if (length(initial.loglik)==0)
            stop("unable to evalutate the likelihood at the intial parameters")

        # Compute the intial penalty too
        if (!is.null(penmat))
            penalty0 <- sum(param *(penmat %*% param))/2
        else penalty0 <- 0

        fit <- do.call(mfun, mpar)
    }
        
    if (mc.cores > 1 & !fork) stopCluster(hmm_cluster)
    time2 <- proc.time() 
    # find the fitted coefs in the output, and the loglik
    nfit <- names(fit)
    indx <- pmatch(c("coef", "par", "log", "value"), nfit, nomatch=0)

    fcoef <- if (indx[1] >0) fit[[indx[1]]] else
                 if (indx[2]>0) fit[[indx[2]]] else {
                     zz <- seq_along(nfit)[-indx]
                     fit[[zz[1]]]
                 }

    if (is.matrix(fcoef)) {
        if (ncol(fcoef) == nparm) {
            param <- fcoef[nrow(fcoef),]  # the last parameters used
            beta[cmap>0] <- param[cmap]
        } else stop("wrong number of columns in coefficient matrix")
    }
    else {
        param <- fcoef
        beta[cmap>0] <- param[cmap]
    }

    flog <-  if (indx[3] >0) fit[[indx[3]]] else
                 if (indx[4]>0) fit[[indx[4]]] else NULL

    # Compute the penalties
    if (!is.null(penmat)) {  #matrix ones later
        penalty <- sum(param * (penmat %*% param))/2
        pderiv  <-  c(param %*% penmat)
    }
    else penalty <- 0

    # Undo any scaling and centering
    if (!is.null(xtrans)) {
        beta <- xtrans %*% beta
        param[cmap[cmap>0]] <- beta[cmap>0]
    }

    # Add nice dimnames
    bcol <- paste0(row(qmatrix)[qmatrix!=0], ":",
                   col(qmatrix)[qmatrix!=0])
    if (bcount[2]>0) {
        temp <- unique(rcoef[, c("response", "lp")])
        lp <- temp$lp + 1 - temp$lp[match(temp$response, temp$response)]
        bcol <- c(bcol, paste0("R", paste(temp$response, lp, sep='.')))
    }
    if (bcount[3]>0) bcol <- c(bcol, paste0("p", 1:bcount[3]))
    dimnames(beta) <- list(dimnames(X)[[2]], bcol)
    dimnames(cmap) <- dimnames(beta)


    time3 <- proc.time()
    compute.time <- rbind(setup= time1-time0,
                          compute= time2- time1,
                          finish = time3 - time2)

    final <- list(coefficients= param, 
                  loglik = c(intial=initial.loglik, final=flog),
                  penalty= c(initial=penalty0, final=penalty),
                  beta=beta,
                  fit=fit, 
                  evals= hmm_count_of_calls, 
                  time = compute.time,
                  scale =  scale,
                  bcount = bcount,
                  cmap=cmap, rmap=rindex,
                  qmatrix = qmatrix,   # the structure and state names
                  rfun= rfun, pfun=pfun,
                  nstate = nstate,
                  n = c(rows=nrow(mf), subjects=nid),
                  na.action = na.action,
                  removed = removed.obs, 
                  call=Call,  xlevels=xlevels,
                  contrasts= attr(X, "contrasts"),
                  terms = Terms
                  )
    if (!is.null(penmat)) fit$pen.deriv <- pderiv

    # add in model stuff, present if there were any factors
    if (length(xlevels) >0)   final$xlevels <- xlevels
    if (length(contrasts) >0) final$contrasts <- contrasts

    class(final) <- "hmm"
    final
}

