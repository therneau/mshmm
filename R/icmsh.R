# Fit an interval censored mult-state hazard model
icmsh <- function(formula, data, subset, weights,
                id, qmatrix, mfun, mpar= list(), iter=20,
                mc.cores = getOption("mc.cores", 2L),
                icoef, scale=TRUE, penalty, constraint,
                statedata, exact ="death",
                control= icmsh.control(), ...) {
    Call <- match.call()
    time0 <- proc.time()

    ## We want to pass any ... args to cmsp.control, but not pass things
    ##  like "dats=mydata" where someone simply made a typo.  The use of ...
    ##  is simply to allow things like "eps=1e6" with easier typing
    extraArgs <- list(...)
    if (length(extraArgs)) {
        controlargs <- names(formals(icmsp.control)) #legal arg names
        indx <- pmatch(names(extraArgs), controlargs, nomatch=0L)
        if (any(indx==0L))
            stop(gettextf("Argument %s not matched", 
                          names(extraArgs)[indx==0L]), domain = NA)
    }

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

    if (!is.numeric(iter) || length(iter) >1 || iter <=0) 
        stop("iter must be a postive integer")
    else iter <- ceiling(iter)

    if (missing(exact) || length(exact) ==0) exact <- NULL
    else {
        exact <- match(exact, statenames)
        if (any(is.na(exact)))
            stop("exact argument contains a state not in qmatrix")
    } 

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

    weights <- model.weights(mf)
    if (length(weights) >0) warning("weights are not yet supported")
    weights <- rep(1.0, nrow(mf))

   
    # mark out the exact states
    iexact <- match(exact, statenames, nomatch=0)
    if (missing(exact)) {
        iexact <- iexact[iexact >0]  #don't complain if our default "death" is
                                     # not present
    } else if (any(iexact==0))
        stop("exact argument not found in the qmatrix: ", exact[iexact==0])

    # Grab the response and validate it
    Y <- model.response(mf)
    if (!inherits(Y, "Surv") || attr(Y, "type") != "mright")     
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
    id <- model.extract(mf, "id")
    id2 <- match(id, unique(id))
    index <- order(id2, Y[,1])
    if (any(diff(index) != 1)) stop("data not sorted by time within each id")

    # Get the first pass of the X matrix, and from that information create
    # cmap and tmap
    X <- model.matrix(Terms, mf)
    Xassign <- attr(X, "assign")
    # the next two are included in the output object, and used by later
    #  model.frame and model.matrix calls
    xlevels <- .getXlevels(delete.response(Terms), mf)
    contrasts <- attr(X, "contrasts")

    parse2 <- parsecovar2(parse1, statedata, dformula, Terms, qmatrix,
                          colnames(X), Xassign)
    cmap <- parse2$cmap # coefficients for the transitions
    tmap <- parse2$tmap # terms for the transitions

    # Now deal with missings, which we couldn't do before
    # Y and id can't be missing,
    # missing rate variables are filled in with lvcf
    idmiss <- is.na(id)
    ymiss <-  is.na(Y)
    first <- !duplicated(id)
    # apply last-value-carried-forward to the rate variables
    # this doesn't work if the first obs for a subject is missing
    xmiss <- apply(is.na(X[first,,drop=FALSE]), 2, any) 

   # we are cruel: anyone with a hole is no longer a valid timeline
    #  Toss all rows for that subject
    tossid <- unique(id[ymiss | idmiss | xmiss | is.na(weights)])
    if (length(tossid) >0 ) {
        keep <- !(id %in% tossid)
        na.action <- which(!keep)
        class(na.action) <- "omit"
        Y <- Y[keep,, drop=FALSE]
        X    <- X[keep,, drop=FALSE]
        id   <- id[keep]
        ytime <- ytime[keep]
        ystate <- ystat[keep]
        # message for printout
        removed <- c(subjects= length(toss), y= sum(ymiss), rate=sum(xmiss),
                     id= sum(idmiss))
    }
    else {
        na.action <- NULL
        removed <- NULL
    }

    # replace missing covariates for the rates using lvcf
    last <- !duplicated(id, fromLast=TRUE)
    for (i in 1:ncol(X)) {
        if (any(is.na(X[,i] & !last))) X[,i] <- lvcf(id, X[,i], ytime)
    }

    # Initialize the coefficients.  We do this before scaling X, since the
    #  user's view of coefficients is always on the original scale.
    # First icoef, then any overrides from options
    #  We allow initial values from a prior model that has fewer terms
    # 
    B <- 0*cmap
    if (!missing(icoef)) {
        if (!missing(intercept)) stop("only one of intercept or icoef is allowed")
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
    else {
        itemp <- simplehaz(Y, id, .2)
        B[1,] <- log(itemp[qmatrix>0])
    }
    # Import any init() or fixed() from the options
    # not yet implemented

    # Standardize the X matrix.  In the extrememly rare case that there
    # is a linear predictor that does not involve the intercept, e.g. a user had
    # factor(group)-1, we can't do so.
    if (TRUE) {
        # use the qr decomp (tentative idea)
        if (scale) {
            qrX <- qr(X)
            X <- qr.Q(qrX)
            btrans <- qr.R(qrX)
            xtrans <- solve(btrans)  #ToDo: will fail for singular X
            B <- btrans %*% B
        }
   } else {
       if (Xassign[1]!=0 || any(cmap[1,] ==0)) {
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
    
    fid <- match(id, unique(id)) # use id of 1, 2, 3, ... in fit routines
    # the otype variable: 1= interval censored = is in a known state at
    #  this time point; 2= exact = known state and we know exactly when it
    #  was entered (e.g. death), 0 = neither of these
    # interval censored or exact + markers = an error
    temp1 <- ystate %in% iexact
    otype <- 1L*(ystate>0) + 1L*temp1

    # set the maximizer function part of the icmfit call
    if (missing(mfun)) {
        if (iter ==0) {
            mfun <- "hmmloglik"
            if (is.null(mpar)) mpar <- list(fn= mfun)
            if (is.null(mpar$fn)) mpar$fn <- mfun
        }
        else { # the default
            mfun <- "hmmscore"
            if (is.null(mpar)) mpar$fn <- "hmmboth" # special for hmmscore
            if (is.null(mpar$fn)) mpar$fn <- "hmmboth"
            mpar$iter <- iter
            if (length(constraint)) mpar$constraint <- constraint
            mpar$debug <- control$debug
        }
    } else {
        if (!inherits(mfun, "function")) 
            stop("mfun argument must be a function")
        if (mfunname == "hmmscore") {
            if (is.null(mpar)) mpar <- list(fn= hmmloglik)
            if (is.null(mpar$fn)) mpar$fn <- hmmloglik
        } else if (mfunname == "hmmscore") {
            if (is.null(mpar)) mpar <- list(fn=hmmboth) # special for hmmscore
            if (is.null(mpar$fn)) mpar$fn <- hmmboth
            mpar$iter <- iter
            if (length(constraint)) mpar$constraint <- constraint
            mpar$debug <- control$debug
        } else if (mfunname== "optim") {
            if (is.null(mpar)) mpar$fn <- hmm1
            if (is.null(mpar$fn)) mpar$fn <- hmm1
            if (is.null(mpar$gr)) mpar$gr <- hmmgrad
            control <- mpar$control
            if (is.null(control$iter)) control$maxit <- iter
            if (is.null(control$fnscale)) control$fnscale <- -1
            mpar$control <- control
        }  else {
            ff <- names(formals(mfun))
            if (any(ff== "fn") && is.null(mpar$fn)) mpar$fn <- hmm1
            if (any(ff== "gr") && is.null(mpar$gr)) mpar$gr <- hmmgrad
            if (any(ff== "iter") && is.null(mpar$iter)) mpar$iter <- iter
        }
        if (!is.null(constraint) && any(ff== "constraint") && 
            is.null(mpar$constraint)) mpar$constraint <- constraint
    }

    fit <- icmfit(ytime, ystat, X, fid, otype, qmatrix, cmap, B, 
                  mfun, mpar, iter, constraint, penalty, mc.cores,      
                  control)
    browser()
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

