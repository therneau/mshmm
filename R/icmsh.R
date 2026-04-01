# The main function
icmsh <- function(formula, data, subset, weights, 
                id, qmatrix, markers, iprob, init, fixed,  
                penalty,  constraint, statedata,
                iter=30, exact= "death", mfun=hmmscore, mpar=list(), mfattr, 
                mc.cores=getOption("mc.cores", 2L), 
                control= hmm.control(), ...) {
    Call <- match.call()
    time0 <- proc.time()

    ## We want to pass any ... args to icmsh.control, but not pass things
    ##  like "dats=mydata", i.e., where someone made a typo.  The use of ...
    ##  is simply to allow things like "eps=1e6" with easier typing
    extraArgs <- list(...)
    if (length(extraArgs) && missing(control)) {
        controlargs <- names(formals(hmm.control)) #legal arg names
        indx <- pmatch(names(extraArgs), controlargs, nomatch=0L)
        if (any(indx==0L))
            stop(gettextf("Argument %s not matched", 
                          names(extraArgs)[indx==0L]), domain = NA)
        control <- do.call(hmm.control, extraArgs)
    } else if (missing(control)) control <- hmm.control()


    # create a call to model.frame() that contains the formula (required)
    #  and any other of the relevant optional arguments
    #  but don't evaluate it just yet
    indx <- match(c("formula", "data", "subset", "weights", 
                    "id"), names(Call), nomatch=0)
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
    # nlp = number of linear predictors used by transitions, markers, and
    #  initial state.  Set the first element, others will be done later.
    nlp <- integer(3)
    nlp[1] <- length(qmap)
    # mark each element of qmatrix with the linear predictor which it maps to
    #  this is used when computing for transitons to death
    qmatrix[qmap] <- 1:(nlp[1])

    # If the user did not provide an exact argument, don't complain if
    #  our default value of "death" is not one of the known states, instead
    #  treat the default as NULL
    temp <- match(exact, statenames)
    if (any(is.na(temp))) {
        if (missing(exact) || length(exact) ==0) exact <- NULL
        else stop("exact argument contains a state not in qmatrix")
    } 

    # a 0 row in qmap = an absorbing state (you never leave)
    absorb <- (rowSums(qmatrix>0) ==0)
    # an exact state that is absorbing: a censored obs is known to not 
    #  be in an 'exactabsorb' state at that time (there is no way to get there)
    exactabsorb <- which(absorb & (statenames %in% exact))

    if (!is.numeric(iter) || length(iter) >1 || iter <0) 
        stop("iter must be a non-negative integer")
    else iter <- ceiling(iter)

    # Is there state data?
    if (!missing(statedata)) { # check that it is okay
        if (!inherits(statedata, "data.frame"))
            stop("statedata must be a data frame")
        if (names(statedata)[1] != "state" || !is.character(statedata$state))
            stop("first variable in statedata must be a character variable named 'state'")
        indx <- match(statenames, statedata$state, nomatch=0)
        if (any(indx==0))
            stop("statedata$state does not contain all the states")
        statedata <- statedata[indx,]  # same row order as the states
        # Statedata might have rows for states that are not in the data set,
        #  for instance if the hmm call had used a subset argument.  Any of
        #  those are eliminated by the above line.
    } else statedata <- data.frame(state=statenames)
     
    # if the formula is a list, do the first level of processing on it,
    #  which is to pick off the list of variable names
    if (is.list(formula)) {
        if (length(formula)==1 && inherits(formula[[1]], "formula")) {
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

    # grab markers for the hidden states
    if (missing(markers)) nmarker <- 0
    else {
        marker1 <- parsemarker1(markers, statedata)
        nmarker <- length(unique(marker1$marker))  # number of markers
        # the result has a separate list of markers (character) and formulas for
        #  the covariates of the markers (most or all of which might be ~1)
    }
        
    # create the master formula, used for model.frame
    # the term.labels + reformulate + environment trio is used in [.terms;
    #  if it's good enough for base R it's good enough for me
    tlab <- attr(delete.response(terms(dformula)), "term.labels") #rhs of dform
    if (!is.null(parse1))
        tlab <- c(tlab, unlist(lapply(parse1$rhs, function(x){
            attr(terms.formula(x), "term.labels")})))
    if (nmarker > 0) {
        if (any(marker1$marker %in% tlab)) {
            stop("a variable can not be both a marker and a predictor")
            # The above test can be fooled: use log(pib) as a marker and pib 
            #  pib for a rate
        }
        mlab <- unlist(lapply(marker1$mterm, function(x) {
            attr(terms.formula(x), "term.labels")}))
        tlab <- c(tlab, mlab, marker1$marker) #markers last
    }
    newform <- reformulate(unique(tlab), dformula[[2]])
    environment(newform) <- environment(dformula)
    formula <- newform  # used for model.frame, not reported to user

    # Evaluate the expanded formula to create the model frame
    tform$formula <- formula
    mf <- eval(tform, parent.frame())
    if (nrow(mf) ==0) stop("data has 0 rows")

    # create a new Terms that doesn't have the marker variables, they don't
    #  become part of the X matrix 
    if (nmarker==0) Terms <- terms(mf)
    else {
        if (length(tlab) == nmarker) {
            # only a ~1 for the state transitions
            Terms <- terms(dformula)
        } else {
            dummy <- tlab[seq(from=1, to=length(tlab)- nmarker)]
            dummyform <- reformulate(unique(dummy), dformula[[2]])
            Terms <- terms(dummyform)
        }
    } 
    termnames <- attr(Terms, 'term.labels')

    # check that the data is sorted by time within subject, all rows for a
    # subject need to be contiguous. 
    Y <- model.response(mf)
    id2 <- model.extract(mf, "id")
    uid <- unique(id2) # the id's themselves need not be ordered
    id  <- match(id2, unique(id2))   # relabel as 1,2, etc for code
    if (is.matrix(Y)) index <- order(id2, Y[,1])
    else index <- order(id, Y)
    if (any(diff(index) != 1)) stop("data not sorted by time within each id")

    # Get the first pass of the X matrix, and from that information create
    # cmap and tmap (remove missings later)
    X <- model.matrix(Terms, mf)
    xassign <- attr(X, "assign")
    # next 2 will be saved as part of the hmm object, and used by later
    # model.frame or model.matrix calls with new data
    xlevels <- .getXlevels(delete.response(Terms), mf)
    contrasts <- attr(X, "contrasts")
    
    # Do the second pass on the formula and markers
    parse2 <- parsecovar2(parse1, statedata, dformula, Terms, qmatrix,
                          colnames(X), xassign)
    cmap <- parse2$cmap # coefficients for the transitions
    tmap <- parse2$tmap # terms for the transitions

    # number of linear predictors (cols of cmap) that are
    #  devoted to rates, markers, and intial state.  The last is filled
    #  in further below in the initial state section.
    # For categorical markers we will also want to know the number of 
    #   categories (used to set up response functions)
    # 
    if (nmarker >0) {
        markerlevels <- lapply(marker1$marker, function(x) 
            levels(mf[[x]]))
        marker2 <- parsemarker2(marker1, statedata, Terms, colnames(X), 
                                xassign, markerlevels)
        nlp[2] <- ncol(marker2$cmap)
        if (nlp[2] ==0) b2 <- 0 else b2 <- nlp[1] + 1:nlp[2]
        cmap <- cbind(cmap,
                      ifelse(marker2$cmap==0, 0, marker2$cmap +max(cmap))) 
    } 

    nparam <- max(cmap) # total estimated parameters
                           
    # mark out the exact states as an integer, used in the fitting portion
    if (length(exact)==0) iexact <- NULL
    else iexact <- match(exact, statenames)

    # The response will normally be a Surv object, with known states as the
    # status
    if (inherits(Y, "Surv")) {
        if (attr(Y, "type") == "right") {
            if (length(exactabsorb)==1) {
                # special case: 0/1 status can be used if there is only one
                # exact absorbing state
                # make it 0/exact rather than 0/1
                ystate <- ifelse(Y[,2]==0, 0L, exactabsorb)
            } else stop("simple Surv() only allowed if there is 1 exact state")
        } else if (attr(Y, "type") == "mright") {
            ystate <- attr(Y, "states")
            # all states in Y must appear in statenames, but are allowed to
            #  be in a different order. ystat will be in statenames order, with
            #  0= censored
            temp <- match(ystate, statenames)
            if (any(is.na(temp)))
                stop("response has a state not found in qmatrix")
            ystate <- c(0L, temp)[1 + Y[,2]]
         }
        ytime <- Y[,1]
    } else {
        # "time" will be the response, all states are latent.  Unusual since
        # death is normally one of the states, and it is not latent.
        # If the user didn't specify an exact argument, ignore our default
        #  of 'death', i.e.,don't give an error that there is value for the 
        #  exact option that is not in the set of states.
        if (!is.null(exact)) # user specified one
            stop("the response must be a Surv object if there are exact states")
        if (!is.numeric(Y)) stop("response must be numeric or Surv")
        ytime <- Y
        ystate <- rep(0L, nrow(mf)) # all censored
    }   

    weights <- model.weights(mf)
    if (length(weights) >0) warning("weights are not yet supported")
    weights <- rep(1.0, nrow(X))

    # Now deal with missings, which we couldn't do before
    # Y, id and weights can't be missing
    # rate variables will be fixed up using lvcf
    # markers can be missing
    idmiss <- is.na(id)
    ymiss <-  is.na(Y)
    first <- !duplicated(id)
    # apply last-value-carried-forward to the rate variables
    # this doesn't work (statistically) if the first obs for a person is NA
    for (i in 1:ncol(X)) {
        if (any(is.na(X[,i]))) X[,i] <- lvcf(id, X[,i])
    }
    xfirstmiss <- apply(is.na(X[first,,drop=FALSE]),1, any)         
    # we are cruel: anyone with a hole is no longer a valid timeline
    #  toss the entire subject
    tossid <- unique(c(id[ymiss | idmiss | is.na(weights)], 
                       (id[first])[xfirstmiss]))
    if (length(tossid) >0 ) {
        keep <- !(id %in% tossid)
        na.action <- which(!keep)
        class(na.action) <- "omit"
        Y <- Y[keep,, drop=FALSE]
        X    <- X[keep,, drop=FALSE]
        id   <- id[keep]
        id   <- match(id, unique(id)) # renumber them as 1,2 3...
        ytime <- ytime[keep]
        ystate <- ystate[keep]
        weights <- weights[keep]
        mf <- mf[keep,]  # the markers have not yet been pulled out
        # message for printout
        removed <- c(subjects= length(tossid), y= sum(ymiss),
                     rate=sum(xfirstmiss),
                     id= sum(idmiss))
    }
    else {
        na.action <- NULL
        removed <- NULL
    }
    
    # Initialize the coefficients.  We do this before scaling X, since the
    #  user's view of coefficients is always on the original scale.
    # First init, then any overrides from options
    #  We allow initial values from a prior model which might have fewer terms
    # 
    B <- 0*cmap
    indx <- match(1:nparam, cmap)  # nparam = number of unique coefs
    beta.names <- paste(rownames(cmap)[indx], colnames(cmap)[indx], sep='.')
    param <- rep(0, nparam)
    if (!missing(init)) {
        if (inherits(init, "hmm")) { # a prior hmm model
            priormod <- init
            init <- coef(priormod, matrix=TRUE, fixed=TRUE)
        }
        if (is.matrix(init)) {
            # allow for partial matching, so that a smaller model can feed a 
            #  larger
            rmatch <- match(rownames(init), rownames(cmap))
            cmatch <- match(colnames(init), colnames(cmap))
            if (any(is.na(rmatch))) 
                stop("init has covariates not in the current model")
            if (any(is.na(cmatch)))
                stop("init has linear predictors not in the current model")
            B[rmatch, cmatch] <- init
        } else if (is.numeric(init)) {
            if (!is.null(names(init))) {
                index <- match(names(init), beta.names)
                if (any(is.na(index)))
                    stop("init has an coefficient not found in the model: ",
                         (names(init)[is.na(index)])[1])
                else param[index] <- init
            } else {
                if (length(init) != nparam) stop("wrong length for init")
                else param <- init
            }
            B <- coef.to.B(param, cmap, B)
        } else stop("init must be a numeric vector, matrix, or prior fit")
    }

    if (!missing(fixed)) {
        if (is.matrix(fixed)) {
            rmatch <- match(rownames(fixed), rownames(cmap))
            cmatch <- match(colnames(fixed), colnames(cmap))
            if (any(is.na(rmatch))) 
                stop("fixed has covariates not in the current model")
            if (any(is.na(cmatch)))
                stop("fixed has linear predictors not in the current model")
            if (is.logical(test)) cmap[test] <- -cmap(test)
            else if (is.numeric(test)) cmap[test>0] <- - cmap[test>0]
            else stop("fixed argument must be logical or numeric")
        } else { # it should be a vector
            fixed <- asLogical(fixed)  # change numeric to T/F
            if (!is.null(names(fixed))) {
                index <- match(names(fixed), param.names)
                if (any(is.na(index)))
                    stop("fixed has an coefficient not found in the model: ",
                         (names(fixed)[is.na(index)])[1])
                i2 <- which(cmap %in% index[fixed])
                cmap[i2] <- -cmap[i2]  # a negative index marks it
            }
            else if (length(fixed) != nparam) stop("wrong length for fixed")
            else {                
                i2 <- (which(cmap>0))[fixed]
                cmap[i2] <- -cmap[i2]
            }
        } 
    }

    # Import any init() or fixed() from the options
    # not yet done

    # Standardize the X matrix.  All linear preditors must include an intercept
    # To do otherwise, e.g., allow "~ group -1" as a formula, makes our formula
    #   processing just too difficult.
    # Markers are not in the X matrix, so don't get scaled
    if (xassign[1]!=0 || any(cmap[1,] ==0)) 
        stop("-1 in formulas not allowed")
    Xmean <-  rep(0, ncol(X)) # don't scale
    Xscale <- rep(1, ncol(X))
    if ((control$scale || control$center) && ncol(X) >1) {
        if (control$scale & !control$center) {
            warning("scale=TRUE implies center=TRUE")
            doscale <- docenter <- TRUE
            control$center <- TRUE
        }
        rvar <- 2:ncol(X) # don't scale the intercept!
        Xmean <-  rep(0, ncol(X))
        Xscale <- rep(1, ncol(X))
        if (control$center) Xmean[rvar] <- colMeans(X[,rvar, drop=FALSE])
        if (control$scale)  Xscale[rvar] <- apply(X[,rvar,drop=FALSE], 2, sd)
        for (i in rvar) X[,i] <- (X[,i]- Xmean[i])/Xscale[i]
        # we have XB = (X T^{-1}) (T B) where T is a transformation matrix
        #  don't forget the markers were exempt, only rvar cols transformed
        # see rescaling in the code vignette for more detail
        btrans <- diag(Xscale)   # T matrix, transforms B
        btrans[1, rvar] <- Xmean[rvar]
        xtrans <- diag(1/Xscale)      # btrans inverse, transforms X
        xtrans[1, rvar] <- -(Xmean/Xscale)[rvar]
        B <- btrans %*% B  # the coefs were in terms of unscaled X
    } else {
        xtrans <- diag(ncol(X)) # no transformation
        btrans <- xtrans
    }
    param <- B.to.coef(B, cmap)  # don't use "coef" as variable name

    # preprocess constraint and penalty
    if (!missing(constraint)) {
        constraint <- hmmconstraint(constraint, Terms, cmap)
        constraint <- constraint %*% xtrans # users write for untransformed X
    } else constraint <- NULL
    if (!missing(penalty)) {
        penalty <- hmmconstraint(penalty, Terms)
        penalty <- penalty %*% xtrans # the are written for untransformed X
        penmat <- crossprod(penalty)
    } else penmat <- NULL
    
    # initial probabilities for each subject
    # Anyone who has a non-censored first obs is assumed to start in
    #  that state
    fstate <- ystate[!duplicated(id)]  # first state for each subject
    nid <- length(fstate)
    if (all(fstate>0)) { 
        # simple interval censored, or known start, iprob ignored if present
        iprob <- matrix(0, nid, nstate)
        iprob[cbind(1:nid), fstate] <- 1
    }
    else if (missing(iprob))
        stop("iprob argument is required")
    else if (inherits(iprob, "formula"))  
        stop("iprob = formula, code not yet completed")
    else if (!is.numeric(iprob))
        stop("iprob must be numeric or a formula")
    else {
        if (any(is.na(iprob))) stop("iprob cannot contain missing values")
        if (any(iprob<0) || any(iprob >1))
            stop("iprob must contain values between 0 and 1")
        if (is.vector(iprob)) {
            if (length(iprob) != nstate) stop("wrong length for iprob")
            iprob <- iprob/sum(iprob)
            # an init per subject makes later code a bit easier
            iprob <- matrix(rep(iprob, each= length(uid)), ncol=nstate)
        } else if (!(is.matrix(iprob) && length(dim(iprob))==2))
            stop ("iprob must be a vector or matrix")
        else {
            if (ncol(iprob) != nstate)
                stop("iprob should have one column per state")
            iname <- rownames(iprob)
            if (!is.null(iname)) {
                indx <- match(uid, iname)
                if (any(is.na(indx))) {
                    badid <- uid[is.na(indx)]
                    # only print the first few
                    if (length(badid) > 4) badid <- badid[1:4]
                    stop("id values not found in iprob:", 
                         paste(badid, collapse=' '))
                    }
                # if anyone was removed due to missing, iprob[indx,] makes sure
                #  the right initial prob is used
                iprob <- iprob/rowSums(iprob)
                iprob <- iprob[indx,]
            } else stop("an iprob matrix must have id as the row names")
        }
    }  

    if (!(all(fstate> 0) || all(fstate==0))) {
        # someone is using exact start for some, and probabalistic start
        # for others. 2:1 odds this was a mistake
        warning("Exact starting state for some subjects but not others. ",
                "Really?") 
        for (i in 1:nid) {
            if (fstate[i] !=0) {
                iprob[i,] <- 0
                iprob[i, fstate[i]] <-1
            }
        }
    }

    # Set up response functions
    if (nmarker > 0) {
        ymarker <- mf[unique(marker1$marker)]  # a data frame
        nmarker <- ncol(ymarker)
        rlist <- marker2$response
        if (b2>0) eta <- X%*% B[,b2]
        for (i in 1:nmarker) {
            # check for a valid response vector
            if (!is.null(rlist[[i]]$checkfun)) rlist[[i]]$checkfun(ymarker[,i])
            tfun <- rlist[[i]]$rfun
            keep <- !is.na(ymarker[,i])
            if (length(rlist$e2map[[i]]) >0 ) {
                # the e2map indices from parsemarker are within the cmap cols
                #  for markers, hmm1 and hmm2 will want overall cmap index
                j <- rlist$e2map[[i]] + nlp[1]
                test <- tfun(ymarker[keep,i], eta[keep, j])
                }
            else test <- tfun(ymarker[keep,i]) #e.g. a fixed missclass matrix
            if (nrow(test) != nstate || ncol(test) != sum(keep)) 
                stop("wrong result from marker function ",i)
        }
    } else {
        rlist <- NULL
        ymarker <- NULL
    }

    if (missing(mfun)) mfunname <- "hmmscore"
    else mfunname <- Call[["mfun"]]
    # if the default for mfun is changed in the hmm call, the above needs to 
    #  change too, since Call does not contain default arguments
    if (mfunname== "hmmloglik") { # no iteration
        iter <- 0
        mfattr <- list(param="par", result=1)
    }
    else if (mfunname == "hmmscore") 
        mfattr<- list(param= "par", fn="fn", deriv=TRUE, gfun=NULL, 
                     iter= "iter", result=2)
    else if (mfunname== "optim")
        mfattr <- list(param="par", fn="fn", deriv= TRUE, gfun="gr", 
                     iter= list(control="maxit"), result=3)
    else if (mfunname== "hmm1" || mfunname== "hmm2") 
        mfattr <- list(param="par", fun="fn", deriv=FALSE, 
                       iter= "iter", result=4)
    else if (!missing(mfattr)) {
        mfname <- c("param", "fun", "deriv", "iter", "result")
        if (any(is.na(match(mfname, names(mfattr)))))
            stop("the mfattr argument is not complete")
    } else stop("user supplied optimizer must include mfattr argument")
    time1 <- proc.time()
    mfit <- msh.fit(id, ytime, ystate, X, iprob, B,
                     cmap, nlp, ymarker, rlist, qmatrix,
                     mc.cores, control, mfun, mfattr, mpar, iter,
                     iexact, conmat, penmat)


    time2 <- proc.time()
    if (mc.cores > 1 & control$makecluster) stopCluster(hmm_cluster)

    # Undo any scaling and centering
    if (control$scale || control$center) {
        B <- coef.to.B(mfit$param, cmap, B)
        Bscale <- xtrans %*% B
        param <- B.to.coef(Bscale, cmap, fixed=TRUE)
    }

    # Add names
    pname <- outer(rownames(cmap), colnames(cmap), paste, sep='_')
    names(param) <- pname[match(unique(cmap[cmap>0]), cmap)]
    compute.time <- rbind(setup= time1-time0,
                          compute= time2- time1)

    final <- list(coefficients= param, 
                  loglik = mfit$loglik,
                  time = compute.time,
                  cmap= cmap, nlp= nlp,
                  qmatrix = qmatrix,   # the structure and state names
                  n = c(observations =nrow(mf), id =nid)
                  )
    if (!is.null(removed)) final$removed <- removed
    if (!is.null(mfit$penalty) && mfit$penalty >0)   
        final$penalty <- c(initial= penalty0, final= mfit$penalty)
    if (!is.null(mfit$fit)) final$fit <- mfit$fit
    if (control$center && ncol(X) >1) final$xmean <- Xmean
    if (control$scale  && ncol(X) >1) final$xscale <- Xscale
    if (control$detail) final <- c(final, 
                                   list(alpha=mfit$alpha, deriv=mfit$deriv))
    final <- c(final, list(call=Call,  xlevels=xlevels,
                  contrasts= attr(X, "contrasts"),
                  terms = Terms))
    class(final) <- "icmsh"
    final
}
