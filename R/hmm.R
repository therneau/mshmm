# The main function
hmm <- function(formula, data, subset, weights, na.action, 
                id, qmatrix, markers, iprob,
                mfun= hmmscore, mpar= list(), mgrad, iter=20,
                mc.cores= getOption("mc.cores", 2L),
                init, fixed, scale=TRUE, penalty, constraint,
                statedata, exact ="death",
                control= icmsh.control(), ...) {
    Call <- match.call()
    time0 <- proc.time()

    ## We want to pass any ... args to icmsh.control, but not pass things
    ##  like "dats=mydata", i.e., where someone made a typo.  The use of ...
    ##  is simply to allow things like "eps=1e6" with easier typing
    extraArgs <- list(...)
    if (length(extraArgs) && missing(control)) {
        controlargs <- names(formals(icmsh.control)) #legal arg names
        indx <- pmatch(names(extraArgs), controlargs, nomatch=0L)
        if (any(indx==0L))
            stop(gettextf("Argument %s not matched", 
                          names(extraArgs)[indx==0L]), domain = NA)
        control <- do.call(icmsh.control, extraArgs)
    } else if (missing(control)) control <- icmsh.control()


    # create a call to model.frame() that contains the formula (required)
    #  and any other of the relevant optional arguments
    #  but don't evaluate it just yet
    indx <- match(c("formula", "data", "subset", "weights", "na.action",
                    "id"), names(Call), nomatch=0)
    if (indx[1] ==0) stop("a formula argument is required")
    if (indx[6] ==0) stop("an id argument is required")
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

    # If the user did not provide an exact argument, don't complain if
    #  our default value of "death" is not one of the known states, instead
    #  treat the default as NULL
    if (missing(exact) || length(exact) ==0) exact <- NULL
    else {
        exact <- match(exact, statenames)
        if (any(is.na(exact)))
            stop("exact argument contains a state not in qmatrix")
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
        indx <- match(statename, statedata$state, nomatch=0)
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

    # nlp[3] = number of linear predictors (cols of cmap) that are
    #  devoted to rates, markers, and intial state.  The last is filled
    #  in further below in the initial state section.
    # For categorical markers we will also want to know the number of 
    #   categories (used to set up response functions)
    if (nmarker >0) {
        markerlevels <- sapply(marker1$marker, function(x) 
            length(levels(mf[[x]])))
        marker2 <- parsemarker2(marker1, stateddata, Terms, colnames(X), 
                                xassign, markerlevels)
        nlp <- c(ncol(cmap), ncol(marker2$cmap), 0)
        cmap <- cbind(cmap,
                      ifelse(marker2$cmap==0, 0, marker2$cmap +max(cmap))) 
    } else nlp <- c(ncol(cmap), 0, 0)

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
        # "time" will be the response, all states are latent
        # if the user didn't specify an exact argument, ignore our default
        if (!missing(exact)) # user specified one
            stop("the response must be a Surv object if there are exact states")
        exact <- NULL
        if (!is.numeric(Y)) stop("response must be numeric or Surv")
        ytime <- Y
        ystate <- rep(0L, nrow(mf))
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
        if (any(is.na(X[,i]))) X[,i] <- lvcf(id, X[,i], ytime)
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
        ytime <- ytime[keep]
        ystate <- ystate[keep]
        mf <- mf[keep,]  # the markers have not yet been pulled out
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
            init <- coef.to.B(priormod$coefficients, priormod$cmap, fixed=TRUE)
        }
        if (is.matrix(init)) {
            # allow for partial matching, so that a smaller model can feed a 
            #  larger
            rmatch <- match(row.names(init), row.names(cmap))
            cmatch <- match(col.names(init), col.names(cmap))
            if (any(is.na(rmatch))) 
                stop("init has covariates not in the current model")
            if (any(is.na(cmatch)))
                stop("init has linear predictors not in the current model")
            B[rmatch, cmatch] <- init
        } else if (is.numeric(init)) {
            if (!is.null(names(init))) {
                index <- match(names(init), param.names)
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
            rmatch <- match(row.names(fixed), row.names(cmap))
            cmatch <- match(col.names(fixed), col.names(cmap))
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

    # Standardize the X matrix.  In the extrememly rare case that there
    # is a linear predictor that does not involve the intercept, e.g. a user had
    # factor(group)-1, we can't do so.
    #  Markers are not in the X matrix, so not scaled
    if (xassign[1]!=0 || any(cmap[1,] ==0)) {
        if (scale)
            warning("not possible to scale the data")
        scale <- FALSE
    }
    if (scale && ncol(X) >1) {
        rvar <- 2:ncol(X) # don't scale the intercept!
        Xmean <-  rep(0, ncol(X))
        Xscale <- rep(1, ncol(X))
        Xmean[rvar] <- colMeans(X[,rvar])
        Xscale[rvar] <- apply(X[,rvar], 2, sd)
        for (i in rvar) X[,i] <- (X[,i]- Xmean[i])/Xscale[i]
        # we have XB = (X T^{-1}) (T B) where T is a transformation matrix
        #  don't forget the markers were exempt, only rvar cols transformed
        # see rescaling in the code vignette for more detail
        btrans <- diag(Xscale)   # T matrix, transforms B
        btrans[1, rvar] <- Xmean[rvar]
        xtrans <- diag(1/Xscale)      # btrans inverse, transforms X
        xtrans[1, rvar] <- -(Xmean/Xscale)[rvar]
        B <- btrans %*% B  # the coefs were in terms of unscaled X
    }

    param <- B.to.coef(B, cmap)  # don't use "coef" as variable name

    # preprocess constraint and penalty
    if (!missing(constraint)) {
        constraint <- hmmconstraint(constraint, Terms, cmap)
        constraint <- constraint %*% xtran # users write for untransformed X
    } else constraint <- NULL
    if (!missing(penalty)) {
        penalty <- hmmconstraint(penalty, Terms)
        penalty <- penalty %*% xtran # the are written for untransformed X
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
        stop("initial probability not available for all subjects")
    else if (is.formula(iprob))  
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
        ymarker <- mf[, unique(marker1$marker)]
        nmarker <- length(ymarker)
        rlist <- marker2$response
        eta <- X%*% B[,b2]
        for (i in 1:nmarker) {
            # check for a valid response vector
            if (!is.null(rlist[[i]]$check)) rlist[[i]]$check(ymarker[,i])
            keep <- !is.na(ymarker[,i])
            test <- rlist$rfun(ymarker[keep,i], eta[keep, rlist$rindex[[i]]])
            if (nrow(test) != nstate || ncol(test) != sum(keep)) 
                stop("wrong result from marker function ",i)
        }
    } else rlist <- NULL

    # the otype variable: 1= interval censored = is in a known state at
    #  this time point; 2= exact = known state and we know exactly when it
    #  was entered (e.g. death), 3= one or more markers, 0 = none of these
    # ymarker is a data.frame, not a matrix, hence sapply
    # What to do with an obs where markers are measured, but the state is
    #  known, i.e., ystate>0?  I think it is data dependent, so we warn the
    #  user. The current hmm1/hmm2 code will ignore the markers.
    temp1 <- ystate %in% iexact
    temp2 <- ystate %in% exactabsorb
    if (nmarker >0 ) {
        temp3 <- rowSums(sapply(ymarker, is.na)) >0
        if (any(temp3 & ystate>0))
            warning("marker variables present for an obs with known state")
    } else temp3 <- 0
    otype <- ifelse(temp1, 2L,
                    ifelse(ystate>0, 1L, 3L*temp3))
    
    # nlp has the number of colums of cmap for the rates, markers, and initial
    #  values. Now make two helpers
    # b1/b2/b3 are vectors containing the column number of cmap for each of
    #  the three, parmcount is the number of estimated coefficients 
    #  associated with each
    temp <- rep(1:3, nlp)
    b1 <- which(temp==1)
    b2 <- which(temp==2)
    b3 <- which(temp==3)
    utemp <- function(zed) 
        if (length(zed)>0) length(unique(zed[zed>0])) else 0L
    parmcount <- c(utemp(cmap[,b1]), utemp(cmap[,b2]), utemp(cmap[,b3]))

    # Creat the mapping matrix from linear predictors eta to the estimated
    #  portion of the coefficients vector (ignore fixed coefs), which is the
    #  part that the maximization function will "see".
    # See "derivatives, eta to beta" in the code vignette for details, along
    #  with a couple of open questions.
    eta.to.beta <- matrix(0, ncol(cmap), sum(parmcount))
    temp <- match(cmap, unique(cmap[cmap>0]), nomatch=0)
    ebindex1 <- (col(cmap)[temp>0] -1L)*sum(parmcount) + temp[temp>0]
    ebindex2 <- row(cmap)[temp>0]
    
    # Set up copies of the hmm1 and hmm2 functions to have the scope
    #  of this function. See "scope" in the code vignette for details.
    # I pass them down the calling chain as "logfun" as a way (I hope)
    #  to make it a little clearer in following routines that they are copies
    hmm1x <- hmm1; hmm2x <- hmm2
    environment(hmm1x) <- environment()
    environment(hmm2x) <- environment()

    # Set up parallel, Windows can't fork, others can
    if (mc.cores >1 && control$makecluster) {
        fork <- FALSE
        hmm_cluster <- makeCluster(mc.cores) #start up parallel
        }
    else fork <- TRUE
    time1 <- proc.time()

    # get the initial loglik and penalty
    rindex <- which(qmatrix >0)
    param <- B.to.coef(B, cmap)
    initial.loglik <- hmmloglik(param, B, cmap, id, mc.cores, fork, 
                                logfun= hmm1x)
    if (length(initial.loglik) ==0) 
        stop("unable to evaluate at the intial parameters")

    if (!is.null(penmat))
        penalty0 <- sum(param *(penmat %*% param))/2
    else penalty0 <- 0

    mfunname <- Call[["mfun"]]
    # if the default for mfun is changed in the hmm call, this needs to 
    #  change too, Call does not contain default arguments
    if (is.null(mfunname)) mfunname <- "hmmscore"

    # Do we need to iterate?
    # setting mfun to hmmloglik was an older way of doing 0 iterations, if
    #  so then we are already done via the initial loglik above
    if (mfunname== "hmmloglik" || (!missing(iter) && iter ==0)) {
        # no iteration
        loglik <- NULL
        penalty <- NULL
        fit <- NULL
    }
    else { # Yes, iterate
        if (mfunname== "hmmscore") {
            # hmmscore is currently the default maximizer
            # the user may have supplied an mpar arg with things specific
            #  to hmmscore.  Add onto it the par, fn, B, etc args
            mpar$par <- param
            mpar$fn <- hmmboth
            mpar <- c(mpar, list(B=B, cmap=cmap, id=id, 
                                 mc.cores= mc.cores, fork= fork, logfun= hmm2x))
            if (!missing(iter)) mpar$iter <- iter
            if (length(constraint)) mpar$constraint <- constraint
            mpar$debug <- control$debug
            hmm_count_of_calls <- c(0,0)  # a debugging line, see hmm2
            fit <- do.call(hmmscore, mpar)
            if (control$debug >0) print(hmm_count_of_calls)
        } else if (mfunname== "optim") {
            mpar$par <- param
            mpar$fn <- hmmloglik
            mpar$gr <- hmmgrad
            mpar <- c(mpar, list(B=B, cmap= cmap, id=id, 
                                 mc.cores= mc.cores, fork= fork))
            mpar$logfun <- hmm1x
            mpar$grfun  <- hmm2x
            if (is.null(control$iter)) control$maxit <- iter
            if (is.null(control$fnscale)) control$fnscale <- -1
            fit <- do.call(optim, mpar)
        } else if (mfunname== "mcmc0" || mfunname=="mcmc1") {
            mpar$par <- param
            mpar$logfun = "hmmloglik"  # no derivatives needed
            mpar <- c(mpar, list(B=B, cmap= cmap, id=id, 
                                 mc.cores= mc.cores, fork= fork, logfun=hmm1x))
            if (length(constraint)) mpar$constraint <- constraint
            fit <- do.call(mfunname, mpar)
        } else {
            # User supplied function.  We assume the first 3 args,
            #  whatever their names, are the starting estimate, the
            #  loglik function, and optionally the gradient. The
            #  grad argument tells us which
            # (par, fn) or (par, fn, gn)
            if (!inherits(mfun, "function"))
                stop("mfun argument must be a function")
           if (missing(mgrad)) 
                stop("mgrad is needed for a user supplied maximizer")
            else if (!(mgrad %in% 0:2))
                stop("valid mgrad arguments are 0-2")
            mfunarg <- formalArgs(mfun)  # the names of their args
            mpar[[mfunarg[1]]] <- param
            if (mgrad == 1) {
                mpar[[mfunarg[2]]] <- hmmboth
                mpar$logfun= hmm2x
            }
            else {
                mpar[[mfunarg[2]]]  <- hmmloglik
                mpar$logfun= hmm1x
            }
            if (mgrad==2) {
                mpar[[mfunarg[3]]] <- hmmgrad
                mpar$grfun <- hmm2x
            } 
            mpar <- c(mpar, list(B=B, cmap= cmap, id=id, 
                                 mc.cores= mc.cores, fork= fork))
            if (mgrad==2) mpar$gr <- hmmgrad
            ff <- names(formals(mfun))
            if (any(ff== "iter") && is.null(mpar$iter)) mpar$iter <- iter
            fit <- do.call(mfun, mpar)
        }
        
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
    }

    time2 <- proc.time()
    if (mc.cores > 1 & control$makecluster) stopCluster(hmm_cluster)

    # Undo any scaling and centering
    if (scale) {
        B <- coef.to.B(param, cmap, B)
        Bscale <- xtrans %*% B
        param <- B.to.coef(Bscale, cmap)
    }

    # Add nice dimnames
    bcol <- paste0(row(qmatrix)[qmatrix!=0], ":",
                   col(qmatrix)[qmatrix!=0])
    if (nlp[2]>0) {
        temp <- unique(rcoef[, c("marker", "lp")])
        lp <- temp$lp + 1 - temp$lp[match(temp$response, temp$response)]
        bcol <- c(bcol, paste0("M", paste(temp$response, lp, sep='.')))
    }
    if (nlp[3]>0) bcol <- c(bcol, paste0("p", 1:nlp[3]))
    dimnames(B) <- list(dimnames(X)[[2]], bcol)
    dimnames(cmap) <- dimnames(beta)


    time3 <- proc.time()
    compute.time <- rbind(setup= time1-time0,
                          compute= time2- time1,
                          finish = time3 - time2)

    final <- list(coefficients= param, 
                  loglik = c(intial=initial.loglik, final=loglik),
                  penalty= c(initial=penalty0, final=penalty),
                  beta=beta,
                  time = compute.time,
                  scale =  scale,
                  nlp = nlp,
                  cmap=cmap, rmap=rindex,
                  qmatrix = qmatrix,   # the structure and state names
                  nstate = nstate,
                  n = c(rows=nrow(mf), subjects=nid),
                  na.action = na.action,
                  removed = removed, 
                  call=Call,  xlevels=xlevels,
                  contrasts= attr(X, "contrasts"),
                  terms = Terms
                  )
    if (!is.null(penmat)) fit$pen.deriv <- pderiv
    if (!is.null(fit)) fit$fit <- fit

    class(final) <- "hmm"
    final
}
