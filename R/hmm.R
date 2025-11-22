# The main function
hmm <- function(formula, data, subset, weights, na.action, 
                id, qmatrix, markers,
                pfun= hmminit, pcoef, entry, istate, 
                mfun, mpar= list(), iter=20,
                mc.cores= getOption("mc.cores", 2L),
                icoef, intercept, scale=TRUE, penalty, constraint,
                statedata, exact= "death", control= cmsh.control(), ...) {
    Call <- match.call()
    time0 <- proc.time()

    ## We want to pass any ... args to cmsp.control, but not pass things
    ##  like "dats=mydata" where someone just made a typo.  The use of ...
    ##  is simply to allow things like "eps=1e6" with easier typing
    extraArgs <- list(...)
    if (length(extraArgs)) {
        controlargs <- names(formals(cmsp.control)) #legal arg names
        indx <- pmatch(names(extraArgs), controlargs, nomatch=0L)
        if (any(indx==0L))
            stop(gettextf("Argument %s not matched", 
                          names(extraArgs)[indx==0L]), domain = NA)
    }

    # create a call to model.frame() that contains the formula (required)
    #  and any other of the relevant optional arguments
    #  but don't evaluate it just yet
    indx <- match(c("formula", "data", "subset", "weights", "na.action",
                    "id", "istate"), names(Call), nomatch=0)
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

    if (!is.numeric(iter) || length(iter) >1 || iter <=0) 
        stop("iter must be a postive integer")
    else iter <- ceiling(iter)

    if (missing(exact)) {
        # don't complain if the default value of "death" is not in the list,
        #  death might not be a state for this model
        if (is.na(match(exact, statenames))) exact <- NULL
    } else {
        exact <- match(exact, statenames)
        if (any(is.na(exact)))
            stop("exact argument contains a state not in qmatrix")
    }

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
    # Deal with an initial formula (not yet done)
    iformula <- NULL
        
    # create the master formula, used for model.frame
    # the term.labels + reformulate + environment trio is used in [.terms;
    #  if it's good enough for base R it's good enough for me
    tlab <- attr(delete.response(terms(dform)), "term.labels") #rhs of dform
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
    dummy <- tlab[seq(from=1, to=length(tlab)- nmarker)]
    dummyform <- reformulate(unique(dummy), dformula[[2]])
    Terms <- terms(dummyform)
    termnames <- attr(Terms, 'term.labels')

    # check that the data is sorted by time within subject, all rows for a
    # subject need to be contiguous. 
    Y <- model.response(mf)
    id <- model.extract(mf, "id")
    if (is.matrix(Y)) index <- order(id, Y[,1])
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

    # for categorical markers we will want to know the number of categories
    if (nmarker >0) {
        markerlevels <- sapply(marker1$marker, function(x) 
            length(levels(mf[[x]])))
        marker2 <- parsemarker2(marker1, stateddata, Terms, colnames(X), 
                                xassign, markerlevels)
        bcount <- c(ncol(cmap), ncol(marker2$cmap), 0)
        cmap <- cbind(cmap,
                      ifelse(marker2$cmap==0, 0, marker2$cmap +max(cmap))) 
    } else bcount <- c(ncol(cmap), 0, 0)
    nparam <- max(cmap) # total estimated parameters
    # For transitions we will want only the first bcount[1] columns of cmap
    #  sometimes the marker columns or initial state cols, other times 
    #  we will want them all. Hence bcount.
                           
    # The response will normally be a Surv object, with known states as the
    # status
    if (inherits(Y, "Surv")) {
        if (attr(y, "type") == "right") {
            if (length(exact)==1 && exact %in% statenames) {
                # special case: 0/1 status can be used if there is 1 exact state
                iexact <- match(exact, statenames)
                ystat <- Y[,2]* iexact  # make it  0/exact instead of 0/1
            } else stop("simple Surv() only allowed if there is 1 exact state")
        } else if (attr(Y, "type") == "mright") {
            ystate <- attr(Y, "states")
            iexact <- (match(ystate, statenames))
            if (any(is.na(iexact)))
                stop("response has a state not found in qmatrix")
            ystat <- Y[,2] -1L  # 0 = censored

            # states with biomarkers should not appear in ystate
            if (nmarker >1) {
                browser()  # to be filled in
            }
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
        ystat <- rep(0, nrow(mf))
    }   

    weights <- model.weights(mf)
    if (length(weights) >0) warning("weights are not yet supported")
    weights <- rep(1.0, nrow(X))

    # Now deal with missings, which we couldn't do before
    # Y, id and weights can't be missing
    # markers can be missing
    # rate variables can only be missing for the last obs of a subject
    idmiss <- is.na(id)
    ymiss <-  is.na(Y)
    first <- !duplicated(id)
    xmiss <- apply(is.na(X[first,,drop=FALSE], 2, any)) #LVCF won't work here

    if (!missing(istate) && any(missing(istate)))
        stop("the istate argument cannot contain missing values")

    # we are cruel: anyone with a hole is no longer a valid timeline
    #  toss the entire subject
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
    temp <- rbind(cmap, mmap)
    indx <- match(1:nparam, temp)
    param.names <- paste(rownames(temp)[indx], colnames(temp)[indx], sep='.')
    param <- rep(0, nparam)
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
            if (!is.null(names(icoef))) {
                index <- match(names(icoef), param.names)
                if (any(is.na(index)))
                    stop("icoef has an coefficient not found in the model: ",
                         (names(icoef)[is.na(index)])[1])
                else param[index] <- icoef
            } else {
                if (length(icoef) != nparam) stop("wrong length for icoef")
                else param <- icoef
            }
            B <- coef.to.B(param, cmap, mmap)
        } else stop("icoef must be a numeric vector or matrix")
    }

    # Import any init() or fixed() from the options
    # not yet done

    # Standardize the X matrix.  In the extrememly rare case that there
    # is a linear predictor that does not involve the intercept, e.g. a user had
    # factor(group)-1, we can't do so.
    #  Markers are not in the X matrix, so not scaled
    if (Xassign[1]!=0 || any(cmap[1,] ==0)) {
        if (!missing(scale) && scale)
            warning("not possible to scale the data")
        scale <- FALSE
    }
    if (scale) {
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
        xtrans <- diag(1/Xscale)      # T inverse, transforms X
        xtrans[1, rvar] <- -(Xmean/Xscale)[rvar]
        B <- btrans %*% B  # the coefs were in terms of unscaled X
    }
    
    if (!missing(intercept)) {
        if (!scale) stop("intercept initialization requires scale=TRUE")
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
            # either the rates, or all, is allowed
            if (length(intercept) == bcount[1]) 
                B[1, 1:bcount[1]] <- intercept
            else if (length(intercept= ncol(cmap)))
                B[1,] <- intercept
            else stop("wrong length for intercept")
        }
    }

    b1 <- 1:bcount[1]
    b2 <- seq(bcount[1]+1, length=bcount[2])  
    b3 <- seq(bcount[1] + bcount[2] +1, length=bcount[3]) #might be none
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
    

    if (missing(entry)) entry <- rep(1.0, nstate)  # so it has no effect
    tempfun <- function(x) length(unique(x[x>0]))
    parmcount <- c(tempfun(cmap[,b1]), tempfun(cmap[,b2]), tempfun(cmap[,b3]))

    # initial probabilities for each subject
    if (missing(iprob)) iprob <- NULL
    else if (nmarker ==0) {
        warning("no markers, iprob ignored")
        iprob <- NULL
    } else {
         if (!is.numeric(iprob) || any(iprob<0) || any(iprob >1))
            stop("iprob must contain values between 0 and 1")
        if (is.vector(iprob)) {
            if (length(iprob) != nstate) stop("wrong length for iprob")
            iprob <- iprob/sum(iprob)
            # an init per subject makes later code a bit easier
            iprob <- matrix(rep(iprob, each= length(uid)), ncol=nstate)
        } else {
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
                # if anyone was removed due to missing, above makes sure that
                #  the right initial prob is used
                iprob <- iprob/rowSums(iprob)
                iprob <- iprob[indx,]
                
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

    # the otype variable: 1= interval censored, 2= exact, 3= has one or
    #  more markers, 0 = none of the above
    temp1 <- ystat %in% iexact
    temp2 <- ystat >0
    if (nmarker >0 ) {
        temp3 <- rowSums(sapply(ymarker, is.na))
        otype <- ifelse(temp1, 2L, ifelse(temp2, 1L, 3L*(temp3>0)))
    } else otype <- ifelse(temp1, 2L, 1L * temp2)


    # initial values
    if (bcount[3] ==0 && is.null(iprob)) p0fixed <- pfun(nstate)
    else p0fixed <- NULL

    # set the argument list for the function call
    if (missing(mfun)) {
        if (iter ==0) mfun <- hmmloglik
        else { # the default
            mfun <- hmmscore
            mpar$fn <- hmmboth # special for hmmscore
            mpar$iter <- iter
            if (length(constraint)) mpar$constraint <- constraint
            mpar$debug <- control$debug
        }
    } else {
        if (!inherits(mfun, "function")) 
            stop("mfun argument must be a function")
        ff <- names(formals(mfun))
        if (any(ff== "fn") && is.null(mpar$fn)) mpar$fn <- hmm1
        if (any(ff== "gr") && is.null(mpar$gr)) mpar$gr <- hmmgrad
        if (any(ff== "iter") && is.null(mpar$iter)) mpar$iter <- iter
        if (!is.null(constraint) && any(ff== "constraint") && 
            is.null(mpar$constraint)) mpar$constraint <- constraint
    }

    hfit <- hmmfit(ytime, ystat, X, id, otype, qmatrix, cmap, rlist,
                    mfun, mpar, iter, penalty, control)

