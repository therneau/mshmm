# The main function
hmm <- function(formula, data, subset, weights, na.action, 
                id, qmatrix, markers,
                pfun= hmminit, pcoef, entry, istate, 
                mfun=hmmtest, mpar= list(), 
                mc.cores= getOption("mc.cores", 2L),
                icoef, intercept, scale=TRUE, penalty, constraint,
                statedata, exact= "death",
                debug=0, fork=.Platform$OS.type=="unix") {
    Call <- match.call()
    time0 <- proc.time()
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
    if (all(qmatrix[row(qmatrix) > col(qmatrix)] == 0)) uppertri <- TRUE
    else  uppertri <- FALSE
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
    if (missing(markers)) stop("hmm model must have markers")
    marker1 <- parsemarker1(markers, statedata)
    nmarker <- length(marker1$marker)  # number of markers
    # the result has a separate list of markers (character) and covariates 
    #  for markers (list of NULL or char))
    
    # Deal with an initial formula (not yet done)
    iformula <- NULL
        
    # create the master formula, used for model.frame
    # the term.labels + reformulate + environment trio is used in [.terms;
    #  if it's good enough for base R it's good enough for me
    tlab <- attr(delete.response(terms(dform)), "term.labels") #rhs of dform
    if (!is.null(parse1))
        tlab <- c(tlab, unlist(lapply(parse1$rhs, function(x){
            attr(terms.formula(x), "term.labels")})))
    if (any(marker1$marker %in% tlab)) {
        stop("a variable can not be both a marker and a predictor")
        # The above test can be fooled: use log(pib) as a marker and pib for 
        #  a rate.
    }
    tlab <- c(tlab, unlist(marker1$mterm), marker1$marker) #markers last
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
    id2 <- match(id, unique(id))
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
    markerlevels <- sapply(marker1$marker, function(x) 
        length(levels(mf[[x]])))
    marker2 <- parsemarker2(marker1, stateddata, Terms, colnames(X), xassign,
                            markerlevels)
    bcount <- c(ncol(cmap), ncol(marker2$cmap), 0)
    cmap <- cbind(cmap,
                   ifelse(marker2$cmap==0, 0, marker2$cmap +max(cmap))#+ markers
    nparam <- max(cmap) # total estimated parameters
    # For transitions we will want only the first bcount[1] columns of cmap
    #  sometimes the marker columns or initial state cols, other times 
    #  we will want them all. Hence bcount.
                           
    # Now deal with missings, which we couldn't do before
    # Y and id can't be missing: toss those rows
    # markers can be missing
    # rate variables can only be missing for the last obs of a subject
    indx <- (rowSums(cmap >0) >0) # this variable is used in a rate
    xmiss <- (rowSums(is.na(X[,indx])) >0) # a missing in this row
    idmiss <- is.na(id)
    last <- !duplicated(id, fromLast=TRUE)
    xmiss <- xmiss & !last     # don't worry about a missing in a last row
    ymiss <- is.na(Y)
    
    # Check for missing values in istate as well, they are fatal
    if (!missing(istate) && any(missing(istate))
        stop("the istate argument cannot contain missing values")
    
    # we are cruel: anyone with a hole is no longer a valid timeline
    #  toss em
    tossid <- unique(id[ymiss | xmiss | idmiss])
    if (length(tossid) >0 ) {
        keep <- !(id %in% tossid)
        na.action <- which(!keep)
        class(na.action) <- "omit"
        Y <- Y[keep,, drop=FALSE]
        X    <- X[keep,, drop=FALSE]
        id   <- id[keep]
        # message for printout
        removed <- c(subjects= length(toss), y= sum(ymiss), rate=sum(xmiss),
                     id= sum(idmiss))
    }
    else {
        na.action <- NULL
        removed <- NULL
    }

    # If there are exact states, e.g., death, they are marked by the status 
    #  portion of a Surv() response
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
            ystat <- c(0, iexact)[Y[,2]] # recode
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
    ytime <- c(diff(Y[,1]), 0)  # time to next obs

    weights <- model.weights(mf)
    if (length(weights) >0) stop("weights are not yet supported")

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
        constraint <- hmmconstraint(constraint, Terms)
        constraint <- constraint %*% xtran # the are written for untransformed X
    } else constraint <- NULL
    if (!missing(penalty)) {
        penalty <- hmmconstraint(penalty, Terms)
        penalty <- penalty %*% xtran # the are written for untransformed X
        penmat <- crossprod(penalty)
    } else penmat <- NULL
    

    rindex <- which(qmatrix > 0)

    if (missing(entry)) entry <- rep(1.0, nstate)  # so it has no effect
    tempfun <- function(x) length(unique(x[x>0]))
    parmcount <- c(tempfun(cmap[,b1]), tempfun(cmap[,b2]), tempfun(cmap[,b3]))

    # initial probabilities for each subject
    if (missing(iprob)) iprob <- NULL
    else if (!is.null(marker1)) {
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


    

    # Do a dummy call to the response function(s), and make sure they
    #  return an object of the right shape.
    eta <- X%*% beta
    for (i in 1:ny) {
        keep <- (doresponse & !is.na(yobs[,i]))
        fit <- try(rfun[[i]](yobs[keep,i], nstate, 
                             eta[keep,b2map[[i]], drop=FALSE], gradient=FALSE))
        if (is.character(fit)) 
            stop("call failed for response function ", i, ": ", fit)
        else if (!is.numeric(fit)) stop("fitting function did not return numeric")
        if (length(fit) ==1) fit <- as.matrix(fit)
        if (!is.matrix(fit) || nrow(fit)!= nstate || ncol(fit)!= sum(keep))
            stop("response function must return a matrix with nstate rows and one column per response")
    }
    if (bcount[3] ==0 && is.null(iprob)) p0fixed <- pfun(nstate)
    else p0fixed <- NULL
    eps <- 1e-4
    hmm1 <- function(who, beta) {
        rows <- which(id ==uid[who])  # the subjects of interest
        eta <- X[rows,] %*% beta
        # starting probability
        if (!is.null(iprob)) alpha <- iprob[who,]
        else if (is.null(p0fixed))
            alpha <- pfun(nstate, eta[1,b3], gradient=FALSE)
        else alpha <- p0fixed
     
        # Now the response functions for this set
        rneed <- (otype[rows]==1 | otype[rows]==3)
        rlist <- vector("list", ny)
    #    cat("in hmm1\n"); browser()
        for (k in 1:ny) {
            j <- b2map[[k]]  #columns of beta for this response
            indx <- rneed & !is.na(yobs[rows, k])
            yy <- yobs[rows[indx], k]
            if (length(yy) >0) {
                if (length(j)==0) rlist[[k]] <-rfun[[k]](yy, nstate, gradient=FALSE)
                else rlist[[k]] <- rfun[[k]](yy, nstate, eta[indx, j, drop=FALSE], 
                                  gradient=FALSE)
    #            browser()
            }
        }
    #    cat("rlist done\n"); browser()
        
        # Compute the collection of matrix exponentials for the subject
        # The upper routine sends back the array of results as a vector
        #  along with the number of times there were tied eigenvalues
        #  we'll send the ties back as an attribute
        # We don't need the last row for each subject.
        # The call to upper uses the Ward approx (nterm=0) rather than
        #  the Higham09.  The former seems to better match my pade routine.
        r2 <- length(rows)  # there should always be at least 2 rows for a subject
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
        
            if (debug >2 & any(Pmat < -eps)) {
                cat ("stop1\n"); browser()}
            if (any(Pmat > (1+eps) | Pmat < -eps)) return("underflow")
            Pmat <- pmax(Pmat, 0)  # we sometimes get tiny negative numbers
        } 
        
        # Now walk through the visits one by one
        offset <- 0   # watch out for underflow
        nc <- integer(ny)  # the number of non-censored & non-missing so far
        rmat <- matrix(0., nstate, nstate)

        for (jj in seq_along(rows)) {
            j <- rows[jj] # j is the index in the original data, jj in our subset
            if (otype[j] == 2) {
                # exact event time (death)
                rmat[rindex] <- exp(eta[jj-1, b1]) #covariate just before this point
                dtemp <- rmat[,death]
                alpha <- alpha * dtemp
                 if (debug > 2) cat("A2: j=", j, "alpha=", alpha, "\n")
            }
            else if (otype[j] == 3) {
                # entry to the study
                temp <- sum(alpha*entry)
                if (temp <= 0) 
                  return(paste("subject", uid[who],"enters in an impossible state"))
                alpha <- alpha*entry /temp
                if (debug > 2) cat("A3: j=", j, "alpha=", alpha, "\n")
            }
        
            if (otype[j]==3 || otype[j]==1) {  # an outcome was observed
                temp <- rep(1, nstate)
                for (k in 1:ny) {
                    if (!is.na(yobs[j,k])) {
                        nc[k] <- nc[k] +1
                        alpha <- alpha* rlist[[k]][,nc[k]]
                        temp <- temp * rlist[[k]][,nc[k]]
                    }
                }
                if (!all(is.finite(alpha)) || sum(alpha) <=0) {
                    if (debug > 1) browser()
                    return("underflow") 
                }
                if (debug>2) cat("B: j=", j, "alpha=", alpha, "\n")
            }

            if (jj< r2) {  # if not the last obs
                # transition matrix
                alpha <- alpha %*% Pmat[,,jj]  # transition to next time point
                if (debug > 2) cat("C: j=", j, "alpha=", alpha, "\n")

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
            }
        }

        loglik <- offset + log(sum(alpha))
        attr(loglik, "counts") <- ucount
        loglik
    }
    makeindex <- function(cmap, all=cmap) {
        nonzero <- (cmap > 0)
        parms <- sort(unique(all[all>0]))  # the parameter numbers for this group
        p <- ncol(cmap)    # number of linear predictors
        k <- match(cmap[nonzero], parms)  #parameter number
        nparm <- length(parms)
        rr <- row(cmap)[nonzero]  # which X to use
        cc <- col(cmap)[nonzero]  #which eta this is
        list(xindex= rr, tindex= cc + (k-1)*p, dim=c(p, nparm))
        }

    Rtrans <- vector("list", ny)  #one element per response function
    if (bcount[2]) { #if there are response parameters
        tfun <- function(dmat, x, map) {
            tmat <- matrix(0., map$dim[1], map$dim[2])
            tmat[map$tindex] <- x[map$xindex]
            dmat %*% tmat
        }
        for (i in 1:ny) {
            if (length(b2map[[i]]) >0 && any(cmap[, b2map[[i]]] > 0)) {
                formals(tfun)[[3]] <- makeindex(cmap[,b2map[[i]], drop=FALSE],
                                                cmap[,b2])
                Rtrans[[i]] <- tfun
            }
        }
    }
    if (bcount[3]) { #if there are initial probability  parameters
        pitrans <- function(dmat, x, map) {
            tmat <- matrix(0., map$dim[1], map$dim[2])
            tmat[map$tindex] <- x[map$xindex]
            dmat %*% tmat
        }
        formals(pitrans)[[3]] <- makeindex(cmap[,b3, drop=FALSE])
    }
    cmap.b1 <- makeindex(cmap[,b1, drop=FALSE])
    Ptrans <- function(alpha, dmat, x, map=cmap.b1) {
        tmat <- matrix(0., map$dim[1], map$dim[2])
        tmat[map$tindex] <- x[map$xindex]
        #treat dmat as though it were a matrix with first dim nstate*nstate
        dim(dmat) <- c(nstate*nstate, map$dim[1])
        dmat2 <- dmat %*% tmat  #transform
        t(rowsum(dmat2 * rep(alpha, nstate*map$dim[2]), rep(1:nstate, each=nstate),
                 reorder=FALSE))
    }    
    if (!missing(death)) {  
        dtemp <- col(qmatrix)[rindex]
        deathcol  <- which(dtemp== death)  # list of linear predictors
        deathtrans <- function(R, x, map=cmap.b1) {
            rows <- which(qmatrix[,death] > 0)  #non-zero elements of column d
            dmat <- matrix(0, nstate, map$dim[1])
            for (i in 1:length(rows)) 
                dmat[rows[i], deathcol[i]] <- R[rows[i], death]

            tmat <- matrix(0., map$dim[1], map$dim[2])
            tmat[map$tindex] <- x[map$xindex]
            dmat %*% tmat
        }
    }
    psetup <- function(rmat, rindex) {
        n.eta <- length(rindex)
        out <- array(0., c(nstate, nstate, n.eta))
        temp <- matrix(0., nstate, nstate)
        rr <- row(temp)[rindex]
        for (i in 1:n.eta) {
            temp2 <- temp
            exp.eta <- rmat[rindex[i]]  # elements of rmat are exp(eta)
            temp2[rindex[i]] <-  exp.eta
            temp2[rr[i], rr[i]] <- -exp.eta
            out[,,i] <- temp2
            }
        out
    }
    hmm2 <- function(who, beta) {
        rows <- which(id ==uid[who])  # the subjects of interest
        eta <- X[rows,] %*% beta
        P.d  <- matrix(0., parmcount[1], nstate)
        R.d  <- matrix(0., nstate, parmcount[2])
        pi.d <- matrix(0., nstate, parmcount[3])

        # starting probability
        if (!is.null(iprob)) alpha <- iprob[who,]
        else if (is.null(p0fixed)) alpha <- pfun(nstate, eta[1,b3], gradient=TRUE)
        else alpha <- p0fixed
        if (bcount[3]) {
            pi.d <- pitrans(attr(alpha, 'gradient'), X[rows[1],])
            attr(alpha, 'gradient') <- NULL  # no longer needed
        }
        
        # Execute the response functions, over the uncensored obs
        rlist <- rgrad <- vector("list", ny)
        rneed <- (otype[rows] ==1 | otype[rows]==3)  # non-censored rows
        for (k in 1:ny) {
            index <- rneed & !is.na(yobs[rows,k])
            j <- b2map[[k]]  #linear predictors for this response
            yy <- yobs[rows[index], k]
            if (length(yy) > 0) {
                temp <- rfun[[k]](yy, nstate, eta[index, j, drop=FALSE], 
                                    gradient= TRUE)
                rlist[[k]] <- temp
                rgrad[[k]] <- attr(temp, "gradient")
            }
        }

        # Walk through the observations one by one
        P.d  <- matrix(0., pcount[1], nstate)
        offset <- 0  # watch out for underflow
        ecount <- c(length(rows), 0)
        nc <- integer(ny)    #number non-censored so far
        r2 <- length(rows)
        rmat <- matrix(0., nstate, nstate)

        for (jj in seq_along(rows)) {
            j <- rows[jj]
            if (otype[j] == 2 & jj> 1) {
                # exact event time (death)
                dtemp <- rmat[,death]  #rate at this point
                dtemp[death] <- 0      # this line should be redundant
                if (pcount[3]) pi.d <- pi.d * rep(dtemp, pcount[3])
                if (pcount[2]) R.d  <- R.d  * rep(dtemp, pcount[2])
                # Why the j-1 below?  A death density has to depend on covariates
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
            else if (otype[j] == 3) {
                # entry to the study
                temp <- sum(alpha*entry)
                if (temp <= 0) return(paste("subject", uid[who],
                                          "enters in an impossible state"))
                if (pcount[3]) pi.d <- (entry/temp)*( pi.d -
                                   alpha %*% (entry %*% pi.d)/temp )
                if (pcount[2]) R.d <- (entry/temp)* (R.d - 
                                   alpha %*% (entry %*% R.d)/temp)
                # remember that P.d is (nparm, nstate)
                if (pcount[1]) P.d  <- (rep(entry, each=pcount[1])/temp) * (P.d -
                                   outer(c(P.d %*% entry), alpha) /temp)
                alpha <- alpha*entry /temp
                if (debug > 2) {
                    cat("\n entry: alpha=", format(alpha), "\n")
                    #if (pcount[3]) print(pi.d)
                    # if (pcount[2]) print(R.d)
                    # print(P.d)
                }
                if (debug > 2) cat("A3: j=", j, "alpha=", alpha, "\n")
            }

            if (otype[j]==3 || otype[j]==1) {  # an outcome was observed
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
            
            if (jj < r2) { # not the last row
                # state matrix transformation P
                rmat[rindex] <- exp(eta[jj,b1])
                if (!all(is.finite(rmat))) {
                    # a horrible beta can overflow
                    if (debug > 1) 
                        save(rmat, beta, file=paste0("rfail", who, ".rda")) 
                    return("underflow")
                }
                
                diag(rmat) <- diag(rmat) -rowSums(rmat)
                tder <- psetup(rmat, rindex)
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
        }
        if (debug>3) browser()
        list(alpha=alpha, deriv= rbind(P.d, t(R.d), t(pi.d)), ecount=ecount,
             offset = offset)
    }
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

