# Parse the list of markers for latent states.  This will be a list of
#  formulas, each left hand side a marker, predicted by a state and possibly
#  other covariates, with the options after a /
# Separate out the states, markers, and options.
# It is possible to have multiples on the left, e.g. "pib + tau ~ A/lognormal"
#
parsemarker1 <- function(flist, statedata) {
    if (inherits(flist, "formula")) flist <- list(flist) # only one marker
    if (any(sapply(flist, function(x) !inherits(x, "formula"))))
        stop("an element of the marker list is not a formula")
    if (any(sapply(flist, length) != 3))
        stop("all formulas must have a left and right side")
    
    rhs2 <- lapply(flist, function(x) rightslash(x[[3]]))
    options <- lapply(rhs2, function(x) x[[2]]) # the list of options
    
    # The state information is saved simple vector, each element points to the
    #  column of statedata associated with the marker.  Markers are saved
    #  as a vector of names.  A list of formulas, appropriate for the the
    #  parent routine to create the model frame + the list of options rounds
    #  out the result.
    #  The reformulate call in the parent hmm routine works best if the response
    # (the marker) is moved to the right hand side of the formula.
    nform <- length(flist)
    sname    <- colnames(statedata)
    statecol <- nleft <- integer(nform)
    extraterm <- vector("list", nform) # extra predictors for the
    marker <- NULL
    for (i in 1:nform) {
        rhs <- rhs2[[i]][[1]]
        lhs <- flist[[i]][1:2] # keep the ~, so it looks like a 1 sided formula
        # grab the marker names
        temp <- attr(terms(lhs), "term.labels")
        nleft[i] <- length(temp)
        marker <- c(marker, temp)
        if (length(rhs)==1) { # usual case, only a state name on the right
            first <- rhs
            if (is.name(first)) jcol <- match(as.character(first), sname)    
            else if (is.character(first)) jcol <- match(first, sname)
            else stop("unrecognized state vector: ", deparse(first))

            if (is.na(jcol)) stop("unrecognized state vector: ", deparse(first))
            statecol[i] <- jcol
            flist[[i]][[3]] <- NULL  # left hand side is now right hand
        } else {
            if (!(is.call(rhs) && rhs[[1]]== as.name("+")))
                stop("unrecognized right hand side of marker formula")
            dummy <- ~ x
            #environment(dummy) <- environment(rhs)
            dummy[[2]] <- rhs[[3]]
            extraterm[[i]] <- dummy # save it as a 1 sided formula

            first <- rhs[[2]]
            if (is.name(first)) jcol <- match(as.character(first), sname)    
            else if (is.character(first)) jcol <- match(first, sname)
            else stop("unrecognized state vector: ", deparse(first))

            if (is.na(jcol)) stop("unrecognized state vector: ", deparse(first))
            statecol[i] <- jcol
            
            # this is a bit sneaky
            #  rhs is state + rest of covariates, element 1 is '+', 2 the state,
            #   3 the rest
            #  replace element 2 with flist[[[i]][[2]], which is the left hand
            #   side of the formula in flist
            #  replace flist[[i]][[2]] with rhs
            #  set flist[[i]][[3]] to NULL, making it a 1 sided formula
            # All this so that flist[[i]] is in the right form for the parent to
            #  use, so as to build a model.frame.
            #
            rhs[[2]] <- flist[[i]][[2]]  # replace state with the lhs
            flist[[i]][[2]] <- rhs   #keep the right
            flist[[i]] <- flist[[i]][1:2]  # make it 1 sided
        }
    }
      
    # nleft = number of markers on the left, in each equation
    # The statecol vector and options list need to replicate as well
    # extraterm contains any extra predictors, besides state

    list(formula = flist, marker= marker, statecol= rep(statecol, nleft), 
         options=options[rep(1:nform, nleft)], 
         extraterm= extraterm[rep(1:nform, nleft)])
}
