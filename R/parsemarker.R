# Parse the list of markers for latent states.  This will be a list of
#  formulas, each left hand side a marker, predicted by a state and possibly
#  other covariates, with the options after a /
# Separate out the states, markers, and options.
parsemarker1 <- function(flist, statedata) {
    if (inherits(flist, "formula")) flist <- list(flist) # only one marker
    if (any(sapply(flist, function(x) !inherits(x, "formula"))))
        stop("an element of the marker list is not a formula")
    if (any(sapply(flist, length) != 3))
        stop("all formulas must have a left and right side")
    
    # some clever soul will use "log(pib) + log(tau) ~ A /gaussian" to specify
    # that both pib and tau are markers for A.  Deal with this by breaking the
    # above into two formulas.
    formbreak <- function(x) {
        if (as.name(x)== '('
    

    temp <- lapply(flist, function(x) rightslash(x[[3]])
    options <- lapply(temp, function(x) x[[2]]) # the list of options
    
    #  Save out the state information, which is the first term on the right hand
    # side of each formula.  It must be a column name in statedata.
    #  The reformulate call in the parent hmm routine works best if the response
    # (the marker) is moved to the right hand side of the formula.
    # is moved to the right hand side of the formula, which is simple.
    nform <- length(flist)
    statecol <- integer(nform)
    scol <- colnames(statedata)
    for (i in 1:nform) {
        rhs <- temp[[i]][[1]]
        if (length(rhs)==1) { # usual case, only a state name on the right
            first <- rhs[[1]]
            if (is.name(first)) jcol <- match(as.character(first), scol)    
            else if (is.character(first)) jcol <- match(first, scol)
            else stop("unrecognized state vector: ", deparse(first))

            if (is.na(jcol)) stop("unrecognized state vector: ", deparse(first))
            statecol[i] <- j
            flist[[i]][[3]] <- NULL  # left hand side is now right hand
        } else {
            # the rhs has state + covariates
            if (!is.call(rhs[[1]]) && rhs[[1]]== as.name("+"))
                stop("unrecognized state vector")
            first <- rhs[[2]]
            if (is.name(first)) jcol <- match(as.character(first), scol)    
            else if (is.character(first)) jcol <- match(first, scol)
            else stop("unrecognized state vector: ", deparse(first))

            if (is.na(jcol)) stop("unrecognized state vector: ", deparse(first))
            statecol[i] <- j
            rhs[[2]] <- flist[[i]][[2]]  # replace state with the lhs
            flist[[i]] <- rhs   #keep the right
        }
    }
        
    list(formula = flist, statecol=statecol, options=options)
}
