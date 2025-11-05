# Parse the list of markers for latent states.  This will be a list of
#  formulas, each left hand side a marker, predicted by a state and possibly
#  other covariates, with the options after a /
# Separate out the states, markers, and options.
#
parsemarker1 <- function(flist, statedata) {
    if (inherits(flist, "formula")) flist <- list(flist) # only one marker
    if (any(sapply(flist, function(x) !inherits(x, "formula"))))
        stop("an element of the marker list is not a formula")
    if (any(sapply(flist, length) != 3))
        stop("all formulas must have a left and right side")
    
    rhs2 <- lapply(flist, function(x) rightslash(x[[3]]))
    options <- lapply(rhs2, function(x) x[[2]]) # the list of options
    
    # Common formula:  log(pib) ~ A / gaussian
    # Harder case: pzmemory + pzglobal ~ N + practice + sex / gaussian + common
    # 1. Return "~ pib" as the formula for case 1, for case 2 return 
    #    "~ pzmemory + pzglobal + practice + sex"; two on the left fails in
    #   the reformulate() code in hmm.
    # 2. For case 2, the returned marker component will have a separate element
    #   for each marker, statecol and options will also have one element
    #   per marker.
    # 3. ptplot(~N + practice + sex) shows that we can't just use a simple
    # "grab left term of the first +" strategy to find the state, "N" is further
    #  down the parse tree than that.  Hence the grableft() function.
    #  It returns a list of length 2, first element will be "N", second
    #  'practice + sex'.  See 'parsing' in the code vignette for 
    #  deeper explanations.
    grableft <- function(x) {
        if (length(x)==1) return(list(x, NULL))  #only a single term
        if (is.call(x) && x[[1]]== as.name('+')) {
            if (is.name(x[[2]])) list(x[[2]], x[[3]])
            else {
                temp <- grableft(x[[2]])
                x[[2]] <- temp[[2]]
                list(temp[[1]], x)
            }
        } else stop("state identifier can only be followed by '+'")
    }
            
    nform <- length(flist)
    nmarker <- integer(nform)
    stateinfo <- vector("list", nform)  # a list of info per formula
    marker <- NULL
    sname <- colnames(statedata)
    gtemp <- lapply(rhs2, function(x) grableft(x[[1]])) # not the options

    for (i in 1:nform) {
        # how many on the left?
        ftemp <- flist[[i]][1:2] # becomes a no-response formula
        temp2 <- attr(terms(ftemp), "term.labels")
        nmarker[i] <- length(temp2)
        marker <- c(marker, temp2)
        
        # what was the state?
        first <- gtemp[[i]][[1]] # the first word right of the ~
        if (is.call(first)) {
            # the user has something like "A(0:2)", which means that any state
            #  with statedata$A==0 will get assigned the first guassian peak,
            #  rows with statedata$A ==1 the second, and statedata$A==2
            #  the third one. States that don't match will get a density of 0.
            # This last is normally the death state, we don't need parameters
            #  for a "death" peak because the marker will never be measured for
            #  a death obs.
            jcol <- match(as.character(first[[1]]), sname)
            if (is.na(jcol)) stop("unrecognized state vector: ", deparse(first))
            first[[1]] <- as.name("c")
            temp <- eval(first)  # A(1:3) becomes the vector 1,2,3
            temp <- unique(temp[!is.na(temp)]) # users do weird things....
            stateinfo[[i]] <- list(sname= names(statedata)[jcol], levels=temp,
                               index=match(statedata[,jcol], temp, nomatch=0))
        } else {
            if (is.name(first)) jcol <- match(as.character(first), sname)    
            else if (is.character(first)) jcol <- match(first, sname)
            else stop("unrecognized state vector: ", deparse(first))

            if (is.na(jcol)) stop("unrecognized state vector: ", deparse(first))
            temp <- unique(statedata[,jcol])
            stateinfo[[i]] <-list(sname= names(statedata)[jcol], levels=temp,  
                              index=match(statedata[,jcol], temp[!is.na(temp)]))
        }

        # Further predictors?
        temp <- lapply(gtemp, function(x) {
            if (is.null(x[[2]])) NULL
            else {
                dummy <- ~x
                dummy[[2]] <- x[[2]]  # right hand side, less state and options
                attr(terms(dummy), "term.labels")
            }
        })
     
    list(marker= marker, statinfo= stateinfo, options=options, 
         mterm = temp, nmarker = nmarker)

# Why not expand the options list above so there is one element per marker?
#   The parsemarker2 routine, which builds tmap and cmap, needs to know
#   if there is a single sex effect for the 2 markers * 3 N states = 6
#   linear predictors (as in the above), or if the user had two separate
#   /common options, which would lead to 2 sex effects in the 6 linear
#   predictors. 
# Worse, there are actually 12 linear predictors, 6 means and 6 std, and
#   common has multiple meanings.

# The term map (tmap) will have one row per state, followed by one row
#  per covariate. Most often there are no additional covariates.
# There will be one column per (marker, parameter of the
#  distribution) pair, i.e., one per linear predictor. 
# The coefficient map cmap is an expansion of tmap.


parsemarker2 <- function(parse1, statedata, Terms, Xname, Xassign,
                         markerlevels) {
    nstate <- nrow(statedata)
    nform  <- length(parse1$nmarker)
    umarker <- unique(parse1$marker)
    nmarker <- length(umarker)

    # grab the distribution name and parameter names of the distribution
    ocheck <- function(x) {
        cmatch <- function(x) {
        if (is.call(x)) eval(x, env= env2)
        else if (length(x)==1) cdist[[as.character(x)]]
        else stop("unrecognized distribution")
        }
            
        # at present only dist, dist + common, or common + dist
        if (is.call(x)  && x[[1]] == as.name("+")) {
            t1 <- cmatch(x[[2]])
            if (!is.null(t1$error)) stop(t1$error)
            t2 <- cmatch(x[[3]])
            if (!is.null(t2$error)) stop(t2$error)

            if (t1$npar > 0 && t2$npar >0) 
                stop("more than one distribution given: ", deparse(x))
            else if (t1$npar >0) t1
            else if (t2$npar >0) t2
            else stop("unrecognized distribution: ", deparse(x))
        } else if (length(x) ==1) cmatch(x)
        else stop("unrecognized distribution: ", deparse(x))
    }
    odist <- lapply(parse1$opt, ocheck)

    # 
    #  First pass: each marker should point to only one state set, and only
    #   one distribution. Columns of tmap will be marker:dist:param.  Any
    #   'common' arguments affect values within a column, a gaussian for 
    #   instance always has two parameters.
    #  A marker with additive affects on states is weird: an equation with
    #   pib ~ A and another with pib ~ N; but someday someone will want it. 
    #   When that day comes we will think about allowing it.
    #  A marker might appear in 2 formula, e.g., a covariate applies to the
    #   mean but not the std, but it will still be parameters of the same dist.
    #  One formula can also have 2 markers.  
    statecol <- rep(parse1$statecol, parse1$nmarker)
    check1 <- sapply(umarker, function(x) {
        length(unique(statecol[parse1$marker==x]))
    })
    if (any(check1 >1)) stop("marker points to two states: ", umarker[check1>1])

    oindex <- rep(1:nform, parse1$nmarker) # the option list for each marker
    # first element in 'marker' points to oindex[1] element of options, etc
    tdist <- sapply(odist[oindex], function(x) x$dist)
    check2 <- sapply(umarker, function(x) {
        length(unique(tdist[parse1$marker==x]))
    })
    if (any(check2 >1)) stop("marker points to two distributions: ",
                            umarker[check2>1])
        
    # set of the columns for tmap and cmap, one for each marker:(parm of dist)
    #  pair. We've already checked for conflicts, so can grab one
    # Indexing below even confuses the author
    #  j <- match(umarker, parse1$marker) = first time each unique marker 
    #    appears in the marker list
    #  oindex[j] = the options formula that goes with it, which is the index of
    #    odist as well
    # 
    uindex <- oindex[match(umarker, parse1$marker)] # first formula for a marker
    dpar  <- lapply(odist[uindex], function(x) x$param)
    ddist <- sapply(odist[uindex], function(x) x$dist)
    npar <- sapply(dpar, length)
    coldata <- data.frame(marker= rep(umarker, npar),
                          dist  = rep(ddist, npar),
                          param = unlist(dpar))
    cname <- with(coldata, paste(marker, param, sep=':'))
    # column names are marker:param, row names are states followed by variables
    # row names for tmap
    rname <- c(statedata$state, unlist(parse1$mterm))
    tmap <- matrix(0L, nrow=length(rname), ncol=length(cname))

    # Walk through the formulas, from first to last. Succeeding ones
    #  trump prior ones
    # dmap is a set of unique integers, so I don't end up reusing an index
    dmap <- matrix(1:length(tmap), nrow=nrow(tmap), ncol=ncol(tmap))
    
    for (i in 1:nform) {
        stemp <- statedata[, parse1$statecol[i]]
        index1 <- m     
    }
}                 
                                     
                        
