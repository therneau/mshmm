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
    
    rhs2 <- lapply(flist, function(x) rightslash(x[[3]])
    options <- lapply(rhs2, function(x) x[[2]]) # the list of options
    
    # Common formula:  log(pib) ~ A / gaussian
    # Harder case: pzmemory + pzglobal ~ N + practice + sex / gaussian + common
    # 1. Return "~ pib" as the formula for case 1, for case 2 return 
    #    "pzmemory + pzglobal ~ practice + sex"; both work for the reformulate()
    #  code in hmm.
    # 2. For case 2, the returned marker component will have a separate element
    #   for each marker, statecol and options will also have one element
    #   per marker.
    # 3. ptplot(~N + practice + sex) shows that we can't just use a simple
    # 'grab left term of first + call" strategy.  Hence the grableft() function.
    #  It returns a list of length 2, first element will be "N", second
    #  'practice + sex'.  Note that we can't use terms + reformulate, if the
    #  right portion had an interaction, it would get lost.
    grableft <- function(x) {
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
    nmarker <- statecol <- integer(nform)
    marker <- NULL
    sname <- colnames(statedata)
    for (i in 1:nform) {
        # how many on the left?
        temp <- flist[[i]][1:2] # becomes a no-response formula
        temp2 <- attr(terms(temp), "term.labels")
        nmarker[i] <- length(temp2)
        marker <- c(marker, temp2)
        
        rhs <- rhs2[[i]][[1]] # current right hand side, without options
        if (length(rhs)==1) { # usual case, only a state name on the right
            first <- rhs[[1]]
            if (is.name(first)) jcol <- match(as.character(first), sname)    
            else if (is.character(first)) jcol <- match(first, sname)
            else stop("unrecognized state vector: ", deparse(first))

            if (is.na(jcol)) stop("unrecognized state vector: ", deparse(first))
            statecol[i] <- jcol
            flist[[i]][[3]] <- NULL  # left hand side is now right hand
        } else {
            gtemp <- grableft(rhs)
            first <- gtemp[[1]]
            if (is.name(first)) jcol <- match(as.character(first, sname)) 
            else if (is.character(first)) jcol <- match(first, sname)
            else stop("unrecognized state vector: ", deparse(first))

            if (is.na(jcol)) stop("unrecognized state vector: ", deparse(first))
            statecol[i] <- j
            flist[[i]][[3]] <- gtemp[[2]]
        }
    }
      
    temp <- sapply(rhs, function(x) length(x[[1]]))
    rval <- list(formula = flist, statecol= statecol, marker = marker, 
                 options=options, nmarker = nleft)
    if (any(temp>0)) {
         rhs= lapply(rhs2, function(x) {
             if (length(x[[1]]) >0) {
                 dummy <- formula[[1]]  # a convenient 3 element formula
                 dummy[[2]] <- x
                 attr(terms(dummy[1:2], "term.labels"))
             } else NULL
         })   
         rval$mvar <- rhs
    }
    rval
}

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

# setup for response functions
env2 <- new.env(parent.frame=parent(2))
assign(gaussian, env=env2,
       value= function(param=c("mean", "std")) {
           if (missing(param)) # assume both
               list(dist="gaussian", param=c("mean", "std"), npar=2)
           else {
               param= match.call(param),
               list(dist="gaussian", param= param, npar=2)
           }
       })  

assign(gamma, env=env2,
       value= function(param=c("mean", "std")) {
           if (missing(param)) # assume both
               list(dist="gamma", param=c("mean", "std"), npar=2)
           else {
               param= match.call(param),
               list(dist="gamma", param= param, npar=2)
           }
       })  

assign(multinomial, env=env2,
       value= function(pattern) {
           nstate <- nrow(pattern)
           npar = sum(pattern >0)
           param= paste0("p", 1:npar)
       })

cdist <- list(gaussian= list(dist="gaussian", param=c("mean", "std"), npar=2),
              gamma = list(dist="gamma", param=c("mean", "std"), npar=2)
              common= list(dist="null", npar=0)
              )

parsemarker2 <- function(parse1, statedata, Terms, Xname, Xassign) {
    nstate <- nrow(statedata)
    nform  <- length(parse1$nmarker)
    umarker <- unique(parse1$marker)
    nmarker <- length(umarker)

    # grab the distribution name and parameter names of the distribution
    ocheck <- function(x) {
        cmatch <- function(x) {
        if (is.function(x)) eval(x, env= env2)
        else if (length(x)==1) cdist[as.character(x)]
        else stop("unrecognized distribution")
        }
            
        # at present only dist, dist + common, or common + dist
        if (is.call(x)  && x[[2]] == as.name("+")) {
            t1 <- cdist(x[[2]])
            t2 <- cdist(x[[3]])
            if (t1$npar==0) t1
            else if (t2$npar==0) t2
            else stop("more than one distribution given: ", deparse(x))
        } else if (length(x) ==1) cdist(x)
        else stop("unrecognized distribution: ", deparse(x))
    }
    
    odist <- lapply(parse1$opt, oparse)

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
        length(unique(statecol[marker==x]))
    })
    if (any(check1 >1)) stop("marker points to two states: ", umarker[check1>1])

    oindex <- rep(1:nform, parse1$nmarker) # the option list for each marker
    # first element in 'marker' points to oindex[1] element of options, etc
    tdist <- sapply(odist, function(x) x$dist)
    check2 <- sapply(umarker, function(x) {
        length(unique(tdist[marker==x]))
    })
    if (any(check2 >1)) stop("marker points to two distributions: ",
                            umarkder[check2>1])
        
    # set of the columns for tmap and cmap
    uindex <- oindex[match(umarker, marker)] # first formula for a marker
    dpar <- lapply(odist[uindex], function(x) x$param)
    ddist <- tdist[unidex]
    npar <- sapply(dpar, length)
    coldata <- data.frame(marker= rep(umarker, each=npar),
                          dist  = rep(ddist, each=npar),
                          param = unlist(dpar))
    cname <- with(coldata, paste(marker, dist, param, collapse=':'))
    
    # And the row names
    rname <- c(statedata$state, unlist(parse1$mvar))
    
}
                   
                                     
                        
