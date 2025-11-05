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
    # Harder case: 
    #    pzmemory + pzglobal ~ N(1:2) + practice + sex /gaussian + common
    # The state information is required to be the first thing on the rhs
    #  of the formula.  However,  ptplot( ~N(1:2) + practice + sex) shows
    #  that we can't simply grab element the first thing off the parse tree:
    #  N is further down than that.
    # The following function splits the non-options part of rhs2 into the
    #  the leftmost part of the formula from our point of view (state info)
    #  and the rest. In practice, for most formulas the "rest" will be NULL.
    grableft <- function(x) {
        if (is.call(x) && x[[1]]== as.name('+')) {
            temp <- grableft(x[[2]])
            if (is.null(temp[[2]])) list(state=temp[[1]], covar=x[[3]])
            else {
                x[[2]] <- temp[[2]]
                list(state= temp[[1]], covar= x)
            }
        }
        else list(state=x, covar=NULL)
    }
    
    nform <- length(flist)    # number of formulas
    nmarker <- integer(nform) # number of markers in each formula
    stateinfo <- vector("list", nform)  # a list of info per formula
    marker <- NULL            # vector of all the markers
    sname <- colnames(statedata)
    gtemp <- lapply(rhs2, function(x) grableft(x[[1]])) # not the options

    # now go through the formlas one by one
    for (i in 1:nform) {
        # how many on the left?
        ftemp <- flist[[i]][1:2] # becomes a no-response formula
        temp2 <- attr(terms(ftemp), "term.labels")
        nmarker[i] <- length(temp2)
        marker <- c(marker, temp2)
        # what was the state?
        first <- gtemp[[i]][[1]] # the first term right of the ~
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
            temp <- statedata[,jcol]
            temp <- unique(temp[!is.na(temp)]) # an NA matches nothing
            stateinfo[[i]] <-list(sname= names(statedata)[jcol], levels=temp,  
                              index=match(statedata[,jcol], temp, nomatch=0))
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
    } 
    list(marker= marker, stateinfo= stateinfo, options=options, 
         mterm = temp, nmarker = nmarker)
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

parsemarker2 <- function(parse1, statedata, Terms, Xname, Xassign,
                         markerlevel) {
    nstate <- nrow(statedata)
    nform  <- length(parse1$nmarker)
    umarker <- unique(parse1$marker)
    nmarker <- length(umarker)

    # separate out the three parts of an options list
    # the hmm.dist vector has the list of legal distributions, see response.R
    osplit <- function(x) {
        # look at a single part of the formula
        if (is.call(x)) {
            if (x[[1]]== as.name("+")){
                # we want a single list, not a list of lists
                c(osplit(x[[2]]), osplit(x[[3]]))
            } else if (x[[1]]== as.name("init")) {
                x[[1]] <- as.name("c")
                list(init = eval(x))
            } else if (as.character(x[[1]]) %in% hmm.dist) 
                list(dist = as.character(x[[1]]), call=x)
            else stop("unrecognized function after / in a marker specification")
        } else {
            if (as.character(x) == "common") list(common=TRUE)
            else if (as.character(x) %in% hmm.dist) list(dist= as.character(x))
            else stop("unrecognized keyword after / in a marker specification")
        }
    }
    odist <- lapply(parse1$options, osplit)
            
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

    # Walk through the formulas from first to last, and call the setup
    #  function of the distribution used for each.
    # E.g. gaussian() is the setup function a gaussian peak, and it will
    #  return a response function along with ancillary information.
    rlist <- lapply(1:nform, function(i) {
        arglist <- list(stateinfo= parse1$stateinfo[[i]],
                        markerlevel = markerlevel[[i]])
        if (!is.null(odist[[i]]$call)) {
            # add on user arguments
            temp <- odist[[i]]$call
            temp[[1]] <- as.name("list")
            arglist <- c(arglist, eval(temp))
        }
        do.call(odist[[i]]$dist, arglist)
    })    
    browser()    
    
    # The total number of linear predictors = col names for tmap and cmap
    temp <- lapply(1:nform, function(i)
        paste(parse1$marker[i], rlist[[oindex[i]]]$pname, sep=':'))
    lpname <- unique(unlist(temp))
    numlp <- length(lpname)
    cmap <- matrix(0L, nrow=length(Xname), ncol= numlp,
                   dimnames=list(Xname, lpname))

    if (length(unlist(parse1$mterm)) >0) {
        # there are extra variables, create tmap
        tmap <- matrix(0L, nrow=length(attr(Terms, "term.labels")), ncol=numlp)
    } else tmap <- NULL

    # Walk through the formulas one at a time and build cmap,
    # dealing with "common"
    #
    
        
        
   uindex <- oindex[match(umarker, parse1$marker)] # first formula for a marker
    ddist <- sapply(odist[uindex], function(x) x$dist)

    # Walk through the formulas, from first to last. Succeeding ones
    #  trump prior ones
    # dmap is a set of unique integers, so I don't end up reusing an index
    dmap <- matrix(1:length(tmap), nrow=nrow(tmap), ncol=ncol(tmap))
    
    for (i in 1:nform) {
        stemp <- statedata[, parse1$statecol[i]]
        index1 <- m     
    }
}             
                                     
                        
