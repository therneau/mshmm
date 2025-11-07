# Parse the list of markers for latent states.  This will be a list of
#  formulas, each left hand side a marker, predicted by a state and possibly
#  other covariates, with the options after a /
# Separate out the states, markers, and options.
#
# Common formula:  A:log(pib) ~ 1 / gaussian
# Harder case: 
#    N(1:3):pzmemory + N(1:3):pzglobal ~ 1 / gaussian,
#    N(1:3):pzmemory + N(1:3):pzglobal ~ practice /
#                     gaussian(param="mean") + common
# The state predicts the marker, not vise-versa, the order of state:marker
#  is intentional. 
# By default each state/marker/parameter of the distribution will be a
#  separate linear predictor. So in the first, if the "A" column of 
#  statedata had values of "neg", "pos" and NA (for death), there will
#  be 4 linear predictors Aneg:log(pib).mean Aneg:log(pib).std, etc.
# The harder case has 12 linear predictors, there are 12 separate coefs
#  for the means. Assume practice is a factor for visits 1, 2:3, and 4+,
#  then there will be 2 practice coefficients, which add to each of
#  the 6 linear predictors for the mean (3 states by 2 markers).
#
parsemarker1 <- function(flist, statedata) {
    if (inherits(flist, "formula")) flist <- list(flist) # only one equation
    if (any(sapply(flist, function(x) !inherits(x, "formula"))))
        stop("an element of the marker list is not a formula")
    if (any(sapply(flist, length) != 3))
        stop("all formulas must have a left and right side")
    
    rhs2 <- lapply(flist, function(x) rightslash(x[[3]]))
    options <- lapply(rhs2, function(x) x[[2]]) # the list of options
    
    nform <- length(flist)    # number of formulas
    lhs <- lapply(flist, function(x) markerpair(x[[2]], statedata))

    # lhs[[k]] will have a stateinfo, marker, stateinfo, marker, ... pairs for
    #  however many markers were in this formula.  I could make multiple
    #  markers on one line illegal but I have one potential use case.
    # remake this to something more simple
    #   nmarker = number of markers in the formula, so I can match with options
    #     which has one element per formula
    #   marker = character vector of markers (might contain duplicates)
    #   stateinfo = list with one element per marker
    nmarker <- sapply(lhs, length)/2
    temp <- unlist(lhs, recursive=FALSE)
    marker <- unname(unlist(temp[names(temp)== "marker"]))
    stateinfo <- unname(temp[names(temp)== "stateinfo"])

    # Further covariates and intercepts.  We will later need to distinguish
    #  between ~ x and ~ 1+x, so pass the equation forward, make it
    #  a 1 sided formula
    mterm <- lapply(rhs2, function(x) {
        dummy <- flist[[1]][1:2]  # keep the proper environment
        dummy[[2]] <- x[[1]]
        dummy
    })      
    
    # the parent needs the markers and the extra variables (mterm) to
    #  build the data frame, stateinfo is passed to parsemarker2
    list(marker=marker, mterm=mterm, options=options, nmarker=nmarker,
         stateinfo= stateinfo)
}

# expand the set of state:marker implied by a left hand side.
markerpair <- function(x, statemap) {
    nstate <- nrow(statemap)
    sname <-  names(statemap)
    if (x[[1]]== as.name("*") || x[[1]]== as.name("/") || 
        x[[1]]== as.name("-") || x[[1]]== as.name("%in%"))
        stop("invalid operator on the left of a marker formula")
    else if (is.call(x) && x[[1]]== as.name("+")) 
        c(markerpair(x[[2]], statemap), 
             markerpair(x[[3]], statemap)) #concatonate two lists
    else if (is.call(x) && x[[1]]== as.name(":")) {
        state <- x[[2]]
        marker<- x[[3]]
        # process the left hand side of a :
        if (is.call(state)) {
            # the user has something like "A(0:2)", which means that any state
            #  with statedata$A==0 will get assigned the first guassian peak,
            #  rows with statedata$A ==1 the second, and statedata$A==2
            #  the third one. States that don't match will get a density of 0.
            # This last case is often the death state, we don't need parameters
            #  for a "death" peak because the marker will never be measured for
            #  a death obs.
            jcol <- match(as.character(state[[1]]), sname)
            if (is.na(jcol)) {
                # a numeric vector like c(1,3,5) is also legal, but very odd
                # if there were 5 states it would mean that 2 and 4 map to
                # nothing
                index <- eval(state)
                if (is.numeric(index) && index== as.integer(index) &&
                    all(index>0) && all(index <= nstate)) {
                    temp <- rep(0L, nstate)
                    temp[index] <- seq(along=length(index))
                    stateinfo <- list(sname="state", 
                                      levels= statedata[1,index],
                                      index= temp)
                } else stop("unrecognized state vector: ", deparse(state))
            } else {
                state[[1]] <- as.name("c")
                temp <- eval(state)  # A(1:3) becomes the vector 1,2,3
                temp <- unique(temp[!is.na(temp)]) # users do weird things....
                stateinfo <- list(sname= names(statedata)[jcol], levels=temp,
                               index=match(statedata[,jcol], temp, nomatch=0))
            }
        } else {
            if (is.name(state)) jcol <- match(as.character(state), sname)    
            else if (is.character(state)) jcol <- match(state, sname)
            else stop("unrecognized state vector: ", deparse(state))

            if (is.na(jcol)) stop("unrecognized state vector: ", deparse(state))
            temp <- statedata[,jcol]
            temp <- unique(temp[!is.na(temp)]) # an NA matches nothing
            stateinfo <-list(sname= names(statedata)[jcol], levels=temp,  
                              index=match(statedata[,jcol], temp, nomatch=0))
        }

        # The right hand side must be a single marker variable
        # However, log(pib) is certainly legal
        if (class(marker)== "(") 
            stop("a parenthesised list of markers is not supported")
        list(stateinfo= stateinfo, marker= deparse(marker))
    } else stop("unrecognized portion of a marker formula", deparse(x))
}

# Why not expand the options list above so there is one element per marker?
#   We will need to know if a common option applies to both. 
#   For the 'hard' equation at the top of this file, 2 markers * 3 N states = 6
#   linear predictors for the mean, all get a single icvol effect.
#
# The term map (tmap) will have one row per state, followed by one row
#  per covariate. Most often there are no additional covariates.
# There will be one column per (marker, parameter of the
#  distribution) pair, i.e., one per linear predictor. 
# The coefficient map cmap is an expansion of tmap.

parsemarker2 <- function(parse1, statedata, Terms, Xname, Xassign,
                         markerlevel) {
    nstate <- nrow(statedata)
    nform  <- length(parse1$nmarker)  # number of formulas = number of options
    umarker <- unique(parse1$marker)
    nmarker <- length(umarker)

    # separate out the three parts of an options list
    # the hmm.dist vector has the list of legal distributions, see response.R
    osplit <- function(x) {
        # look at a single part of the formula
        if (is.call(x)) {
            if (x[[1]]== as.name("+")){
                # we want a single list, not a list of lists, hence c()
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
    #   one distribution. Columns of tmap will be dist:marker.param.  Using
    #   'common' across parameters would be odd, but is legal; a gaussian for 
    #   instance always has two parameters.
    #  A marker with additive affects on states is weird: an equation with
    #   pib ~ A and another with pib ~ N; but someday someone will want it. 
    #   When that day comes we will think about allowing it.
    #  A marker might appear in 2 formula, e.g., a covariate applies to the
    #   mean but not the std, but it will still be parameters of the same dist.
    #  One formula can also have 2 markers.  
    statecol <- sapply(parse1$stateinfo, function(x) x$sname)
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

    # Walk through the markers from first to last, and call the setup
    #  function of the distribution used for each.  
    # E.g. gaussian() is the setup function a gaussian peak, and it will
    #  return a response function along with ancillary information.
    mindex <- match(parse1$marker, umarker)  # markerlevel will be in this order
    rlist <- lapply(seq(along=mindex), function(i) {
        arglist <- list(stateinfo= parse1$stateinfo[[i]],
                        markerlevel = markerlevel[[mindex[i]]])
        tdist <- odist[[oindex[i]]] #options for this marker
        if (!is.null(tdist$call)) {
            # add on user arguments
            temp <- tdist$call
            temp[[1]] <- as.name("list")
            arglist <- c(arglist, eval(temp))
        }
        do.call(tdist$dist, arglist)
    })    
    
    # Create col names for tmap and cmap, one per linear predictor
    temp <- lapply(1:nform, function(i)
        paste(parse1$marker[i], rlist[[i]]$pname, sep=':'))
    lpname <- unique(unlist(temp))
    numlp <- length(lpname)
    cmap <- matrix(0L, nrow=length(Xname), ncol= numlp,
                   dimnames=list(Xname, lpname))

    if (any(sapply(parse1$mterm, function(x) x != ~1))) {
        # at least one marker formula is more than "~1", we will need tmap
        tmap <- matrix(0L, nrow=length(attr(Terms, "term.labels")), ncol=numlp)
        intercept <- sapply(parse1$mterm, hasintercept)
    } else tmap <- NULL
    browser()
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
                                                           
# Run down the formula parse tree and see if there is an explicit
#  intercept term, i.e. a "1".  An intercept won't be part of an interaction
#  so no need to chase *, : or (.
hasintercept <- function(x) {
    if (class(x) == "formula") hasintercept(x[[2]])
    else if (is.numeric(x) && x==1) TRUE
    else if(is.call(x) && x[[1]] == as.name("+")) 
       (hasintercept(x[[2]]) || hasintercept(x[[3]]))
    else FALSE
}
