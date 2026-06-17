#
# Various methods for probability in state estimates, created by
#  predict.cmsh.  They act a lot like survfit objects in survival
#
dim.pstate.cmsh <- function(x) {
    # to a user, the pstate.cmsh object will always have a data dimension of
    #  at least 1, and state dimension of at least 1
    c(data= nrow(x$data), states=length(x$states))
}
    
# subscripting is most often used to plot a subset of curves
#  the cumulative hazard is lost if we subset on states  
# if dim(x)= (1,1) then, depending on whether prior [] calls used drop or not,
#  x$pstate could be a vector of length n, where n= number of time points, 
#  a matrix of dimension (n,1), or an array of dimension(n,1,1)
# similar issues for dim of (1,p) and (k,1). The goal is that the user need
#  not know or worry about it
"[.pstate.cmsh" <- function(x, ..., drop=TRUE) {
    if (!inherits(x, "pstate.cmsh")) 
        stop("[.pstate.cmsh called on non-pstate object")

    ndots <- ...length()  # not avail in R 3.4
    if (ndots==0) return(x) # return the object as is
    if (ndots>0 && !missing(..1)) i <- ..1 else i <- NULL
    if (ndots==2 && !missing(..2)) j <- ..2 else j <- NULL
    if (ndots >2) stop("incorrect number of subscripts for pstate.cmsh object")
    
    if (is.null(i) & is.null(j)) return(x)

    # if the pstate object has data or state dimension of 1, we allow a
    #  single subscript.  Two subscripts is always legal
    dd <- dim(x)
    if (ndots==1) {
        if (all(dd>1)) 
            stop("incorrect number of subscripts for pstate.cmsh object")

        x$pstate <- (as.matrix(x$pstate))[,i, drop=drop]
        if (dd[1]>1) {
            x$data <- x$data[i,,drop=FALSE] # retain it as a data frame
            x$cumhaz <- x$cumhaz[,i,, drop=drop]
        }
        else {
            x$states <- x$states[i]
            x$cumhaz <- NULL
        }
    } else { # two subscripts, one of which might be missing
        if (!is.null(j)) { #subscript on the states
            temp <- x$states[j]  # an invalid subscript will fail here
            if (length(temp)== length(x$states)) return(x) #kept em all
            else if (length(temp) ==0) return (NULL) # tossed em all
            else x$states <- temp

            # dd[2] >1, or I won't get here
            x$cumhaz <- NULL
            if (is.null(i)) { # only on states
                if (length(dim(x$pstate))==3)
                    x$pstate <- x$pstate[,,j, drop=drop]
                else x$pstate <- (as.matrix(x$pstate))[,j, drop=drop]
            } else {
                # both i and j
                x$data <- x$data[i,, drop=FALSE] # retain it as a data frame
                if (length(dim(x$pstate))==3) # if dd[1]==1 it might be a matrix
                    x$pstate <- x$pstate[,i,j, drop=drop]
                else x$pstate <- x$pstate[,j, drop=drop]
            }
        } else {
            # j is null, only i
            browser()
            x$data <- x$data[i,, drop=FALSE]
            x$cumhaz <- x$cumhaz[,i,, drop=drop]
            if (length(dim(x$pstate))==3) 
                 x$pstate <- x$pstate[,i,,drop=drop]
            else x$pstate <- (as.matrix(x$pstate))[,i, drop=drop]
        }
    }       
    x
}
   
# arguments are modeled on those in survival::plot.survfit 
plot.pstate.cmsh <- function(x, cumhaz=FALSE, cumprob=FALSE, 
                             xlab, ylab, 
                             xscale=1, yscale=1, xmax, xlim, ylim, ...) {
    # we don't yet have CI for the curves
    if (!inherits(x, "pstate.cmsh")) stop("x is not a pstate.cmsh object")
    xtime <- if(is.vector(x$time)) x$time else x$time[,1]
    ntime <- length(xtime)
    if (missing(xlab)) xlab <- x$tname[1]
    nstate <- dim(x)[2]

    if (!missing(cumhaz) && !(is.logical(cumhaz) && all(!cumhaz))) {
        # user entered a cumhaz argument, and it's not all "FALSE"
        if (missing(ylab)) ylab <- "Cumulative hazard"
        if (is.null(x$cumhaz)) stop("x does not have a cumhaz component")
        cdim <- dim(x$cumhaz)
        # a check for users who created their own pstate.cmsh object
        if (cdim[1] != ntime) stop("invalid cumhaz component")
        nhaz <- cdim[length(cdim)] # number of unique hazards, within data dim
        if (is.logical(cumhaz)) {
            if (length(cumhaz) ==1) yhat <- matrix(x$cumhaz, nrow=ntime)
            else if (length(cumhaz) != nhaz)
                stop("wrong length for cumhaz argument")
            else {
                if (cdim==3) yhat <- matrix(x$cumhaz[,,cumhaz], nrow= ntime) 
                else yhat <- x$cumhaz[,cumhaz, drop=FALSE]
            }
        } else if(is.numeric(cumhaz)) {
            if (!all(cumhaz== floor(cumhaz))) 
                stop("cumhaz argument is not integer")
            if (any(cumhaz < 1 | cumhaz > nhaz)) 
                stop("cumhaz subscript out of range") 
            if (cdim==3) yhat <- matrix(x$cumhaz[,,cumhaz], nrow= ntime)
            else yhat <- x$cumhaz[,cumhaz, drop=FALSE]
        } else stop("invalid cumhaz argument")
    }
    else if (!missing(cumprob) && !(is.logical(cumprob) && all(!cumprob))) {
        # user has a cumprob argument and it is not "FALSE"
        if (missing(ylab)) ylab <- "Cumulative probability in state"
        if (is.logical(cumprob)) {
            if (length(cumprob) ==1) jj <- 1:nstate
            else if (length(cumprob) != nstate)
                stop("wrong length for cumprob argument")
            else jj <- (1:nstate)[cumprob]
        } else if (is.numeric(cumprob)) {
            if (!all(cumprob== floor(cumprob))) 
                stop("cumprob argument is not integer")
            if (any(cumprob < 1 | cumprob > nstate)) 
                stop("cumprob subscript out of range") 
            jj <- cumprob
        } else stop("invalid cumprob argument")
        
        # Now create the cumulative curves, jj is the order
        # if nstate=1 then jj will be 1
        if (nstate==1) yhat <- matrix(x$pstate, nrow=ntime) # nothing to do
        else {
            temp <- drop(x$pstate)  # just in case the data dimension is 1
            if (length(jj) ==1) {  # a use could say "cumprob=2"
                if (is.matrix(temp)) yhat <- temp[,jj, drop= FALSE]
                else yhat <- temp[,,jj]
            } else {
                if (is.matrix(temp)) 
                    yhat <- t(apply(temp[,jj, drop=FALSE], 1, cumsum))
                else {
                    temp <- apply(x$pstate[,,jj, drop=FALSE], 1:2, cumsum)
                    # the state dimension is now first
                    yhat <- matrix(aperm(temp, c(2,3,1)), nrow= ntime)
                }
            }
        }
    } else {
        # neither a cumprob or cumhaz argument
        yhat <- matrix(x$pstate, nrow= ntime) # fold it into a matrix
        if (missing(ylab)) ylab <- "Probability in state"
    }

    # check consistency of range arguments
    if (!missing(xlim)) {
        if (!(is.numeric(xlim) && length(xlim)==2)) {
            warning("invalid xlim, value ignored")
            xlim <- NULL
        }
    } else xlim <- NULL
    if (!missing(ylim)) {
        if (!(is.numeric(ylim) && length(ylim)==2)) {
            warning("invalid ylim, value ignored")
            ylim <- NULL
        }
    } else ylim <- NULL
    
    if (!missing(xmax)) {
        if (!is.null(xlim)) {
            warning("cannot have both xlim and xmax arguments, xmax ignored")
            options(plot.survfit = NULL)
        }
        else {
            if (!is.numeric(xmax) || length(xmax) !=1)
                stop("invalid xmax value")
            keep <- (xtime <= xmax)
            xtime <- xtime[keep]
            yhat <- yhat[keep,, drop=FALSE]
            options(plot.survfit = list(xmax=xmax))
        }
    } else options(plot.survfit = NULL)

    #
    # Draw the basic box
    #
    plot(range(xtime, finite=TRUE, na.rm=TRUE)/xscale, 
         range(yhat, finite=TRUE, na.rm=TRUE)*yscale, 
         type='n', xlab=xlab, ylab=ylab, ...)
    if(yscale != 1) par(usr =par("usr")/c(1, 1, yscale, yscale))   
    if (xscale !=1) par(usr =par("usr")*c(xscale, xscale, 1, 1))   
    # The use of [[par(usr)]] just above is a bit sneaky.  I want the
    # lines and points routines to be able to add to the plot, *without*
    # passing them a global parameter that determines the y-scale or forcing
    # the user to repeat it.

    matlines(xtime, yhat, type='l', ...)
}

lines.pstate.cmsh <- function(x, cumhaz=FALSE, cumprob=FALSE, ...) { 
    # we don't yet have CI for the curves
    if (!inherits(x, "pstate.cmsh")) stop("x is not a pstate.cmsh object")
    xtime <- if(is.vector(x$time)) x$time else x$time[,1]
    ntime <- length(xtime)
    nstate <- dim(x)[2]

    if (!missing(cumhaz) && !(is.logical(cumhaz) && all(!cumhaz))) {
        # user entered a cumhaz argument, and it's not "FALSE"
        if (is.null(x$cumhaz)) stop("x does not have a cumhaz component")
        cdim <- dim(x$cumhaz)
        # a check for users who created their own pstate.cmsh object
        if (cdim[1] != ntime) stop("invalid cumhaz component")
        nhaz <- cdim[length(cdim)] # number of unique hazards, within data dim
        if (is.logical(cumhaz)) {
            if (length(cumhaz) ==1) yhat <- matrix(x$cumhaz, nrow=ntime)
            else if (length(cumhaz) != nhaz)
                stop("wrong length for cumhaz argument")
            else {
                if (cdim==3) yhat <- matrix(x$cumhaz[,,cumhaz], nrow= ntime) 
                else yhat <- x$cumhaz[,cumhaz, drop=FALSE]
            }
        } else if(is.numeric(cumhaz)) {
            if (!all(cumhaz== floor(cumhaz))) 
                stop("cumhaz argument is not integer")
            if (any(cumhaz < 1 | cumhaz > nhaz)) 
                stop("cumhaz subscript out of range") 
            if (cdim==3) yhat <- matrix(x$cumhaz[,,cumhaz], nrow= ntime)
            else yhat <- x$cumhaz[,cumhaz, drop=FALSE]
        } else stop("invalid cumhaz argument")
    }
    else if (!missing(cumprob) && !(is.logical(cumprob) && all(!cumprob))) {
        # user has a cumprob argument and it is not "FALSE"
        if (is.logical(cumprob)) {
            if (length(cumprob) ==1) jj <- 1:nstate
            else if (length(cumprob) != nstate)
                stop("wrong length for cumprob argument")
            else jj <- (1:nstate)[cumprob]
        } else if (is.numeric(cumprob)) {
            if (!all(cumprob== floor(cumprob))) 
                stop("cumprob argument is not integer")
            if (any(cumprob < 1 | cumprob > nstate)) 
                stop("cumprob subscript out of range") 
            jj <- cumprob
        } else stop("invalid cumprob argument")
        
        # Now create the cumulative curves, jj is the order
        # if nstate=1 then jj will be 1
        if (nstate==1) yhat <- matrix(x$pstate, nrow=ntime) # nothing to do
        else {
            temp <- drop(x$pstate)  # just in case the data dimension is 1
            if (length(jj) ==1) {  # a use could say "cumprob=2"
                if (is.matrix(temp)) yhat <- temp[,jj, drop= FALSE]
                else yhat <- temp[,,jj]
            } else {
                if (is.matrix(temp)) 
                    yhat <- t(apply(temp[,jj, drop=FALSE], 1, cumsum))
                else {
                    temp <- apply(x$pstate[,,jj, drop=FALSE], 1:2, cumsum)
                    # the state dimension is now first
                    yhat <- matrix(aperm(temp, c(2,3,1)), nrow= ntime)
                }
            }
        }
    } else {
        # neither a cumprob or cumhaz argument
        yhat <- matrix(x$pstate, nrow= ntime) # fold it into a matrix
    }


    matlines(xtime, yhat, ...)
}

points.pstate.cmsh <- function(x, cumhaz=FALSE, cumprob=FALSE, ...) {
    # we don't yet have CI for the curves
    if (!inherits(x, "pstate.cmsh")) stop("x is not a pstate.cmsh object")
    xtime <- if(is.vector(x$time)) x$time else x$time[,1]
    ntime <- length(xtime)
    nstate <- dim(x)[2]

    if (!missing(cumhaz) && !(is.logical(cumhaz) && all(!cumhaz))) {
        # user entered a cumhaz argument, and it's not "FALSE"
        if (is.null(x$cumhaz)) stop("x does not have a cumhaz component")
        cdim <- dim(x$cumhaz)
        # a check for users who created their own pstate.cmsh object
        if (cdim[1] != ntime) stop("invalid cumhaz component")
        nhaz <- cdim[length(cdim)] # number of unique hazards, within data dim
        if (is.logical(cumhaz)) {
            if (length(cumhaz) ==1) yhat <- matrix(x$cumhaz, nrow=ntime)
            else if (length(cumhaz) != nhaz)
                stop("wrong length for cumhaz argument")
            else {
                if (cdim==3) yhat <- matrix(x$cumhaz[,,cumhaz], nrow= ntime) 
                else yhat <- x$cumhaz[,cumhaz, drop=FALSE]
            }
        } else if(is.numeric(cumhaz)) {
            if (!all(cumhaz== floor(cumhaz))) 
                stop("cumhaz argument is not integer")
            if (any(cumhaz < 1 | cumhaz > nhaz)) 
                stop("cumhaz subscript out of range") 
            if (cdim==3) yhat <- matrix(x$cumhaz[,,cumhaz], nrow= ntime)
            else yhat <- x$cumhaz[,cumhaz, drop=FALSE]
        } else stop("invalid cumhaz argument")
    }
    else if (!missing(cumprob) && !(is.logical(cumprob) && all(!cumprob))) {
        # user has a cumprob argument and it is not "FALSE"
        if (is.logical(cumprob)) {
            if (length(cumprob) ==1) jj <- 1:nstate
            else if (length(cumprob) != nstate)
                stop("wrong length for cumprob argument")
            else jj <- (1:nstate)[cumprob]
        } else if (is.numeric(cumprob)) {
            if (!all(cumprob== floor(cumprob))) 
                stop("cumprob argument is not integer")
            if (any(cumprob < 1 | cumprob > nstate)) 
                stop("cumprob subscript out of range") 
            jj <- cumprob
        } else stop("invalid cumprob argument")
        
        # Now create the cumulative curves, jj is the order
        # if nstate=1 then jj will be 1
        if (nstate==1) yhat <- matrix(x$pstate, nrow=ntime) # nothing to do
        else {
            temp <- drop(x$pstate)  # just in case the data dimension is 1
            if (length(jj) ==1) {  # a use could say "cumprob=2"
                if (is.matrix(temp)) yhat <- temp[,jj, drop= FALSE]
                else yhat <- temp[,,jj]
            } else {
                if (is.matrix(temp)) 
                    yhat <- t(apply(temp[,jj, drop=FALSE], 1, cumsum))
                else {
                    temp <- apply(x$pstate[,,jj, drop=FALSE], 1:2, cumsum)
                    # the state dimension is now first
                    yhat <- matrix(aperm(temp, c(2,3,1)), nrow= ntime)
                }
            }
        }
    } else {
        # neither a cumprob or cumhaz argument
        yhat <- matrix(x$pstate, nrow= ntime) # fold it into a matrix
    }

    matpoints(xtime, yhat, ...)
}
