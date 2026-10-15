# Marker distribution functions
# hmm(...  marker= list(log(pib) ~ gaussian, 
#                             dx ~ multinomial(pattern=pmat)),..
# The parsemarker2 routine will create a call of
#           gaussian(stateinfo, ...) where 'gaussian' is the function below that
#              sets and returns a gaussian response function + other info
# likewise  multinomial(stateinfo, nclass, pattern)
#
# hmm.dist is used by parsemarker to check for legal names
hmm.dist <- c("gaussian", "logistic", "beta", "gamma", "discrete")

# stateinfo: name and level information that is used to create labels
#  for the linear predictors, and to match linear predictors to states.
#  A list with one element per marker, with elements of
#    name: used for creating a label
#    levels= a vector of labels, e.g. if of length 2 there would be 2 gaussian
#       densities
#    index: a vector of length nstate with values of 1,2, ... which density
#       goes to each state; 0= no density
# static:  used to indicate no linear predictors.
# mlevel: factor levels, if the marker was a factor, otherwise NULL
# param: optional params, eg."mean" or "std", not present= both
#
# return a list with 
#   rfun: the response function
#   pname: labels for the parameters
#   subset: which subset are referred to by the param argument
#   checkfun, optional
#   format: a function that turns the marker's linear predictors (at
#     covariate values of 0) into the parameters on the scale shown to the
#     user: the misclassification matrix for discrete, (mean, std) for
#     gaussian and logistic, the shapes for beta, (shape, rate) for gamma.
#     print(), summary() and fit$misclass / fit$markers use it, so each
#     distribution controls its own display (see "Response functions" in the
#     code vignette).
# The functions below are called by parsemarker2, and the presence or
#  absence of levels is sometimes enough for an error, .e.g. a factor can't
#  be a gaussian marker.  But the full set of y values that will be presented
#  is not yet available. The returned checkfun, if present, is called a
#  bit later by hmm when the response list is generated, and can be more
#  complete wrt a valid y.   For instance state:zed/gaussian where zed was
#  a list will be caught in the later check.
#  
# ---------------------------------------------------------------------
# Continuous markers: gaussian, logistic, beta and gamma.
# Each group of states (stateinfo$levels) has its own two parameters, held
#  in two linear predictors per group, in the order
#     (group 1: par 1, par 2), (group 2: par 1, par 2), ...
#   gaussian, logistic:  mean, log(std)
#   beta              :  log(shape1), log(shape2)
#   gamma             :  log(shape), log(rate)
# rfun(y, eta, gradient) returns an nstate x m matrix of densities
#  f(y | state) (0 for states that cannot produce the marker); with
#  gradient=TRUE it has attribute "gradient", an nstate x m x nlp array of
#  derivatives with respect to each linear predictor.
# Options, given in the marker formula:
#   param = "mean" or "std" (...): which parameters this formula line applies to
#   init  = matrix of starting values on the natural scale, one row per
#           group of states, columns (mean, std), (shape1, shape2) or
#           (shape, rate).
#           With ~ 0 (no estimated parameters) init gives the fixed values.
#  With no init, cmsh() starts every group at the overall mean and standard
#  deviation of the marker; give init when the states are fully hidden.
# ---------------------------------------------------------------------
contmarker <- function(name, parnames, stateinfo, mlevel, param, init,
                       static, logdens, link, inverse, check) {
    if (length(mlevel) > 0)
        stop("a factor variable cannot be a ", name, " marker")
    if (missing(static) || is.null(static)) static <- FALSE
    npeak  <- length(stateinfo$levels)
    index  <- stateinfo$index
    nstate <- length(index)
    pname  <- makedistlabels(stateinfo, parnames)
    if (missing(param) || is.null(param)) subset <- 1:ncol(pname)
    else {
        ii <- match(param, parnames)
        if (any(is.na(ii))) 
            stop("unrecognized ", name, " parameter: ", param[is.na(ii)])
        if (all(1:2 %in% ii)) subset <- 1:ncol(pname)
        else if (ii[1] == 1) subset <- 1:npeak * 2L - 1L
        else subset <- 2L * 1:npeak
    }

    initval <- NULL
    if (!missing(init) && !is.null(init)) {
        if (!is.matrix(init)) init <- matrix(init, nrow=npeak, byrow=TRUE)
        if (nrow(init) != npeak || ncol(init) != 2)
            stop("init for a ", name, " marker must have one row per group of ",
                 "states (", npeak, ") and columns ", 
                 paste(parnames, collapse=", "))
        if (!is.null(colnames(init)) && all(parnames %in% colnames(init)))
            init <- init[, parnames, drop=FALSE]
        if (!is.null(rownames(init)) && 
            all(as.character(stateinfo$levels) %in% rownames(init)))
            init <- init[as.character(stateinfo$levels), , drop=FALSE]
        initval <- c(t(link(init)))    # (par1, par2) for group 1, group 2, ...
    }
    if (static && is.null(initval))
        stop("init is needed for ~ 0 / ", name)

    # format: linear predictors (par1, par2 for group 1, group 2, ...) to the
    #  natural-scale parameters, one row per group of states
    format <- function(eta) {
        m <- inverse(matrix(eta, ncol=2, byrow=TRUE))
        colnames(m) <- parnames
        m
    }

    # eta: m rows, 2*npeak columns
    core <- function(y, eta, gradient) {
        ny <- length(y)
        if (!is.matrix(eta)) eta <- matrix(eta, nrow=ny, ncol=2L*npeak, 
                                           byrow= (length(eta) == 2L*npeak))
        yprob <- matrix(0., nstate, ny)
        if (gradient) ygrad <- array(0., dim=c(nstate, ny, 2L*npeak))
        for (i in seq_len(npeak)) {
            sg <- which(index == i)
            if (length(sg) == 0) next
            j <- 2L*i - 1L
            ld <- logdens(y, eta[, j], eta[, j+1L], gradient)
            f <- exp(ld$logf)
            yprob[sg, ] <- rep(f, each=length(sg))
            if (gradient) {
                ygrad[sg, , j]     <- rep(f * ld$d1, each=length(sg))
                ygrad[sg, , j+1L]  <- rep(f * ld$d2, each=length(sg))
            }
        }
        if (gradient) attr(yprob, "gradient") <- ygrad
        yprob
    }

    if (static) {  # fixed parameters, no linear predictors
        fixed <- initval
        rfun <- function(y, eta, gradient=FALSE)
            core(y, matrix(fixed, nrow=length(y), ncol=2L*npeak, byrow=TRUE),
                 FALSE)
        return(list(name=name, rfun=rfun, pname=NULL, static=TRUE,
                    checkfun=check, index=index, glevels=stateinfo$levels,
                    parnames=parnames, fixedpar=fixed, format=format,
                    type="continuous"))
    }
    rfun <- function(y, eta, gradient=FALSE) core(y, eta, gradient)
    list(name=name, rfun=rfun, pname=pname, subset=subset, checkfun=check,
         init=initval, static=FALSE, index=index, glevels=stateinfo$levels,
         parnames=parnames, format=format, type="continuous")
}

gaussian <- function(stateinfo, mlevel, param, init, static, ...) {
    logdens <- function(y, mu, lsd, gradient) {
        sd <- exp(lsd); z <- (y - mu)/sd
        out <- list(logf = dnorm(y, mu, sd, log=TRUE))
        if (gradient) { out$d1 <- z/sd; out$d2 <- z^2 - 1 }
        out
    }
    link <- function(m) {
        if (any(m[,2] <= 0)) stop("std must be positive")
        cbind(m[,1], log(m[,2]))
    }
    inverse <- function(m) cbind(m[,1], exp(m[,2]))
    check <- function(y)
        if (!is.numeric(y)) stop("a gaussian marker must be numeric")
    contmarker("gaussian", c("mean", "std"), stateinfo, mlevel, 
               if (missing(param)) NULL else param,
               if (missing(init)) NULL else init,
               if (missing(static)) FALSE else static, logdens, link, inverse,
               check)
}

makedistlabels <- function(stateinfo, parms) {
    nlev <- length(stateinfo$levels)
    nparm <-length(parms)
    # I don't need to say "state(dementia)", "dementia" will do
    if (stateinfo$sname == "state") 
        rbind(rep(stateinfo$levels, each=nparm), rep(parms, nlev))
    else {
        temp <- paste0(stateinfo$sname, '(', stateinfo$levels, ')')
        rbind(rep(temp, each=nparm), rep(parms, nlev))
    }
}

# logistic, a bit fatter tails.  Parameters: mean and log(std); the scale
#  of the logistic distribution is std * sqrt(3)/pi
logistic <- function(stateinfo, mlevel, param, init, static, ...) {
    logdens <- function(y, mu, lsd, gradient) {
        sc <- exp(lsd) * sqrt(3)/pi; z <- (y - mu)/sc
        out <- list(logf = dlogis(y, mu, sc, log=TRUE))
        if (gradient) {
            th <- tanh(z/2)
            out$d1 <- th/sc; out$d2 <- z*th - 1
        }
        out
    }
    link <- function(m) {
        if (any(m[,2] <= 0)) stop("std must be positive")
        cbind(m[,1], log(m[,2]))
    }
    inverse <- function(m) cbind(m[,1], exp(m[,2]))
    check <- function(y)
        if (!is.numeric(y)) stop("a logistic marker must be numeric")
    contmarker("logistic", c("mean", "std"), stateinfo, mlevel, 
               if (missing(param)) NULL else param,
               if (missing(init)) NULL else init,
               if (missing(static)) FALSE else static, logdens, link, inverse,
               check)
}

# beta distribution for a marker in (0, 1).  Parameters: log(shape1),
#  log(shape2)
beta <- function(stateinfo, mlevel, param, init, static, ...) {
    logdens <- function(y, la, lb, gradient) {
        a <- exp(la); b <- exp(lb)
        out <- list(logf = dbeta(y, a, b, log=TRUE))
        if (gradient) {
            out$d1 <- a * (digamma(a + b) - digamma(a) + log(y))
            out$d2 <- b * (digamma(a + b) - digamma(b) + log1p(-y))
        }
        out
    }
    link <- function(m) {
        if (any(m <= 0)) stop("shape parameters must be positive")
        log(m)
    }
    check <- function(y) {
        if (!is.numeric(y) || any(y <= 0 | y >= 1, na.rm=TRUE)) 
            stop("a beta marker must lie strictly between 0 and 1")
    }
    contmarker("beta", c("shape1", "shape2"), stateinfo, mlevel, 
               if (missing(param)) NULL else param,
               if (missing(init)) NULL else init,
               if (missing(static)) FALSE else static, logdens, link, 
               function(m) exp(m), check)
}

# gamma distribution for a positive marker.  Parameters: log(shape), log(rate)
#  log f = a log(b) + (a-1) log(y) - b y - lgamma(a), a = shape, b = rate
#  d log f / d log(a) = a (log(b) + log(y) - digamma(a))
#  d log f / d log(b) = a - b y
# The mean is shape/rate and the variance shape/rate^2.
# cmsh() has no automatic starting values for a gamma marker; give init
#  (e.g. shape = mean^2/var, rate = mean/var) when the states are hidden.
gamma <- function(stateinfo, mlevel, param, init, static, ...) {
    logdens <- function(y, la, lb, gradient) {
        a <- exp(la); b <- exp(lb)
        out <- list(logf = dgamma(y, shape=a, rate=b, log=TRUE))
        if (gradient) {
            out$d1 <- a * (lb + log(y) - digamma(a))
            out$d2 <- a - b*y
        }
        out
    }
    link <- function(m) {
        if (any(m <= 0)) stop("shape and rate must be positive")
        log(m)
    }
    check <- function(y) {
        if (!is.numeric(y) || any(y <= 0, na.rm=TRUE)) 
            stop("a gamma marker must be positive")
    }
    contmarker("gamma", c("shape", "rate"), stateinfo, mlevel, 
               if (missing(param)) NULL else param,
               if (missing(init)) NULL else init,
               if (missing(static)) FALSE else static, logdens, link,
               function(m) exp(m), check)
}

# multivariate logit, first category is the reference. If there are k
#  groups, eta will have k-1 columns
# This is used by other functions
mlogit <- function(eta, gradient=FALSE) {
    m <- ncol(eta)
    denom <- 1 + rowSums(exp(eta))
    pi <- cbind(1, exp(eta)) /denom

    if (gradient) {
        if (nrow(eta) ==1) {  # only one subject
            pi2 <- drop(pi)
            #temp <- diag(pi) - outer(pi, pi) # next line is a touch faster
            temp <- diag(pi2) -  rep(pi2, m+1) * rep(pi2, each= m+1)
            dmat <- temp[,-1,drop=FALSE] 
            dim(dmat) <- c(1, dim(dmat))  # make it a "row" per subject
        }
        else {
            dmat <- array(0., dim=c(nrow(eta), m+1, m))
            dmat[,1,] <- -pi[,-1]/denom
            for (j in 1:m) { # for each column of eta
                dmat[,j+1,] <- -pi[,j+1]* pi[,-1]
                dmat[,j+1, j] <- pi[,j+1]*(1-pi[,j+1])
            }
        }
        attr(pi, "gradient") <- dmat
    }
    pi
}
# 
#  The rfun for discrete has the usual y and eta arguments, along with a
# matrix etamap. The response y is a factor, each level maps to a column
# of etamap, rows of etamap correspond to states, and elements of etamap index
# to the linear predictors.  The details of this corresondence are controled
# by a pattern argument.  If there are no linear predictors and instead
# an init matrix, that is used directly: a row per state, col per response.
#  This function sets this all up and creates rfun
#
discrete <- function(stateinfo, mlevel, init, pattern, static) {
    # A discrete (categorical) marker: P(observed level | true state).
    # stateinfo$levels: the groups of states that share one distribution
    # stateinfo$index : for each state, which group (0 = marker impossible,
    #                   e.g. death)
    # mlevel          : the levels of the marker (a factor)
    # static=TRUE     : no parameters, init is a fixed probability matrix
    # static=FALSE    : multinomial-logit parameters are estimated
    ngroup <- length(stateinfo$levels)
    nstate <- length(stateinfo$index)
    index  <- stateinfo$index
    if (length(mlevel) ==0)
        stop("a discrete marker must be a factor")
    nlevel <- length(mlevel)

    # helper: put a user matrix into (group x marker level) order
    #  rows may be given per group (ngroup rows) or per state (nstate rows)
    to_group <- function(m, what) {
        if (!is.matrix(m)) stop(what, " must be a matrix")
        if (ncol(m) != nlevel)
            stop(what, " must have one column per level of the marker")
        if (!is.null(colnames(m))) {
            k <- match(mlevel, colnames(m))
            if (any(is.na(k))) stop(what, " column names must match the marker levels")
            m <- m[, k, drop=FALSE]
        }
        if (nrow(m) == ngroup && 
            (is.null(rownames(m)) || all(rownames(m) %in% stateinfo$levels))) {
            if (!is.null(rownames(m))) m <- m[match(stateinfo$levels, rownames(m)),,drop=FALSE]
            return(m)
        }
        if (nrow(m) == nstate) {   # one row per state: take a state from each group
            return(m[match(seq_len(ngroup), index), , drop=FALSE])
        }
        stop(what, " should have one row per ", stateinfo$sname, " level (",
             ngroup, ") or one row per state (", nstate, ")")
    }

    if (static) { # no parameters
        if (missing(init)) stop("init argument is needed for ~0 / discrete")
        if (nrow(init) == nstate) {
            missmat <- init
            if (!is.null(colnames(init))) 
                missmat <- init[, match(mlevel, colnames(init)), drop=FALSE]
        } else {
            g <- to_group(init, "init")
            missmat <- matrix(0, nrow=nstate, ncol= nlevel)
            missmat[index>0,] <- g[index[index>0],]
        }
        missmat <- unname(as.matrix(missmat))
        rfun <- function(y, error= missmat, gradient=FALSE) 
            error[,as.integer(y), drop=FALSE]  
        # format: the fixed matrix, one row per state
        format <- function(eta) missmat
        return(list(name="discrete", rfun=rfun, pname=NULL, static=TRUE,
                    missmat=missmat, mlevel=mlevel, index=index,
                    format=format, type="discrete"))
    }

    # ---- estimated case ----
    # emap: ngroup x nlevel.  -1 = reference cell, 0 = impossible cell,
    #  k > 0 = cell uses linear predictor k
    if (!missing(pattern)) {
        pattern <- to_group(pattern, "pattern")
        if (any(pattern != floor(pattern)) || any(pattern < 0))
            stop("a pattern matrix for discrete must contain integers >= 0")
        emap <- matrix(0L, ngroup, nlevel)
        for (i in 1:ngroup) {
            nz <- which(pattern[i,] > 0)
            if (length(nz) == 0) stop("pattern row ", i, " has no possible levels")
            ref <- nz[which.min(pattern[i, nz])]   # smallest code = reference
            emap[i, ref] <- -1L
            others <- setdiff(nz, ref)
            emap[i, others] <- pattern[i, others]
        }
        # number the linear predictors 1, 2, ... (shared codes stay shared)
        # (numbered row by row, so lp 1 belongs to the first group, etc.)
        te <- t(emap)
        uval <- unique(te[te > 0])
        emap[emap > 0] <- match(emap[emap > 0], uval)
    } else {
        # fully parameterized: first marker level is the reference in each row
        emap <- cbind(-1L, matrix(seq_len(ngroup*(nlevel-1)), nrow=ngroup,
                                 byrow=TRUE))
    }
    nlp <- max(emap)
    if (nlp == 0) stop("discrete marker has no free parameters; use ~0")

    # labels for each linear predictor: (group, marker level)
    lp_cell <- match(seq_len(nlp), emap)            # first cell using each lp
    glab <- if (stateinfo$sname == "state") as.character(stateinfo$levels) else
               paste0(stateinfo$sname, "(", stateinfo$levels, ")")
    pname <- rbind(glab[row(emap)[lp_cell]], as.character(mlevel)[col(emap)[lp_cell]])

    # starting values from an init probability matrix, if given
    initval <- NULL
    if (!missing(init)) {
        pinit <- to_group(init, "init")
        pinit <- pinit / rowSums(pinit)
        initval <- numeric(nlp)
        for (k in seq_len(nlp)) {
            cells <- which(emap == k, arr.ind=TRUE)
            r <- cells[1,1]; cc <- cells[1,2]
            ref <- which(emap[r,] == -1L)
            initval[k] <- log(pmax(pinit[r, cc], 1e-10) / pmax(pinit[r, ref], 1e-10))
        }
    }

    rfun <- function(y, eta, gradient=FALSE) {
        # y   : marker values (factor or integer codes), length m
        # eta : m x nlp matrix of linear predictors
        # value: nstate x m matrix of P(y | state); attr "gradient" is
        #        an nstate x m x nlp array of derivatives wrt eta
        y  <- as.integer(y)
        ny <- length(y)
        if (!is.matrix(eta)) eta <- matrix(eta, ncol=nlp)
        phat <- matrix(0., nrow=nstate, ncol=ny)
        if (gradient) gmat <- array(0., dim=c(nstate, ny, nlp))
        for (g in seq_len(ngroup)) {
            sg <- which(index == g)
            if (length(sg) == 0) next
            cells <- emap[g,]
            E <- matrix(-Inf, ny, nlevel)
            E[, cells == -1L] <- 0
            for (l in which(cells > 0)) E[,l] <- eta[, cells[l]]
            emax <- apply(E, 1, max)
            P <- exp(E - emax)
            P <- P / rowSums(P)
            py <- P[cbind(seq_len(ny), y)]
            phat[sg,] <- rep(py, each=length(sg))
            if (gradient) {
                for (k in unique(cells[cells > 0])) {
                    lk <- which(cells == k)
                    dk <- py * ((y %in% lk) - rowSums(P[, lk, drop=FALSE]))
                    gmat[sg,,k] <- rep(dk, each=length(sg))
                }
            }
        }
        if (gradient) attr(phat, "gradient") <- gmat
        phat
    }

    check <- function(y) {
        if (!is.factor(y))
            stop("marker for a discrete distribution must be categorical")
    }

    # format: the nlp linear predictors to the misclassification matrix,
    #  one row per group of states, one column per marker level
    format <- function(eta) {
        P <- matrix(0, ngroup, nlevel)
        for (g in seq_len(ngroup)) {
            e <- rep(-Inf, nlevel)
            e[emap[g,] == -1L] <- 0
            k <- which(emap[g,] > 0)
            e[k] <- eta[emap[g, k]]
            e <- exp(e - max(e))
            P[g,] <- e/sum(e)
        }
        P
    }
            
    list(name="discrete", rfun=rfun, pname=pname, subset=seq_len(nlp),
         init=initval, emap=emap, checkfun= check, static=FALSE,
         mlevel=mlevel, index=index, glevels=stateinfo$levels,
         format=format, type="discrete")
}
        
