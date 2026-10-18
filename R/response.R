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
# The functions below are called by parsemarker2, and the presence or
#  absence of levels is sometimes enough for an error, .e.g. a factor can't
#  be a gaussian marker.  But the full set of y values that will be presented
#  is not yet available. The returned checkfun, if present, is called a
#  bit later by hmm when the response list is generated, and can be more
#  complete wrt a valid y.   For instance state:zed/gaussian where zed was
#  a list will be caught in the later check.
#  
.continuous_marker_init <- function(init, npeak, family) {
    if (is.null(init)) return(NULL)
    pnames <- switch(family,
                     beta=c("shape1", "shape2"),
                     gamma=c("shape", "rate"),
                     c("mean", "std"))
    if (is.matrix(init)) {
        if (!all(dim(init) == c(npeak, 2L)))
            stop(family, " init must have one row per latent marker level and two columns")
        if (!is.null(colnames(init))) {
            indx <- match(pnames, colnames(init))
            if (anyNA(indx)) stop(family, " init columns must be ", paste(pnames, collapse=", "))
            init <- init[, indx, drop=FALSE]
        }
        init <- c(t(init))
    } else {
        if (!is.numeric(init) || length(init) != 2L * npeak)
            stop(family, " init must contain two values per latent marker level")
        init <- as.numeric(init)
    }
    if (family %in% c("beta", "gamma")) {
        if (any(!is.finite(init) | init <= 0))
            stop("beta init shape parameters must be finite and positive")
        log(init)
    } else {
        scale_index <- 2L * seq_len(npeak)
        if (any(!is.finite(init)) || any(init[scale_index] <= 0))
            stop(family, " init standard deviations must be finite and positive")
        init[scale_index] <- log(init[scale_index])
        init
    }
}

.continuous_marker_density <- function(y, eta, gradient=FALSE, npeak, map,
                                       family) {
    if (!is.matrix(eta)) eta <- matrix(eta, nrow=length(y))
    if (nrow(eta) != length(y) || ncol(eta) != 2L * npeak)
        stop("marker linear predictor has the wrong dimension")
    nstate <- length(map)
    yprob <- matrix(0, nstate, length(y))
    if (gradient)
        ygrad <- array(0, dim=c(nstate, length(y), ncol(eta)))

    for (i in seq_len(npeak)) {
        states <- which(map == i)
        if (!length(states)) next
        j <- 2L * i - 1L
        if (family == "gaussian") {
            mu <- eta[, j]
            std <- exp(eta[, j + 1L])
            dens <- dnorm(y, mu, std)
            d1 <- (y - mu) / std^2
            d2 <- ((y - mu) / std)^2 - 1
        } else if (family == "logistic") {
            mu <- eta[, j]
            std <- exp(eta[, j + 1L])
            scale <- std * sqrt(3) / pi
            z <- (y - mu) / scale
            dens <- dlogis(y, mu, scale)
            centered <- 2 * plogis(z) - 1
            d1 <- centered / scale
            d2 <- -1 + z * centered
        } else if (family == "beta") {
            shape1 <- exp(eta[, j])
            shape2 <- exp(eta[, j + 1L])
            dens <- dbeta(y, shape1, shape2)
            common <- digamma(shape1 + shape2)
            d1 <- shape1 * (log(y) - digamma(shape1) + common)
            d2 <- shape2 * (log1p(-y) - digamma(shape2) + common)
        } else {
            shape <- exp(eta[, j])
            rate <- exp(eta[, j + 1L])
            dens <- dgamma(y, shape=shape, rate=rate)
            d1 <- shape * (eta[, j + 1L] + log(y) - digamma(shape))
            d2 <- shape - rate * y
        }
        yprob[states, ] <- rep(dens, each=length(states))
        if (gradient) {
            ygrad[states, , j] <- rep(dens * d1, each=length(states))
            ygrad[states, , j + 1L] <- rep(dens * d2, each=length(states))
        }
    }
    if (gradient) attr(yprob, "gradient") <- ygrad
    yprob
}

.continuous_marker <- function(stateinfo, mlevel, param=NULL, init=NULL,
                               static=FALSE, family) {
    if (length(mlevel)) stop("a factor variable cannot be a ", family, " marker")
    npeak <- length(stateinfo$levels)
    pnames <- switch(family,
                     beta=c("shape1", "shape2"),
                     gamma=c("shape", "rate"),
                     c("mean", "std"))
    initial <- .continuous_marker_init(init, npeak, family)
    if (is.null(param)) subset <- seq_len(2L * npeak)
    else {
        index <- match(param, pnames)
        if (anyNA(index)) stop("unrecognized ", family, " parameter: ", param[is.na(index)][1])
        subset <- sort(unique(unlist(lapply(index, function(k)
            2L * seq_len(npeak) - (2L - k)))))
    }
    state_index <- stateinfo$index
    marker_family <- family
    marker_npeak <- npeak
    rfun <- function(y, eta, gradient=FALSE) {
        .continuous_marker_density(
            y, eta, gradient=gradient, npeak=marker_npeak,
            map=state_index, family=marker_family
        )
    }

    check <- if (family == "beta") {
        function(y) {
            observed <- y[!is.na(y)]
            if (!is.numeric(y) || any(!is.finite(observed)) ||
                any(observed <= 0 | observed >= 1))
                stop("a beta marker must be numeric and strictly between 0 and 1")
        }
    } else if (family == "gamma") {
        function(y) {
            observed <- y[!is.na(y)]
            if (!is.numeric(y) || any(!is.finite(observed)) || any(observed <= 0))
                stop("a gamma marker must be finite, numeric, and strictly positive")
        }
    } else {
        function(y) {
            observed <- y[!is.na(y)]
            if (!is.numeric(y) || any(!is.finite(observed)))
                stop("a ", family, " marker must be finite and numeric")
        }
    }

    if (static) {
        if (is.null(initial)) stop("init is required for a fixed ", family, " marker")
        fixed_initial <- initial
        fixed_basefun <- rfun
        fixed_rfun <- function(y, gradient=FALSE) {
            eta <- matrix(fixed_initial, nrow=length(y),
                          ncol=length(fixed_initial), byrow=TRUE)
            fixed_basefun(y, eta, gradient=FALSE)
        }
        return(list(name=family, rfun=fixed_rfun, pname=NULL,
                    subset=integer(), checkfun=check, initial=initial,
                    initial_user=TRUE, fixed=TRUE, family=family,
                    npeak=npeak, latent_levels=stateinfo$levels,
                    state_index=stateinfo$index))
    }

    list(name=family, rfun=rfun,
         pname=makedistlabels(stateinfo, pnames), subset=subset,
         checkfun=check, initial=initial, initial_user=!is.null(initial),
         fixed=FALSE, family=family, npeak=npeak,
         latent_levels=stateinfo$levels, state_index=stateinfo$index)
}

.marker_default_initial <- function(response, y) {
    y <- y[!is.na(y)]
    if (!length(y) || response$family == "discrete" || response$fixed)
        return(response$initial)
    npeak <- response$npeak
    if (response$family %in% c("gaussian", "logistic")) {
        center <- mean(y)
        spread <- sd(y)
        if (!is.finite(spread) || spread <= 0) spread <- 1
        rep(c(center, log(spread)), npeak)
    } else if (response$family == "beta") {
        center <- mean(y)
        variance <- var(y)
        common <- center * (1 - center) / variance - 1
        if (!is.finite(common) || common <= 0) {
            shape1 <- shape2 <- 2
        } else {
            shape1 <- center * common
            shape2 <- (1 - center) * common
        }
        rep(log(c(shape1, shape2)), npeak)
    } else {
        center <- mean(y)
        variance <- var(y)
        if (!is.finite(variance) || variance <= 0 || center <= 0) {
            shape <- rate <- 1
        } else {
            shape <- center^2 / variance
            rate <- center / variance
        }
        rep(log(c(shape, rate)), npeak)
    }
}

gaussian <- function(stateinfo, mlevel, param, init, static=FALSE, ...) {
    .continuous_marker(stateinfo, mlevel,
                       if (missing(param)) NULL else param,
                       if (missing(init)) NULL else init,
                       static, "gaussian")
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

# logistic, a bit fatter tails
logistic <- function(stateinfo, mlevel, param, init, static=FALSE, ...) {
    .continuous_marker(stateinfo, mlevel,
                       if (missing(param)) NULL else param,
                       if (missing(init)) NULL else init,
                       static, "logistic")
}

# beta distribution
beta <- function(stateinfo, mlevel, param, init, static=FALSE, ...) {
    .continuous_marker(stateinfo, mlevel,
                       if (missing(param)) NULL else param,
                       if (missing(init)) NULL else init,
                       static, "beta")
}

# gamma distribution for a positive marker, parameterized by log(shape) and
# log(rate).  The natural-scale initialization uses shape and rate.
gamma <- function(stateinfo, mlevel, param, init, static=FALSE, ...) {
    .continuous_marker(stateinfo, mlevel,
                       if (missing(param)) NULL else param,
                       if (missing(init)) NULL else init,
                       static, "gamma")
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
.discrete_init_matrix <- function(init, stateinfo, mlevel) {
    if (!is.matrix(init)) stop("discrete init must be a probability matrix")
    ngroup <- length(stateinfo$levels)
    nstate <- length(stateinfo$index)
    if (ncol(init) != length(mlevel))
        stop("discrete init columns must match marker levels")
    if (!is.null(colnames(init))) {
        indx <- match(mlevel, colnames(init))
        if (anyNA(indx) || anyDuplicated(indx))
            stop("discrete init columns must match marker levels")
        init <- init[, indx, drop=FALSE]
    }
    if (nrow(init) == ngroup) {
        if (!is.null(rownames(init))) {
            indx <- match(as.character(stateinfo$levels), rownames(init))
            if (anyNA(indx) || anyDuplicated(indx))
                stop("discrete init rows must match latent marker levels")
            init <- init[indx, , drop=FALSE]
        }
    } else if (nrow(init) == nstate) {
        grouped <- matrix(NA_real_, ngroup, ncol(init))
        for (i in seq_len(ngroup)) {
            rows <- which(stateinfo$index == i)
            grouped[i, ] <- colMeans(init[rows, , drop=FALSE])
        }
        init <- grouped
    } else stop("discrete init must have one row per latent marker level or state")
    if (any(!is.finite(init)) || any(init < 0) || any(rowSums(init) <= 0))
        stop("discrete init must contain finite nonnegative probabilities")
    init / rowSums(init)
}

.static_discrete_init <- function(init, stateinfo, mlevel) {
    if (!is.matrix(init)) stop("discrete init must be a probability matrix")
    ngroup <- length(stateinfo$levels)
    nstate <- length(stateinfo$index)
    if (ncol(init) != length(mlevel))
        stop("discrete init columns must match marker levels")
    if (!is.null(colnames(init))) {
        indx <- match(mlevel, colnames(init))
        if (anyNA(indx) || anyDuplicated(indx))
            stop("discrete init columns must match marker levels")
        init <- init[, indx, drop=FALSE]
    }
    if (any(!is.finite(init)) || any(init < 0))
        stop("discrete init must contain finite nonnegative probabilities")

    if (nrow(init) == nstate) {
        return(list(state=unname(init), display=init))
    }
    grouped <- matrix(0, ngroup, ncol(init),
                      dimnames=list(as.character(stateinfo$levels), mlevel))
    if (nrow(init) > ngroup) {
        rows <- if (is.null(rownames(init))) seq_len(ngroup) else
            match(as.character(stateinfo$levels), rownames(init))
        if (anyNA(rows)) rows <- seq_len(ngroup)
        if (anyDuplicated(rows))
            stop("discrete init rows must uniquely match latent marker levels")
        grouped[,] <- init[rows, , drop=FALSE]
    } else if (is.null(rownames(init))) {
        grouped[seq_len(nrow(init)), ] <- init
    } else {
        rows <- match(rownames(init), as.character(stateinfo$levels))
        if (anyNA(rows) || anyDuplicated(rows))
            stop("discrete init rows must match latent marker levels")
        grouped[rows, ] <- init
    }
    missmat <- matrix(0, nstate, ncol(init), dimnames=list(NULL, mlevel))
    keep <- stateinfo$index > 0L
    missmat[keep, ] <- grouped[stateinfo$index[keep], , drop=FALSE]
    list(state=unname(missmat), display=grouped)
}

.discrete_initial_eta <- function(prob, emap) {
    ncoef <- max(c(0L, emap))
    if (!ncoef) return(numeric())
    values <- vector("list", ncoef)
    floorp <- sqrt(.Machine$double.eps)
    for (i in seq_len(nrow(emap))) {
        ref <- which(emap[i, ] == -1L)
        if (!length(ref)) next
        refprob <- exp(mean(log(pmax(prob[i, ref], floorp))))
        for (k in unique(emap[i, emap[i, ] > 0L])) {
            cells <- which(emap[i, ] == k)
            values[[k]] <- c(values[[k]], log(pmax(prob[i, cells], floorp) / refprob))
        }
    }
    vapply(values, function(x) if (length(x)) mean(x) else 0, numeric(1))
}

discrete <- function(stateinfo, mlevel, init, pattern, static=FALSE, ...) {
    ngroup <- length(stateinfo$levels)
    nstate <- length(stateinfo$index)
    if (!length(mlevel)) stop("a discrete marker must be a factor")
    nlevel <- length(mlevel)

    if (!missing(pattern)) {
        if (!is.matrix(pattern) || nrow(pattern) != ngroup || ncol(pattern) != nlevel)
            stop("discrete pattern must have one row per latent marker level and one column per marker level")
        if (!is.null(colnames(pattern))) {
            indx <- match(mlevel, colnames(pattern))
            if (anyNA(indx) || anyDuplicated(indx))
                stop("discrete pattern columns must match marker levels")
            pattern <- pattern[, indx, drop=FALSE]
        }
        if (!is.null(rownames(pattern))) {
            indx <- match(as.character(stateinfo$levels), rownames(pattern))
            if (anyNA(indx) || anyDuplicated(indx))
                stop("discrete pattern rows must match latent marker levels")
            pattern <- pattern[indx, , drop=FALSE]
        }
        if (any(!is.finite(pattern)) || any(pattern != floor(pattern)) || any(pattern < 0))
            stop("discrete pattern entries must be nonnegative integers")
        if (any(rowSums(pattern > 0) == 0))
            stop("each discrete pattern row must contain a positive entry")
        emap <- matrix(match(pattern, unique(c(0L, as.integer(pattern)))) - 1L,
                       nrow=ngroup, ncol=nlevel)
        nz <- emap > 0L
        maxp <- max(emap)
        shared <- integer()
        if (maxp > 0L) {
            locations <- lapply(seq_len(maxp), function(k)
                unique(row(emap)[emap == k]))
            shared <- which(lengths(locations) > 1L)
        }
        if (length(shared)) {
            for (k in shared) emap[emap == k] <- emap[emap == k] + maxp
        }
        for (i in seq_len(ngroup)) {
            candidates <- unique(emap[i, emap[i, ] > 0L])
            reference <- min(candidates)
            emap[i, emap[i, ] == reference] <- -1L
        }
        old <- unique(emap[emap > 0L])
        if (length(old)) emap[emap > 0L] <- match(emap[emap > 0L], old)
    } else if (!static) {
        if (nlevel == 1L) emap <- matrix(-1L, ngroup, 1L)
        else {
            emap <- cbind(-1L,
                          matrix(seq_len(ngroup * (nlevel - 1L)),
                                 nrow=ngroup, byrow=TRUE))
        }
    } else emap <- matrix(-1L, ngroup, nlevel)

    dimnames(emap) <- list(paste0(stateinfo$sname, stateinfo$levels),
                           marker=mlevel)
    ncoef <- max(c(0L, emap))
    if (ncoef) {
        labels <- outer(rownames(emap), colnames(emap), paste, sep="_")
        pname <- labels[match(seq_len(ncoef), emap)]
    } else pname <- NULL

    init_prob <- if (missing(init)) NULL else if (static) NULL else
        .discrete_init_matrix(init, stateinfo, mlevel)
    if (static) {
        if (missing(init)) stop("init is required for a fixed discrete marker")
        static_init <- .static_discrete_init(init, stateinfo, mlevel)
        fixed_probability <- static_init$state
        rfun <- function(y, gradient=FALSE) {
            if (is.factor(y)) y <- as.integer(y)
            fixed_probability[, y, drop=FALSE]
        }
        return(list(name="discrete", family="discrete", rfun=rfun,
                    pname=NULL, subset=integer(), initial=numeric(),
                    initial_user=TRUE, fixed=TRUE,
                    probability=static_init$display,
                    marker_levels=mlevel, latent_levels=stateinfo$levels,
                    state_index=stateinfo$index))
    }

    initial <- if (is.null(init_prob)) NULL else
        .discrete_initial_eta(init_prob, emap)
    response_map <- emap
    state_index <- stateinfo$index
    rfun <- function(y, eta, gradient=FALSE) {
        ny <- length(y)
        if (is.factor(y)) y <- as.integer(y)
        if (anyNA(y) || any(y < 1L) || any(y > ncol(response_map)))
            stop("invalid discrete marker level")
        if (!is.matrix(eta)) eta <- matrix(eta, nrow=ny)
        if (nrow(eta) != ny || ncol(eta) < max(c(0L, response_map)))
            stop("discrete marker linear predictor has the wrong dimension")
        phat <- matrix(0, nrow=length(state_index), ncol=ny)
        if (gradient)
            gmat <- array(0, dim=c(length(state_index), ny, ncol(eta)))
        obs <- seq_len(ny)
        for (i in seq_len(nrow(response_map))) {
            states <- which(state_index == i)
            if (!length(states)) next
            rowmap <- response_map[i, ]
            logits <- matrix(-Inf, ny, ncol(response_map))
            logits[, rowmap == -1L] <- 0
            for (k in unique(rowmap[rowmap > 0L]))
                logits[, rowmap == k] <- eta[, k]
            rowmax <- apply(logits, 1L, max)
            weights <- exp(logits - rowmax)
            prob <- weights / rowSums(weights)
            observed <- prob[cbind(obs, y)]
            phat[states, ] <- rep(observed, each=length(states))
            if (gradient) {
                for (k in unique(rowmap[rowmap > 0L])) {
                    total <- rowSums(prob[, rowmap == k, drop=FALSE])
                    deriv <- observed * (as.numeric(rowmap[y] == k) - total)
                    gmat[states, , k] <- rep(deriv, each=length(states))
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
    if (ncoef) {
        pname <- vapply(seq_len(ncoef), function(k) {
            cell <- which(emap == k, arr.ind=TRUE)
            paste(rownames(emap)[unique(cell[, 1L])], collapse="/")
        }, character(1L))
        plevel <- vapply(seq_len(ncoef), function(k) {
            cell <- which(emap == k, arr.ind=TRUE)
            paste(colnames(emap)[unique(cell[, 2L])], collapse="/")
        }, character(1L))
        pname <- rbind(pname, plevel)
    }
    list(name="discrete", family="discrete", rfun=rfun, pname=pname,
         subset=seq_len(ncoef), checkfun=check, initial=initial,
         initial_user=!is.null(initial), fixed=FALSE, emap=emap,
         probability=init_prob, marker_levels=mlevel,
         latent_levels=stateinfo$levels, state_index=stateinfo$index)
}
        
