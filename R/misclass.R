# The misclassification (emission) matrix of a discrete marker is stored
#  in the fit as its own component, fit$misclass, one element per marker:
#    matrix   : P(observed level | true state); rows = states that can
#               produce the marker, columns = marker levels
#    se, lower, upper : standard error and confidence limits of each
#               probability (NA for a fixed matrix or a reference-only row)
#    coefficients : the estimated parameters on the logit scale
#               (log odds against the row's reference level), with se
#    fixed    : TRUE if the matrix was supplied (~ 0 / discrete(init=))
# The helper below is called once, at the end of cmsh().  The matrix itself
#  comes from the format() function that discrete() returns, so the
#  distribution controls its own display (see response.R).

build_misclass <- function(rlist, mnames, B, cmap, nlp, vcov, statenames,
                           conf.int = 0.95) {
    out <- list(); ecols <- integer(0)
    z <- qnorm(1 - (1 - conf.int)/2)
    for (i in seq_along(rlist)) {
        r <- rlist[[i]]
        if (!identical(r$type, "discrete")) next
        has <- which(r$index > 0)          # states that can show the marker
        if (isTRUE(r$static)) {
            has <- which(r$index > 0 & rowSums(abs(r$format())) > 0)
            P <- r$format()[has, , drop=FALSE]
            dimnames(P) <- list(state = statenames[has], observed = r$mlevel)
            NAm <- P; NAm[,] <- NA_real_
            out[[mnames[i]]] <- list(matrix = P, se = NAm, lower = NAm,
                                     upper = NAm, coefficients = NULL,
                                     fixed = TRUE, covariates = FALSE,
                                     conf.int = conf.int)
            next
        }
        emap <- r$emap
        cols <- r$e2map + nlp[1]           # this marker's columns of cmap
        ecols <- c(ecols, cols)
        eta0 <- B[1, cols]                 # intercepts (logit scale)
        pfun <- r$format                  # group x level probabilities
        Pg <- pfun(eta0)
        # delta method for the probabilities
        idx <- cmap[1, cols]               # parameter numbers (<=0: fixed)
        nk <- length(cols)
        V <- matrix(0, nk, nk)
        ok <- !is.null(vcov) && all(idx[idx > 0] <= nrow(vcov))
        if (ok && any(idx > 0))
            V[idx > 0, idx > 0] <- vcov[idx[idx > 0], idx[idx > 0]]
        G <- array(0, dim = c(dim(Pg), nk))
        h <- 1e-5
        for (k in seq_len(nk)) {
            ek <- rep(0, nk); ek[k] <- h
            G[,,k] <- (pfun(eta0 + ek) - pfun(eta0 - ek))/(2*h)
        }
        seg <- Pg
        for (a in seq_len(nrow(Pg))) for (b in seq_len(ncol(Pg))) {
            g <- G[a, b, ]
            seg[a, b] <- sqrt(max(0, drop(g %*% V %*% g)))
        }
        if (!ok) seg[,] <- NA_real_
        # confidence limits on the logit(p) scale
        sel <- seg/(Pg*(1 - Pg))
        lo <- plogis(qlogis(Pg) - z*sel); hi <- plogis(qlogis(Pg) + z*sel)
        bad <- !is.finite(sel)
        lo[bad] <- NA_real_; hi[bad] <- NA_real_
        expand <- function(m) {
            m <- m[r$index[has], , drop=FALSE]
            dimnames(m) <- list(state = statenames[has], observed = r$mlevel)
            m
        }
        # the parameters themselves
        grp <- sapply(seq_len(nrow(emap)), function(g)
                      paste(statenames[r$index == g], collapse = "/"))
        lab <- sapply(seq_len(nk), function(k) {
            cell <- which(emap == k, arr.ind = TRUE)
            paste0("P(obs ", paste(r$mlevel[unique(cell[, 2])], collapse = "/"),
                   " | true ", paste(grp[unique(cell[, 1])], collapse = ", "), ")")
        })
        sek <- if (ok) sqrt(diag(V)) else rep(NA_real_, nk)
        sek[idx <= 0] <- NA_real_
        ctab <- cbind(logit = eta0, se = sek, prob = NA_real_,
                      lower = NA_real_, upper = NA_real_)
        for (k in seq_len(nk)) {
            cell <- which(emap == k, arr.ind = TRUE)[1, ]
            ctab[k, "prob"]  <- Pg[cell[1], cell[2]]
            ctab[k, "lower"] <- lo[cell[1], cell[2]]
            ctab[k, "upper"] <- hi[cell[1], cell[2]]
        }
        rownames(ctab) <- lab
        out[[mnames[i]]] <- list(matrix = expand(Pg), se = expand(seg),
                                 lower = expand(lo), upper = expand(hi),
                                 coefficients = ctab, fixed = FALSE,
                                 covariates = any(cmap[-1, cols] != 0),
                                 conf.int = conf.int)
    }
    epar <- if (length(ecols)) sort(unique(abs(cmap[, ecols][cmap[, ecols] != 0])))
            else integer(0)
    list(misclass = out, ecols = ecols, epar = epar)
}

# Continuous markers.  fit$markers has one element per gaussian, logistic,
#  beta or gamma marker:
#    distribution : "gaussian", "logistic", "beta" or "gamma"
#    parameters   : one row per group of states, the two parameters on
#                   their natural scale, at covariate values of 0, as given
#                   by the distribution's format() function (response.R)
#    se, lower, upper : standard errors and confidence limits of those
#    coefficients : every coefficient of the marker's linear predictors
#                   (mean, log(std) / log(shape); plus any covariates)
#    fixed, covariates
build_markers <- function(rlist, mnames, B, cmap, nlp, vcov, statenames,
                          conf.int = 0.95) {
    out <- list(); ecols <- integer(0)
    z <- qnorm(1 - (1 - conf.int)/2)
    for (i in seq_along(rlist)) {
        r <- rlist[[i]]
        if (!identical(r$type, "continuous")) next
        ng <- length(r$glevels)
        grp <- sapply(seq_len(ng), function(g)
                      paste(statenames[r$index == g], collapse = "/"))
        nat <- function(eta) {      # natural scale, ng x 2, from format()
            m <- r$format(eta)
            dimnames(m) <- list(state = grp, parameter = r$parnames)
            m
        }
        if (isTRUE(r$static)) {
            P <- nat(r$fixedpar); NAm <- P; NAm[,] <- NA_real_
            out[[mnames[i]]] <- list(distribution = r$name, parameters = P,
                                     se = NAm, lower = NAm, upper = NAm,
                                     coefficients = NULL, fixed = TRUE,
                                     covariates = FALSE, conf.int = conf.int)
            next
        }
        cols <- r$e2map + nlp[1]
        ecols <- c(ecols, cols)
        eta0 <- B[1, cols]
        idx <- cmap[1, cols]
        ok <- !is.null(vcov)
        se0 <- rep(NA_real_, length(cols))
        if (ok) {
            use <- idx > 0 & idx <= nrow(vcov)
            se0[use] <- sqrt(pmax(0, diag(vcov)[idx[use]]))
        }
        P <- nat(eta0)
        # each natural parameter is an increasing function of its own linear
        #  predictor: the delta method for the se, and the limits are the
        #  transformed limits of the linear predictor
        h <- 1e-6
        dP <- (nat(eta0 + h) - nat(eta0 - h))/(2*h)
        sem <- dP * matrix(se0, ncol = 2, byrow = TRUE)
        lo <- nat(eta0 - z*se0); hi <- nat(eta0 + z*se0)
        # all coefficients of these linear predictors
        cm <- cmap[, cols, drop = FALSE]
        w <- which(cm != 0, arr.ind = TRUE)
        ctab <- NULL
        if (nrow(w)) {
            pn <- abs(cm[w])
            sec <- rep(NA_real_, length(pn))
            if (ok) {
                use <- cm[w] > 0 & pn <= nrow(vcov)
                sec[use] <- sqrt(pmax(0, diag(vcov)[pn[use]]))
            }
            ctab <- cbind(coef = B[, cols, drop = FALSE][w], se = sec)
            rownames(ctab) <- paste0(colnames(cmap)[cols][w[, 2]], " ",
                                     rownames(cmap)[w[, 1]])
            ctab <- ctab[order(w[, 2], w[, 1]), , drop = FALSE]
        }
        out[[mnames[i]]] <- list(distribution = r$name, parameters = P,
                                 se = sem, lower = lo, upper = hi,
                                 coefficients = ctab, fixed = FALSE,
                                 covariates = any(cmap[-1, cols] != 0),
                                 conf.int = conf.int)
    }
    epar <- if (length(ecols)) sort(unique(abs(cmap[, ecols][cmap[, ecols] != 0])))
            else integer(0)
    list(markers = out, ecols = ecols, epar = epar)
}

print_markers <- function(mk, digits = 4, ci = FALSE) {
    for (nm in names(mk)) {
        m <- mk[[nm]]
        cat("Marker ", nm, " (", m$distribution,
            if (m$fixed) ", fixed" else ", estimated",
            "): parameters by true state\n", sep = "")
        if (ci && !m$fixed) {
            tab <- NULL
            for (k in 1:2) {
                t1 <- cbind(m$parameters[, k], m$se[, k], m$lower[, k],
                            m$upper[, k])
                rownames(t1) <- paste(rownames(m$parameters),
                                      colnames(m$parameters)[k])
                tab <- rbind(tab, t1)
            }
            colnames(tab) <- c("estimate", "se",
                               paste(c("lower", "upper"), format(m$conf.int)))
            print(round(tab, digits))
        } else print(round(m$parameters, digits))
        if (isTRUE(m$covariates)) {
            cat("  (the marker has covariates; parameters shown at covariate",
                "values of 0)\n  coefficients of its linear predictors:\n")
            print(round(m$coefficients, digits))
        }
        cat("\n")
    }
    invisible(mk)
}

# Extract the misclassification matrix (or the whole component)
misclass <- function(fit, marker, what = c("matrix", "se", "lower", "upper",
                                           "coefficients", "all")) {
    if (!inherits(fit, "cmsh")) stop("fit must be a cmsh object")
    what <- match.arg(what)
    mc <- fit$misclass
    if (is.null(mc)) return(NULL)
    if (missing(marker)) {
        if (length(mc) > 1 && what != "all")
            return(lapply(mc, function(m) m[[what]]))
        marker <- 1
    }
    if (what == "all") mc[[marker]] else mc[[marker]][[what]]
}

print_misclass <- function(mc, digits = 4, ci = FALSE) {
    for (nm in names(mc)) {
        m <- mc[[nm]]
        cat("Misclassification matrix for ", nm,
            if (m$fixed) " (fixed)" else " (estimated)",
            ": P(observed | true state)\n", sep = "")
        print(round(m$matrix, digits))
        if (isTRUE(m$covariates))
            cat("  (the marker has covariates; matrix shown at covariate values of 0)\n")
        if (!m$fixed && ci) {
            cat("\n")
            tab <- m$coefficients
            colnames(tab) <- c("logit", "se(logit)", "prob",
                               paste(c("lower", "upper"), format(m$conf.int)))
            print(round(tab, digits))
        }
        cat("\n")
    }
    invisible(mc)
}

# summary: transition coefficients with hazard ratios and confidence
#  limits, one block per transition; the misclassification matrix separately
summary.cmsh <- function(object, conf.int = 0.95, ...) {
    z <- qnorm(1 - (1 - conf.int)/2)
    cf <- object$coefficients
    se <- rep(NA_real_, length(cf))
    if (!is.null(object$var) && nrow(object$var) == length(cf))
        se <- sqrt(diag(object$var))
    cmap <- object$cmap
    cols <- seq_len(ncol(cmap))
    if (length(object$ecols)) cols <- setdiff(cols, object$ecols)
    tabs <- list()
    for (j in cols) {
        k <- cmap[, j]
        if (all(k == 0)) next
        use <- which(k != 0)
        b <- cf[abs(k[use])]; s <- se[abs(k[use])]
        s[k[use] < 0] <- NA_real_                  # fixed coefficients
        tab <- cbind(coef = b, "exp(coef)" = exp(b), "se(coef)" = s,
                     lower = exp(b - z*s), upper = exp(b + z*s),
                     z = b/s, p = pchisq((b/s)^2, 1, lower.tail = FALSE))
        colnames(tab)[4:5] <- paste(c("lower", "upper"), format(conf.int))
        rownames(tab) <- rownames(cmap)[use]
        tabs[[colnames(cmap)[j]]] <- tab
    }
    mc <- object$misclass
    if (!is.null(mc) && conf.int != 0.95) {       # recompute the limits
        for (i in seq_along(mc)) if (!mc[[i]]$fixed) {
            P <- mc[[i]]$matrix; sel <- mc[[i]]$se/(P*(1 - P))
            lo <- plogis(qlogis(P) - z*sel); hi <- plogis(qlogis(P) + z*sel)
            lo[!is.finite(sel)] <- NA; hi[!is.finite(sel)] <- NA
            mc[[i]]$lower <- lo; mc[[i]]$upper <- hi
            ct <- mc[[i]]$coefficients
            for (k in seq_len(nrow(ct))) {         # match each parameter's cell
                cell <- which(abs(P - ct[k, "prob"]) < 1e-12, arr.ind = TRUE)[1, ]
                ct[k, "lower"] <- lo[cell[1], cell[2]]
                ct[k, "upper"] <- hi[cell[1], cell[2]]
            }
            mc[[i]]$coefficients <- ct
            mc[[i]]$conf.int <- conf.int
        }
    }
    npar <- length(unique(cmap[cmap > 0]))
    ll <- object$loglik
    out <- list(call = object$call, coefficients = tabs, misclass = mc,
                markers = object$markers,
                loglik = ll, df = npar,
                AIC = -2*ll[length(ll)] + 2*npar,
                n = object$n, states = object$states,
                converged = object$fit$converged, conf.int = conf.int,
                na.action = object$na.action)
    class(out) <- "summary.cmsh"
    out
}

print.summary.cmsh <- function(x, digits = max(1L, getOption("digits") - 3L),
                               signif.stars = FALSE, ...) {
    if (!is.null(x$call)) { cat("Call:\n"); dput(x$call); cat("\n") }
    cat(" States: ", paste(paste(seq_along(x$states), x$states, sep = "= "),
                           collapse = ", "), "\n\n")
    for (nm in names(x$coefficients)) {
        cat("Transition ", nm, "\n", sep = "")
        printCoefmat(x$coefficients[[nm]], digits = digits, P.values = TRUE,
                     has.Pvalue = TRUE, signif.stars = signif.stars,
                     cs.ind = c(1, 3), tst.ind = 6, na.print = "", ...)
        cat("\n")
    }
    if (!is.null(x$misclass)) print_misclass(x$misclass, digits = 4, ci = TRUE)
    if (!is.null(x$markers)) print_markers(x$markers, digits = 4, ci = TRUE)
    ll <- x$loglik
    if (length(ll) == 1) cat("Log-likelihood: ", format(ll), "\n")
    else cat("Log-likelihood: (", paste(format(ll), collapse = ", "),
             "), 2*diff: ", round(diff(ll), 1), " on ", x$df, " df\n", sep = "")
    cat("AIC: ", format(unname(x$AIC)), "\n", sep = "")
    if (!is.null(x$converged)) cat("Converged: ", x$converged, "\n", sep = "")
    cat(x$n[1], "observations,", x$n[2], "subjects\n")
    invisible(x)
}
