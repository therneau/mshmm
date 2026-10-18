# User-facing marker summaries.  Marker distributions are fit on linear-
# predictor scales, but reported on their natural scales at the reference
# covariate pattern (the intercept row of the model matrix).

.marker_coefficient_table <- function(columns, B, cmap, vcov,
                                      conf.int=0.95) {
    cells <- which(cmap[, columns, drop=FALSE] != 0L, arr.ind=TRUE)
    if (!nrow(cells)) return(NULL)
    cells <- cells[order(cells[, 2L], cells[, 1L]), , drop=FALSE]
    ids <- abs(cmap[cbind(cells[, 1L], columns[cells[, 2L]])])
    free <- cmap[cbind(cells[, 1L], columns[cells[, 2L]])] > 0L
    estimate <- B[cbind(cells[, 1L], columns[cells[, 2L]])]
    se <- rep(NA_real_, length(ids))
    if (!is.null(vcov)) {
        usable <- free & ids <= nrow(vcov)
        se[usable] <- sqrt(pmax(diag(vcov)[ids[usable]], 0))
    }
    z <- qnorm(1 - (1 - conf.int) / 2)
    label <- paste(colnames(cmap)[columns[cells[, 2L]]],
                   rownames(cmap)[cells[, 1L]])
    data.frame(estimate=estimate, std.error=se,
               lower=estimate - z * se, upper=estimate + z * se,
               row.names=make.unique(label), check.names=FALSE)
}

.discrete_probability <- function(response, eta) {
    emap <- response$emap
    probability <- matrix(0, nrow(emap), ncol(emap),
                          dimnames=list(response$latent_levels,
                                        response$marker_levels))
    for (i in seq_len(nrow(emap))) {
        logits <- rep(-Inf, ncol(emap))
        logits[emap[i, ] == -1L] <- 0
        positive <- emap[i, ] > 0L
        logits[positive] <- eta[emap[i, positive]]
        logits <- logits - max(logits)
        value <- exp(logits)
        probability[i, ] <- value / sum(value)
    }
    probability
}

.eta_vcov <- function(ids, free, vcov) {
    value <- matrix(0, length(ids), length(ids))
    if (is.null(vcov)) return(value + NA_real_)
    for (i in seq_along(ids)) for (j in seq_along(ids)) {
        if (free[i] && free[j] && ids[i] <= nrow(vcov) && ids[j] <= ncol(vcov))
            value[i, j] <- vcov[ids[i], ids[j]]
    }
    value
}

marker.output <- function(response, marker, B, cmap, nlp, coefficients, vcov,
                          conf.int=0.95) {
    columns <- nlp[1] + response$e2map
    ids <- sort(unique(abs(cmap[, columns, drop=FALSE])))
    ids <- ids[ids > 0L]
    ctable <- .marker_coefficient_table(columns, B, cmap, vcov, conf.int)
    z <- qnorm(1 - (1 - conf.int) / 2)

    if (response$family == "discrete") {
        if (response$fixed) {
            probability <- response$probability
            standard_error <- lower <- upper <- probability + NA_real_
        } else {
            reference <- unname(B[1L, columns])
            probability <- .discrete_probability(response, reference)
            parameter_ids <- abs(cmap[1L, columns])
            free <- cmap[1L, columns] > 0L
            eta_var <- .eta_vcov(parameter_ids, free, vcov)
            gradient <- array(0, c(dim(probability), length(reference)))
            h <- 1e-6
            for (j in seq_along(reference)) {
                delta <- numeric(length(reference)); delta[j] <- h
                gradient[, , j] <-
                    (.discrete_probability(response, reference + delta) -
                     .discrete_probability(response, reference - delta)) / (2 * h)
            }
            standard_error <- probability
            for (i in seq_len(nrow(probability)))
                for (j in seq_len(ncol(probability))) {
                    g <- gradient[i, j, ]
                    standard_error[i, j] <-
                        sqrt(pmax(drop(g %*% eta_var %*% g), 0))
                }
            if (is.null(vcov)) standard_error[,] <- NA_real_
            logit_se <- standard_error / (probability * (1 - probability))
            lower <- plogis(qlogis(probability) - z * logit_se)
            upper <- plogis(qlogis(probability) + z * logit_se)
            invalid <- !is.finite(logit_se)
            lower[invalid] <- upper[invalid] <- NA_real_
        }
        return(list(marker=marker, family="discrete", fixed=response$fixed,
                    matrix=probability, probabilities=probability,
                    se=standard_error, lower=lower, upper=upper,
                    coefficients=ctable, coefficient_ids=ids,
                    covariates=length(columns) > 0L &&
                        any(cmap[-1L, columns, drop=FALSE] != 0L),
                    conf.int=conf.int))
    }

    reference <- if (response$fixed) response$initial else unname(B[1L, columns])
    parameters <- matrix(reference, nrow=response$npeak, ncol=2L, byrow=TRUE)
    transformed <- response$family %in% c("beta", "gamma")
    if (transformed) {
        parameters <- exp(parameters)
        colnames(parameters) <- if (response$family == "beta")
            c("shape1", "shape2") else c("shape", "rate")
    } else {
        parameters[, 2L] <- exp(parameters[, 2L])
        colnames(parameters) <- c("mean", "std")
    }
    rownames(parameters) <- response$latent_levels

    standard_error <- lower <- upper <- parameters + NA_real_
    if (!response$fixed && length(columns) && !is.null(vcov)) {
        parameter_ids <- abs(cmap[1L, columns])
        free <- cmap[1L, columns] > 0L
        linear_se <- rep(NA_real_, length(columns))
        usable <- free & parameter_ids <= nrow(vcov)
        linear_se[usable] <- sqrt(pmax(diag(vcov)[parameter_ids[usable]], 0))
        deriv <- if (transformed) c(t(parameters)) else
            c(t(cbind(rep(1, nrow(parameters)), parameters[, 2L])))
        standard_error[,] <- matrix(deriv * linear_se, ncol=2L, byrow=TRUE)
        lo_lp <- reference - z * linear_se
        hi_lp <- reference + z * linear_se
        if (transformed) {
            lower[,] <- matrix(exp(lo_lp), ncol=2L, byrow=TRUE)
            upper[,] <- matrix(exp(hi_lp), ncol=2L, byrow=TRUE)
        } else {
            lower[,] <- matrix(lo_lp, ncol=2L, byrow=TRUE)
            upper[,] <- matrix(hi_lp, ncol=2L, byrow=TRUE)
            lower[, 2L] <- exp(lower[, 2L])
            upper[, 2L] <- exp(upper[, 2L])
        }
    }
    list(marker=marker, family=response$family,
         distribution=response$family, fixed=response$fixed,
         parameters=parameters, se=standard_error,
         lower=lower, upper=upper, coefficients=ctable,
         coefficient_ids=ids,
         covariates=length(columns) > 0L &&
             any(cmap[-1L, columns, drop=FALSE] != 0L),
         conf.int=conf.int)
}

# Extract a fitted categorical emission (misclassification) matrix.
misclass <- function(fit, marker,
                     what=c("matrix", "se", "lower", "upper",
                            "coefficients", "all")) {
    if (!inherits(fit, "cmsh")) stop("fit must be a cmsh object")
    what <- match.arg(what)
    value <- fit$misclass
    if (!length(value)) return(NULL)
    if (missing(marker)) {
        if (length(value) > 1L && what != "all")
            return(lapply(value, `[[`, what))
        marker <- 1L
    }
    if (what == "all") value[[marker]] else value[[marker]][[what]]
}

print_markers <- function(markers, digits=4L, ci=FALSE) {
    for (name in names(markers)) {
        item <- markers[[name]]
        cat("Marker ", name, " (", item$family,
            if (item$fixed) ", fixed" else ", estimated",
            "): parameters by latent level\n", sep="")
        if (!ci || item$fixed) print(round(item$parameters, digits))
        else {
            table <- do.call(rbind, lapply(seq_len(ncol(item$parameters)), function(j) {
                value <- cbind(item$parameters[, j], item$se[, j],
                               item$lower[, j], item$upper[, j])
                rownames(value) <- paste(rownames(item$parameters),
                                         colnames(item$parameters)[j])
                value
            }))
            colnames(table) <- c("estimate", "se", "lower", "upper")
            print(round(table, digits))
        }
        cat("\n")
    }
    invisible(markers)
}

print_misclass <- function(items, digits=4L, ci=FALSE) {
    for (name in names(items)) {
        item <- items[[name]]
        cat("Misclassification matrix for ", name,
            if (item$fixed) " (fixed)" else " (estimated)",
            ": P(observed | latent level)\n", sep="")
        print(round(item$matrix, digits))
        if (ci && !item$fixed && !is.null(item$coefficients)) {
            cat("\nLinear-predictor coefficients:\n")
            print(round(item$coefficients, digits))
        }
        cat("\n")
    }
    invisible(items)
}
