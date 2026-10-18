# The coefficients/coef and print methods

coef.cmsh <- function(object, matrix=FALSE, markers=TRUE,
                     fixed=TRUE, matrix1= matrix, matrix2= matrix, 
                     matrix3= matrix, ...) {
    cmap <- object$cmap
    if (missing(matrix)) matrix <- (matrix1 | matrix2 | matrix3)
    if (matrix) {
        B <- coef.to.B(object$coefficients, cmap, fixed=fixed)
        dimnames(B) <- dimnames(cmap)
        keep <- rep(c(matrix1, matrix2, matrix3), object$nlp)
        if (!markers && object$nlp[2] > 0L)
            keep[object$nlp[1] + seq_len(object$nlp[2])] <- FALSE
        B[,keep, drop=FALSE]
    }
    else {
        ids <- if (fixed) seq_along(object$coefficients)
               else sort(unique(cmap[cmap>0]))
        if (!markers && length(object$marker.ids))
            ids <- setdiff(ids, object$marker.ids)
        object$coefficients[ids]
    }
}

summary.cmsh <- function(object, conf.int=0.95, markers=TRUE, ...) {
    estimate <- object$coefficients
    se <- if (is.null(object$var)) rep(NA_real_, length(estimate))
          else sqrt(pmax(diag(object$var), 0))
    fixed_ids <- unique(abs(object$cmap[object$cmap < 0L]))
    if (length(fixed_ids)) se[fixed_ids] <- NA_real_
    zcrit <- qnorm(1 - (1 - conf.int) / 2)
    table <- cbind(
        estimate=estimate,
        `exp(estimate)`=exp(estimate),
        std.error=se,
        lower=estimate - zcrit * se,
        upper=estimate + zcrit * se,
        z=estimate / se,
        p=2 * pnorm(-abs(estimate / se))
    )
    if (length(object$marker.ids))
        table <- table[-object$marker.ids, , drop=FALSE]
    value <- list(call=object$call, coefficients=table,
                  loglik=object$loglik, n=object$n,
                  markers=if (markers) object$markers else NULL,
                  misclass=if (markers) object$misclass else NULL,
                  AIC=-2 * object$loglik[length(object$loglik)] +
                      2 * length(coef(object, fixed=FALSE, markers=TRUE)),
                  converged=object$fit$converged,
                  conf.int=conf.int)
    class(value) <- "summary.cmsh"
    value
}

print.summary.cmsh <- function(x, digits=max(1L, getOption("digits") - 3L), ...) {
    if (!is.null(x$call)) {
        cat("Call:\n")
        dput(x$call)
        cat("\n")
    }
    printCoefmat(x$coefficients, digits=digits, P.values=TRUE,
                 has.Pvalue=TRUE, ...)
    if (length(x$markers)) print_markers(x$markers, digits=digits, ci=TRUE)
    if (length(x$misclass)) print_misclass(x$misclass, digits=digits, ci=TRUE)
    cat("AIC: ", format(unname(x$AIC)), "\n", sep="")
    if (!is.null(x$converged)) cat("Converged: ", x$converged, "\n", sep="")
    invisible(x)
}

# functions to go from B to coef, and coef to B
# B is a matrix such X%*% B gives all the linear predictors
#   coef is a vector of the parameters
# Within the iteration, fixed coeffficients are saved in B and coef only
#  contains the paramters to be maximized; as returned to the user, the
#  coefficients vector contains both fixed and non-fixed values.
# If fixed= F and B is present, only replace the non-fixed portions in B.
# I purposefully don't copy dimnames of coef over to B, it would just slow
#  X%*%B. The coef/coefficients method, which is user facing, does otherwise.
coef.to.B <- function(coef, cmap, B, fixed=FALSE) {
    if (fixed || missing(B)) B <- matrix(0, nrow=nrow(cmap), ncol=ncol(cmap))
    if (fixed) B[cmap!=0] <- coef[abs(cmap)]
    else B[cmap>0] <- coef[pmax(cmap,0)]                   
    # The B matrix might contain fixed coefficients, leave those in
    #  peace, cmap will have negative values there
    B[cmap>0] <- coef[pmax(cmap,0)]
    B
}

B.to.coef <- function(B, cmap, fixed=FALSE) {
    if (fixed) c2 <- abs(cmap) # copy all coefs
    else c2 <- ifelse(cmap>0, cmap, 0L) # ignore negatives in cmap

    # cmap might have the same integer twice
    ncoef <- max(c2)
    new <- double(ncoef)
    new <- double(ncoef)
    i1 <- unique(c2[c2>0])
    new[i1] <- B[match(i1, c2)]
    new
}       

print.cmsh <- function(x, digits=max(1L, getOption("digits")- 3L), 
                       signif.stars= FALSE, ...) {
     if (!is.null(cl<- x$call)) {
	cat("Call:\n")
	dput(cl)
	cat("\n")
	}
   
     extrarows <- NULL
     if (is.null(x$var)) {
         # when iter=0 is used derivatives and variance are not computed
         print(coef(x, matrix=TRUE, markers=FALSE), digits=digits)
         cat('\n')
     }
     else {
         coef <- x$coefficients
         se <- sqrt(diag(x$var))
         fixed_ids <- unique(abs(x$cmap[x$cmap < 0L]))
         if (length(fixed_ids)) se[fixed_ids] <- NA_real_
         tmp <- cbind(coef, exp(coef), se, coef/se,
                      pchisq((coef/se)^2, 1, lower.tail=FALSE))
         dimnames(tmp) <- list(names(coef), c("coef", "exp(coef)",
                                              "se(coef)", "z", "p"))
         # print it group by group
         # lazy: I don't want to type x$cmap many times
         #  remove transitions with no covariates
         cmap <- x$cmap
         if (length(x$ecols)) cmap <- cmap[, -x$ecols, drop=FALSE]
         cmap <- cmap[, colSums(cmap != 0) > 0, drop=FALSE]
         cname <- colnames(cmap)
         printed <- rep(FALSE, length(cname))
         rname <- rownames(cmap)
         for (i in 1:length(cname)) {
             # if multiple colums of tmat are identical, only print that
             #  set of coefficients once
             if (!printed[i]) {
                 j <- apply(cmap, 2, function(x) all(x == cmap[,i])) 
                 printed[j] <- TRUE
                 tlab <- c(paste(cname[j], collapse=", "), "") #label for block
                 if (nchar(tlab[1]) > 20) {     
                     k <- length(extrarows) +1L
                     extrarows <- c(extrarows, 
                                    paste(LETTERS[k], tlab, sep=': '))
                     tlab <- c(LETTERS[k], "") # new block label
                 }
                 ids <- abs(cmap[, i])
                 use <- which(ids > 0L)
                 tmp2 <- tmp[ids[use],, drop=FALSE]
                 names(dimnames(tmp2)) <- tlab
                 # restore character row names
                 rownames(tmp2) <- rname[use]
                 printCoefmat(tmp2, digits=digits, P.values=TRUE, 
                              has.Pvalue=TRUE, signif.stars=signif.stars, ...)
                 cat("\n")   
             }
         }
     }

     cat(" States: ", paste(paste(seq(along.with=x$states), x$states, 
                                  sep='= '), collapse=", "), '\n')
     if (!is.null(extrarows)) cat(extrarows, sep='\n')
     cat("\n")
     if (length(x$markers)) print_markers(x$markers, digits=digits)
     if (length(x$misclass)) print_misclass(x$misclass, digits=digits)
     loglik <- round(x$loglik, 2)
     if (length(loglik)==1)
         cat("Log-likelihood: ", format(x$loglik[1]), "\n")
     else {
         ldiff <- diff(x$loglik)
         df <- length(unique(x$cmap[x$cmap!=0]))
         cat("Log-likelihood: (", paste(format(x$loglik), collapse=", "),
             "), 2*diff:",
             round(ldiff,1), "on", df, "df\n")
     }

     cat(x$n[1], "observations,", x$n[2], "subjects\n")
     if (length(x$na.action)) {
         cat("   ", length(x$na.action), "observations removed")
         temp <- paste(x$removed[x$removed>0], names(x$removed[x$removed>0]))
         cat(" (", paste(temp, collapse=", "), ")\n", sep='')
     }  
     else cat("\n")
     invisible(x)
}
