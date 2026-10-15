# The coefficients/coef and print methods

coef.cmsh <- function(object, matrix=FALSE, 
                     fixed=TRUE, matrix1= matrix, matrix2= matrix, 
                     matrix3= matrix, markers=FALSE, ...) {
    # markers=FALSE (default): the parameters of the markers are not part of
    #  the result; they live in object$misclass (discrete markers) and
    #  object$markers (gaussian, logistic, beta, gamma), on the scale set by
    #  each distribution's format() function in response.R.
    #  markers=TRUE returns every parameter, on the linear predictor scale.
    cmap <- object$cmap
    ecols <- if (markers) NULL else object$ecols
    epar  <- if (markers) NULL else object$epar
    if (missing(matrix)) matrix <- (matrix1 | matrix2 | matrix3)
    if (matrix) {
        B <- coef.to.B(object$coefficients, cmap, fixed=fixed)
        dimnames(B) <- dimnames(cmap)
        keep <- rep(c(matrix1, matrix2, matrix3), object$nlp)
        if (length(ecols)) keep[ecols] <- FALSE
        B[,keep, drop=FALSE]
    }
    else {
        cf <- object$coefficients
        if (fixed) use <- seq_along(cf)
        else use <- sort(unique(cmap[cmap>0]))
        cf[setdiff(use, epar)]
    }
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
   
     coef <- x$coefficients
     se <- sqrt(diag(x$var))
     tmp <- cbind(coef, exp(coef), se, coef/se,
                  pchisq((coef/se)^2, 1, lower.tail=FALSE))
     dimnames(tmp) <- list(names(coef), c("coef", "exp(coef)",
                                          "se(coef)", "z", "p"))
     # print it group by group
     # lazy: I don't want to type x$cmap many times
     #  remove transitions with no covariates
     cmap <- x$cmap
     if (length(x$ecols)) cmap <- cmap[, -x$ecols, drop=FALSE]  # not the marker lp
     cmap <- cmap[, colSums(cmap) > 0, drop=FALSE]
     cname <- colnames(cmap)
     printed <- rep(FALSE, length(cname))
     rname <- rownames(cmap)
     extrarows <- NULL
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
                                paste(LETTERS[k], tlab[1], sep=": "))
                 tlab <- c(LETTERS[k], "") # new block label
             }
             tmp2 <- tmp[cmap[,i],, drop=FALSE]
             names(dimnames(tmp2)) <- tlab
             # restore character row names
             rownames(tmp2) <- rname[cmap[,i]>0]
             printCoefmat(tmp2, digits=digits, P.values=TRUE, 
                          has.Pvalue=TRUE, signif.stars=signif.stars, ...)
             cat("\n")   
         }
     }

     cat(" States: ", paste(paste(seq(along.with=x$states), x$states, 
                                  sep='= '), collapse=", "), '\n')
     if (!is.null(extrarows)) cat(extrarows, sep='\n')
     cat("\n")
     # marker results, each on the scale given by its distribution's format()
     if (!is.null(x$misclass)) print_misclass(x$misclass, digits=digits)
     if (!is.null(x$markers)) print_markers(x$markers, digits=digits)
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
