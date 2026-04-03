# The coefficients/coef and print methods

coef.cmsh <- function(object, matrix=FALSE, 
                     fixed=TRUE, matrix1= matrix, matrix2= matrix, 
                     matrix3= matrix, ...) {
    cmap <- object$cmap
    if (missing(matrix)) matrix <- (matrix1 | matrix2 | matrix3)
    if (matrix) {
        B <- coef.to.B(object$coefficients, cmap, fixed=fixed)
        dimnames(B) <- dimnames(cmap)
        keep <- rep(c(matrix1, matrix2, matrix3), object$nlp)
        B[,keep, drop=FALSE]
    }
    else if (fixed) object$coefficients
    else object$coefficients[sort(unique(cmap[cmap>0]))]
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

print.cmsh <- function(x, digits=max(options()$digits - 4, 3), ...) {
     if (!is.null(cl<- x$call)) {
	cat("Call:\n")
	dput(cl)
	cat("\n")
	}
   
     B <- coef(x, matrix=TRUE, fixed=TRUE)
     printCoefmat(B, has.Pvalue=FALSE)
     cat(" States: ", paste(paste(seq(along.with=x$states), x$states, sep='= '),
                            collapse=", "), '\n')
     cat("\n")
     loglik <- round(x$loglik, 2)
     if (length(loglik)==1)
         cat("Log-likelihood: ", format(x$loglik[1]), "\n")
     else cat("Log-likelihood: initial=", format(x$loglik[1]), 
              " final=", format(x$loglik[2]), "\n")

     cat(x$n[1], "observations", x$n[2], "subjects\n")
     if (length(x$na.action)) {
         cat("   ", length(x$na.action), "observations removed")
         temp <- paste(x$removed[x$removed>0], names(x$removed[x$removed>0]))
         cat(" (", paste(temp, collapse=", "), ")\n", sep='')
     }  
     else cat("\n")
     invisible(x)
}
