# The coefficients/coef and print methods

coef.hmm <- function(object, matrix=FALSE, fixed=TRUE, ...) {
    cmap <- object$cmap
    if (matrix) {
        B <- coef.to.B(object$coefficients, cmap, fixed=fixed)
        dimnames(B) <- dimnames(cmap)
    }
    else if (fixed) object$coefficients
    else object$coefficients[cmap[cmap>0]]
}

print.hmm <- function(x, digits=max(options()$digits - 4, 3), ...) {
     if (!is.null(cl<- x$call)) {
	cat("Call:\n")
	dput(cl)
	cat("\n")
	}
   
     B <- coef(x, matrix=TRUE, fixed=TRUE)
     printCoefmat(B, has.Pvalue=FALSE)
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
