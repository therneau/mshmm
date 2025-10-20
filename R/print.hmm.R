print.hmm <- function(x, digits=max(options()$digits - 4, 3), ...) {
     if (!is.null(cl<- x$call)) {
	cat("Call:\n")
	dput(cl)
	cat("\n")
	}
   
     printCoefmat(x$beta, has.Pvalue=FALSE)
     cat("\n")
     loglik <- round(x$loglik, 2)
     cat("Log-likelihood: initial=", format(x$loglik[1]), 
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
