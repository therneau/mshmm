vcov.cmsh <- function (object, complete = TRUE, markers = FALSE, ...) 
{
    # markers=FALSE (default): rows/columns for the parameters of a
    #  misclassification matrix are left out, to match coef(object)
    if (is.null(object$var)) return(NULL)
    vmat <- object$var
    cf <- object$coefficients
    if (nrow(vmat) == length(cf)) {
        vname <- names(cf)
        dimnames(vmat) <- list(vname, vname)
    }
    keep <- rep(TRUE, nrow(vmat))
    if (!markers && length(object$epar) && nrow(vmat) == length(cf))
        keep[object$epar] <- FALSE
    if (!complete && nrow(vmat) == length(cf)) keep <- keep & !is.na(cf)
    vmat[keep, keep, drop = FALSE]
}
