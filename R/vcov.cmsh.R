vcov.cmsh <- function (object, complete = TRUE, ...) 
{
    if (is.null(object$var)) return(NULL)
    vmat <- object$var
    vname <- names(object$coefficients)
    dimnames(vmat) <- list(vname, vname)
    if (!complete && any(is.na(coef(object)))) {
        keep <- !is.na(coef(object))
        vmat[keep, keep, drop = FALSE]
    }
    else vmat
}
