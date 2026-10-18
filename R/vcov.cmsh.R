vcov.cmsh <- function (object, complete = TRUE, markers=TRUE, ...)
{
    if (is.null(object$var)) return(NULL)
    vmat <- object$var
    vname <- names(object$coefficients)
    dimnames(vmat) <- list(vname, vname)
    keep <- seq_len(nrow(vmat))
    if (!markers && length(object$marker.ids))
        keep <- setdiff(keep, object$marker.ids)
    if (!complete) keep <- keep[!is.na(object$coefficients[keep])]
    vmat[keep, keep, drop=FALSE]
}
