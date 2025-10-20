#
# Given an error matrix, compute the coefficients that would
# lead to it.  Assume they are listed in usual R matrix order
#
error_to_coef <- function(emat) {
    if (nrow(emat) != ncol(emat)) stop("error matrix must be square")
    if (any(emat <0 | emat >1)) stop("elements must be between 0 and 1")
    if (any(rowSums(emat) > 1)) stop("row sum is > 1")

    nstate <- ncol(emat)
    nonzero <- which(emat >0 & diag(nstate)==0)  # coefs for each of these
    diag(emat) <- 1 + diag(emat) - rowSums(emat)  #make rows sum to 1
    
    lmat <- log(ifelse(emat>0, emat, 1))
    coef <- lmat - rep(diag(lmat), nstate)
    coef[nonzero]
}
