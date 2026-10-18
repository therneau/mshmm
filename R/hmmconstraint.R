# The user provides constraints as special data sets, each row specifies
#  a linear predictor _lp_, a constraint id _cid_ and a constant _cwt_
#
hmmconstraint <- function(cdata, Terms, cmap, xtrans, B) {
    if (!inherits(cdata, "data.frame")) 
        stop("constraint or penalty must be a data frame containing ",
             " _lp__, _cid_, _cwt_ and data set variables")
    indx <- match(c("_lp__", "_cid_", "_cwt_"), names(cdata))
    if (any(is.na(indx))) 
        stop("constraint or penalty must be a data frame containing ",
             " _lp__, _cid_, and _cwt_ variables")
    cX <- model.matrix(Terms, data=cdata)  # this won't have _lp_, _cid_, _cwt_

    etaid <- match(cdata[,"_lp__"], colnames(cmap))
    if (any(is.na(etaid))) {
        bad <- unique(cdata[,"_lp__"][is.na(etaid)])
        stop("linear predictor not found in model: ", bad[1])
    }

    cid <- cdata[,"_cid_"]
    if (!is.numeric(cid) || any(is.na(cid)) || any(cid != floor(cid)) || 
        any(cid<1))
        stop("the constrast id _cid_ must be a postive integer")
    idcount <- table(cid)
    # I'm not so sure about this check. Perhaps someone does want a certain,
    #   single predicted value to be >=0 (?)
    #if (any(idcount) ==1)
    #    stop("contrasts must be between at least 2 predicted values")

    cwt <- cdata[, "_cwt_"]
    if (!is.numeric(cwt) || any(is.na(cwt)))
        stop("contrast weights _cwt_ must be numeric, and not missing")

    cid <- match(cid, unique(cid)) # make them 1, 2, ... for convenience
    if (missing(xtrans)) xtrans <- diag(nrow(cmap))
    if (!is.matrix(xtrans) || any(dim(xtrans) != nrow(cmap)))
        stop("xtrans has the wrong dimension")
    if (missing(B)) B <- matrix(0, nrow(cmap), ncol(cmap))

    cX <- cX %*% xtrans
    nparam <- length(B.to.coef(B, cmap))
    const <- matrix(0, max(cid), nparam + 1L)
    for (i in seq_len(max(cid))) {
        for (j in which(cid==i)) {
            pmap <- cmap[,etaid[j]]
            z <- cwt[j] * cX[j,]
            for (p in unique(pmap[pmap > 0]))
                const[i, 1L + p] <- const[i, 1L + p] + sum(z[pmap == p])
            fixed <- which(pmap < 0)
            if (length(fixed))
                const[i, 1L] <- const[i, 1L] +
                    sum(z[fixed] * B[fixed, etaid[j]])
        }
    }
    const
}    
