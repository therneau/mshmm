# The user provides constraints as special data sets, each row specifies
#  a linear predictor _eta_, a constraint id _cid_ and a constant _cval_
#
hmmconstraint <- function(cdata, Terms, cmap) {
    if (!inherits(cdata, "data.frame")) 
        stop("constraint or penalty must be a data frame containing ",
             " _eta_, _cid_, and _cval_ variables")
    indx <- match(c("_eta_", "_cid_", "_cval_"), names(cdata))
    if (any(is.na(indx)) 
        stop("constraint or penalty must be a data frame containing ",
             " _eta_, _cid_, and _cval_ variables")
    X <- model.matrix(Terms, data=cdata)  # this won't have _eta_, cid, cval

    etaid <- match(cdata[,"_eta_"], colnames(cmap))
    if (any(is.na(etaid))) {
        bad <- unique(cdata[,"_eta_"][is.na(etaid)])
        stop("linear predictor not found in model: ", bad[1])
    }

    cid <- cdata(,"_cid_")
    if (!is.numeric(cid) || any(is.na(cid)) || any(cid != floor(cid)) || 
        any(cid<1))
        stop("the constrast id _cid_ must be a postive integer")
    idcount <- table(cid)
    if (any(idcount) ==1)
        stop("contrasts must be beteen at least 2 predicted values")

    cval <- cdata[, "_cval_"]
    if (!is.numeric(cval) || any(is.na(cval)))
        stop("contrast weights _cval_ must be numeric, not missing")

    cid <- match(cid, unique(cid)) # make them 1, 2, ... for convenience
    cmat <- matrix(0, max(cid), nrow(cmat))
    for (i in 1:max(cid)) {
        for (j in which(cid==i)) {
            xvar <- which(cmap[,eta[j]] > 0) # variables for this lp
            cmat[i,xvar] <- cmat[i, xvar] <- cval[j]* X[j, xvar]
    }
    cmat
}    
