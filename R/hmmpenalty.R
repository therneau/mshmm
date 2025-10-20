# Create a constraint matrix
#  These will be contrasts across linear predictors, i.e., the beta
#  component of a fit, which are then translated into effects on the
#  fitted parameters.
#  
hmmpenalty <- function(fit, data, lp, contrast) {
    if (!inherits(fit, "hmm")) 
        stop("the fit argument must be the result of an hmm call")
    cmap <- fit$cmap
    if (is.numeric(lp)) {
        if (any(lp != floor(lp) | lp < 1)) stop ("lp must be integers > 0")
        if (any(lp > ncol(cmap))) stop("lp is out of range")
        indx1 <- as.integer(lp)
    }
    else {
        indx1 <- match(lp, colnames(cmap))
        if (any(is.na(indx1))) 
            stop("element ", which(is.na(indx1)), 
                 " of lp does not match a linear predictor name")
    }

    # create the starting coefficients
    X <- model.matrix(delete.response(terms(fit)), data, xlev= fit$xlevels, 
                         contrast.arg = fit$contrast)

    if (length(indx1) != nrow(X))
        stop("number of rows in the data set does not match length of lp")
    
    # Now create the contrast array, which will be an array with
    #  first two dimensions equal to the beta matrix within hmm, and one
    #  slice per contrast.
    ntest <- nrow(contrast)
    if (!is.matrix(contrast)) 
        stop("contrast argument must be a matrix")
    cmat <- array(0., dim=c(dim(fit$cmap), ntest))
                  
    if (ncol(contrast) != length(lp)) {
        if (ncol(contrast) != 2)
            stop("contrast matrix has the wrong number of columns")
        # alternate form, turn it into a matrix of contrasts
        temp <- contrast
        contrast <- matrix(0., ntest, length(lp))
        for (i in 1:ntest) 
            contrast[i, temp[i,]] <- c(1, -1)
    }
    
    for (i in 1:nrow(contrast)) {
        for (j in 1:length(lp)) 
            cmat[,indx1[j],i] <- cmat[,indx1[j],i] + contrast[i,j]*X[j,]
    }

    # translate this to a matrix with one row per constraint and one
    #  column per parameter
    if (any(fit$cmap ==0)) {
        dmat <- matrix(0., nrow(contrast), 1 + length(fit$coefficients))
        for (i in 1:nrow(contrast)) {
            dmat[i,] <- tapply(c(cmat[,,i]), c(fit$cmap), sum)
        }
        if (!all.equal(dmat[,1], rep(0, nrow(contrast))))
            warning("some contrasts involved fixed parameters")
        dmat[,-1]
    }
    else {  # this is unlikely to ever happen
        dmat <- matrix(0., nrow(contrast), length(fit$coefficients))
        for (i in 1:nrow(contrast)) {
            dmat[i,] <- tapply(c(cmat[,,i]), c(fit$cmap), sum)
        }
        dmat
    }
}

        
# The function to create tests of two rates against one another.
# It has the exact same arguments as hmmpenalty, and 90% the same code.
# The variance portion is similar to predict.hmm
# gsolve and qform are stolen from the survival package
gsolve <- function(mat, y, eps=sqrt(.Machine$double.eps)) {
    # solve using a generalized inverse
    # this is very similar to the ginv function of MASS
    temp <- svd(mat, nv=0)
    dpos <- (temp$d > max(temp$d[1]*eps, 0))
    dd <- ifelse(dpos, 1/temp$d, 0)
    # all the parentheses save a tiny bit of time if y is a vector
    if (all(dpos)) x <- drop(temp$u %*% (dd*(t(temp$u) %*% y)))
    else if (!any(dpos)) x <- drop(temp$y %*% (0*y)) # extremely rare
    else x <-drop(temp$u[,dpos] %*%(dd[dpos] * (t(temp$u[,dpos, drop=FALSE]) %*% y)))
    attr(x, "df") <- sum(dpos)
    x
}

qform <- function(var, beta) { # quadratic form b' (V-inverse) b
    temp <- gsolve(var, beta)
    list(test= sum(beta * temp), df=attr(temp, "df"))
}

hmmcontrast <- function(fit, data, lp, contrast) {    
    if (!inherits(fit, "hmm")) 
        stop("the fit argument must be the result of an hmm call")
    cmap <- fit$cmap
    if (!is.null(fit$fit$hessian)) hessian <- fit$fit$hessian
    else stop("hessian not found")

    if (is.numeric(lp)) {
        if (any(lp != floor(lp) | lp < 1)) stop ("lp must be integers > 0")
        if (any(lp > ncol(cmap))) stop("lp is out of range")
        indx1 <- as.integer(lp)
    }
    else {
        indx1 <- match(lp, colnames(cmap))
        if (any(is.na(indx1))) 
            stop("element ", which(is.na(indx1)), 
                 " of lp does not match a linear predictor name")
    }

    # create the starting coefficients
    X <- model.matrix(delete.response(terms(fit)), data, xlev= fit$xlevels, 
                         contrast.arg = fit$contrast)

    if (length(indx1) != nrow(X))
        stop("number of rows in the data set does not match length of lp")
    
    # Create or verify the contrast matrix between linear predictors
    ntest <- nrow(contrast)
    if (!is.matrix(contrast)) 
        stop("contrast argument must be a matrix")
                  
    if (ncol(contrast) != length(lp)) {
        if (ncol(contrast) != 2)
            stop("contrast matrix has the wrong number of columns")
        # alternate form, turn it into a matrix of contrasts
        temp <- contrast
        contrast <- matrix(0., ntest, length(lp))
        for (i in 1:ntest) 
            contrast[i, temp[i,]] <- c(1, -1)
        dimnames(contrast) <- list(paste(lp[temp[,1]], lp[temp[,2]], sep=' - '),
                                   NULL)
    }

    # The tests and variance involve fit$coef, any fixed coefficients
    #  found in fit$beta do not play a role.  
    # If a column of cmap has two elements the same that's just wierd (two
    #  coefs with the same coefficient), but the double loop allows for it.
    Z <- matrix(0, nrow(X), length(fit$coef))
    for (i in 1:nrow(X)) {
        for (k in unique(fit$cmap[,indx1[i]])) {
            if (k>0) {
                indx <- which(fit$cmap[,indx1[i]] ==k)
                # length(indx) >1 is very odd -- two covariates forced to
                #   have the same coef.  It will likely never happen, but...
                if (length(indx) ==1) Z[i,k] <- X[i,indx]
                else Z[i,k] <- sum(X[i, indx])
            }
        }
    }
    rownames(Z) <- lp

    # The estimated linear predictors and their variance matrix
    estimate  <- Z %*% fit$coef
    varest    <- Z %*% gsolve(hessian, t(Z))
    tests <- cbind(test = drop(contrast %*% estimate),
                   std = drop(sqrt(diag(contrast %*% varest %*% t(contrast)))))
    tests <- cbind(tests, p= 2*pnorm(-abs(tests[,1]/tests[,2])))
    list(estimate= drop(estimate), var= varest, test=tests)
}
