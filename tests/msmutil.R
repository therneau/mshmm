#
# The msm package stores its coefficients in a somewha non-obvious (to me) order
#
# Create a cmap object as is found in coxph.
#
coefmat <- function(fit) {
    qq <- fit$Qmatrices
    nvar <- length(qq)-1
    indx <- which(qq[[1]] !=0)  # logbaseline = intercept
    vname <- names(qq)[-length(qq)]
    vname[1] <- "Intercept"
    beta <- matrix(0, nvar, length(indx))
    for (i in 1:(length(qq)-1)) { 
        beta[i,] <- (qq[[i]])[indx]
    }
    trans <- paste(row(qq[[1]])[indx], col(qq[[1]])[indx], sep=':')
    dimnames(beta) <- list(vname, trans)
    beta
}

# The fit$estimates object has a lot of zeros in it, as does the
#  fit$covmat object.  Chris leaves a placeholder for every rate*covariate pair
#  even if a given covariate is not used for that transtition.

# The msm object does not retain a terms object or a model.frame. 
# (There is a data$mf component, but it has class data.frame.  It may have
#  started out as a model frame, since it still retains a terms attribute.)
#   
# For my work, I need a model.matrix function. It requires a saved model with
#  terms and contrasts components. We know that the model matrix for msm
#  will be the same as a model.matrix from an lm model, so do a "fake" call
#  to lm.
model.matrix.msm <-  function(object, data=NULL, center=TRUE, ...) {
    # I don't need the response, but lm does
    tform <- `(time)` ~ x1   # dummy formula, x1 will be replaced
    tform[[3]] <-  object$covariates[[2]]
    lfit <- lm(tform, data=object$data$mf)
    if (is.null(data)) X <- model.matrix(lfit)
    else X <- model.matrix(delete.response(terms(lfit)), data=data, 
                      contrast= lfit$contrasts, xlev= lfit$xlevels)
    if (center) {
        means <- c("(Intercept)"= 0, attr(object$data$mm.cov, "means"))
        scale(X, center=means, scale=FALSE)
    } else X
}

anova.msm <- function(fit1, fit2) {
    loglik <- c(fit1$minus2loglik , fit2$minus2loglik)
    nvar1 <- sum(unlist(fit1$Qmatrices[-1]) !=0)
    nvar2 <- sum(unlist(fit2$Qmatrices[-1]) !=0)
    df <- nvar2 - nvar1
    c(loglik= -diff(loglik), df=df, 
      pchi= pchisq(-diff(loglik), df, lower.tail=FALSE))
}

# This function is used to get contrasts of coefficients.
# A missing part of R.  Used for survival data.
contrast <- function(fit, data, global=FALSE, weight, 
                     tol=sqrt(.Machine$double.eps), transition) {
    newx <- model.matrix(fit, data= data)
    if (global)  # test for all differences =0
        newx <- scale(newx[-1,], center=newx[1,], scale=FALSE)
    else {    
        if (missing(weight))  # subtract first row
            newx <- scale(newx, center=newx[1,], scale=FALSE) 
        else  { # use a weighted average
            wt <- weight/sum(weight)# use a weighted average
            if (length(wt) != nrow(newx)) stop("wrong length for weight")
            newx <- scale(newx, center= wt %*% newx , scale=FALSE)
        }
    }

    if (!missing(transition)) { # this is a multistate model
        if (inherits(fit, "coxphms")) {
            indx <- fit$cmap[, transition]
            beta <- c(0, coef(fit))[1L+ indx]
            vmat <- cbind(0, rbind(0, vcov(fit)))[indx+1L, indx +1L]
        }
        else if (inherits(fit, "msm")) {
            indx <- makemap(fit)[,transition]
            beta <- c(0, fit$estimates)[1L+ indx]
            vmat <- cbind(0, rbind(0, fit$covmat))[indx+1L, indx+1L]  
        }
        else if (inherits(fit, "hmm")) {
            indx <- fit$cmap[,transition]
            beta <- c(0, fit$coef)[1L+ indx]
            if (is.null(fit$fit$hessian)) vtemp <- solve(fit$fit$S2)
            else vtemp <- solve(fit$fit$hessian)
            vmat <- cbind(0, rbind(0, vtemp))[indx+1L, indx+1L] 
        }   
        else stop("transition not applicable")
    } else {
        beta <- coef(fit)
        vmat <- vcov(fit)
    }

    test <- drop(newx %*% beta)
    V    <- newx %*% vmat %*% t(newx)
    if (global) {
        stemp <- svd(V)
        nonzero <- (stemp$d > tol)
        ctemp <- test %*% stemp$u[,nonzero]
        chi <- ctemp %*% diag(1/stemp$d[nonzero]) %*% c(ctemp)
        c(chisq= drop(chi), df= sum(nonzero))
    }
    else {
        std <- sqrt(diag(V))
        z <- ifelse(std==0, 0, test/std)
        cbind(estimate=test, std.err=std, z= z)
    }
}


# older version of the above, specific to msm
msm.eta <- function(object, data=NULL, intercept=TRUE, weight) {   
    cmap <- makemap(object)     
    if (!intercept) cmap[1,] <- 0

    newX <- model.matrix.msm(object, data)
    ntran <- ncol(cmap)
    eta <- matrix(0, nrow(newX), ncol(cmap),
                  dimnames=list(NULL, colnames(cmap)))
    std <- eta
    zbeta <- c(0, object$estimates)
    zvar <- cbind(0, rbind(0, object$covmat))
    for (j in 1:ntran) {
        index <- 1L + cmap[,j]
        eta[,j] <- newX %*% zbeta[index]
        temp  <- (newX %*% zvar[index, index])*newX
        std[,j] <- sqrt(rowSums(temp))
    }

    if (!missing(weight)) {
        # recenter them
        wt <- c(weight/sum(weight))
        for (j in 1:ntran) {
            wmean <- sum(wt*eta[,j])
            eta[,j] <- eta[,j] - wmean
        }
    }                  
    list(fit=eta, se.fit=std)
}

lvcf <- function(id, x, time) {
    if (!missing(time)) indx <- order(id, time)
    else indx <- order(id)   

    for (i in seq(along=x)) {
        j <- indx[i]
        if (!is.na(x[j]) || i==1 || id[j]!= id[jlag]) current <- x[j]
        else x[j] <- current
        jlag <- j
    }
    x
}

nostutter <- function(id, x, censor=0) {
    # censor is the code to use for censoring, 
    # the output will have censor as the first code
    if (is.character(x) | is.numeric(x)) x <- as.factor(x)
    if (is.factor(x)) {
        newlev <- unique(c(censor, levels(x)))
        iscensor <-( x== censor)  # already marked as censored
        x <- as.integer(x)
        x[iscensor] <- 0
    } else stop("invalid variable type")        
    
    n <- length(id)
    if (length(x) != n) stop("wrong length for x or id")
    for (i in 1:n) {
        if (i==1 || id[i] != id[i-1]) {
            if (is.na(x[i])) current <- 0 else current <- x[i]
        } else if (!is.na(x[i])) {
            if (x[i]== current) x[i] <- 0
            else if (x[i] >0) current <- x[i]
        }
    }
    # if censor were level 3 of 5 in input, then the unique x at
    #  this point would be 0, 1, 2, 4, 5
    factor(x, sort(unique(x)), newlev)  
}
