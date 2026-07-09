contrast <- function(fit, data, global=FALSE, weight, transition,
                     tol=sqrt(.Machine$double.eps)) {
    newx <- model.matrix(fit, data= data)
    if (!missing(transition)) {
        cmap <- fit$cmap
        if (is.numeric(transition)) {
            if (any(transition<1 | transition > ncol(cmap))) 
                stop("invalid transition number")
            k <- as.integer(transition)
        } else {
            k <- match(transition, colnames(cmap))
            if (any(is.na(k))) stop("unrecognized transition")
        }
        if (length(k)==1) j <- cmap[,k]  else stop("not done")
        vmat <- vcov(fit)[j,j, drop=FALSE]
        beta <- coef(fit)[j]
        if (any(cmap[,k]==0)) newx <- newx[, cmap[k]>0, drop=FALSE]
    }
    else {
        vmat <- vcov(fit)
        beta <- coef(fit)
    }
    
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
