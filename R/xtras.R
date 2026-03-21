#  The transition matrix and response functions return derivatives with respect
# to the linear predictors eta, which then become derivatives wrt the parameter
# vector beta via an eta.beta function. The routine below creates that function,
# given a portion of cmap. 
#  The set of linear predictors in cmap is partitioned into transition matrix,
# then response functions, then initial state predictors; the three do not
# share parameters (coefficients) thus the separate functions.
#  If cmap[i,j] =k for k>0, then (partial eta_j)/(partial beta_k) = x_i, so
# there is a new transition matrix for each row of X. The pattern of nonzero
# values remains the same however.
# Fixed parameters are invisible to the maximizer, and don't appear in derivative
#  matrices
derivfun <- function(cmap) {
    # here is the function
    tfun <- function(x, neta, nbeta, index1, index2) {
        dmat <- matrix(0, neta, nbeta)
        dmat[index1] <- x[index2]
        dmat
    }
    # Now set default values for the last 4 arguments
    args <- as.list(formals(tfun))
    cx <- which(cmap>0)
    betas <- unique(cmap[cx])  # the betas of interest
    args$neta  <- ncol(cmap)
    args$nbeta <- length(betas)
    cx <- which(cmap>0)
    i <- col(cmap)[cx]
    j <- match(cmap[cx], betas)
    args$index1 <- (i-1)* ncol(cmap) + j
    args$index2 <- row(cmap)[cx]
    formals(tfun) <- args
    tfun
}
    
hmm.control <- function(smallpos= 1e-3, debug= 0, center=TRUE, scale=FALSE,
                        makecluster= .Platform$OS.type=="windows",
                        detail= FALSE) {
    if (!is.numeric(smallpos) || length(smallpos) >1 || smallpos <=0)
        stop("smallpos must be a single value >0")
    if (!is.logical(scale)) stop("scale must be TRUE/FALSE")
    if (!is.logical(makecluster)) stop("makecluster must be TRUE/FALSE")
    if (is.logical(debug)) debug <- as.integer(debug)
    if (!is.numeric(debug) || debug != floor(debug)) 
        stop("debug must be TRUE/FALSE or an integer")
    if (!(is.logical(detail))) stop("detail option must be TRUE/FALSE")

    list(smallpos= smallpos, debug=debug, center=center, scale=scale, 
         makecluster=makecluster, detail=detail)
}
                         
