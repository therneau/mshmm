#  The transition matrix and response functions return derivatives with respect
# to the linear predictors eta, which then become derivatives wrt the parameter
# vector beta via a $Z$ matrix (see the code vignette). The function below 
# creates the functions (eta.beta1 etc) that create those matrices, based on the
# relevant portion of cmap.
#  The set of linear predictors in cmap is partitioned into transition matrix,
# then response functions, then initial state predictors; the three do not
# share parameters (coefficients) thus the separate functions.
#  If cmap[i,j] =k for k>0, then (partial eta_j)/(partial beta_k) = x_i, i.e.,
#  the j,k element of dmat is x[i], dmat = the Z of my vignette.
# There is a new transition matrix for each row of X. The pattern of nonzero
# values remains the same, however.
# Fixed parameters are invisible to the maximizer, and don't appear in 
#  derivative matrices
derivfun <- function(cmap) {
    # here is the function
    # x is a row of the X matrix (given patient and time)
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
    i <- row(cmap)[cx]
    j <- col(cmap)[cx]
    k <- match(cmap[cx], betas)
    args$index1 <- (k-1)* ncol(cmap) + j #points to [j,k] element
    args$index2 <- i
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
                         
