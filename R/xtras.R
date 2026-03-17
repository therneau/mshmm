# functions to go from B to coef, and coef to B
# B is a matrix such X%*% B gives all the linear predictors
#   coef is a vector of the parameters
# Within the iteration, fixed coeffficients are saved in B and coef only
#  contains the paramters to be maximized; as returned to the user, the
#  coefficients vector contains both fixed and non-fixed values.
# If fixed= F and B is present, only replace the non-fixed portions in B.
# I purposefully don't copy dimnames of coef over to B, it would just slow
#  it down. The coef/coefficients method, which is user facing, does otherwise.
coef.to.B <- function(coef, cmap, B, fixed=FALSE) {
    if (fixed || missing(B)) B <- matrix(0, nrow=nrow(cmap), ncol=ncol(cmap))
    if (fixed) B[cmap!=0] <- coef(abs(cmap))
    else B[cmap>0] <- coef[pmax(cmap,0)]                   
    # The B matrix might contain fixed coefficients, leave those in
    #  peace, cmap will have negative values there
    B[cmap>0] <- coef[pmax(cmap,0)]
    B
}

B.to.coef <- function(B, cmap, fixed=FALSE) {
    if (fixed) c2 <- abs(cmap) # copy all coefs
    else c2 <- ifelse(cmap>0, cmap, 0L) # ignore negatives in cmap

    # cmap might have the same integer twice
    ncoef <- max(c2)
    new <- double(ncoef)
    new <- double(ncoef)
    i1 <- unique(c2[c2>0])
    new[i1] <- B[match(i1, c2)]
    new
}       

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
    
hmm.control <- function(smallpos= 1e-3, debug= 0, scale=TRUE,
                         makecluster= .Platform$OS.type=="windows") {
    if (!is.numeric(smallpos) || length(smallpos) >1 || smallpos <=0)
        stop("smallpos must be a single value >0")
    if (!is.logical(scale)) stop("scale must be TRUE/FALSE")
    if (!is.logical(makecluster)) stop("makecluster must be TRUE/FALSE")
    if (is.logical(debug)) debug <- as.integer(debug)
    if (!is.numeric(debug) || debug != floor(debug)) 
        stop("debug must be TRUE/FALSE or an integer")

    list(smallpos= smallpos, debug=debug, scale=scale, makecluster=makecluster)
}
                         
