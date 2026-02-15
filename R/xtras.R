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

icmsh.control <- function(smallpos= 1e-3, debug= FALSE,
                         makecluster= .Platform$OS.type=="windows") {
    if (!is.numeric(smallpos) || length(smallpos) >1 || smallpos <=0)
        stop("smallpos must be a single value >0")
     if (!is.logical(makecluster)) stop("makecluster must be TRUE/FALSE")
    if (is.logical(debug)) debug <- as.integer(debug)
    if (!is.numeric(debug) || debug < 0) 
        stop("debug must be TRUE/FALSE or an integer >=0")

    list(smallpos= smallpos, debug=debug, makecluster=makecluster)
}
                         
