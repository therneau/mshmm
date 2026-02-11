# functions to go from B to coef, and coef to B
# B is a matrix st X%*% B gives all the linear predictors
# coef is a vector of all the estimated parameters, it is what the maximizer
#  works on
coef.to.B <- function(coef, cmap, B) {
    # The B matrix might contain fixed coefficients, leave those in
    #  peace, cmap will have negative values there
    B[cmap>0] <- coef[pmax(cmap,0)]
    B
}

B.to.coef <- function(B, cmap) {
    # fixed coefs don't transmit to coef
    ncoef <- max(cmap)  # negative for fixed ones
    new <- double(ncoef)
    i1 <- unique(cmap[cmap>0])
    new[i1] <- B[match(i1, cmap)]
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
                         
