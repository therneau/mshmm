# functions to go from B to coef, and coef to B
# B is a matrix st X%*% B gives all the linear predictors
# coef is a vector of all the estimated parameters, it is what the maximizer
#  works on
coef.to.B <- function(coef, cmap, mmap, B) {
    # The B matrix might contain fixed coefficients, leave those in
    #  peace
    if (is.null(mmap)) {
        # B and cmaps are the same size
        B[cmap>0] <- coef(cmap)
        B
    } else {
        temp1 <- B[,1:ncol(cmap)]
        temp1[cmap>0] <- coef(cmap)
        temp2 <- B[,-(1:ncol(cmap))]
        temp2[mmap>0] <- coef[mmap]
        cbind(temp1, temp2)
    }
}
B.to.coef <- function(B, cmap, mmap) {
    # fixed coefs don't transmit to coef
    ncoef <- max(c(cmap, mmap))
    new <- double(ncoef)
    i1 <- unique(cmap[cmap>0])
    new[i1] <- B[match(i1, cmap)]
    if (!is.null(mmap)) {
        i2 <- unique(mmap[mmap>0])
        ncoef[i2] <- B[match(i2, mmap)]
    }
    ncoef
}       
