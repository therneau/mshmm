#
# Test the Ptrans routine
#
library(mshmm)
alpha <- 1:4/10   # 4 states
dP <- array(runif(48, -.2, 2), dim=c(4,4,3))  # 3 linear predictors
x <- c(1, 2.2, 0, -1.1, 3.3) # 5 covariates for this subject
cmap <- cbind(c(1,2,0,0,3), c(4,0, 5,6,0), c(7,8,0,0, 9))

test1 <- mshmm:::Ptrans(alpha, dP, cmap, x)  # not exported

ptest <- function(alpha, dP, cmap, x) {
    dd <- dim(dP)  
    # dP will have dim(nstate, nstate, number of etas)
    #  it contains the deriviative wrt each eta, for each element of P
    # The line defining Z is arcane, see "derivatives, eta to beta" in
    #  the code vignette, Z transforms from d/deta to d/dbeta (derivatives)
    nbeta <- max(cmap)
    Z <- matrix(0,dd[3], nbeta)
    cz <- which(cmap>0)
    Z[cbind(col(cmap)[cz], cmap[cz])] <- x[row(cmap)[cz]]

    # transform each derivative of P
    nstate <- dd[1]
    newd <- array(0, dim=c(nstate, nstate, nbeta))
    for (i in 1:nstate) {
        for (j in 1:nstate) newd[i,j,] <- dP[i,j,] %*% Z
    }
   
    # create the product alpha %*% deriv, one row for each coefficient
    aderiv <- matrix(0, nbeta, nstate)
    for (i in 1:nbeta) aderiv[i,] <- alpha %*% newd[,,i]
    aderiv
}
 
test2 <- ptest(alpha, dP, cmap, x)
all.equal(test1, test2)

if (FALSE) {
    # did the cleverness of Ptrans pay off?
    # do a bigger problem as a test
    cmap <- matrix(1:45, 5, 9)  # 5 predictors, 9 linear predictors
    x <- 1:5
    dP <- array(runif(6*6*9), dim=c(6,6,9)) # 6 states, 9 transtions
    alpha <- rep(1/6,6)

    nrep <- 1e5 # 5000 obs x 20 iterations
    t1 <- system.time({
        for (i in 1:nrep) temp <- Ptrans(alpha, dP, cmap, x)})
    t2 <- system.time({
        for (i in 1:nrep) temp <- ptest(alpha, dP, cmap, x)})
    rbind(t1, t2)
    # For the example above the time was 4.5 vs 25.5 seconds, a factor of
    #  of 5.7.  On the other hand 25 seconds is a small part of overall
    #  compute time for many of our fits.  So maybe worth the complexity
}
