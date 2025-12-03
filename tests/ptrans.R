#
# Test the Ptrans routine
#
alpha <- 1:4/10   # 4 states
dP <- array(runif(48, -.2, 2), dim=c(4,4,3))  # 3 linear predictors
x <- c((-4):4 # 5 covariates for this subject
cmap <- cbind(c(1,2,0,0,3), c(4,0, 5,6,0), c(7,8,0,0, 9))
#test1 <- mshmm:::Ptrans(alpha, dP, cmap, x)
test1 <- Ptrans(alpha, dP, cmap, x)

ptest <- function(alpha, dP, x) {
    # dP will have dim(nstate, nstate, number of etas)
    # A simple transform is
    dd <- c(dim(dP), length(x))
    nstate <- dd[1]
    newd <- array(0, dim=c(nstate, nstate, dd[4]))
    for (i in 1:nstate) {
        for (j 1:nstate) newd[i,j,] <- dP[i,j,] * ifelse(x==0, 0, 1/x)
    }
    dim(newd) <- c(nstate, nstate, dd[4])
    
    temp <- matrix(0, dd[4], nstate)
    for (i in 1:dd[4]) temp[i,] <- alpha %*% newd[,,i]
    temp
}
