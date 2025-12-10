# Checks of the survexpm routine
#
library(survival)
library(Matrix)

q1 <- matrix(0, 5, 5)  # the simple model of the NAFLD data
q1[1,2] <- q1[2,3] <- q1[3,4] <- 1
q1[1:4,5] <- 1

set.seed(1960)
rmat <- q1
rmat[rmat>0] <- exp(runif(7, -1, 1))
diag(rmat) <- -rowSums(rmat)

s1 <- survexpmsetup(rmat)
e1 <- survexpm(rmat, 2, s1)  # use the decomposition method
e2 <- survexpm(rmat, 2)      # use my Pade
e3 <- as.matrix(expm(2*rmat)) # use Matrix

all.equal(e1, e2)
all.equal(e1, e3)

#
# Compute derivatives
#
d1 <- survexpm(rmat, 2, s1, deriv=2)  # use the eigen decomp
d2 <- survexpm(rmat, 2, deriv=2)      # use the Pade approach
all.equal(d1$P, e1)
all.equal(d1$deriv, d2$deriv)


eps <- 1e-6
indx <- cbind(row(rmat)[rmat>0], col(rmat)[rmat>0])
for (i in 1:7) {
    rtemp <- rmat
    # the rows of rmat have to sum to zero, which is what defines the
    #  diagonal. So we must perturb both
    rtemp[indx[i,1], indx[i,]] <- rtemp[indx[i,1], indx[i,]] + c(-eps, eps)
    ptemp <- as.matrix(expm(2*rtemp))
    delta <- (ptemp- e1)/eps
    print(all.equal(d1$deriv[,,i], delta, tol=eps))
}

if (FALSE) {
    # is one method faster than another?  The value of 2 pushes the pade into
    #  a longer calculation (more terms), so use something smaller. 
    n <- 4e5
    tt <- .5
    t1a <- system.time({for (i in 1:n) survexpm(rmat, tt, s1, deriv=FALSE)})    
    t1b <- system.time({for (i in 1:n) survexpm(rmat, tt, deriv=FALSE)})
    t1c <- system.time({for (i in 1:n) expm(tt*rmat)})
    t2a <- system.time({for (i in 1:n) survexpm(rmat, tt, s1, deriv=TRUE)})    
    t2b <- system.time({for (i in 1:n) survexpm(rmat, tt, deriv=TRUE)})
    temp <- rbind(eigen= t1a, pade=t1b, expm=t1c, 
                  "eigen/deriv"= t2a, "pade/deriv"=t2b)
}
