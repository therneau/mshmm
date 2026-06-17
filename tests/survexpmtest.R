# Checks of the survexpm routine
#
library(mshmm)
aeq <- function(x,y) all.equal(as.vector(x), as.vector(y))

q1 <- matrix(0, 5, 5)  # the simple model of the NAFLD data
q1[1,2] <- q1[2,3] <- q1[3,4] <- 1
q1[1:4,5] <- 1

set.seed(1960)
rmat <- q1
rmat[rmat>0] <- exp(runif(7, -1, 1))
diag(rmat) <- -rowSums(rmat)

e1 <- survexpm(rmat, 2.1)  # use the decomposition method
e2 <- survexpm(rmat, 2.1, method="pade")      # use my Pade
e3 <- as.matrix(expm(2.1*rmat)) # use Matrix

all.equal(e1, e2)
all.equal(e1, e3)

#
# Compute derivatives
#
d1 <- survexpm(rmat, 2.1, deriv=1)  # use the eigen decomp
d2 <- survexpm(rmat, 2.1, deriv=1, method="pade") # use the Pade approach
all.equal(d1$P, e1)
all.equal(d1$deriv, d2$deriv)

d3 <- survexpm(rmat, 2.1, deriv=2)
d.eta <- rmat[rmat>0]  # deriv of each element wrt eta
all.equal(d3$deriv, d1$deriv * rep(d.eta, each=25))

# brute force derivatives
eps <- 1e-6
indx <- cbind(row(rmat)[rmat>0], col(rmat)[rmat>0])
for (i in 1:7) {
    rtemp <- rmat
    # the rows of rmat have to sum to zero, which is what defines the
    #  diagonal. So we must perturb both
    rtemp[indx[i,1], indx[i,]] <- rtemp[indx[i,1], indx[i,]] + c(-eps, eps)
    ptemp <- as.matrix(expm(2.1*rtemp))
    delta <- (ptemp- e1)/eps
    print(all.equal(d1$deriv[,,i], delta, tol=eps))
}

if (FALSE) {
    # is one method faster than another?  The value of 2 pushes the pade into
    #  a longer calculation (more terms), so use something smaller. 
    n <- 1e5
    tt <- .5
    t1a <- system.time({for (i in 1:n) survexpm(rmat, tt, deriv=FALSE)})    
    t1b <- system.time({for (i in 1:n) survexpm(rmat, tt, deriv=FALSE, 
                                                method="pade")})
    t1c <- system.time({for (i in 1:n) expm(tt*rmat)})
    t2a <- system.time({for (i in 1:n) survexpm(rmat, tt, deriv=TRUE)})    
    t2b <- system.time({for (i in 1:n) survexpm(rmat, tt, deriv=TRUE, 
                                                method="pade")})
    temp1 <- rbind(eigen= t1a, pade=t1b, expm=t1c, 
                  "eigen/deriv"= t2a, "pade/deriv"=t2b)

    # try a bigger matrix, corresponds to A1-A4 x C1-C4, + death
    sdata <- data.frame(state= c(paste0("S", 1:16), 'death'),
                        A= c(rep(1:4, 4), 5),
                        C= c(rep(1:4, each=4), 5),
                        D= rep(0:1, c(16,1)))
    r2 <- qplus(sdata, death='death')
    r2[r2>0] <- exp(runif(sum(r2>0), -1,0))
    diag(r2) <- -rowSums(r2)

    tt <- .5
    x1a <- system.time({for (i in 1:n) survexpm(r2, tt, deriv=FALSE)})    
    x1b <- system.time({for (i in 1:n) survexpm(r2, tt, deriv=FALSE, 
                                                method="pade")})
    x1c <- system.time({for (i in 1:n) expm(tt*r2)})
    x2a <- system.time({for (i in 1:n) survexpm(r2, tt, deriv=TRUE)})    
    x2b <- system.time({for (i in 1:n) survexpm(r2, tt, deriv=TRUE, 
                                                method="pade")})
    temp2 <- rbind(eigen= x1a, pade=x1b, expm=x1c, 
                  "eigen/deriv"= x2a, "pade/deriv"=x2b)
}

# Non-triangular
set.seed(1953)
type <- logical(100)
good <- logical(100)
for (i in 1:100) {
    temp <- matrix(runif(49, .05, .15), 7, 7)
    diag(temp) <- diag(temp) - rowSums(temp)
    test <- survexpm(temp, deriv=TRUE)
    type[i] <- (test$method=="eigen")
    good[i] <- all.equal(test$P, expm(temp))
}
table(type,good)  

#
# A case where Pade occurs, tied eigenvalues, upper triangular, 
#  eigenvectors and values are real
#
tied <- rbind(c(-.1, .02, .03,   0, .05),
              c( 0 , -.2,   0, .04, .16),
              c( 0,   0,  -.3,  .2,  .1),
              c( 0,   0,    0, -.2,  .2), 0)

test <- survexpm(tied, deriv=T)
test$method
all.equal(test$P, as.matrix(expm(tied)))

etest <- mshmm:::hmmeigen(tied)  # this fcn is not exported
rinv <- solve(etest$right)
scale <- colSums(etest$right * etest$left)
norm(etest$right)* norm(rinv)  # huge condition number
max(1/scale)                   # and approx condition number as well

# Both of these fail, showing that survexp was correct to avoid the eigen
test2 <- etest$right %*% diag(etest$values) %*% rinv
all.equal(test2, tied)
p1 <- etest$right %*% diag(exp(etest$values)) %*% rinv
all.equal(p1, test$P)

# The expm eigen has a more forgiving cutoff before switching to Pade
# than my survexpm function
all.equal(expm(tied), expm(tied, method="R_Eigen"))

#
# Another check of Pade, with another matrix that has tied eigenvalues
#
R <- matrix(0, 6, 6)
R[1,2:3] <- 1
R[2:3, 4] <- 1
R[3:4, 5] <- 1
R[-6,6] <- 1
R[R>0] <- c(.05, .05, .06, .07, .05, .15, rep(c(.05,.07, .15), c(3,1,1)))
diag(R) <- -rowSums(R)

q1 <- expm(R* 1.3)
q2 <- pade(R* 1.3)
q3 <- survexpm(R, 1.3, deriv=FALSE)
aeq(q1, q2$P)
aeq(q1, q3)
