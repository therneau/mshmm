# Checks of the survexpm routine
#
library(Matrix)

q1 <- matrix(0, 5, 5)  # the simple model of the NAFLD data
q1[1,2] <- q1[2,3] <- q1[3,4] <- 1
q1[1:4,5] <- 1

set.seed(1960)
rmat <- q1
rmat[rmat>0] <- exp(runif(7, -1, 1))
diag(rmat) <- -rowSums(rmat)

e1 <- survexpm(rmat, 2)  # use the decomposition method
e2 <- survexpm(rmat, 2, method="pade")      # use my Pade
e3 <- as.matrix(expm(2*rmat)) # use Matrix

all.equal(e1, e2)
all.equal(e1, e3)

#
# Compute derivatives
#
d1 <- survexpm(rmat, 2, deriv=TRUE)  # use the eigen decomp
d2 <- survexpm(rmat, 2, deriv=TRUE, method="pade") # use the Pade approach
all.equal(d1$P, e1)
all.equal(d1$deriv, d2$deriv)

# brute force derivatives
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
for (i in 1:100) {
    temp <- matrix(runif(49, .05, .15), 7, 7)
    diag(temp) <- diag(temp) - rowSums(temp)
    test <- survexpm(temp, deriv=TRUE)
    type[i] <- (test$method=="eigen")
}
table(type) 

# something is weird
e1 <- eigen(temp)
e2 <- eigen(t(temp))
e3 <- hmmeigen(temp)
e4 <- hmmeigen(t(temp))

scale3 <- colSums(e3$right*e3$left)
scale4 <- colSums(e4$right*e4$left)

rinv3 <- solve(e3$right)


ord3 <- order(Mod(e3$values), decreasing=TRUE)  # eigen sorts them post dgeev
ord4 <- order(Mod(e4$values), decreasing=TRUE)

all.equal(e1$value,  e2$values)
all.equal(e1$values, e3$values[ord3])
all.equal(e1$values, e4$values[ord4])
all.equal(e1$vectors, e3$right[,ord3)
all.equal(e2$vectors, e3$left[,ord3]

