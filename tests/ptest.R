#
# Test the pade function
#
library(hmm)

qmat <- matrix(0, 6, 6)
qmat[1,2] <- qmat[1,3] <- .01
qmat[2,4] <- .03
qmat[3,4] <- .04
qmat[3,5] <- .03
qmat[4,5] <- .1
qmat[,6]  <- c(.02, .02, .02, .02, .06, 0)

t1 <- expm(qmat)
t2 <- pade(qmat)
all.equal(t1, t2$S)

rmat <- qmat
diag(rmat) <- diag(rmat) - rowSums(rmat)  # a transformation matrix

all.equal(expm(rmat), pade(rmat)$S)

# The value of 50 forces the squaring portion of the routine
for (i in c(2, 5, 20, 50)) {
    test <- pade(i*rmat)
 #   print(test$nterm)
    print(all.equal(test$S, expm(i*rmat)))
    }

# Now test derivatives
#  Assume that elements of R are exp(eta)
#   for element 13 eta = -4 + 1.1 beta, for element 33 -4 + 2.1 beta
eps <- 1e-7
coef <- c(1.1, 2.1)
beta <- pi/3

rmat1 <- rmat2 <- rmat
rmat1[c(13,33)] <- exp(-4 + coef*beta)
diag(rmat1) <- diag(rmat1) - rowSums(rmat1)

rmat2[c(13,33)] <- exp(-4 + coef*(beta + eps))
diag(rmat2) <- diag(rmat2) - rowSums(rmat2)

d1 <- (expm(rmat2) - expm(rmat1))/ eps  # the result by brute force

# derivative of rmat wrt beta
drmat <- matrix(0., 6, 6)
eta <- (-4 + coef*beta)
drmat[c(13, 33)] <- coef * exp(eta)
diag(drmat) <- diag(drmat) - rowSums(drmat)

all.equal((rmat2-rmat1)/eps, drmat, tol=eps)  #dR/ dbeta


test <- pade(rmat1, array(drmat, c(6,6,1)))  #derivative from pade
all.equal(d1, test$dmat[,,1], tol=eps)

#Force it through the other loop of pade
#  the value of 10 gives no squaring step (test2$nterm =6)
d2 <- (expm(10*rmat2) - expm(10*rmat1))/eps
test2 <- pade(10*rmat1, 10*array(drmat, c(6,6,1)))
all.equal(d2, test2$dmat[,,1], tol=eps)

# Now force it to square a couple of times
d3 <- (expm(50*rmat2) - expm(50*rmat1))/eps
test3 <- pade(50*rmat1, 50*array(drmat, c(6,6,1)))
all.equal(d3, test3$dmat[,,1], tol=eps)


# Do the derivatives wrt eta and then combine them
rr <- (row(rmat))[c(13,33)]
temp1 <- temp2 <- matrix(0,6,6)
temp1[rr[1], rr[1]] <- -exp(eta[1])
temp1[13] <- exp(eta[1])
temp2[rr[2], rr[2]] <- -exp(eta[2])
temp2[33] <- exp(eta[2])

p2 <- pade(rmat1, array(c(temp1, temp2), c(6,6,2)))

all.equal(p2$dmat[,,1]*coef[1] + p2$dmat[,,2]*coef[2], test$dmat[,,1])

