# Test the upper call
library(hmm)

beta <- log(c(.01, .01, .03, .04, .03, .1, .02, .02, .02, .02, .06))
index <- as.integer(c(7, 13, 20, 21, 27, 28, 31:35))

eta <- outer(1:5, beta)
times <-  c(1.1, 1.2, .8, .3, .2)

# compute the true values
true <- array(0, c(6,6,5))
for (i in 1:5) {
    rmat <- matrix(0, 6,6)
    rmat[index] <- exp(eta[i,])
    diag(rmat) <- diag(rmat) - rowSums(rmat)
    true[,,i] <- expm(rmat * times[i])
}

# A call with "0" terms for the H09 algorithm uses the Ward approach instead
test1 <- .Call("upper", 6L, eta, times, index, 1e-6, 0)
temp <- array(test1$P, c(6,6,5))
all.equal(temp, true)

# A matrix which caused issues
index2 <- as.integer(c(5, 10, 13:15))
eta <- c(-2.769494, -1.181372, -3.433431, -2.839075, -1.287341)
time2 <- seq(.003, 16, length=50)

test2 <- .Call("upper", 4L, outer(rep(1,50), eta), time2, index2, 1e-7, 0)
true2 <- array(0, c(4,4,50))
rmat <- matrix(0, 4,4)
rmat[index2] <- exp(eta)
diag(rmat) <- diag(rmat) - rowSums(rmat)
for (i in 1:50) 
    true2[,,i] <- expm(rmat * time2[i])
all.equal(test2$P, c(true2))
    
