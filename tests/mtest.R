library(hmm)
#
# Test the mlogit function
#
# First a simple vector
eta <- matrix(c(-1, -2, .1), nrow=1)

eps <- 1e-7
test <- mlogit(eta, gradient=TRUE)

temp <- c(1, exp(eta))
p <- temp/sum(temp)
all.equal(as.vector(test), p)

deriv <- matrix(0,4,3)
for (i in 1:3) {
    eta2 <- eta
    eta2[1,i] <- eta2[1,i] + eps
    temp2 <- c(1, exp(eta2))/(1 + sum(exp(eta2)))
    deriv[,i] <- (temp2- p)/eps
}
all.equal(deriv, attr(test, "gradient")[1,,], tol=sqrt(eps))

# Now with multiple rows
eta <- matrix(-(1:15)/15, nrow=5)
test <- mlogit(eta, gradient=TRUE)

phat <-  cbind(1,exp(eta))/ (1+ rowSums(exp(eta)))
all.equal(test, phat, check.attributes=FALSE)

dmat <- attr(test, "gradient")
for (i in 1:nrow(phat)) {
    temp <- diag(phat[i,]) - outer(phat[i,], phat[i,])
    print(c(all.equal(dmat[i,,], temp[,-1]), 
          all.equal(dmat[i,,],  
                    attr(mlogit(eta[i,,drop=FALSE], TRUE), "gradient")[1,,])))
}
