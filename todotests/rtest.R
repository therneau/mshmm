library(hmm)
#
# Test the continuous response functions
#
f1 <- function(y, eta, gradient=FALSE) {
    cuts <- matrix(c(-Inf, 1.3, 1.5, 1.3, 1.5, Inf), 3, 2)
    hmmncut(y, nstate=4, eta=eta, gradient=gradient, cuts, c(1,2,3,0))
}

# Check out the function
yy <- seq(1, 2, by=.25)
test1 <- f1(yy, -2)
ny <- length(yy)
tmat <- matrix(0, 4, ny)
estd <- exp(-2)

tmat[1,] <- pnorm(1.3-yy, 0, estd)
tmat[2,] <- pnorm(1.5-yy, 0, estd) - pnorm(1.3-yy, 0, estd)
tmat[3,] <- 1 - pnorm(1.5-yy, 0, estd)
all.equal(tmat, test1)

# Now the derivatives wit respect to eta
tfun <- function(y, eta, c1, c2) {
    s <- exp(eta)
    pnorm(c2-y, 0, s) - pnorm(c1-y, 0, s)
}
eps= 1e-8
test2 <- attr(f1(yy, -2, gradient=TRUE), "gradient")

tmat[1,] <- (tfun(yy, eps-2, -Inf, 1.3) - tfun(yy, -2, -Inf, 1.3))/eps
tmat[2,] <- (tfun(yy, eps-2,  1.3, 1.5) - tfun(yy, -2,  1.3, 1.5))/eps
tmat[3,] <- (tfun(yy, eps-2,  1.5, Inf) - tfun(yy, -2,  1.5, Inf))/eps

all.equal(tmat, test2[,,1], tolerance=sqrt(eps))

# Repeat with the logistic
f2 <- function(y, eta, gradient=FALSE) {
    cuts <- matrix(c(-Inf, 1.3, 1.5, 1.3, 1.5, Inf), 3, 2)
    hmmlcut(y, nstate=4, eta=eta, gradient=gradient, cuts, c(1,2,3,0))
}

test1 <- f2(yy, -2)
estd <- exp(-2)*sqrt(3)/pi

tmat[1,] <- plogis(1.3-yy, 0, estd)
tmat[2,] <- plogis(1.5-yy, 0, estd) - plogis(1.3-yy, 0, estd)
tmat[3,] <- 1 - plogis(1.5-yy, 0, estd)
all.equal(tmat, test1)

# Now the derivatives wit respect to eta
tfun <- function(y, eta, c1, c2) {
    s <- exp(eta) * sqrt(3)/pi
    plogis(c2-y, 0, s) - plogis(c1-y, 0, s)
}
eps= 1e-8
test2 <- attr(f2(yy, -2, gradient=TRUE), "gradient")

tmat[1,] <- (tfun(yy, eps-2, -Inf, 1.3) - tfun(yy, -2, -Inf, 1.3))/eps
tmat[2,] <- (tfun(yy, eps-2,  1.3, 1.5) - tfun(yy, -2,  1.3, 1.5))/eps
tmat[3,] <- (tfun(yy, eps-2,  1.5, Inf) - tfun(yy, -2,  1.5, Inf))/eps

all.equal(tmat, test2[,,1], tolerance=sqrt(eps))


#
#  Gaussian response
#
f3 <- function(y, eta, gradient=FALSE) {
    hmmgauss(y, nstate=5, eta=eta, gradient=gradient, c(1,3,2,2, 0))
}

# the most common use is with the same eta for everyone
# this is for a 3 class problem, common std of exp(-1),
#  means of 0.1, 1.1, 2.2
eta3 <- outer(c(1,1,1,1,1), c(0.1, 1.1, 2.2, -1), '*')
test3 <- f3(yy, eta3, gradient=TRUE)
temp <- dnorm(outer(yy, c(0.1, 2.2, 1.1, 1.1), '-'), sd=exp(-1)) 
all.equal(t(temp), test3[1:4,], check.attributes=FALSE)
all.equal(test3[5,], rep(0,5))  # expect to fail

#derivatives wrt eta
zz <- test3
attr(zz, 'gradient') <- NULL
for (i in 1:4) {
    etax <- eta3
    etax[,i] <- etax[,i] + eps
    temp <- f3(yy, etax)
    dd <- (temp - zz)/eps
    print(all.equal(dd, attr(test3, "gradient")[,,i], tol=eps*10))
}   

# Use 3 means and 3 std values
eta4 <- cbind(eta3, -1.3, -2)
test4 <- f3(yy, eta4, gradient=TRUE)
temp <- dnorm(outer(yy, c(0.1, 2.2, 1.1, 1.1), '-'),
              sd= exp(outer(rep(1,length(yy)), c(-1, -2, -1.3, -1.3))))
all.equal(t(temp), test4[1:4,], check.attributes=FALSE)

#derivatives
zz <- test4
attr(zz, 'gradient') <- NULL
for (i in 1:6) {
    etax <- eta4
    etax[,i] <- etax[,i] + eps
    temp <- f3(yy, etax)
    dd <- (temp - zz)/eps
    print(all.equal(dd, attr(test4, "gradient")[,,i], tol=eps*10))
}   


# beta distribution
yy <- seq(0, 1, length=7)[2:6]  
f4 <- function(y, eta, gradient=FALSE) {
    hmmbeta(y, nstate=4, eta=eta, gradient=gradient, c(1,3,2,2))
}
test4 <- f4(yy, eta4, gradient=TRUE)
temp <- rbind(dbeta(yy, exp(eta4[,1]), exp(eta4[,2])),
              dbeta(yy, exp(eta4[,3]), exp(eta4[,4])),
              dbeta(yy, exp(eta4[,5]), exp(eta4[,6])))
                    
all.equal(test4, temp[c(1,3,2,2),], check.attributes=FALSE)

#derivatives
zz <- test4
attr(zz, 'gradient') <- NULL
for (i in 1:6) {
    etax <- eta4
    etax[,i] <- etax[,i] + eps
    temp <- f4(yy, etax)
    dd <- (temp - zz)/eps
    print(all.equal(dd, attr(test4, "gradient")[,,i], tol=sqrt(eps)))
}   

