#
# A test case for the Prentice algorithm
#
library(mshmm)
library(Matrix)

qfun <- function(parm) {
    # simple 3 by 3
    p <- exp(parm)
    matrix(c(-(p[1] + p[2]), 0, 0, p[1], -p[3], 0, p[2], p[3], 0), 3,3)
    }

par <- c(1, 2.2, 1.8)  # linear predictors for the 3 off diagonal values
qq <- qfun(par)  # a transition rates matrix

# Obtain the eigenvectors and eigenvalues using the eigen function
temp <- eigen(qq)
a <- temp$vectors
ainv <- solve(a)
d <- temp$values

P1 <- as.matrix(expm(pi*qq/10)) # directly use expm
P2 <- a %*% diag(exp(d * pi/10)) %*% ainv # Prentice formula
all.equal(P1, P2, check.attributes=FALSE)

# Now for derivative wrt parms
# d1 = deriv of R wrt par[1]
d1 <- matrix(0,3,3)
d1[1,1] <- -par[1]*exp(par[1])
d1[1,2] <- par[1] * exp(par[1])

gfun <- function(d, time, eps=1e-8) {
    outer(d, d, function(a, b)
        ifelse(abs(a-b)< eps, time*exp(time * (a+b)/2),
               (exp(a*time) - exp(b*time))/(a-b)))
}
V <- (ainv %*% d1 %*% a) * gfun(d, pi/10)
D1 <- a %*% V %*% ainv  # deriv of P wrt p[1], Prentice formula

# brute force derivative
eps <- 1e-6
D1x <- (expm((pi/10)*qfun(par + c(eps, 0, 0))) - expm((pi/10)*qq)) / eps
all.equal(D1, as.matrix(D1x), check.attributes=FALSE, tol=eps)

#
# use survexpm to access internal version
#
test2 <- survexpm(qq, pi/10, deriv=1)
all.equal(D1, test2$deriv[,,1])
