#
# Try a test case a la Prentice
#
library(hmm)

qfun <- function(parm) {
    # simple 3 by 3
    p <- exp(parm)
    matrix(c(-(p[1] + p[2]), 0, 0, p[1], -p[3], 0, p[2], p[3], 0), 3,3)
    }

par <- c(1, 2.2, 1.8)
qq <- qfun(par)

temp <- eigen(qq)
a <- temp$vectors
ainv <- solve(a)
d <- temp$values

# exp works
all.equal(a %*% diag(exp(d * pi/10)) %*% ainv, as.matrix(expm(pi*qq/10)),
          check.attributes=FALSE)

# Now for derivative wrt parms
# p1
d1 <- matrix(0,3,3)
d1[1,1] <- -par[1]*exp(par[1])
d1[1,2] <- par[1] * exp(par[1])

gfun <- function(d, time, eps=1e-8) {
    outer(d, d, function(a, b)
        ifelse(abs(a-b)< eps, time*exp(time * (a+b)/2),
               (exp(a*time) - exp(b*time))/(a-b)))
}
V <- (ainv %*% d1 %*% a) * gfun(d, pi/10)
D1 <- a %*% V %*% ainv

eps <- 1e-6
D1x <- (expm((pi/10)*qfun(par + c(eps, 0, 0))) - expm((pi/10)*qq)) / eps
all.equal(D1, as.matrix(D1x), check.attributes=FALSE, tol=eps)

#
# use the test function to get the same thing
#
qmatrix <- matrix(0,3,3)
qmatrix[c(4, 7, 8)] <- 1:3

test2 <- derivtest(pi/10, 1, matrix(par,1), qmatrix, matrix(1:3, 1))
all.equal(D1, test2$dmat[,,1])
