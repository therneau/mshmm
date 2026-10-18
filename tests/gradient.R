# Test gradients
library(mshmm)
aeq <- function(x, y, ...) all.equal(as.vector(x), as.vector(y), ...)

sname <- levels(test1$state)[-1]
qmat <- matrix(0, 6, 6, dimnames=list(from= sname, to=sname))
qmat[1,2] <- qmat[1,3] <- .01
qmat[2,4] <- .03
qmat[3,4] <- .04
qmat[3,5] <- .03
qmat[4,5] <- .1
qmat[,6]  <- c(.02, .02, .02, .02, .06, 0)
icoef <- log(c(.01, .01, .03, .04, .03, .1, .02, .02, .02, .02, .06))

# First, a non-hmm
hfit1 <- cmsh(Surv(age,state) ~1, 
             data=test1, id=id, qmatrix=qmat, init=icoef, iter=0)
hfit1b <- cmsh(Surv(age,state) ~1, mc.cores=1,
             data=test1, id=id, qmatrix=qmat, init=icoef, iter=0, detail=TRUE)
   
# Do a numeric derivative for all parameters
eps <- 1e-8
pp <- which(qmat >0)
dvec <- double(length(pp))
for (j in 1:length(pp)) {
    itemp <- icoef
    itemp[j] <- itemp[j] + eps
    tfit <-  cmsh(Surv(age,state) ~1, 
             data=test1, id=id, qmatrix=qmat, init=itemp, iter=0)
    dvec[j] <- tfit$loglik
}
d1 <- (dvec- hfit1$loglik)/eps

# hfit1b$derv has the derivatives of alpha wrt parameters
# The loglik is log(sum(alpha)) for each subject, summing over states
t1 <- apply(hfit1b$deriv, c(1,3), sum)
alpha <- colSums(hfit1b$alpha)
t1 <- t1%*% diag(1/alpha)   # derivative of log, per subject
stopifnot(isTRUE(aeq(rowSums(t1), d1, tolerance=sqrt(eps))))

# Round 2: add a covariate
icoef2 <- rbind('(Intercept)'=icoef, male= 0:10/20)
hfit2 <-   cmsh(Surv(age, state) ~ male, data=test1, mc.cores=1,
               id =  id, qmatrix = qmat, init=icoef2, iter=0,
               control=cmsh.control(center=FALSE, scale=FALSE, detail=TRUE))

dvec <- double(length(icoef2))
for (j in 1:length(icoef2)) {
    temp <- icoef2
    temp[j] <- temp[j]+ eps
    tfit <-  cmsh(Surv(age, state) ~ male, data=test1, mc.cores=1,
               id =  id, qmatrix = qmat, init=temp, iter=0,
               control=cmsh.control(center=FALSE, scale=FALSE))
    dvec[j] <- tfit$loglik
}
d2 <- (dvec- hfit2$loglik)/eps

t2 <- apply(hfit2$deriv, c(1,3), sum)  # 22 parameters, 4 subjects
alpha <- colSums(hfit2$alpha)
t2 <- t2%*% diag(1/alpha)   # derivative of log, per subject
stopifnot(isTRUE(aeq(rowSums(t2), d2, tolerance=sqrt(eps))))


 
