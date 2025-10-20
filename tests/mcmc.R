#
# a short test of the mcmc routine to check the calling sequence
#   use the msm data since it is fast
#
library(msm)  # for the cav data
library(hmm)

Qm <- rbind(c(0, .148, 0, .0171),
            c(0,  0,  .202, .081),
            c(0,  0,   0,  .126),
            c(0,  0,   0,   0))  #page 38, msm manual

cinit <- function(nstate, ...) {
    init <- rep(0.0, nstate)
    init[1] <- 1   # everyone starts in state 1
    init
}

# the first state is exact, so no response function.
otype <- ifelse(cav$firstobs, 0, ifelse(cav$state==4, 2,1))

errfun <- function(y, nstate, eta, gradient) {
    # true states 1, 2, and 3 have separate linear predictors
    temp1 <- hmulti(y, nstate, eta[,1], gradient,
                    statemap= rbind(1:2, 0,0,0))
    temp2 <- hmulti(y, nstate, eta[,2:3], gradient,
                    statemap= rbind(0, c(2,1,3), 0, 0))
    temp3 <- hmulti(y, nstate, eta[,4], gradient,
                    statemap= rbind(0, 0, c(3,2), 0))
    pmat <- rbind(temp1[1,], temp2[2,], temp3[3,], ifelse(y==4,1,0))
    if (gradient) {
        gmat <- array(c(attr(temp1, 'gradient'), attr(temp2, 'gradient'),
                        attr(temp3, 'gradient')), dim=c(nstate, length(y), 4))
        attr(pmat, "gradient") <- gmat
    }
    pmat
}

rcoef <- data.frame(lp=1:4, term=0, coef=1:4, 
                    init= log(c(2/8, 1/8, 1/8, 3/7)))
qcoef <- data.frame(state1=c(1,2,1), state2=c(2,3,4), term=1,
                    coef=1:3)

hfit1 <- hmm(cbind(years, state) ~ sex, data=cav, mc.cores=6,
             id = PTNUM, qmatrix=Qm, rfun=errfun, pfun=cinit,
             death=4, otype=otype,  rcoef=rcoef, qcoef=qcoef,
             scale=FALSE,
             mfun=hmmscore, mpar=list(gr="hmmboth", hessian=TRUE))

vmat <- solve(hfit1$fit$hessian)
groups <- rep(c(1,2), c(8,4))  # error parms in group 2
mpar <- list(iter=100, burnin=100, 
             prior.means = hfit1$coef,  # centers of distribution
             prior.std = rep(10, length(hfit1$coef)), # fairly vague prior
             pvar= vmat,   # starting variance for trial values
             group=groups, # bundle the parameters
             reset = 50,   # how often to reset during burnin
             trace = 0,    # level of chattiness (every trace say something)
             bfit  = TRUE) # return the burnin iterations

set.seed(1953)  # force common results
hfit2 <- hmm(cbind(years, state) ~ sex, data=cav, mc.cores=6,
              id = PTNUM, qmatrix=Qm, rfun=errfun, pfun=cinit,
              death=4, otype=otype,  rcoef=rcoef, qcoef=qcoef,
              scale=FALSE, icoef=hfit1$coef,
              mfun=mcmc1, mpar=mpar)

v2 <- var(hfit2$fit$par)  # empirical variance of the mcmc process
signif(rbind(hessian = sqrt(diag(vmat)), mcmc=sqrt(diag(v2))), 2)

# with only 100 iterations the mcmc will have large simulation error
# A better check on both it and the hessian would be 2-3 thousand iters,
#  but we don't want such a slow test in the default suite
