#
# Try using Fisher scoring on the cav data
#
library(hmm)
library(msm)


# Per the msm manual, it is important to start iteration with a 
#  sensible intercepts.  For this model the valid transitions are
#  1-2, 2-3, and (1:3) to 4
valid <- rbind(c(0, 1, 0, 1),
               c(0, 0, 1, 1),
               c(0, 0, 0, 1), 0)

# Here is a so-so estimate
qmat1 <- crudeinits.msm(state ~ years, PTNUM, data=cav, qmatrix=valid)

# msm fit
emat <- matrix(0.0, 4,4)
emat[1,2] <- emat[2,1] <- emat[2,3] <- emat[3,2] <- .1  # 10% error rate
mfit1 <- msm(state ~ years, data=cav, subject=PTNUM, qmatrix=qmat1,
             ematrix=emat, death=4, obstrue=firstobs, method="BFGS")

# hmm fits
cinit <- function(nstate, ...) {
    init <- rep(0.0, nstate)
    init[1] <- 1   # everyone starts in state 1
    init
}
 
# the first state is exact, so no response function.
otype <- ifelse(cav$firstobs, 0, ifelse(cav$state==4, 2,1))

rcoef <- data.frame(lp=1:4, term=0, coef=1:4, 
                    init= log(c(2/8, 1/8, 1/8, 3/7)))
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

hfit1 <- hmm(cbind(years, state) ~ sex, data=cav, mc.cores=6,
            id = PTNUM, qmatrix=qmat1, rfun=errfun, pfun=cinit,
            death=4, otype=otype,  rcoef=rcoef, 
            mfun= optim,  
            mpar=list(control=list(fnscale= -1), method='BFGS', gr= "hmmgrad",
                      hessian=TRUE))

# Now try out Fisher scoring.  The optim run took many steps
hfit2 <- hmm(cbind(years, state) ~ sex, data=cav, mc.cores=6,
            id = PTNUM, qmatrix=qmat1, rfun=errfun, pfun=cinit,
            death=4, otype=otype,  rcoef=rcoef, 
            mfun=hmmscore, mpar=list(gr="hmmboth", eps=1e-7))

hfit2$fit$trace$loglik
all.equal(hfit1$loglik, hfit2$loglik, tolerance=1e-7)
all.equal(unname(hfit1$loglik[2]), mfit1$minus2loglik/-2, tolerance=1e-7)


#plot(-hfit$fit$hessian, hfit2$fit$S)
#abline(0,1)   #not too bad

