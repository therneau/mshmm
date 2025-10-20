# This is an MCMC version of mfit4 in ../tests/jackson.R

library(msm)
#library(hmm)
source("../trial/loadall.R") # for testing without loading
Qm <- rbind(c(0, .148, 0, .0171),
            c(0,  0,  .202, .081),
            c(0,  0,   0,  .126),
            c(0,  0,   0,   0))  #page 38, msm manual

ematrix <- rbind(c(.8, .2, 0, 0),
                 c(.1, .8, .1, 0),
                 c(0, 0.3, .7, 0),
                 c(0, 0,    0, 1))

cinit <- function(nstate, ...) {
    init <- rep(0.0, nstate)
    init[1] <- 1   # everyone starts in state 1
    init
}

otype <- with(cav, ifelse(state==4, 2, 1))
first <- !duplicated(cav$PTNUM)
otype[first] <- 0   # the first state is exact, so no response function.

efun <- function(y, nstate, eta, ...) {
    ptemp <- exp(eta[1,])
    emat <- diag(4)
    emat[1,2] <- ptemp[2]    # fill in by column
    emat[2,1] <- ptemp[1]
    emat[2,3] <- ptemp[4]
    emat[3,2] <- ptemp[3]
    emat <- emat/rowSums(emat)
    emat[,y, drop=FALSE]
}

rcoef <- data.frame(lp=1:4, term=0, coef=1:4, 
                    init= log(c(1/8, 2/8, 3/7, 1/8)))

# Verify that I have it set up correctly
# On the linear predictor scale each element is log(e[i,j]/e[i,i])
all.equal(ematrix, efun(1:4, 4, matrix(rcoef$init, 1)))

# With iteration and a covariate
mfit4  <- msm(state ~ years, subject=PTNUM, data=cav,
            qmatrix=Qm, ematrix=ematrix, death=4,
            obstrue= firstobs, covariates=list("1-2"= ~sex, "1-4"= ~sex))

q4 <- data.frame(state1= c(1, 1),
                 state2= c(2, 4),
                 term  = c(1, 1),
                 coef  = c(1, 2))

# Use Nelder-Mead to give it a head start
hfit4a <- hmm(cbind(years, state) ~ sex, data=cav, mc.cores=detectCores()-1,
            id = PTNUM, qmatrix=Qm, rfun=efun, pfun=cinit,
            death=4, otype=otype,  rcoef=rcoef, qcoef= q4, 
            mfun= optim, mpar=list(control=list(fnscale=-1, maxit=1000)))

temp <- hfit4a$beta[,1:5]
q4a <- data.frame(state1 = c(1,1,2,1,1,2,3),
                  state2 = c(2,2,3,4,4,4,4),
                  term   = c(0,1,0,0,1,0,0),
                  coef   = 1:7,
                  init = temp[temp!=0])
r4a <- rcoef
r4a$init <- hfit4a$beta[1, 6:9]

hfit4 <- hmm(cbind(years, state) ~ sex, data=cav, mc.cores=detectCores()-1,
            id = PTNUM, qmatrix=Qm, rfun=efun, pfun=cinit,
            death=4, otype=otype,  rcoef=r4a, qcoef= q4a, 
            mfun= mcmc0, mpar=list(burnin= 1000, iter=3000, print=250, 
                                   reset=250, prior.std=5)
save(hfit4, hfit4a, mfit4, file='jackson2.rda')
