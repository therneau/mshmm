#
# Fit the same model with msm and hmm, to verify that I have the
#  right likelihood.
library(msm)
library(hmm)
aeq <- function(x, y, ...) all.equal(as.vector(x), as.vector(y), ...)

# test1 has 4 subjects, 17 rows
# The age intervals are about a year, so make transitions around 5-15% per
# year.  This makes the loglik far from 1.
sname <- levels(test1$state)
qmat <- matrix(0, 6, 6, dimnames=list(from= sname, to=sname))
qmat[1,2:3] <- 1
qmat[2:3, 4] <- 1
qmat[3:4, 5] <- 1
qmat[-6,6] <- 1
# statefig(c(1,2,2,1), qmat)
icoef <- log(c(.05, .05, .06, .07, .05, .15, rep(c(.05,.07, .15), c(3,1,1))))

# make it a survival endpoint
test1$state <- factor(test1$istate, 0:6, c("censor", sname))

# The simplest model
hfit1 <- hmm(Surv(age,state) ~1, data=test1, id=id, qmatrix=qmat,
             init=icoef, iter=0)

# do the computation by hand
byhand <- function(data, eta, q=qmat, missmat) {
    idlist <- unique(data$id)
    nid <- length(idlist)
    phat <- matrix(0, nid, 6)  # n subjects, 4 states
    if (missing(missmat)) missmat <- diag(6)  # no errors
    for (i in 1:nid) {
        tdata <- subset(data, id== idlist[i])
        n <- nrow(tdata)
        e2 <- eta[data$id==idlist[i],]
        phat[i, tdata$istate[1]] <- 1  # initial state
        delta <- diff(tdata$age)
        for (j in 1:(n-1)) {
            # go forward
            rmat <- q
            rmat[rmat>0] <- exp(e2[j,])
            diag(rmat) <- diag(rmat) - rowSums(rmat)
            phat[i,] <- phat[i,] %*% as.matrix(expm(rmat*delta[j]))
            # multiply by D
            k <- tdata$istate[j+1]
            if (tdata$state[j+1]== "death") {
                temp <- rmat
                temp[,-6] <- 0
                temp[6,]  <- 0
                phat[i,] <-phat[i,] %*% temp
            } else  phat[i,] <- phat[i,]* missmat[k,]
        }
    }       
    phat
}
eta1 <- outer(rep(1, nrow(test1)), icoef)
true1 <- byhand(test1, eta1)
truelog <- sum(log(rowSums(true1)))
aeq(hfit1$loglik, truelog)

hfit1b <- hmm(Surv(age,state) ~1, data=test1, id=id, qmatrix=qmat,
             init=icoef, detail=TRUE, mc.cores=1)
aeq(hfit1b$loglik, truelog)
# derivatives
eps <- 1e-7
deriv <- double(11)
for (i in 1:11) {
    i2 <- icoef
    i2[i] <- i2[i]+ eps
    tfit <- hmm(Surv(age,state) ~1, data=test1, id=id, qmatrix=qmat,
             init=i2, iter=0)
    deriv[i] <- (tfit$loglik - hfit1$loglik)/eps
}
aeq(deriv, apply(hfit1b$deriv,1,sum))


# msm wants the intital rates in qmat, hmm only uses 0 vs >0
qmat[qmat>0] <- exp(icoef)
mfit1 <- msm(istate ~ age, data=test1, subject=id, death=6,
             qmatrix = qmat, fixedpar=TRUE)
# Below verifies that Chris Jackson and I agree wrt formulas
aeq(mfit1$minus2loglik, -2*truelog)  

# Add select covariates
i2 <- rbind(icoef,0,0)
dimnames(i2) <- list(c("(Intercept)", "educ", "male"), colnames(hfit1$cmap))
i2["educ", "1:3"] <- .1
i2["male", c("1:6", "2:6")] <- c(.2, .3)
# init arg expectes an object that looks look like coef(hfit2, matrix=TRUE)

hfit2 <- hmm(list(Surv(age,state) ~1, 
                  1:3 ~ educ, 1:6+ 2:6 ~ male),
             data=test1, id=id, qmatrix=qmat, init=i2, iter=0)
eta2 <- model.matrix(hfit2) %*% i2
true2 <- byhand(test1, eta2)
aeq(hfit2$log, sum(log(rowSums(true2))))


# hmm centered the covariates internally *and* also transformed the 
#  coefficients, i.e., the user never sees the change.  msm on the
#  other hand presents coefs wrt recentered data
mfit2 <- msm(istate ~ age, data=test1, subject=id, 
              qmatrix = qmat, fixedpar=TRUE, death=6, 
              covariates= list("1-3"= ~educ,  "1-6"= ~male, "2-6"= ~male),
              covinits= list(educ=.1, male=c(.2, .3)))

# To match msm we need to precenter our data to match it
test1b <- test1
#center <- attr(mfit1b$data$mm.cov, "means")
center <- c("educ"= 13.1538462,  "male"= 0.230769)
test1b$educ <- test1b$educ - center["educ"]
test1b$male <- test1b$male - center["male"]

hfit2b <- hmm(list(Surv(age,state) ~1, 
                  1:3 ~ educ, 1:6+ 2:6 ~ male),
             data=test1b, id=id, qmatrix=qmat, init=i2, iter=0)
aeq(-2*hfit2b$loglik[2], mfit2$minus2loglik)  


# Repeat with misclassification probabilities
# Here is a miss function for 6 states and fixed probs
e1 <- .12  # an A- as A+ or vice versa
e2 <- .2   # an N- as N+ or vice versa
temp <- outer(c("Acorrect"= (1-e1), "Afalse"= e1), 
              c("Ncorrect"=(1-e2), "Nfalse"= e2), '*')

missmat <- rbind(temp[c(1,2,3,4)], temp[c(2, 1, 4,3)],
                 temp[c(3,4,1,2)], temp[c(4, 3, 2, 1)])
missmat <- cbind(missmat, 0, 0)
missmat <- rbind(missmat, c(0,0,0,.1,.9,0), c(0,0,0,0,0,1))

# Treat the entry state as known (and death of course), others
#  hidde
temp <- with(test1, ifelse(duplicated(id)& state!='death', 0, istate))
test1$state3 <- factor(temp, 0:6, levels(test1$state))
hfit3 <- hmm(list(Surv(age,state3) ~1, 
                  1:3 ~ educ, 1:6+ 2:6 ~ male),
             marker = 1:5 ~1/ multinomial + fixed=missmat),
             data=test1, id=id, qmatrix=qmat, init=i2, iter=0)

mfit3 <- msm(istate ~ age, data=test1, subject= id, 
             qmatrix = qmat, fixedpar=TRUE, death=6,
             ematrix=missmat, initprob=c(1,1,1,1,0,0)/4)

init6 <- function(nstate, ...) {
    c(1,1,1,1,0,0)/4
}
hfit3 <-  hmm(hbind(age, state) ~ 1, data=test1, mc.cores=3,
               id = id, qmatrix = qmat, rfun=hmiss,
               pfun=init6, mfun=hmmtest, mpar=list(fn="hmmloglik"),
               otype= otype, death=6)
aeq(-2*hfit3$loglik[2], mfit3$minus2loglik)
 

# These models use the cav data set from the msm package
Qm <- rbind(c(0, .148, 0, .0171),
            c(0,  0,  .202, .081),
            c(0,  0,   0,  .126),
            c(0,  0,   0,   0))  #page 38, msm manual

ematrix <- rbind(c(.8, .2, 0, 0),
                 c(.1, .8, .1, 0),
                 c(0, 0.3, .7, 0),
                 c(0, 0,    0, 1))


mfit4a <- msm(state ~ years, subject=PTNUM, data=cav,
            qmatrix=Qm, ematrix=ematrix, death=4,
            obstrue= firstobs, fixedpar=TRUE)

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


hfit4a <- hmm(cbind(years, state) ~ 1, data=cav, mc.cores=1,
            id = PTNUM, qmatrix=Qm, rfun=efun, pfun=cinit,
            death=4, otype=otype,  rcoef=rcoef, mfun= hmmtest)
aeq(-2*hfit4a$loglik[2], mfit4a$minus2loglik)


# Now with covariates
mfit4b  <- msm(state ~ years, subject=PTNUM, data=cav,
               qmatrix=Qm, ematrix=ematrix, death=4, 
               obstrue= firstobs, covariates=list("1-2"= ~sex, "1-4"= ~sex),
               fixedpar=TRUE, covinits=list(sex=c(1.1, 2.1)))

qcoef <- data.frame(state1=c(1,1), state2=c(2,4), 
                    term =1, coef=1:2, init=c(1.1, 2.1))
mcenter <- attr(mfit4b$data$mm.cov, "means")

hfit4b <- hmm(cbind(years, state) ~ I(sex-mcenter), data=cav, mc.cores=1,
              id = PTNUM, qmatrix=Qm, rfun=efun, pfun=cinit,
              death=4, otype=otype,  rcoef=rcoef, mfun= hmmtest,
              qcoef=qcoef, scale=FALSE)
aeq(-2*hfit4b$loglik[2], mfit4b$minus2loglik)

