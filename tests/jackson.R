#
# Fit the same model with msm and hmm, to verify that I have the
#  right likelihood.
# Getting the initial values in the correct order is the bugbear here, as msm
#  and icmsh use very different styles, plus msm thinks of the state space 
#  matrix in row major order, and icmsh in column major, ie. standard R.
#
library(msm)
library(mshmm)
aeq <- function(x, y, ...) all.equal(as.vector(x), as.vector(y), ...)

# test1 has 4 subjects, 17 rows
# The age intervals are about a year, so make transitions around 5-15% per
# year.  This makes the loglik far from 1.
sname <- levels(test1$state)[-1]   #censor is not a state
qmat <- matrix(0, 6, 6, dimnames=list(from= sname, to=sname))
qmat[1,2:3] <- 1
qmat[2:3, 4] <- 1
qmat[3:4, 5] <- 1
qmat[-6,6] <- 1
# statefig(c(1,2,2,1), qmat)  # draw it: 11 transitions!

icoef <- rbind(c(.05, .05, .06, .07, .05, .15, rep(c(.05,.07, .15), c(3,1,1))),
               c(rep(0,6), 1:5/10))
row.names(icoef) <- c("(Intercept)", "male")
icoef[1,] <- log(icoef[1,])


# icoef has the 11 intercepts for the 11 transitions, and the 5 age
#   coefs for death.  The ordering of the vector form, is, I have to admit,
#   not the most obvious; but we will more often use the matrix form
#   

# Simple models
hfit0 <- cmsh(Surv(age, state) ~1, data= test1, id=id, qmatrix=qmat,
               init= icoef[1,], iter =0)
hfit1 <- cmsh(list(Surv(age,state) ~1, 0:6 ~ male), 
               data=test1, id=id, qmatrix=qmat, init=icoef, iter=0)

# do the computation by hand
byhand <- function(data, eta, q=qmat, missmat, p0, firstobs=FALSE,debug=FALSE){
    idlist <- unique(data$id)
    nid <- length(idlist)
    phat <- matrix(0, nid, 6)  # n subjects, 4 states
    if (missing(missmat)) missmat <- diag(6)  # no errors
    for (i in 1:nid) {
        tdata <- subset(data, id== idlist[i])
        n <- nrow(tdata)
        e2 <- eta[data$id==idlist[i],]
        if (missing(p0)) phat[i, tdata$istate[1]] <- 1 #observed state
        else phat[i,] <- p0

        if (debug) browser()
        if (!firstobs && !missing(missmat)) {
            # if the initial state was estimated and there is a missmat
            #   apply that missmat at time 0
            phat[i,] <- phat[i,]* missmat[, tdata$istate[1]]
            }
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
            } else  phat[i,] <- phat[i,]* missmat[,k]
        }
    }       
    phat
}

eta1 <- model.matrix(hfit1) %*% coef(hfit1, matrix=TRUE)
true1 <- byhand(test1, eta1)
truelog <- sum(log(rowSums(true1)))
aeq(hfit1$loglik, truelog)

# detail=TRUE forces no iteration, and returns extra info
hfit1b <- cmsh(list(Surv(age,state) ~1, 0:6 ~ male), center=FALSE, mc.cores=1,
              data=test1, id=id, qmatrix=qmat, init=icoef, detail=TRUE)

# derivatives
eps <- 1e-8
nbeta <- sum(hfit1$cmap >0)
deriv <- double(nbeta)
for (i in 1:nbeta) {
    i2 <- icoef[icoef!=0] #treat init as a vector
    i2[i] <- i2[i]+ eps
    tfit <- cmsh(list(Surv(age,state) ~1, 0:6 ~ male), center=FALSE,
                data=test1, id=id, qmatrix=qmat, init=i2, iter=0)
    deriv[i] <- (tfit$loglik - hfit1$loglik)/eps
}

# detail returns derivatives of alpha wrt the parameters, hmm maximizes
#  sum(log(alpha)), one alpha per subject
alpha <- colSums(hfit1b$alpha)  # sum over states
tder  <- apply(hfit1b$deriv, c(1,3), sum) # ditto
logder <- unname(rowSums(tder %*% diag(1/alpha)))
aeq(deriv, logder, tol= sqrt(eps))

# msm wants the intital rates in qmat, cmsh only uses 0 vs >0
mqmat <- qmat
mqmat[qmat>0] <- exp(icoef[1,])
mfit0 <-  msm(istate ~ age, data=test1, subject=id, death=6,
             qmatrix = mqmat, fixedpar=TRUE)
aeq(-2*hfit0$loglik, mfit0$minus2loglik)  # we agree

# now with covariates
minit <- 1:5/10
mfit1 <- msm(istate ~ age, data=test1, subject=id, death=6,
             qmatrix = mqmat, fixedpar=TRUE,
             covariates= list("1-6"=~male, "2-6"=~male, "3-6"=~male, 
                              "4-6"= ~male, "5-6"= ~male), 
             covinits= list(male=minit))

# The above loglik doesn't agree with hfit1 or truelog:
# msm centers each x, and reports results for a model using centered variables
#     Where to find that centering value is not particularly obvious,
#     attr(model.matrix(mfit1), "means") is one method
#   (they appear to use mean(x[duplicated(id)])
#  To match msm, create a recentered data set and invoke hmm with center=FALSE,
#   or check against the byhand function
eta1m <- model.matrix(mfit1) %*% coef(hfit1, matrix=TRUE)
true1m<- byhand(test1, eta1m, qmat)
truelogm <- sum(log(rowSums(true1m)))
aeq(-2*truelogm, mfit1$minus2loglik)

test1b <- test1
test1b$male <- test1b$male - attr(mfit1$data$mm.cov, "means")
hfit1c <- cmsh(list(Surv(age, state) ~ 1, 0:6 ~ male), data = test1b,
                id = id, qmatrix = qmat, init = icoef, iter = 0, center=FALSE)
aeq(mfit1$minus2loglik, -2*hfit1c$loglik)
# So Chris Jackson and I agree


# Add a second covariate, this was used in jackson.R in my earlier hmm package
i2 <- rbind("(Intercept)" = coef(hfit1, matrix=TRUE)[1,], educ=0, male=0)
dimnames(i2) <- list(c("(Intercept)", "educ", "male"), colnames(hfit1$cmap))
i2["educ", "1:3"] <- .1
i2["male", c("1:6", "2:6")] <- c(.2, .3)
# init arg expects an object that looks look like coef(hfit2, matrix=TRUE)

hfit2 <- cmsh(list(Surv(age,state) ~1, 
                  1:3 ~ educ, 1:6+ 2:6 ~ male),
             data=test1, id=id, qmatrix=qmat, init=i2, iter=0)
eta2 <- model.matrix(hfit2) %*% i2
true2 <- byhand(test1, eta2)
aeq(hfit2$log, sum(log(rowSums(true2))))

# cmsh centered the covariates internally *and* also transformed the 
#  coefficients, i.e., the user never sees the change.  msm on the
#  other hand returns coefs wrt recentered data
mfit2 <- msm(istate ~ age, data=test1, subject=id, 
              qmatrix = mqmat, fixedpar=TRUE, death=6, 
              covariates= list("1-3"= ~educ,  "1-6"= ~male, "2-6"= ~male),
              covinits= list(educ=.1, male=c(.2, .3)))
eta2m <- model.matrix(mfit2) %*% coef(hfit2, matrix=T)
true2m <- byhand(test1, eta2m)
aeq(mfit2$minus2loglik, -2*sum(log(rowSums(true2m))))

# To match msm we need to precenter our data to match it
test1b <- test1
center <- attr(mfit2$data$mm.cov, "means")
test1b$educ <- test1b$educ - center["educ"]
test1b$male <- test1b$male - center["male"]

hfit2b <- cmsh(list(Surv(age,state) ~1, 
                  1:3 ~ educ, 1:6+ 2:6 ~ male), center=FALSE,
             data=test1b, id=id, qmatrix=qmat, init=i2, iter=0)
aeq(-2*hfit2b$loglik, mfit2$minus2loglik)  


# Repeat with misclassification probabilities
# Here is a miss function for 6 states and fixed probs
e1 <- .12  # an A- as A+ or vice versa
e2 <- .2   # an N- as N+ or vice versa
temp <- outer(c("Acorrect"= (1-e1), "Afalse"= e1), 
              c("Ncorrect"=(1-e2), "Nfalse"= e2), '*')

missmat <- rbind(temp[c(1,2,3,4)], temp[c(2, 1, 4,3)],
                 temp[c(3,4,1,2)], temp[c(4, 3, 2, 1)])
missmat <- cbind(missmat, 0)
missmat <- rbind(missmat, c(0,0,0,.1,.9))
missmat <- cbind(rbind(missmat,0),0)
missmat[6,6] <- 1
dimnames(missmat) <- list(true= sname, obs=sname)

# For cmsh death is always part of Surv, the marker variable(s) for
#  other states treat it as missing.
test1b$death <- 1*(test1b$state=="death")
dx <- with(test1b, ifelse(state=="death", NA, state))
test1b$dx <- factor(dx, 2:6, sname[1:5])

# treat initial state as random from 1-4
iprob <- c(1,1,1,1,0,0)/4
hfit3 <- cmsh(list(Surv(age,death) ~1, 
                  1:3 ~ educ, 1:6+ 2:6 ~ male), scale=FALSE,
               mc.cores=1,
               data= test1b, id= id, qmatrix= qmat, iprob= iprob,
               center=FALSE, iter=0, init=i2,
               marker= state:dx ~0 /discrete(init=missmat[1:5,1:5]))

eta3 <- model.matrix(hfit3) %*% coef(hfit3, matrix=TRUE)
true3 <- byhand(test1b, eta3, missmat=missmat, p0=iprob)
aeq(hfit3$log, sum(log(rowSums(true3))))

mfit3 <- msm(istate ~ age, data=test1, subject= id, 
             qmatrix = mqmat, fixedpar=TRUE, death=6,
             ematrix=missmat, initprob=c(1,1,1,1,0,0)/4,
             covariates= list("1-3"= ~educ,  "1-6"= ~male, "2-6"= ~male),
             covinits= list(educ=.1, male=c(.2, .3)))
eta3m <- model.matrix(mfit3) %*% coef(hfit3, matrix=TRUE)
aeq(eta3m, eta3)  # verifies that test1b is properly centered
aeq(hfit3$log, mfit3$minus2loglik/ -2)


# Now treat the initial state as known, initprobs still has to be set
#  though it isn't used
test1$firstobs <- 1*(!duplicated(test1$id))
mfit4 <- msm(istate ~ age, data=test1, subject= id, 
             qmatrix = mqmat, fixedpar=TRUE, death=6,
             ematrix=missmat, obstrue=firstobs,
             initprob=c(1,1,1,1,0,0)/4, 
             covariates= list("1-3"= ~educ,  "1-6"= ~male, "2-6"= ~male),
             covinits= list(educ=.1, male=c(.2, .3)))
test <- byhand(test1b, eta3, missmat=missmat)

# create the necessary status for cmsh, first obs and death have a state,
#  others unknown
keep <- test1$state=="death" | !duplicated(test1$id)
temp <- ifelse(keep, as.numeric(test1$state), 1)
test1b$fstate <- factor(temp, 1:7, levels(test1$state))

hfit4 <- cmsh(list(Surv(age, fstate) ~1, 
                  1:3 ~ educ, 1:6+ 2:6 ~ male), scale=FALSE,
               mc.cores=1,
               data= test1b, id= id, qmatrix= qmat,
               center=FALSE, iter=0, init=i2,
               marker= state:dx ~0 /discrete(init=missmat[1:5,1:5]))
true4 <- byhand(test1b, eta3, missmat=missmat, firstobs=TRUE)
aeq(hfit4$log, sum(log(rowSums(true4))))
aeq(-2*hfit4$loglik, mfit4$minus2loglik)
# for the above, I haven't yet figured out exactly what msm is doing

if (FALSE) {
    # if we do a subset, the x means change for msm, so drop covariates
    mfit4x <- msm(istate ~ age, data=test1[1:2,], subject= id, 
             qmatrix = mqmat, fixedpar=TRUE, death=6,
             ematrix=missmat, obstrue=firstobs,
             initprob=c(1,1,1,1,0,0)/4)

    hfit4x <- cmsh(Surv(age, fstate) ~1, scale=FALSE,
               mc.cores=1, 
               data= test1b[1:2,], id= id, qmatrix= qmat,
               center=FALSE, iter=0, init=i2[1,],
               marker= state:dx ~0 /discrete(init=missmat[1:5,1:5]))

    # double check that coefs are the same
    t1 <- qmat
    t1[qmat>0] <- coef(hfit4x)
    aeq(t1, mfit4x$Qmatrices$logbaseline)

    Rmat <- mqmat
    diag(Rmat) <- -rowSums(Rmat)
    P1 <- expm(Rmat * diff(test1$age)[1])
    alpha1 <- c(0,0,0,1,0,0) %*% P1
    alpha2 <- alpha1 * missmat[,4]
    all.equal(log(sum(alpha2)), hfit4x$log)

    jprob <- exp(mfit4x$minus2loglik/ -2)
    c("alpha from cmsh"= sum(alpha2), "alpha from msm" = jprob)
    # I've sent an email to Chris Jackson
}
if (FALSE){
# the code for Chris
if (FALSE){
library(msm)
library(survival)  # for statefig
library(Matrix)    # for expm

test <- data.frame(id= c(1,1), age=c(86.1, 87.4), state=rep("A+N+",2),
                   istate=c(4,4), firstobs= c(4,0))
states <- c("A-N-", "A+N-", "A-N+", "A+N+", "dementia", "death")
qmat <- matrix(0, 6,6, dimnames= list(from= states, to= states))
qmat[1,2:3] <- 1
qmat[2:3, 4] <- 1
qmat[3:4, 5] <- 1
qmat[-6,6] <- 1
statefig(c(1,2,2,1), qmat)  # draw it: 11 transitions!

icoef <- c(.05, .05, .06, .07, .05, .15, rep(c(.05,.07, .15), c(3,1,1)))
q2 <- qmat
q2[q2>0] <- icoef

e1 <- .12  # an A- as A+ or vice versa
e2 <- .2   # an N- as N+ or vice versa
temp <- outer(c("Acorrect"= (1-e1), "Afalse"= e1), 
              c("Ncorrect"=(1-e2), "Nfalse"= e2), '*')

missmat <- rbind(temp[c(1,2,3,4)], temp[c(2, 1, 4,3)],
                 temp[c(3,4,1,2)], temp[c(4, 3, 2, 1)])
missmat <- cbind(missmat, 0)
missmat <- rbind(missmat, c(0,0,0,.1,.9))
missmat <- cbind(rbind(missmat,0),0)
missmat[6,6] <- 1
dimnames(missmat) <- list(true= states, obs= states)

mfit4 <- msm(istate ~ age, data=test, subject= id, 
             qmatrix = q2, fixedpar=TRUE, death=6,
             ematrix=missmat, obstrue=firstobs,
             initprob=c(1,1,1,1,0,0)/4)
-.5* mfit4$minus2loglik

# now do it by hand
R <- q2
diag(R) <- -rowSums(R)
alpha1 <- c(0,0,0,1,0,0) %*% expm(R* 1.3)
alpha2 <- alpha1 * missmat[,4]
loglik <- log(sum(alpha2))
loglik
}
