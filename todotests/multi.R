library(msm)
library(icmsh)
# A test of multinomial response

# test1 has 4 subjects, 17 rows
# The age intervals are about a year, so make transitions around 5-15% per
# year.  This makes the loglik far from 1.
sname <- levels(test1$state)[-1]   #censor is not a state
qmat <- matrix(0, 6, 6, dimnames=list(from= sname, to=sname))
qmat[1,2:3] <- 1
qmat[2:3, 4] <- 1
qmat[3:4, 5] <- 1
qmat[-6,6] <- 1
# statefig(c(1,2,2,1), qmat)
icoef <- log(c(.05, .05, .06, .07, .05, .15, rep(c(.05,.07, .15), c(3,1,1)),
              1:5/10))


# A simple fixed error function, dementia and death are never mistaken
#
e1 <- .12  # an A- as A+ or vice versa
e2 <- .2   # an N- as N+ or vice versa
temp <- outer(c("Acorrect"= (1-e1), "Afalse"= e1), 
              c("Ncorrect"=(1-e2), "Nfalse"= e2), '*')

missmat <- rbind(temp[c(1,2,3,4)], temp[c(2, 1, 4,3)],
                 temp[c(3,4,1,2)], temp[c(4, 3, 2, 1)])
dimnames(missmat) <- list(true=sname[1:4], marker=1:4)


# Fit this with msm first, it bundles death and dementia into the hidden marker
miss2 <- diag(6)
miss2[1:4, 1:4] <- missmat
dimnames(miss2) <- list(sname, sname)

qmat2 <- qmat
qmat2[qmat2!=0] <- exp(icoef) # initial estimates are in qmatrix

mfit3 <- msm(istate ~ age, data=test1, subject= id, 
             qmatrix = qmat2, fixedpar=TRUE, death=6,
             ematrix=miss2, initprob=c(1,1,1,1,0,0)/4)

# For hmm we must separate out death. We can code dementia either way, as
#  an observed IC state, or as a biomarker with no error
#  First do it with dementia in the interval censored part
# Create separate state and marker variables
temp <- c(1,1,1,1,2,3)[test1$istate]
test1$state2 <- factor(temp, 1:3, c("none", "dementia", "death"))
test1$marker1 <- c(1,2,3,4,NA, NA)[test1$istate]

alias <- data.frame(state= sname, 
                        AN = c(1,2,3,4,0, 0),
                        A  = c(0,1,0,1,2, 2),
                        N  = c(0,0,1,1,NA,NA))

hfit3 <- hmm(Surv(age, state2) ~ 1, data=test1, id=id, mc.cores=1,
             qmatrix=qmat, init=icoef, statedata =alias,
             markers= list(AN(1:4):marker1 ~0 /multinomial(init=missmat)))
             


hfit3 <-  hmm(hbind(age, state) ~ 1, data=test1, mc.cores=3,
               id = id, qmatrix = qmat, rfun=hmiss,
               pfun=init6, mfun=hmmtest, mpar=list(fn="hmmloglik"),
               otype= otype, death=6)
all.equal(unname(-2*hfit3$loglik[2]), mfit3$minus2loglik)

# Now do this as a set of 3 response
# The combined states are 1=A-N-, 2=A+N-, 3=A-N+, 4=A+N+, 5=dementia, 6=death
# It is not necessary to ignore A and/or N for a demented subject, but we do
#  so as to exactly match the prior analysis.
astate <- c(1,2,1,2, NA, NA)[test1$state]
nstate <- c(1,1,2,2, NA, NA)[test1$state]
cstate <- c(NA, NA, NA, NA, 1,2)[test1$state]

# error1 will be called with a vector of 1s and 2s, NA values are not
#  passed to it.  It returns a column for each y value with nstate rows
#  (which I know is 6).  Each colum contains prob(y | true state).
# The mlogit function returns a matrix with n rows and 2 columns, the first
#  has 1/(1+exp(eta) and the second exp(eta)/(1+exp(eta)).  Per how
#  we set up rcoef the second element is the smaller one, i.e., the error rate.
# Looking at error1 below, suppose someone has y=1 = A-.  The states are set
#  up as A-N-, A+N-, A-N+, A+N+, dementia, death.  So we would want that
#  column to be (large probability, small, large, small, 0, 0).  The last two
#  elements are 0 because someone who is demented/dead will never have y=1.
#
error1 <- function(y, nstate, eta, gradient=FALSE) {
    if (is.vector(eta)) eta <- matrix(eta, nrow=1)
    temp <- mlogit(eta, gradient)
    indx <- cbind(1:length(y), y)
    ptemp <- temp[indx]  # pick col 1 or 2 from each row
    rbind(ptemp, 1-ptemp, ptemp, 1-ptemp, 0, 0)
}

error2 <- function(y, nstate, eta, gradient) {
    if (is.vector(eta)) eta <- matrix(eta, nrow=1)
    temp <- mlogit(eta, FALSE)
    ptemp <- temp[cbind(1:length(y), y)]  # pick col 1 or 2 from each row
    rbind(ptemp, ptemp, 1-ptemp, 1-ptemp, 0, 0)
} 
   
error3 <- function(y, nstate, ...) {
     temp <- matrix(0L, nrow=nstate, ncol=length(y))
     temp[5,] <- ifelse(y==1, 1, 0)
     temp[6,] <- ifelse(y==2, 1, 0)
     temp
}

rcoef <- data.frame(response=1:2, lp=1:2, term=0, coef=1:2,
                    init=log(c(.12/.88, .2/.8)))

#double check the error matrices
temp1 <- error1(astate, 6, matrix(rcoef$init[1], length(astate),1))
temp2 <- error2(nstate, 6, matrix(rcoef$init[2], length(nstate),1))
temp3 <- error3(cstate, 6)
temp4 <- ifelse(is.na(temp1), 1, temp1) * ifelse(is.na(temp2), 1, temp2)
temp5 <- hmiss(test1$state)
all.equal(temp5[,test1$istate<5], temp4[,test1$istate < 5],
          check.attributes=FALSE)
all.equal(temp5[,test1$istate>4], temp3[,test1$istate > 4])

# Now compute
hfit4 <-  hmm(hbind(age, astate, nstate, cstate) ~ 1, data=test1, mc.cores=1,
              id = id, qmatrix = qmat, 
              rfun = list(error1, error2, error3), rcoef=rcoef,
              pfun=init6, mfun=hmmtest, mpar=list(fn="hmmloglik"),
              otype= otype, death=6)
all.equal(hfit3$loglik, hfit4$loglik)


# The gradient version uses the hmulti function.
#  There are no parameters for error3 and hence no derivative is needed.
err1 <- function(y, nstate, eta, gradient) {
    statemap <- matrix(c(1,2,1,2,0,0,2,1,2,1,0,0), ncol=2)
    hmulti(y, nstate, eta, gradient, statemap)
}

err2 <- function(y, nstate, eta, gradient) {
    statemap <- matrix(c(1,1,2,2,0,0,2,2,1,1,0,0), ncol=2)
    hmulti(y, nstate, eta, gradient, statemap)
}

etest <- err1(1:2, 6, matrix(-(1:2), ncol=1), gradient=TRUE)

hfit5 <-  hmm(hbind(age, astate, nstate, cstate) ~ 1, data=test1, mc.cores=1,
              id = id, qmatrix = qmat, 
              rfun = list(err1, err2, error3), rcoef=rcoef,
              pfun=init6, mfun=hmmtest, mpar=list(fn="hmmboth"),
              otype= otype, death=6)
all.equal(hfit5$loglik, hfit4$loglik)

