#
# Check the loglik for a final censoring, i.e., that they are counted as
#  "not dead", not as likelihood of 1
#
library(mshmm)

# semi-competing risks data
ddata <- data.frame(id=c(1,1,2,2,3,3), 
                    time=c(0,1, 0,2, 0,3), state=c(1,2,1,3,1,0), x=1)
ddata$state <- factor(ddata$state, 0:3, c("censor", "entry", "recur", "death"))
ddata$istate <- c(4,1,2,3)[as.numeric(ddata$state)] #integer version

qmat <- rbind(c(0,1,1), c(0,0,1), c(0,0,0))
qmatm <- rbind(c(0, .1, .2), c(0,0,.3), c(0,0,0))
states <- c("entry", "recur", "death")
dimnames(qmat) <- list(states, states)
dimnames(qmatm) <- list(states, states)
# The qmat for msm has the intial values in the state

# First, by hand, note the msm has the initial rates in qmat, not coefs
hmat <- rbind(c(-.3, .1, .2),
              c(0, -.3,  .3), 0)
p1 <- c(1,0,0) %*% as.matrix(expm(hmat))
p2 <- p1 %*% as.matrix(expm(hmat))
p3 <- p2 %*% as.matrix(expm(hmat))

alpha1 <- p1*c(0,1,0)
alpha2 <- c(0, 0, p2 %*% hmat[,3])  # death is exact
alpha3 <- p3*c(1,1,0)  # not dead
ptot <- c(sum(alpha1), sum(alpha2), sum(alpha3))
loglik1 <- sum(log(ptot))  # the true loglik (by hand)

icoef <- matrix(log(c(.1, .2, .3)), nrow=1, 
                dimnames=list("(Intercept)", c("1:2", "1:3", "2:3")))
hfit <- hmm(Surv(time, state) ~ 1, id= id, qmatrix=qmat,
            data=ddata, init=icoef, iter=0, mc.cores=1)
all.equal(hfit$loglik, loglik1)
