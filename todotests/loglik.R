library(mhsmm)
aeq <- function(x, y, ...) all.equal(as.vector(x), as.vector(y), ...)

# Do the loglik by hand for a small data set that has per-subject initial
#  state.
#
mdata <- structure(list(ptnum = c(1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 2L, 
2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 3L, 
3L, 3L, 3L, 3L, 3L, 3L, 4L, 4L, 4L, 4L, 4L, 4L, 4L, 4L, 4L, 5L, 
5L, 5L, 5L, 5L, 5L, 5L, 5L, 5L, 5L, 5L, 5L, 5L, 5L, 5L, 5L, 5L, 
5L, 5L, 5L, 5L, 6L, 6L, 6L, 6L, 6L, 6L, 6L, 6L, 6L, 6L, 6L, 6L, 
6L, 6L, 6L), age = c(87.67, 88, 88.98, 89, 90, 90.27, 90.34, 
90.97, 86.8, 87, 88, 88.07, 89, 90, 91, 92, 93, 94, 95, 96, 97, 
98, 99, 99.54, 90.08, 91, 91.66, 92, 92.98, 93, 93.39, 86.34, 
87, 87.59, 88, 89, 89.27, 90, 91, 91.01, 87.09, 88, 88.32, 89, 
89.56, 90, 90.9, 91, 92, 92.25, 93, 93.51, 94, 94.92, 95, 96, 
97, 98, 99, 99.73, 99.77, 83.69, 84, 84.91, 85, 86, 86.14, 87, 
87.44, 88, 88.66, 89, 89.9, 90, 90.86, 90.89), male = c(0, 0, 
0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 
1, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 
1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 
1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1), dx = structure(c(1L, 6L, 2L, 
6L, 6L, 2L, 5L, 4L, 2L, 6L, 6L, 3L, 6L, 6L, 6L, 6L, 6L, 6L, 6L, 
6L, 6L, 6L, 6L, 4L, 2L, 6L, 3L, 6L, 3L, 6L, 4L, 1L, 6L, 1L, 6L, 
6L, 2L, 6L, 6L, 4L, 1L, 6L, 1L, 6L, 1L, 6L, 1L, 6L, 6L, 1L, 6L, 
1L, 6L, 1L, 6L, 6L, 6L, 6L, 6L, 5L, 4L, 1L, 6L, 1L, 6L, 6L, 1L, 
6L, 1L, 6L, 1L, 6L, 1L, 6L, 5L, 4L), levels = c("CU", "MCI", 
"dementia", "death", "CU/MCI", "censor"), class = "factor")), row.names = c(NA, 
-76L), class = "data.frame")

mdata$iage <- floor(mdata$age)
mdata$death <- 1L*(mdata$dx == "death")

# Minnesota death rates, used as initial values
dummy <- data.frame(y= c(survival::survexp.mn[60:95, 1:2, "2010"]*365.25),
               age= rep(60:95, 2), male=rep(1:0, each=36))
dfit <- lm(log(y) ~ age + male, dummy)

# The X matrix and coefficients matrix; the dementia rates have a higher slope
#  than death, but are still a bit below death rates at age 90. Male effect
#  is is 1/10 as large.  The values below don't have to be correct, just
#  sensible enough to keep exp(eta) sensible

X <- cbind(1, mdata$iage, mdata$male)
dcoef <- coef(dfit)
beta <- rbind(dcoef[1]+ c(-7.5, -16, -.1, 0, .1),
              dcoef[2]+ c(.1, .2, 0,0,0),
              dcoef[3] * c(0,0, 1, 1, 1))
dimnames(beta) <- list(c("(Intercept)", "iage", "male"),
                       c("1:2", "2:3", "1:4", "2:4", "3:4"))
eta <- X %*% beta

# The latent state is cognitive loss of 0, 1, 2
states <- c("C0", "C1", "C2", "death")
Q <- matrix(0, 4, 4, dimnames=list(from=states, to=states))
qindex <- cbind(c(1,2, 1,2,3), c(2,3,4,4,4)) # nonzero transtions 
Q[qindex] <- 1

# Error matrix
emat <- rbind(c(.8, .2, 0, 0, 1),
              c(.15, .65, .2, 0, .8),
              c( 0,  .1,  .9, 0,  .1),
              c( 0,  0 ,   0, 1, 0))
dimnames(emat) <- list(true=states, 
                    observed=c("CU", "MCI", "demented", "death", "CU/MCI"))

# Prob of entering in the CU, MCI, dementia, death states
entry <- matrix(0, 6, 4)
entry[,2] <- .20 + .015*(mdata$age[!duplicated(mdata$ptnum)]-80)
entry[,1] <- 1-entry[,2]

# The matrix of state probabilities for a subject
phat <- function(id, nstate=4) {
    eta2 <- eta[mdata$ptnum== id,, drop=FALSE]  # eta for this set of subjects
    pdata <- subset(mdata, ptnum==id)
    n <- nrow(pdata)
    alpha <- matrix(0,n, 4, dimnames=list(NULL, states))
    y <- as.numeric(pdata$dx)
    delta <- diff(pdata$age)  # time intervals
    alpha[1,] <- entry[id,] * emat[, y[1]]

    for (i in 2:n) {
        Hmat <- matrix(0, 4, 4)
        eta <- q1$init  # no covariates
        Hmat[qindex] <- exp(eta2[i-1,])
        diag(Hmat) <- -rowSums(Hmat)
        Pmat <- expm(Hmat* delta[i-1]) # probability transition matrix
        
        if (pdata$otype[i] == "usual") D <- diag(emat[,y[i]])
        else if (pdata$otype[i] =="exact") D <- diag(c(exp(eta2[i,3:5]),0))
        else D <- diag(nstate)  # censored

        alpha[i,] <- alpha[i-1,] %*% (Pmat %*% D)
    }
    alpha
}

# Look at the first few steps for id 1:
# The dx at entry is CU.  The chance of (true CO, saw CU) is .685 * .8= .549,
#  of (true MCI, saw CU) is .315 * .15= .047, and 0 for dem or death since
#  they can't have enrolled in those states.
# So alpha= (.55, .047, 0, 0) and the loglik is log(.55 + .047).

alpha1 <- entry[1,] * emat[,1]

# Next they have .33 years of time at age 87 followed by .98 years at age 88
#  before the next event, which is an MCI; column 2 of the error map
pmat <- function(t, eta) {
    H <- matrix(0, 4,4)
    H[qindex] <- exp(eta)
    diag(H) <- -rowSums(H)
    expm(H*t)
}
P1 <- pmat(.33, eta[1,])
P2 <- pmat(.98, eta[2,])
alpha2 <- alpha1 %*% (P1 %*% P2) * emat[,2]

# Then .02 years at 88, 1 year at 89 and .27 at age 90, followed by MCI
P3 <- pmat(.02, eta[3,])
P4 <- pmat(1, eta[4,])
P5 <- pmat(.27, eta[5,])
alpha3 <- (alpha2 %*% P3 %*% P4 %*% P5) *emat[,2]

# Next .07 years to CU/MCI
P6 <- pmat(.07, eta[6,])
alpha4 <- (alpha3 %*% P6) * emat[,5]

# And last, .63 years to death
P7 <- pmat(.63, eta[7,])
density <- c(exp(eta[8,3:5]),0)
alpha5 <- (alpha4 %*% P7)*density

# Fit an hmm to subject 1
# levels of the marker need to match the error matrix
mdata$dx2 <- factor(as.numeric(mdata$dx), 1:5, levels(mdata$dx)[1:5])
htest0 <- cmsh(list(Surv(age, death)~ iage,
                    0:"death" ~ male),  mdata, subset= (ptnum==1),
               mdata, id=ptnum, qmatrix=Q, iprob= entry, init= beta,
               marker= state:dx2 ~0 /discrete(init=emat))
aeq(htest0$loglik, log(sum(alpha5)))


