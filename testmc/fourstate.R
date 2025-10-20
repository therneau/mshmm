#library(hmm)
source("../trial/loadall.R")

witse <- "/projects/bsi/neuro/s101846.adir/Jack/R01-AG11378-2017/hmm/data/"
#load(paste0(witse, "/mcsa_agefill_20170105.rda"))
load(paste0(witse, "/mcsa_agefill_20170123.rda"))  # updated for fewer NA

# The above load the data set mcsa.agefill
# The state variable is multifold
#
#  Make a smaller one that will print as one line per subject
tdata <- subset(mcsa.agefill, ,c(clinic, agevis, clinical.status, male, apoepos,
                                 abnormal.thickness, pib.ratio,
                                 abnormal.tau.ratio, obs.type, iage, state))
names(tdata) <- c("clinic", "age", "cstate", "male", "apoe","thickness", "pib",
                  "tau", "obs.type", "iage", "state")
#tdata$cstate <- ifelse(is.na(tdata$cstate), "", tdata$cstate)

#  Make the "time since enrollment" variable
temp <- subset(tdata, obs.type==3)  #enrollment date
indx <- match(tdata$clinic, temp$clinic)
tdata$etime  <- ifelse(tdata$age < temp$age[indx], 0,  
                       pmax(0, 5+ round(temp$age[indx]- tdata$age)))

# Estimate starting values using the non-censored data
# CN to death and dementia to death
ndata <- subset(tdata, obs.type > 0)
first <- !duplicated(ndata$clinic)
pstate <- c("", ndata$cstate[-nrow(ndata)])  #prior state
pthick <- c(0, ndata$thickness[-nrow(ndata)])
table2(prior= pstate[!first], current= ndata$cstate[!first])

# This tells me that transitions from dementia to CN are possible
# Get some rates
delta <- c(0, diff(ndata$age))
ps2 <- as.numeric(factor(pstate, c("CN", "MCI", "Dementia","Dead")))

pfit1.2 <- glm(thickness ~ iage + offset(log(delta)), poisson,
               data=ndata, subset=(pthick==0 & ps2 <3 & !first))
pfit2.3 <- glm((cstate=="Dementia") ~ iage + offset(log(delta)), poisson,
               data=ndata, subset=( ps2 <3 & !first))
pfit12.4 <- glm((cstate=="Dead") ~ iage + male + etime +
                offset(log(delta)), poisson,
                data=ndata, subset=(ps2 <3 & !first))
pfit3.4  <- glm((cstate=="Dead") ~ iage + offset(log(delta)), poisson,
                data=ndata, subset=(ps2==3 & !first))


# Check that detath rates make sense
# etime is an "enrollment credit": your death rate is lower in the first
#  years after enrollment, with a linear decrease out to 5 years.
if (FALSE) {  #skip this for batch mode
ptest <- predict(pfit12.4, newdata=data.frame(iage=50:95, etime=0, delta=1, 
                                           male=0), type="response")
matplot(50:95, cbind(ptest, survexp.mn[as.character(50:95),,"2000"]*365.25),
        log='y', lty=1, col=1:3,
        xlab="Age", ylab="Death rates per year")
# In an ideal world 1 and 3 would overlap on the graph.
}

sname <- c("N-", "N+", "Dementia", "Death")
qmat <- matrix(0, 4,4,
               dimnames=list(from=sname, to=sname))
qmat[,4] <- qmat[1,2] <- qmat[2,3] <- 1

# The model will be 1 + iage + male + etime
qcoef <- data.frame(state1 = c(1,1, 2,2, 1,1,1,1, 2,2,2,2, 3,3),
                    state2 = c(2,2, 3,3, 4,4,4,4, 4,4,4,4, 4,4),
                    term   = c(0,1, 0,1, 0,1,2,3, 0,1,2,3, 0,1),
                    coef   = c(1,5, 2,6, 3,7,8,9, 3,7,8,9, 4,10),
                    init   = c(coef(pfit1.2), coef(pfit2.3), 
                               coef(pfit12.4), coef(pfit12.4), coef(pfit3.4)))

# Two response functions: we see N- or N+
r.thick <- function(y, nstate, eta, ...) {
    # the rows are N-, N+, dementia, death,
    # the columns are normal or abnormal thickness
    p1 <- exp(eta[1,1])/(1+ exp(eta[1,1])) # An N- as abnormal
    p2 <- exp(eta[1,2])/(1+ exp(eta[1,2])) # An N+ as normal
    p3 <- exp(eta[1,2])/(1+ exp(eta[1,2])) # dementia as abnormal
    emat <- rbind(c(1-p1,  p1), 
                  c(p2,    1-p2),
                  c(1-p3,   p3),
                  c(0,      0))
browser()
    emat[,y,drop=FALSE]
}

# We see clinical state
r.clin <- function(y, nstate, eta, ...) {
    # the rows are N-, N+, dementia, death,
    # the columns are CN, MCI, dementia, death
    temp <- exp(eta[1,])
    temp[6] <- temp[5]/5 # if they are really demented, chance of "CN" is low
    emat <-rbind(c(1, temp[1], temp[2], 0)/(1 + temp[1] + temp[2]),
                 c(temp[3], 1, temp[4], 0)/(1 + temp[1] + temp[3]),
                 c(temp[6], temp[5], 1, 0)/(1 + temp[2] + temp[4]),
                 c(0,0,0,1))
    emat[,y, drop=FALSE]
}

rcoef <- data.frame(response=c(1,1,1,2,2,2,2,2),
                    lp = c(1,2,3,1,2,3,4,5),
                    term =0,
                    coef=1:8,
                    init = -3)  # start at 5% errors

e4 <- function(nstate, ...) c(.95, .04, 0, 0)

c2 <- with(tdata, factor(cstate, levels=c("CN", "MCI", "Dementia", "Death")))

if (FALSE) {
    #test case
    hfit0 <- hmm(hbind(age, thickness+1, c2) ~ iage + male + etime,
                 data=tdata, qmatrix=qmat, qcoef=qcoef, id=clinic, death=4,
                 otype= obs.type, entry=c(1,1,.1, 0), mc.cores=1,
                 rfun=list(r.thick, r.clin), rcoef=rcoef, pfun=e4,
                 mfun= hmmtest, debug=1)
}
hfit4a <- hmm(hbind(age, thickness+1, c2) ~ iage + male + etime,
              data=tdata, qmatrix=qmat, qcoef=qcoef, id=clinic, death=4,
              otype= obs.type, entry=c(1,1,.1, 0), mc.cores=20,
              rfun=list(r.thick, r.clin), rcoef=rcoef, pfun=e4,
              mfun= optim, 
              mpar=list(control=list(fnscale= -1, maxit=2000)))

save(tdata, hfit4a, file="hfit4.rda")                         

#
# Use this as starting estimates for the mcmc
#
qcoef2 <- qcoef
temp <- hfit4a$beta[,1:5]
qcoef2$init <- temp[temp!=0]

rcoef2 <- rcoef
rcoef2$init <- hfit4a$beta[1, 6:12]

hfit4 <- hmm(hbind(age, thickness+1, c2) ~ iage + male + etime,
             data=tdata, qmatrix=qmat, qcoef=qcoef, id=clinic, death=4,
             otype= obs.type, entry=c(1,1,.1, 0), mc.cores=20,
             rfun=list(r.thick, r.clin), rcoef=rcoef, pfun=e4,
             mfun=mcmc2,
             mpar=list(burnin=500, iter=2000, prior.std=4,reset=100, 
                       print=500, cmap=quote(cmap)))

save(tdata,hfit4a, hfit4, file="hfit4.rda")
