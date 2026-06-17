# A test using the NAFLD data, which is large
library(mshmm)
library(msm)
options(mc.cores= floor(parallel::detectCores()*.8))

# Build the timeline data set
temp1 <- subset(nafld1, select=c(id, age, male, case.id))  #baseline
temp1$days <- 0
temp1$nafld <- with(temp1, 1*(!is.na(case.id) & id==case.id)) #nafld at baseline
temp1$case.id <- NULL  # no longer needed

# times <=0 in nafld3 are conditions found at or before enrollment date
# create the "number of comborbidities so far" variable
temp2 <- subset(nafld3, event %in% c("diabetes", "dyslipidemia", "htn"),
                c(id, days))      # only events of interest
temp2$days <- pmax(0, temp2$days) # set prior events to day 0

#cumulative count, including entry row
temp3 <- rbind(temp1[,c("id", "days")], temp2)
temp3 <- temp3[order(temp3$id, temp3$days),]
temp3$mc <- 1:nrow(temp3) - (1:nrow(temp3))[match(temp3$id, temp3$id)]
# keep the last of any rows that are tied on id and days (2 MC on same visit) 
temp3 <- temp3[!duplicated(temp3[,1:2], fromLast=TRUE),]

lfu <- subset(nafld1, select= c(id, futime, status)) # death or last fu
names(lfu)[2] <- "days"

ndata1 <- merge(temp1, temp3, by=c("id", "days"), all=TRUE) # other variables
ndata1 <- merge(ndata1, lfu,   by=c("id", "days"), all=TRUE)# add futime
#create the state variable, a factor 
temp <- with(ndata1, ifelse(!is.na(mc), mc+1,
                    ifelse(!is.na(status) & status==1, 5, 0)))
ndata1$state <-factor(temp, 0:5, 
                      c("alive", "0MC", "1MC", "2MC", "3MC", "death"))
ndata1$mc <- ndata1$status <- NULL  # no longer needed

# Add in the office visits. There will occassionally be labs on sequential
#  days for the same visit, but it isn't worth removing these "duplicate" visits.
# But do remove two tests on the same day.
office <- subset(nafld2, days>0, c(id,days))
office <- office[!duplicated(office),]   # tests on the same day
ndata2 <- merge(ndata1, office, by=c("id", "days"), all=TRUE)
ndata2$state <- lvcf(ndata2$id, ndata2$state)

# Add in a time-dependent age
# create a dummy data set with integer ages
# (pretend they had shown up each birthday, but we only recorded "alive")
acount <- floor(nafld1$futime/365.25) # number of extra ages to add
itemp <- unlist(lapply(acount, function(i) 0:i))
dummy <- data.frame(id= rep(nafld1$id, 1+ acount), days=round(itemp*365.25),
                    iage = itemp) 
# merge it in
ndata3 <- merge(ndata2, dummy, by=c("id", "days"), all=TRUE)
ndata3$state[is.na(ndata3$state)] <- "alive"
ndata3$age <- lvcf(ndata3$id, ndata3$age)
ndata3$iage<- lvcf(ndata3$id, ndata3$age + ndata3$iage)
ndata3$years <- ndata3$days/365.25
ndata1$years <- ndata1$days/365.25

survcheck(Surv2(years, state) ~ 1, ndata1, id=id)$transitions
survcheck(Surv2(years, state) ~1,  ndata3, id=id)$transitions

# Do a coxph fit, on age scale
# ndata1 needs to have some bits filled in
indx <- match(ndata1$id, ndata1$id)
ndata1$age <- ndata1$age[indx] + ndata1$years # current age
ndata1$male <- ndata1$male[indx]
ndata1$nafld <- ndata1$nafld[indx]
coxfit <- coxph(list(Surv2(age, state) ~ male + nafld,
                      0:5 ~ male /common,
                      1:2 + 1:3 + 1:4 ~ male + nafld/common,
                      2:3 + 2:4 ~ male + nafld/common),
                ndata1, id=id)

# Do the hazards appear linear in time?
dummy <- expand.grid(male=0:1, nafld=0:1)
csurv <- survfit(coxfit, newdata=dummy, start.time=40, p0=c(1,0,0,0,0))

# Now for interval censored
# Transition matrix 
sname <- levels(ndata3$state)[-1]  # leave out 'alive'
qmat <- matrix(0,5,5, dimnames= list(sname, sname))
qmat[1,2] <- qmat[2,3] <- qmat[3,4] <- 1
qmat[1:4,5] <- 1

# This is known to work as an intial, but it is slow
icoef <- rep(log(.01), 7)
cfit0 <- cmsh(Surv(years, state) ~ 1, ndata3, id=id,
              qmatrix=qmat, init= icoef)

# just a bit more complex
cfit1 <- cmsh(Surv(years, state) ~ iage, ndata3, id=id,
              qmatrix=qmat, init = cfit0)

# and a model of interest
cfit2 <- cmsh(list(Surv(years, state) ~ iage + male + nafld,
                   0:5 ~ male / common),  ndata3, id=id,
              qmatrix=qmat, init=cfit1)

#compare coxfit and cfit2 coefficients
t1 <- coef(coxfit, matrix=TRUE)[,c(1,2,6:10)]
t2 <- coef(cfit2, matrix=TRUE)[3:4,]
t3 <- rbind(t1[1,], t2[1,], t1[2,], t2[2,])
rownames(t3) <- c("coxph, male", "mshmm, male", "coxph, nafld", "mshmm, nafld")
round(t3, 3)



# Now fit msm
qmat2 <- qmat
qmat2[qmat2>0] <- exp(coef(cfit0))
istate <- as.integer(ndata3$state) -1L
mfit0 <- msm(istate ~ years, data=ndata3, subject= id,
            qmatrix=qmat2, deathexact=5,
            censor= 0, censor.states= 1:4)

mfit1 <- msm(istate ~ years, data=ndata3, subject= id,
            qmatrix=qmat2, deathexact=5,
            censor= 0, censor.states= 1:4,
            covariates= ~iage)

mfit2 <- msm(istate ~ years, data=ndata3, subject= id,
            qmatrix=qmat2, deathexact=5,
            censor= 0, censor.states= 1:4,
            covariates= ~iage + male + nafld,
            constraint= list(male=c(1,2,3,2,5,2,2)))
