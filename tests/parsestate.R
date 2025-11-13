#
# Test the subroutines for parsing set of states
#  parsecovar2 is one of more subtle routines

#library(mshmm)
#parsecovar1 <- mshmm:::parsecovar1
#parsecovar2 <- mshmm:::parsecovar2
#parsemarker1 <- mshmm:::parsemarker1
#parsemarker2 <- mshmm:::parsemerker2
library(survival)
 source("ptplot.R")
 source("parsecovar.R")
 source("parsemarker.R")
 source("response.R")

states <- c("A0N0", "A1N0", "A0N1", "A1N1", "A0N2","A1N2", "death")
qmat <- matrix(0, 7,7, dimnames=list(states, states))
qmat[1,2] <- qmat[3,4] <- qmat[5,6] <- 1  # A0 to A1
qmat[1,3] <- qmat[2,4] <- 1 # N0 to N1
qmat[3,5] <- qmat[4,6] <- 1 # N1 to N2
qmat[1:6,7] <- 1
# statefig(cbind(2,2,2,1), qmat, alty= rep(1:2, c(7,6)))

statedata <- data.frame(state= states,
                        A=c(0,1,0,1,0,1,2),
                        N=c(0,0,1,1,2,2,NA),
                        D=c(0,0,0,0,0,0,1))

# a formula with lots of bits
tform <- list(Surv(time, status) ~ 1+ age,
              0:"death" ~ sex,
              A0N0:D(0) ~ x1 / common,
              A(0):A(1) ~ -age + x2,
              N(0):N(1,2) + N(1):N(2) ~ x3/ init(2),
              A0N0: c('A0N1', "A1N0") ~ x4)

dform <- tform[[1]]  # the default formula
test1 <- parsecovar1(tform[-1])

# Build the new right hand side, in the way that hmm does, icvol will come
# from the markers, below
tlab <- lapply(test1$rhs, function(x) {
              attr(terms(x), "term.labels")
          })
newform <- reformulate(c(attr(terms(dform), 'term.labels'), unlist(tlab),
                         "icvol"))

# For testing, assume that x1 was a factor with 4 levels of A, B, C, D
#  This is what the column names and assign values would be. The icvol term
#  shows up in the markder
Xname <- c("(Intercept)", "age", "sex", "x1B", "x1C", "x1D", "x2", "x3", "x4",
           "icvol")
Xassign <- c(0,1,2,3,3,3,4,5,6,7)

test2 <- parsecovar2(test1, statedata, dform, terms(newform), qmat,
                     Xname, Xassign)

vars <-  c("(Intercept)", "age", "sex", "x1", "x2", "x3", "x4", "icvol")
tran <- paste(row(qmat), col(qmat), sep=':')[qmat>0]
check2 <- matrix(0, length(vars), length(tran), dimnames=list(vars, tran))
indx <- cbind(c(1,4,5,7, 1,2,6,7, 1,2,6, 1,5, 1,2,6, 1,2,6, 1,5, rep(1:3, 6)),
              c(1,1,1,1, 2,2,2,2, 3,3,3, 4,4, 5,5,5, 6,6,6, 7,7, 
                rep(8:13, each=3))) 
check2[indx] <- 1:39
check2[4,2] <- 2
all.equal(check2, test2$tmap)

# expand x1 to 3 rows, i.e., expand tmap to cmap
check3 <-check2[c(1,2,3,4,4,4,5,6,7,8),]
check3[5,1:2] <- check3[5, 1:2] + .1
check3[6,1:2] <- check3[6, 1:2] + .2
rownames(check3) <- Xname
# change to integer
check3[,] <- match(c(check3), sort(unique(c(0, check3)))) -1L
all.equal(check3, test2$cmap)


# Now look at a marker list
mlist <- list(A(0:1):log(pib) + A(0:1):log(tau) ~ 1/ gaussian,
              A(0:1):sqrt(p.tau181) ~ 1 /gaussian,
              N:log(wmh) ~ 1/ logistic,
              N:log(wmh) ~ icvol/ logistic(param= "mean") + common)

test1 <- parsemarker1(mlist, statedata)

# the markers, needed to create the model frame
all.equal(test1$marker, c("log(pib)", "log(tau)", "sqrt(p.tau181)", 
                         "log(wmh)", "log(wmh)"))
all.equal(test1$nmarker, c(2,1,1,1))

# stateinfo has a summary of the column of statedata that was used,
# one element per marker
# A(0:1) states that for statedata$A, the '0', and '1' elements will map
#  to unique Gaussian peaks, A2 = death has no pib or tau distribution
#  levels need not be numeric, index is of length nstate and shows which
#  peak each state maps onto.
#  
temp2 <- list(sname="A", levels=0:1, index=c(1,2,1,2,1,2,0))
temp3 <- list(sname="N", levels=0:2, index=c(1,1,2,2,3,3,0))
all.equal(test1$stateinfo, list(temp2, temp2, temp2, temp3, temp3))

all.equal(test1$options, list(as.name("gaussian"), as.name("gaussian"),
                              as.name("logistic"),
                              expression(logistic(param="mean") +common)[[1]])) 

all.equal(test1$mterm, list(~1, ~1, ~1, ~icvol))


test2 <- parsemarker2(test1, statedata, terms(newform), Xname, Xassign,
                      markerlevel =NULL)
# linear predictors for each: 2 peaks * mean/std =4  for pib, tau, and ptau
all.equal(test2$rindex, list('log(pib)'=1:4, 'log(tau)'=5:8, 
                             'sqrt(p.tau181)' = 9:12, 'log(wmh)'= 13:18))
all(test2$cmap[2:9,] ==0)
all(test2$cmap[1,] == c(1:13, 15:19))
all(test2$cmap[10, c(13, 15, 17)] ==14) # shared icvol coef
all(test2$cmap[10,-c(13, 15, 17)] ==0 )

all.equal(sapply(test2$response, function(x) x$name),
          rep(c("gaussian", "logistic"), c(3,1)))


          
