#
# Test the subroutines for parsing set of states
#  parsecovar2 is one of more subtle routines

#library(mshmm)
#parsecovar1 <- mshmm:::parsecovar1
#parsecovar2 <- mshmm:::parsecovar2
#parsemarker1 <- mshmm:::parsemarker1
states <- c("A0N0", "A1N0", "A0N1", "A1N1", "A0N2","A1N2", "death")
qmat <- matrix(0, 7,7, dimnames=list(states, states))
qmat[1,2] <- qmat[3,4] <- qmat[5,6] <- 1  # A0 to A1
qmat[1,3] <- qmat[2,4] <- 1 # N0 to N1
qmat[3,5] <- qmat[4,6] <- 1 # N1 to N2
qmat[1:6,7] <- 1
# statefig(cbind(2,2,2,1), qmat)

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

# Build the new right hand side, in the way that hmm does
tlab <- lapply(test1$rhs, function(x) {
              attr(terms(x), "term.labels")
          })
newform <- reformulate(c(attr(terms(dform), 'term.labels'), unlist(tlab)))

# For testing, assume that x1 was a factor with 4 levels of A, B, C, D
#  This is what the column names and assign values would be
Xname <- c("(Intercept)", "age", "sex", "x1B", "x1C", "x1D", "x2", "x3", "x4")
Xassign <- c(0,1,2,3,3,3,4,5,6)

test2 <- parsecovar2(test1, statedata, dform, terms(newform), qmat, states,
                     Xname, Xassign)

vars <-  c("(Intercept)", "age", "sex", "x1", "x2", "x3", "x4")
tran <- paste(row(qmat), col(qmat), sep=':')[qmat>0]
check2 <- matrix(0, length(vars), length(tran), dimnames=list(vars, tran))
indx <- cbind(c(1,4,5,7, 1,2,6,7, 1,2,6, 1,5, 1,2,6, 1,2,6, 1,5, rep(1:3, 6)),
              c(1,1,1,1, 2,2,2,2, 3,3,3, 4,4, 5,5,5, 6,6,6, 7,7, 
                rep(8:13, each=3))) 
check2[indx] <- 1:39
check2[4,2] <- 2
all.equal(check2, test2$tmap)

# expand x1 to 3 rows
check3 <-check2[c(1,2,3,4,4,4,5,6,7),]
check3[5,1:2] <- check3[5, 1:2] + .1
check3[6,1:2] <- check3[6, 1:2] + .2
rownames(check3) <- Xname
# change to integer
check3[,] <- match(c(check3), sort(unique(c(0, check3)))) -1L
all.equal(check3, test2$cmap)


# Now look at the marker list
mlist <- list(log(pib) + log(tau) ~ A /gaussian,
              sqrt(p.tau181) ~ A/gaussian(cstd),
              whm ~ N + icvol+ sex/ gamma)

test <- parsemarker1(mlist, statedata)

# formula = the bits that need to pasted up so as to all variable in the
#  model frame
all.equal(test$formula, list(~log(pib) + log(tau), ~sqrt(p.tau181), 
                             ~whm + icvol))
all.equal(test$marker, c("log(pib)", "log(tau)", "sqrt(p.tau181)", "whm"))
all.equal(test$statecol, c(2,2,2,3))
all.equal(test$options, list(as.name("gaussian"), as.name("gaussian"),
                          expression(gaussian(cstd))[[1]], as.name("gamma")))

