library(hmm)

# test the error matrix
# 
tmat <- cbind(c(-1, 1,0,0), c(1, -1, 5, 0), c(2,3, -1, 0), c(0,4,6,0),
              c(0,0,0,-1))
dimnames(tmat) = list(c("A", "B", "C", "dead"),
                      c(letters[1:4], 'dead'))
# The above sets up an error matrix where the first linear predictor
#  maps to the A:b and B:a (truth:observed) errors, and A:a, B:b, C:d and
#  dead:dead are the reference mappings. Observed states c and d both
#  occur with true state C.
# There are 6 linear predictors.

temp <- hmmesetup(tmat)
eta <- -rbind(1:6, 2:7, c(3,3,0,0,1,1), c(2,2,2,2,0,1), c(0,1,0,1,0,3))

# Start with the simplest statemap
err1 <- hmmemat(y=1:5, nstate= 4, eta, FALSE, statemap=1:4, 
                 setup=temp)
# expand it a bit to see if statemap works
err2 <-  hmmemat(y=1:5, nstate=6, eta, FALSE, statemap=c(1,2,2,3,3,4), temp)
all.equal(err2, err1[c(1,2,2,3,3,4),])

# but are the numbers right?
# true state 1
err1 <- hmmemat(y=1:5, nstate= 4, eta, gradient=TRUE, statemap=1:4, 
                 setup=temp)
e1 <- tmat[1, tmat[1,]>0]
e2 <- mlogit(eta[, e1, drop=FALSE], gradient=TRUE)
all.equal(err1[1,], c(e2[1,1], e2[2,2], e2[3,3], 0, 0))
g1 <- attr(err1, 'gradient')[1,,]  #state 1
g2 <- attr(e2, 'gradient')
all.equal(g1[,1:2], rbind(g2[1,1,], g2[2,2,], g2[3,3,], 0, 0))

# true state 2
e1 <- tmat[2, tmat[2,] >0]
e2 <- mlogit(eta[, e1], TRUE)
all.equal(err1[2,], c(e2[1,2], e2[2,1], e2[3,3], e2[4,4], 0))
g1 <- attr(err1, 'gradient')[2,,]  #state 2
g2 <- attr(e2, 'gradient')
all.equal(g1[,c(1,3,4)], rbind(g2[1,2,], g2[2,1,], g2[3,3,], g2[4,4,], 0))

# true state 3
e1 <- tmat[3, tmat[3,] >0]
e2 <- mlogit(eta[, e1], TRUE)
all.equal(err1[3,], c(0, e2[2,2], e2[3,1], e2[4,3], 0))
g1 <- attr(err1, 'gradient')[3,,]  #state 3
g2 <- attr(e2, 'gradient')
all.equal(g1[, 5:6], rbind(0, g2[2,2,], g2[3,1,], g2[4,3,], 0))

# True state 4
all.equal(err1[4,], c(0,0,0,0,1))
all(attr(err1, "gradient")[4,,] ==0)
