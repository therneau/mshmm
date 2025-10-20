library(hmm)
aeq <- function(x,y,...) all.equal(as.vector(x), as.vector(y), ...)
#
# Test out the hmmesetup and hmmematrix functions, including errors
#

resp <- c("CU", "MCI", "dementia", "death", "CU/MCI")
states <- c(paste0("A", 0:2), "death")

pmat <- rbind(c(-1, 1, 0, 0),
              c(2, -1, 3, 0),
              c(0, 4, -1, 0),
              c(0, 0, 0,  -1))
dimnames(pmat) <- list(states, resp[1:4])

partial <- list("CU/MCI" = c("CU", "MCI"))

test1 <- hmmesetup(pmat, states, resp, partial=partial)
all.equal(test1$lpmap,
          list(A0= 1, A1=2:3, A2=4, death=numeric(0)))
aeq(dimnames(test1$emap)[[1]], paste0("p", 1:5))
aeq(test1$emap[,,1][c(1,7,21,22)], rep(1,4))       # A0
aeq(test1$emap[,,2][c(3,6,14, 21, 23)], rep(1, 5)) # A1
aeq(test1$emap[,,3][c(10, 11, 25)], rep(1, 3))     # A2
aeq(test1$emap[1,4,4], 1)  #death
aeq(table(test1$emap), c(87,13))   # 13 total 1's in the table
aeq(test1$statemap, 1:4)

eta <- matrix(seq(1:20)/10, ncol=4)
y <- factor(c(1,1,2,3,5), 1:5, resp)

test2 <- hmmemat(as.integer(y), 4, eta, FALSE, setup=test1)


