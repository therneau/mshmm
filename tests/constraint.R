library(mshmm)

# Translation from user-scale predictor contrasts to optimizer-scale affine
# constraints.  Each returned row represents constant + A %*% beta >= 0.
Terms <- terms(~ x)
cmap <- matrix(1:4, nrow=2L,
               dimnames=list(c("(Intercept)", "x"), c("1:2", "2:3")))
B <- matrix(0, 2L, 2L)
cdata <- data.frame(
    x=c(0, 0, 2, 2),
    `_lp__`=c("1:2", "2:3", "1:2", "2:3"),
    `_cid_`=c(1L, 1L, 2L, 2L),
    `_cwt_`=c(-1, 1, -1, 1),
    check.names=FALSE
)
constraint <- mshmm:::hmmconstraint(cdata, Terms, cmap, diag(2L), B)
expected <- rbind(c(0, -1, 0, 1, 0), c(0, -1, -2, 1, 2))
stopifnot(isTRUE(all.equal(constraint, expected, check.attributes=FALSE)))

# Fixed coefficients move into the affine constant rather than disappearing.
fixed_map <- matrix(c(-1L, 1L, 2L, 3L), nrow=2L,
                    dimnames=dimnames(cmap))
fixed_B <- B
fixed_B[1L, 1L] <- 0.7
fixed_constraint <- mshmm:::hmmconstraint(
    cdata[1:2, ], Terms, fixed_map, diag(2L), fixed_B
)
stopifnot(isTRUE(all.equal(
    fixed_constraint, matrix(c(-0.7, 0, 1, 0), nrow=1L),
    check.attributes=FALSE
)))

# The optimizer enforces the same affine convention.  The unconstrained
# optimum is (-1, 2); beta[1] >= 0 puts the solution at (0, 2).
quadratic <- function(par, ...) {
    target <- c(-1, 2)
    list(loglik=-sum((par - target)^2) / 2,
         deriv=target - par, S=diag(2L))
}
fit <- hmmscore(c(0, 0), fn=quadratic, gr=NULL, iter=50,
                constraint=matrix(c(0, 1, 0), nrow=1L))
stopifnot(isTRUE(all.equal(fit$coef, c(0, 2), tolerance=1e-6)))
stopifnot(all(drop(matrix(c(0, 1, 0), nrow=1L) %*%
                  c(1, fit$coef)) >= -1e-8))

bad <- try(hmmscore(c(-1, 0), fn=quadratic, gr=NULL, iter=2,
                    constraint=matrix(c(0, 1, 0), nrow=1L)), silent=TRUE)
stopifnot(inherits(bad, "try-error"))

