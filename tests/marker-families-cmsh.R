library(mshmm)

# End-to-end parser and likelihood smoke test with every supported marker
# family in one model.  Unit tests in marker-families.R separately verify each
# response gradient and the mixed-marker hmm2 score.
nsubject <- 8L
id <- rep(seq_len(nsubject), each=2L)
time <- rep(c(0, 1), nsubject)
visit <- rep(c(FALSE, TRUE), nsubject)
group <- rep(rep(c(0, 1), each=2L), length.out=length(id))
set.seed(20261016)
gmark <- ifelse(visit, rnorm(length(id), -0.5 + 1.5 * group, 0.8), NA)
lmark <- ifelse(visit, rlogis(length(id), -0.4 + 1.3 * group,
                              0.7 * sqrt(3) / pi), NA)
bmark <- ifelse(visit, rbeta(length(id), 2 + 3 * group, 5 - 2 * group), NA)
gamark <- ifelse(visit, rgamma(length(id), shape=2 + 3 * group,
                               rate=2), NA)
dmark <- factor(
    ifelse(!visit, NA_character_,
           ifelse(runif(length(id)) < ifelse(group == 0, 0.15, 0.80),
                  "abnormal", "normal")),
    levels=c("normal", "abnormal")
)
dtest <- data.frame(id, time, gmark, lmark, bmark, gamark, dmark)

states <- c("A0", "A1")
qmatrix <- matrix(0, 2L, 2L, dimnames=list(states, states))
qmatrix[1L, 2L] <- 0.08
statedata <- data.frame(state=states, A=states)
ginit <- matrix(c(-0.5, 0.8, 1.0, 0.9), 2L, 2L, byrow=TRUE)
linit <- matrix(c(-0.4, 0.7, 0.9, 0.8), 2L, 2L, byrow=TRUE)
binit <- matrix(c(2, 5, 5, 2), 2L, 2L, byrow=TRUE)
gainit <- matrix(c(2, 2, 5, 2), 2L, 2L, byrow=TRUE)
dinit <- matrix(c(0.85, 0.15, 0.20, 0.80), 2L, 2L, byrow=TRUE,
                dimnames=list(states, levels(dmark)))

fit <- cmsh(
    time ~ 1, data=dtest, id=id, qmatrix=qmatrix, statedata=statedata,
    iprob=c(0.8, 0.2), iter=0, mc.cores=1,
    markers=list(
        A:gmark ~ 1 / gaussian(init=ginit),
        A:lmark ~ 1 / logistic(init=linit),
        A:bmark ~ 1 / beta(init=binit),
        A:gamark ~ 1 / gamma(init=gainit),
        A:dmark ~ 1 / discrete(init=dinit)
    ),
    control=cmsh.control(center=FALSE, scale=FALSE, detail=TRUE)
)

stopifnot(inherits(fit, "cmsh"))
stopifnot(length(fit$markers) == 4L)
stopifnot(length(fit$misclass) == 1L)
stopifnot(identical(vapply(fit$markers, `[[`, "", "family"),
                    c(gmark="gaussian", lmark="logistic", bmark="beta",
                      gamark="gamma")))
stopifnot(all(is.finite(fit$alpha)))
stopifnot(all(is.finite(fit$deriv)))
stopifnot(isTRUE(all.equal(
    unname(fit$misclass$dmark$probabilities), unname(dinit),
    tolerance=1e-10, check.attributes=FALSE
)))
stopifnot(isTRUE(all.equal(
    unname(misclass(fit, "dmark")), unname(dinit),
    tolerance=1e-10, check.attributes=FALSE
)))
stopifnot(length(coef(fit)) > length(coef(fit, markers=FALSE)))
stopifnot(!any(names(coef(fit, markers=TRUE))[fit$marker.ids] %in%
               names(coef(fit, markers=FALSE))))
stopifnot(identical(colnames(fit$markers$gamark$parameters),
                    c("shape", "rate")))
