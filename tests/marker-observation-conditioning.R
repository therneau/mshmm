library(mshmm)

make_engine <- function(otype, ystate, ymarker) {
    n <- length(otype)
    qmatrix <- matrix(0L, 3L, 3L)
    qmatrix[1L, 2L] <- 1L
    qmatrix[1L, 3L] <- 2L
    qmatrix[2L, 3L] <- 3L
    response <- mshmm:::gaussian(
        list(levels=c("A0", "A1"), index=c(1L, 2L, 0L), sname="A"),
        mlevel=NULL,
        init=matrix(c(-0.4, 0.8, 1.0, 1.1), 2L, 2L, byrow=TRUE)
    )
    response$e2map <- 1:4

    env <- new.env(parent=asNamespace("mshmm"))
    env$control <- list(debug=0, smallpos=1e-10)
    env$id <- rep(1L, n)
    env$X <- matrix(1, n, 1L)
    env$dtime <- if (n > 1L) rep(0.8, n - 1L) else numeric()
    env$nparm <- c(3L, 4L, 0L)
    env$nlp <- c(3L, 4L, 0L)
    env$nstate <- 3L
    env$iprob <- matrix(c(0.8, 0.2, 0), nrow=1L)
    env$p0fixed <- NULL
    env$nmarker <- 1L
    env$otype <- otype
    env$ystate <- ystate
    env$ymarker <- data.frame(marker=ymarker)
    env$rlist <- list(response)
    env$eta.beta1 <- function(x) diag(3L)
    env$eta.beta2 <- function(x) diag(4L)
    env$qmatrix <- qmatrix
    env$exactabsorb <- 3L
    env$e1 <- 1:3
    env$e3 <- integer()

    h1 <- mshmm:::hmm1
    h2 <- mshmm:::hmm2
    environment(h1) <- env
    environment(h2) <- env
    list(hmm1=h1, hmm2=h2)
}

check_engine <- function(engine, B, tolerance=3e-5) {
    analytic <- engine$hmm2(1L, B)
    score <- rowSums(analytic$deriv) / sum(analytic$alpha)
    eps <- 1e-6
    numeric <- vapply(seq_along(B), function(j) {
        upper <- lower <- B
        upper[j] <- upper[j] + eps
        lower[j] <- lower[j] - eps
        (engine$hmm1(1L, upper) - engine$hmm1(1L, lower)) / (2 * eps)
    }, numeric(1))
    stopifnot(isTRUE(all.equal(
        score, numeric, tolerance=tolerance, check.attributes=FALSE
    )))
}

B <- matrix(c(log(0.08), log(0.02), log(0.05),
              -0.4, log(0.8), 1.0, log(1.1)), nrow=1L)

# A marker history followed by exact death requires all derivative mass to be
# collected into the observed absorbing state.
check_engine(make_engine(
    otype=c(3L, 3L, 2L),
    ystate=c(0L, 0L, 3L),
    ymarker=c(-0.2, 0.6, NA_real_)
), B)

# An alive-censoring row rules out the exact absorbing death state and must
# remove both transition and marker derivative mass carried by that state.
check_engine(make_engine(
    otype=c(3L, 0L),
    ystate=c(0L, 0L),
    ymarker=c(-0.2, NA_real_)
), B)

# Conditioning on a known living state must mask marker derivative mass in all
# other states.
check_engine(make_engine(
    otype=c(3L, 1L),
    ystate=c(0L, 2L),
    ymarker=c(-0.2, NA_real_)
), B)
