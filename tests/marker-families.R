library(mshmm)

# Mushfiq's branch calls this component `init`; the integration branch uses
# `initial`.  The numerical test is intentionally neutral to that metadata
# naming difference.
response_initial <- function(response) {
    if (!is.null(response$initial)) response$initial else response$init
}

check_response_gradient <- function(response, y, eta, tolerance=2e-6) {
    analytic <- response$rfun(y, eta, gradient=TRUE)
    agrad <- attr(analytic, "gradient")
    stopifnot(identical(dim(agrad), c(nrow(analytic), length(y), ncol(eta))))
    eps <- 1e-6
    for (j in seq_len(ncol(eta))) {
        upper <- lower <- eta
        upper[, j] <- upper[, j] + eps
        lower[, j] <- lower[, j] - eps
        numeric <- (response$rfun(y, upper) - response$rfun(y, lower)) /
            (2 * eps)
        stopifnot(isTRUE(all.equal(
            agrad[, , j], numeric, tolerance=tolerance,
            check.attributes=FALSE
        )))
    }
    invisible(TRUE)
}

# Each latent amyloid level is repeated across three cognition levels.  This
# exercises the product-state mapping rather than a simple one-state-per-peak
# response.
stateinfo <- list(
    levels=c("A0", "A1"),
    index=c(1L, 2L, 1L, 2L, 1L, 2L, 0L),
    sname="A"
)
y <- c(-0.4, 0.7, 1.1)

gauss <- mshmm:::gaussian(
    stateinfo, mlevel=NULL,
    init=matrix(c(-0.5, 0.8, 1.2, 1.1), 2L, 2L, byrow=TRUE)
)
eta_gauss <- rbind(
    c(-0.5, log(0.8), 1.2, log(1.1)),
    c(-0.4, log(0.9), 1.1, log(1.2)),
    c(-0.3, log(1.0), 1.0, log(1.3))
)
check_response_gradient(gauss, y, eta_gauss)

logis <- mshmm:::logistic(
    stateinfo, mlevel=NULL,
    init=matrix(c(-0.5, 0.8, 1.2, 1.1), 2L, 2L, byrow=TRUE)
)
check_response_gradient(logis, y, eta_gauss)

beta_response <- mshmm:::beta(
    stateinfo, mlevel=NULL,
    init=matrix(c(2.2, 5.1, 5.4, 2.4), 2L, 2L, byrow=TRUE)
)
eta_beta <- rbind(
    log(c(2.2, 5.1, 5.4, 2.4)),
    log(c(2.4, 4.8, 5.0, 2.7)),
    log(c(2.6, 4.5, 4.7, 3.0))
)
check_response_gradient(beta_response, c(0.2, 0.5, 0.8), eta_beta)

gamma_response <- mshmm:::gamma(
    stateinfo, mlevel=NULL,
    init=matrix(c(2.2, 1.1, 5.4, 2.4), 2L, 2L, byrow=TRUE)
)
eta_gamma <- rbind(
    log(c(2.2, 1.1, 5.4, 2.4)),
    log(c(2.4, 1.2, 5.0, 2.2)),
    log(c(2.6, 1.3, 4.7, 2.0))
)
check_response_gradient(gamma_response, c(0.4, 1.5, 3.2), eta_gamma)

# A fixed continuous marker uses the same density implementation and natural-
# scale init values, but contributes no estimated marker parameters.
fixed_gauss <- mshmm:::gaussian(
    stateinfo, mlevel=NULL,
    init=matrix(c(-0.5, 0.8, 1.2, 1.1), 2L, 2L, byrow=TRUE),
    static=TRUE
)
estimated_at_init <- gauss$rfun(
    y,
    matrix(response_initial(gauss), length(y),
           length(response_initial(gauss)), byrow=TRUE)
)
stopifnot(isTRUE(all.equal(
    fixed_gauss$rfun(y), estimated_at_init,
    tolerance=1e-12, check.attributes=FALSE
)))

# Estimated categorical misclassification probabilities, including a pattern
# with equal-probability cells, have analytic probability gradients.
marker_levels <- c("normal", "borderline", "abnormal")
pattern <- matrix(
    c(1L, 2L, 2L,
      3L, 4L, 4L),
    nrow=2L, byrow=TRUE,
    dimnames=list(c("A0", "A1"), marker_levels)
)
discrete_init <- matrix(
    c(0.2, 0.4, 0.4,
      0.6, 0.2, 0.2),
    nrow=2L, byrow=TRUE,
    dimnames=list(c("A0", "A1"), marker_levels)
)
disc <- mshmm:::discrete(
    stateinfo, mlevel=marker_levels,
    pattern=pattern, init=discrete_init, static=FALSE
)
dy <- factor(c("normal", "borderline", "abnormal"), levels=marker_levels)
eta_disc <- rbind(response_initial(disc),
                  response_initial(disc) + c(0.1, -0.1),
                  response_initial(disc) + c(-0.15, 0.05))
check_response_gradient(disc, dy, eta_disc)

# A code repeated across rows represents a genuinely shared log-odds
# parameter.  Non-shared cells supply each row's reference category.
shared_pattern <- matrix(
    c(1L, 2L, 3L,
      4L, 2L, 5L),
    nrow=2L, byrow=TRUE,
    dimnames=list(c("A0", "A1"), marker_levels)
)
shared_init <- matrix(
    c(0.50, 0.25, 0.25,
      0.40, 0.20, 0.40),
    nrow=2L, byrow=TRUE,
    dimnames=list(c("A0", "A1"), marker_levels)
)
shared_disc <- mshmm:::discrete(
    stateinfo, mlevel=marker_levels,
    pattern=shared_pattern, init=shared_init, static=FALSE
)
shared_eta <- rbind(response_initial(shared_disc),
                    response_initial(shared_disc) + c(0.1, -0.1, 0.05))
check_response_gradient(
    shared_disc,
    factor(c("normal", "abnormal"), levels=marker_levels),
    shared_eta
)

fixed_disc <- mshmm:::discrete(
    stateinfo, mlevel=marker_levels,
    init=discrete_init, static=TRUE
)
fixed_prob <- fixed_disc$rfun(dy)
stopifnot(all(dim(fixed_prob) == c(7L, 3L)))
stopifnot(all(fixed_prob[7L, ] == 0))

# Mixed continuous and estimated categorical markers must propagate their
# scores through the same hmm2() update without corrupting marker indices.
simple_stateinfo <- list(levels=c("A0", "A1"), index=1:2, sname="A")
gauss2 <- mshmm:::gaussian(
    simple_stateinfo, mlevel=NULL,
    init=matrix(c(-0.5, 0.8, 1.2, 1.1), 2L, 2L, byrow=TRUE)
)
gauss2$e2map <- 1:4
disc2 <- mshmm:::discrete(
    simple_stateinfo,
    mlevel=c("normal", "abnormal"),
    init=matrix(c(0.85, 0.15, 0.20, 0.80), 2L, 2L, byrow=TRUE),
    static=FALSE
)
disc2$e2map <- 5:6

test_env <- new.env(parent=asNamespace("mshmm"))
test_env$control <- list(debug=0, smallpos=1e-10)
test_env$id <- rep(1L, 2L)
test_env$X <- matrix(1, 2L, 1L)
test_env$dtime <- 1
test_env$nparm <- c(1L, 6L, 0L)
test_env$nlp <- c(1L, 6L, 0L)
test_env$nstate <- 2L
test_env$iprob <- matrix(c(0.65, 0.35), nrow=1L)
test_env$p0fixed <- NULL
test_env$nmarker <- 2L
test_env$otype <- rep(3L, 2L)
test_env$ystate <- integer(2L)
test_env$ymarker <- data.frame(
    continuous=c(-0.2, 0.8),
    categorical=factor(c("normal", "abnormal"),
                       levels=c("normal", "abnormal"))
)
test_env$rlist <- list(gauss2, disc2)
test_env$eta.beta1 <- function(x) matrix(1, 1L, 1L)
test_env$eta.beta2 <- function(x) diag(6L)
test_env$qmatrix <- matrix(c(0L, 1L, 0L, 0L), 2L, 2L, byrow=TRUE)
test_env$exactabsorb <- integer()
test_env$e1 <- 1L
test_env$e3 <- integer()

hmm1_test <- mshmm:::hmm1
hmm2_test <- mshmm:::hmm2
environment(hmm1_test) <- test_env
environment(hmm2_test) <- test_env
B <- matrix(c(log(0.08), -0.5, log(0.8), 1.2, log(1.1),
              response_initial(disc2)), nrow=1L)
analytic <- hmm2_test(1L, B)
analytic_score <- rowSums(analytic$deriv) / sum(analytic$alpha)
eps <- 1e-6
numeric_score <- vapply(seq_along(B), function(j) {
    upper <- lower <- B
    upper[j] <- upper[j] + eps
    lower[j] <- lower[j] - eps
    (hmm1_test(1L, upper) - hmm1_test(1L, lower)) / (2 * eps)
}, numeric(1))
stopifnot(isTRUE(all.equal(
    analytic_score, numeric_score,
    tolerance=2e-5, check.attributes=FALSE
)))
