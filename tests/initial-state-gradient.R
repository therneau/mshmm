library(mshmm)

eta <- c(-0.4, 0.8)
probability <- mshmm:::initfun(4L, eta, gradient=TRUE)
analytic <- attr(probability, "gradient")
eps <- 1e-6
numeric <- vapply(seq_along(eta), function(j) {
    upper <- lower <- eta
    upper[j] <- upper[j] + eps
    lower[j] <- lower[j] - eps
    (mshmm:::initfun(4L, upper) - mshmm:::initfun(4L, lower)) / (2 * eps)
}, numeric(4L))
stopifnot(isTRUE(all.equal(analytic, numeric, tolerance=1e-7)))
stopifnot(abs(sum(probability) - 1) < 1e-12)
stopifnot(probability[4L] == 0)

