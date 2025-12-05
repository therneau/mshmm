#
# Simple initial estimate of rates
#  y:    multistate Surv object
# id:    subject id
# extra: an additional amount that prevents estimates of zero, adding
#   .5 to all the numerator cells and 1 to the denominator is historical
#
simplehaz <- function(y, id, extra= .2) {
    v1 <- which(duplicated(id, fromLast=TRUE))
    v2 <- which(duplicated(id))
    nstate <- length(attr(y, 'states'))
    state <- factor(y[,2], 0:nstate, c("censor", attr(y, "states")))

    n <- length(y)
    # ignore censored obs "in the middle" of a subjects' data
    ignore <- (state[-n]== "censor" & id[-1]==id[-n])
    if (any(ignore)) {
        y <- y[!ignore]
        id <- id[!ignore]
        state <- state[!ignore]
    }

    v1 <- which(duplicated(id, fromLast=TRUE))
    v2 <- which(duplicated(id))
    # (v1,v2) are pairs of visits for a subject
    tmat <- table(state[v1], state[v2])  # matrix of transitions
    py <- tapply(y[v2,1]- y[v1,1], state[v1], sum)  # time at risk
    rate <- (tmat+ extra)/c(ifelse(is.na(py), 2*extra, py + 2*extra))
    rate[-1, -1] # don't include censor
}
