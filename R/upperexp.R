upperexp <- function(R, eps=1e-8) {
    # this .Call routine only works if diagonal elements are unique
    #  Murphy's law says that 1-2 will occur during the iteration, so on
    #  those we punt to the slower but more robust expm routine.
    dd <- diag(R)
    if (any(diff(sort(dd)) < eps)) return(expm(R))
    
    .Call('upperexp', R, PACKAGE = "hmm")
}
