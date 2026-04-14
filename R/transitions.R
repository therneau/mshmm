#
# Return the transitions table for a state variable
#
transitions <- function(id, time, state) {
    idx <- match(id, unique(id))
    oo <- order(idx, time)
    if (all(diff(oo) ==1) { # ordered data
        id1 <- duplicated(id, fromLast=TRUE)
        id2 <- duplicated(id)
        table(from=state[id1], to=state[id2], useNA="ifany")
    } else {
        id1 <- duplicated(id[oo], fromLast=TRUE)
        id2 <- duplicated(id[oo])
        table(from=(state[oo])[id1], to=(state[oo])[id2], useNA="ifany")
    }
}
