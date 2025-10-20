qplus <- function(data, death=NULL) {
    # only one state can change at a time, and it can only get worse
    # Anything can go to death
    states <- data$state
    if (is.null(states) || any(duplicated(states))) 
        stop("data must have a column named 'state', with unique values")
    if (!is.null(death) && all(states != death)) 
        stop("death state not found")
    nstate <- length(states)
    data$state <- NULL
    
    qmat <- matrix(0., nstate, nstate,
                   dimnames=list(from=states, to=states))

    # the as.numeric allows for use of a factor
    dfun <- function(a, b) {
        if (is.na(a) | is.na(b)) 0
        else as.numeric(b) - as.numeric(a)
    }
    dtemp <- integer(ncol(data))
    for (i in 1:nrow(data)) {
        for (j in 1:nrow(data)) {
            for (k in 1:ncol(data)) dtemp[k] = dfun(data[i,k], data[j,k])
            if (sum(dtemp)==1 && all(dtemp >=0)) qmat[i,j] <- 1
        }
    }
    if (!is.null(death)) {
        qmat[,death] <- 1
        qmat[death, death] <- 0
    }
    qmat
}

noworse <- function(qmat, sdata) {
    # Find all the legal transitions
    valid <- cbind(from=row(qmat)[qmat>0], to=col(qmat)[qmat>0])
    # create a matrix of which terms changed, for each pair
    changes <- matrix(0L, nrow(valid), ncol(sdata))
    for (i in 1:ncol(sdata)) {
        temp1 <- sdata[valid[,1], i]  # starting state
        temp2 <- sdata[valid[,2], i]  # ending state
        changes[,i] <- 1L*(is.na(temp1) !=is.na(temp2) | temp1 != temp2)
    }
    
    # transitions that change a single attribute are of interest
    #  (usually they are the only ones present in the matrix)
    keep <- rowSums(changes) ==1
    valid <- valid[keep,, drop=FALSE]
    changes <- changes[keep,, drop=FALSE]
    
    #  All the transitions with a 1 in column k of changes made that change
    #   wrt the kth column of sdata, and etc for other columns
    #  For all those with a 1 in col 1, form all pairs, and count the
    #   changes for the pair.  For instance A0/N0/CN to A1/N0/CN and A0/N1/CN 
    #   and A1/N1/CN are a pair where the "noworse" constraint applies, since
    #   it is the same A0 to A1 transtion, and N1/CN is worse than N0/CN.
    #   We can't rank N0/MCI and N1/CN though.
    # Apply transitivity as well. If we have an A0 to A1 constraint for N0 vs N1
    #   and another for N1 vs N2, we don't need the N0 vs N2 one.  This slightly
    #   reduces the number of saved constraints.
    #   
    pairs <- lapply(1:ncol(changes), function(i) {
        indx <- which(changes[,i] ==1)  # starting transitions of interest
        temp <- cbind(rep(indx, each=length(indx)),
                      rep(indx, length(indx)))
        # since sdata is ordered, transition j can be "worse" than k only if k<j
        keep <- (temp[,1] < temp[,2])   # only keep the valid pairs
        temp <- temp[keep,, drop=FALSE] # pairs where the second could be worse
        # to be comparable a pair needs to start in the same value for i 
        keep <- (sdata[valid[temp[,1],1], i] == sdata[valid[temp[,2],1], i])
        temp <- temp[keep,,drop=FALSE]
        c2 <- matrix(0L, nrow(temp), ncol(sdata))
        for (j in 1:ncol(sdata))
            c2[,j] <- as.numeric(sdata[valid[temp[,2],1], j]) -
                      as.numeric(sdata[valid[temp[,1],1], j])
        keep <- rowSums(c2)==1 & rowSums(c2<0) ==0
        temp[keep,,drop=FALSE] # pairs whose start state differs by 1
    })
    pairs <- do.call("rbind", pairs)
    data.frame(t1 = paste(valid[pairs[,1], 1], valid[pairs[,1], 2], sep=':'),
               t2 = paste(valid[pairs[,2], 1], valid[pairs[,2], 2], sep=':'))
}
