# The routine allows a list of formulas. The first is the default, such as
#    Surv(time, death) ~ 1+ x
# Later ones are state:state ~ covariates
# See the section on parsing formulas in the code vignette for an more
#  complete discussion of how this is all done. It shows the graph of a
#  parse tree.
#
# We use formulas, but with some changes in how we interpret things.
# An advantage is that we gain all the knowlege of the R parser,
#  it understands nested parenthesis for instance. A  disadvantage is that our 
#  formula has to look 'legal' to the parser.
#     "Surv(time,death) ~ x1 / init=5" for instance won't fly since
#  lm/glm formulas don't have equals signs. 
#     "Surv(time, death) ~ x1 + x2 / init(c(3, 2.1)" is okay, it looks like
#  a function call
#     "1:3 + 2:3 ~ x" is okay and "c(1,2):3 ~ x" but not "(1,2):3 ~ x".
#
# The first pass splits out the left side (states), the formula for the 
#  variables, and the options.  The second of these is used to get the model 
#  frame.  The rest of the processing is deferred until after the model frame
#  has been built.
#
parsecovar1 <- function(flist) {
    # flist = all the formulas except the default
    if (any(sapply(flist, function(x) !inherits(x, "formula"))))
        stop("an element of the formula list is not a formula")
    if (any(sapply(flist, length) != 3))
        stop("all formulas must have a left and right side")
    
    # split the formulas into a right hand and left hand side
    lhs <- lapply(flist, function(x) x[[2]])  
    rhs <- lapply(flist, function(x) x[[3]])  
    
    temp <- lapply(rhs, rightslash)
    options <- lapply(temp, function(x) if (is.list(x)) x[[2]] else NULL)
    rightformula <- lapply(temp, function(x) {
        tfun <- ~ z  #dummy function
        if (is.list(x)) tfun[[2]] <- x[[1]] else tfun[[2]] <- x
        tfun
        }) # make sure each element of rightformula is a one sided formula

    list(lhs = lhs, rhs= rightformula, options=options)
}

# The following function splits a formula at the rightmost slash, ignoring
# the inside of any function or parenthesised phrase.
# Recursive functions like this are almost impossible to read, but luckily 
# it is short.
# The function recurrs on the left and right side of +,*,:, and \%in\%, and on 
#  binary - (but not on unary -).
# If there are options the result will be a 2 element list: the formula without
#  options, and the options (unevaluated), if no options x is returned
rightslash <- function(x) {
    if (!inherits(x, 'call')) return(x)
    else {
        if (x[[1]] == as.name('/')) return(list(x[[2]], x[[3]]))
        else if (x[[1]]==as.name('+') || 
                 (x[[1]]==as.name('-') && length(x)==3)  ||
                 x[[1]]==as.name('*') || x[[1]]==as.name(':')  ||
                 x[[1]]==as.name('%in%')) {
                     temp <- rightslash(x[[3]])
                     if (is.list(temp)) {
                         x[[3]] <- temp[[1]]
                         return(list(x, temp[[2]]))
                     } else {
                         temp <- rightslash(x[[2]])
                         if (is.list(temp)) {
                             x[[2]] <- temp[[2]]
                             return(list(temp[[1]], x))
                         } else return(x)
                     }
                 }
        else return(x)
    }
}
   
# expand the set of transitions implied by a left hand side.
statepair <- function(lhs, statemap) {
    # create a dummy function for each column of statemap
    # state() is a created function such that state("s1", "s2") will first
    #  check that s1 and s2 are present in the statemap[,"state"], give an
    #  error if not, and return which rows of statemap[,"state"] have
    #  a value of "s1" or "s2"
    # If statemap has another column "N" say, we want to also create an N()
    #  function of the same type. Do that by making a copy of state(),
    #  then changing the default values of target and cname in situ by
    #  using formals.
    # All of these are put into a separate environment because these short new
    #  function names might well conflict with a prior variable name.
    # The entire reason for this is so that a user can type N(1) as a 
    #  shorthand for "all the states for which variable N is 1"
    env1 <- new.env(parent= parent.frame(2))
    assign("state", env= env1,
           value = function(..., target=statemap[,1], cname= "state") {
               j <- c(...)
               check <- match(j, target)
               if (any(is.na(check))) stop("value ", j[is.na(check)], 
                                           " not found in ", cname)
               which(target %in% j)
           })

    if (ncol(statemap) > 1) {
        cname <- colnames(statemap)
        for (i in 2:ncol(statemap)) {
            temp <- get("state", env= env1)
            ftemp <- formals(temp)
            ftemp$target <- statemap[,i]
            ftemp$cname <- cname[i]
            formals(temp) <- ftemp
            assign(cname[i], temp, env= env1)
        }
    }

    # seup done, the real work is done by a recursive function
    statewalk <- function(x, nstate) {
        # simple ones first
        if (is.character(x) || (length(x)==1 & is.name(x))) {
            z <- match(as.character(x), statemap$state)
            if (any(is.na(z))) stop("invalid state: ", 
                                    (as.character(x))[is.na(z)])
            else return(z)
        } else if (is.numeric(x)) {
            if (x==0) return(1:nstate)
            else {
                if (x != as.integer(x)) stop("non integer state: ", x)
                else if (x<1 || x > nstate) stop("state out of range: ", x)
                else return(x)
            }
        } else if (length(x)==1) stop("unrecognized symbol in statepair: ", x)
 
        if (x[[1]]== as.name("*") || x[[1]]== as.name("/") || 
            x[[1]]== as.name("-") || x[[1]]== as.name("%in%"))
            stop("invalid operator on the left of a formula")

        if (x[[1]]== as.name(":")) { # the heart of the function
            # the left and right should 2 row matrices
            from <- statewalk(x[[2]], nstate)
            to   <- statewalk(x[[3]], nstate)
            pairs <- rbind(from= rep(from, length(to)), 
                           to= rep(to, each=length(from)))
            pairs[,pairs[1,] != pairs[2,]]
        } else if (x[[1]] == as.name("+")) {
            # both left and right should be from:to sets, or neither
            if (is.matrix(x[[2]]) && nrow(x[[2]])==2) {
                if (is.matrix(x[[3]]) && nrow(x[[3]])==2) cbind(x[[2]], x[[3]])
                else stop("statewalk error 1")
            } else {
                if (is.matrix(x[[3]]) && nrow(x[[3]]) ==2) 
                    stop("statewalk error 2")
                else cbind(statewalk(x[[2]], nstate), statewalk(x[[3]], nstate))
            }
        } else if (x[[1]] == as.name("(")) statewalk(x[[2]], nstate)
        else if (x[[1]] == as.name("c")) statewalk(eval(x), nstate)
        else { #match to a row in statemap
            # a user might write 3:"death" or 3:death
            if (is.character(x) || is.name(x)) 
                z <- which(statemap$state == as.character(x)) 
            else z <- eval(x, env= env1)

            if (!is.numeric(z)) stop("non-numeric state: ", deparse(x))
            if (any(z != as.integer(z))) {
                j <- min(which(z!= as.integer(z)))
                stop("non integer state: ", z[j])
            }
            if (any(z<1 | z > nstate)) {
                j <- min(which(z<1 | z>nstate))
                stop("state out of range: ", z[j])
            }
            z
        }
     }
 
    statewalk(lhs, nrow(statemap))
}

        
# A key trick for the rhs (variables) of a formula is to call
#  terms() on it rather than parse it myself.  One glitch is that an
#  interaction might appear as "x1:x2" in term.labels attribute of the master
#  formula (Terms argument) but as "x2:x1" in the terms of a the sub-formula
#  for a particular transition; a simple match() call on term.labels won't
#  work. The function below works around this.  The factors attribute for 
#  each term is a matrix with a column for x1:x2 in one and x2:x1 in the other
#  and row names that include x1 and x2. 
#  Sort rows into the same order and find matching columns.

termmatch <- function(f1, f2) {
    # f1 = attr(terms, 'factors') of smaller formula, f2 =  master formula
    if (length(f1)==0) return(NULL)   # a formula with only ~1
    irow <- match(rownames(f1), rownames(f2))
    if (any(is.na(irow))) stop ("termmatch failure 1") # should never happen
    hashfun <- function(j) sum(ifelse(j==0, 0, 2^(seq(along.with=j))))
    # make a variant of f1 that has the same number of rows as f2. If we instead
    # subset f2, then x1:x2 could match the hashfun of x1:x2:x3
    dummy <- matrix(0, nrow(f2), ncol(f1))
    dummy[irow,] <- f1

    hash1 <- apply(dummy, 2, hashfun)
    hash2 <- apply(f2, 2, hashfun)
    index <- match(hash1, hash2)
    if (any(is.na(index))) stop("termmatch failure 2")
    index
}

# do the options contain "common"?
hascommon <- function(options) {
    if (is.null(options)) return(FALSE)
    tform <- ~ x # dummy formula
    tform[[2]] <- options
    "common" %in% attr(terms(tform), "term.labels")
}
        
# The second part of parsing the formula
# The goal is to create a matrix with a row for each term in the model,
#  e.g., sex or ns(age, 3), and a column for each transition. Elements of
#  the matrix will be integers identifing unique terms, with 0= not present
# A formula with a /common option will lead to repeated values. 
# The code for coxph is more complex because there the "intercept" is the 
#  baseline hazard function; keeping track of shared baselines is more subtle
#
#  parse2 function
#   parse1 = return from parsecovar1
#   statedata= alias data frame,
#   dformula = default formula 
#   Terms = the terms structure from the master formula that was used to 
#    build the data frame.
#   qmatrix = the Q matrix of valid transitions
#   qmatrix = vector of valid states
# The version of this code in survival is more complex, because it also
#  needs to sort out strata, which are not a concept here
parsecovar2 <- function(parse1, statedata, dformula, Terms, qmatrix,
                        Xname, Xassign) {
    nterm <- 1L + length(attr(Terms, "term.labels")) # +1 for (Intercept)
    states <- statedata[,1]  # will always contain the state names
    nstate <- length(states)
    from <- row(qmatrix)[qmatrix>0]
    to   <- col(qmatrix)[qmatrix>0]
    tran.id <- paste(from, to, sep=':') # col labels for tmap and cmap
    ntran <- length(from)
    if (Xassign[1] != 0) {
        # I am not sure that this can every happen, nevetheless
        # make sure that our versions of tmap and cmap have an intercept row
        # If it turns out to be all zeros, we can remove it at the end
        # 
        Xname <- c(Xname, "(Intercept)")
        Xassign <- c(0, Xassign)
    }
    ncoef <- length(Xname)  # Xattr is the same length
    Tfac <- attr(Terms, "factors") # save some typing later

    # Create tmap: a row for each term and a column for each transition.
    #  value of 0 = this term isn't used for this transition
    #  1, 2, etc = marks unique sets of coefficients
    #  dmap = a matrix of unique integers to draw from, so that we don't reuse 
    #  an index
    # cmap: row for each coefficient, col for each transition.  A term like
    #  ns(age, df=4) will map to 3  coefficients.
    # Xattr points to the term (0 = intercept) each column of X descends from
    # 
    tmap <- matrix(0L, nterm, ntran)
    dmap <- matrix(seq_len(ncoef*ntran), ncol= ntran) # term numbers
    cmap <- matrix(0L, ncoef, ntran)    
    init <- matrix(NA, ncoef, ntran) # fill in user supplied initial values
    # init is not yet implemented

    # these two functions deal with the fact that tmap has a row per term
    # and cmap and dmap a row per covariate, in both cases i is a term number
    tindx <- function(i) match(i, 1+Xassign)
    cindx <- function(i) which(Xassign %in% (i-1))

    # initialize every column with the default formula, which cannot have a
    #  /common option
    temp <- delete.response(terms(dformula))
    dterm <- 1L + termmatch(attr(temp, "factors"), Tfac)
    if (attr(temp, "intercept") ==1) dterm <- c(1L, dterm)
    for (i in 1:ntran) {
        tmap[dterm,i] <- dmap[tindx(dterm),i]
        cmap[cindx(dterm),i] <- dmap[cindx(dterm), i]
    }
    # Create a list of working formulas, one per transition. Each starts as the
    #  default formula. 
    formlist <- lapply(1:ntran, function(i) dformula[-2])  

    # Discover whether each formula contains an explicit +1.  The way that I
    #  have found is to paste an explicit "-1 +" to the front, and see if
    #  if the resulting intercept attribute is 0 or 1
    has1 <- sapply(parse1$rhs, function(x) {
        tform <- ~ -1 +zed     # dummy formula
        tform[[2]][[3]] <- x   # replace 'zed' with the element of parse1$rhs
        attr(terms.formula(tform), "intercept")
    })

    # The transitions targeted by each formula in the list of model formulas,
    #  These are the formulas in the model statement,
    #  contained in parse1$lhs.  (If there is only a default formula, lhs is
    #  NULL and translist will be a list with 0 elements.)
    # Each element of translist gives the columns of tmap affected by that
    #  formula line, ie. the set of transitions.The result will be a list,
    translist <- lapply(parse1$lhs, function(x) {
        temp <- statepair(x, statedata)
        if (is.matrix(temp)) id <- paste(temp[1,], temp[2,], sep=':')
        else id <- paste(temp[1], temp[2], sep=':')  #single transition
        indx <- match(id, tran.id)
        # A specification like A(0):A(1) will generate spurious illegal 
        #  combinations, ignore them silently.
        # This means we won't catch all user errors. But if all of indx are
        #  missing, it is almost certainly a direct mistake of an impossible
        #  A:B transtion
        if (all(is.na(indx))) 
            stop("invalid transtion(s): ", paste(id, collapse=', '))
        indx[!is.na(indx)]
    })

    # Process each formula in turn, there might be none (rare)
    for (k in seq(along.with= parse1$lhs)) {
        kform <- (parse1$rhs[[k]]) # the formula for this update,
        kcommon <- hascommon(parse1$options[[k]])
        j <- translist[[k]] # cols of tmat that this kform applies to

        add <- 1L + termmatch(attr(terms(kform), "factors"), Tfac)
        if (has1[k]==1) add <- c(1L, add)
        if (length(add) >0) { # false for a formula with only "-" terms
            if (kcommon) {
                tmap[add,j] <- dmap[tindx(add), j[1]]
                cmap[cindx(add),j] <- dmap[cindx(add), j[1]]
            } else {
                tmap[add,j] <- dmap[tindx(add), j]
                cmap[cindx(add),j] <- dmap[cindx(add), j]
            }
        }
        
        # Update the running formula for each transition using update.formula
        # To do this temporarily paste "~ ." on the front if the formula starts
        #  with a leading minus or minus sign, or "~. +" otherwise.
        # Manipulating strings seems to be the only way to do this.
        # uform = arg 2 of update.formula, that I want to apply.
        if (substring(deparse(kform[[2]]), 1,1) %in% c("-", "+")) 
            utemp <- "~." else utemp <- "~.+"
        uform <-parse(text= paste(utemp, deparse(kform[[2]]), collapse=""))[[1]]
        for (jj in j) { # for each affected term separately
            # update this transition's formula
            had.intercept <- attr(terms(formlist[[jj]]), "intercept")
            formlist[[jj]] <- update.formula(formlist[[jj]], uform)
            fterm <- terms.formula(formlist[[jj]])
            keep <- 1L + termmatch(attr(fterm, "factors"), Tfac)
            if (attr(fterm, "intercept") == 0 && had.intercept) 
                warning("-1 term in a formula was ignored")
            keep <- c(1L, keep)

            # anything not in the formula should be 0 in tmap/cmap
            toss <- (1:nterm)[-keep]
            if (length(toss) >0 && any(tmap[toss,jj] >0)) {
                tmap[toss, jj] <- 0L
                cmap[cindx(toss),jj] <- 0L
            }
        }
    }
    
    # reset the values in tmap to 0,1,2,3,...,  ditto for cmap
    tmap[,] <- match(c(tmap), sort(unique(c(0, tmap)))) -1L
    cmap[,] <- match(c(cmap), sort(unique(c(0, cmap)))) -1L
    dimnames(tmap) <- list(c("(Intercept)", attr(Terms, "term.labels")),
                           tran.id)
    dimnames(cmap) <- list(Xname, tran.id)
    mapid <- rbind(from, to)
    colnames(mapid) <- tran.id
    list(tmap= tmap, cmap= cmap, mapid= mapid)
}
 
