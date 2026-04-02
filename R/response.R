# Marker distribution functions
# hmm(...  marker= list(log(pib) ~ gaussian, 
#                             dx ~ multinomial(pattern=pmat)),..
# The parsemarker2 routine will create a call of
#           gaussian(stateinfo, ...) where 'gaussian' is the function below that
#              sets and returns a gaussian response function + other info
# likewise  multinomial(stateinfo, nclass, pattern)
#
# hmm.dist is used by parsemarker to check for legal names
hmm.dist <- c("gaussian", "logistic", "beta", "discrete")

# stateinfo: name and level information that is used to create labels
#  for the linear predictors, and to match linear predictors to states.
#  A list with one element per marker, with elements of
#    name: used for creating a label
#    levels= a vector of labels, e.g. if of length 2 there would be 2 gaussian
#       densities
#    index: a vector of length nstate with values of 1,2, ... which density
#       goes to each state; 0= no density
# static:  used to indicate no linear predictors.
# mlevel: factor levels, if the marker was a factor, otherwise NULL
# param: optional params, eg."mean" or "std", not present= both
#
# return a list with 
#   rfun: the response function
#   pname: labels for the parameters
#   subset: which subset are referred to by the param argument
#   checkfun, optional
# The functions below are called by parsemarker2, and the presence or
#  absence of levels is sometimes enough for an error, .e.g. a factor can't
#  be a gaussian marker.  But the full set of y values that will be presented
#  is not yet available. The returned checkfun, if present, is called a
#  bit later by hmm when the response list is generated, and can be more
#  complete wrt a valid y.   For instance state:zed/gaussian where zed was
#  a list will be caught in the later check.
#  
gaussian <- function(stateinfo, levels, param, ...) {
    if (length(levels) >0) 
        stop("a factor variable cannot be a guassian marker")
    npeak <- length(stateinfo$levels)
    pname <- makedistlabels(stateinfo, c("mean", "std"))
    if (missing(param)) subset= 1:ncol(pname)
    else {
        index <- match(param, c("mean", "std"))
        if (any(is.na(index))) 
            stop("unrecognized gaussian parameter: ", param[is.na(index)])
        if (all(sort(index) == 1:2)) subset<- 1:ncol(pname)
        else if (index==1) subset <- 1:npeak * 2L - 1L
        else subset <- 2L* 1:npeak
    }
    
    rfun <- function(y, eta, gradient=FALSE,  npeak, map) {
        # y a vector of m values, m= number of measurements of this biomarker
        # eta = matrix of values, m rows by (2*npeak) columns of mean, log(std),
        #    mean, log(std), etc.
        # map= a map of peak to state
        # value: nstate rows by m columns, row j= f(y| state=j)
        #
        nstate <- length(map)
        mcount <- table(map[map!=0])
        yprob <- matrix(0., nstate, length(y))
        if (gradient) ygrad <- array(0, dim=c(nstate, length(y), ncol(eta)))
   
        for (i in 1:npeak) {
            j <- i*2L - 1L  # npeak (mean, log(std)) pairs
            std <- exp(eta[,j+1])
            yprob[map==i,] <- rep(dnorm(y, eta[,i], std, log=TRUE),
                                      each= mcount[i])
            if (gradient) { #(gradient is easier on log scale)
                gradient[map==i,,j] = (y- eta[,j])/std^2
                gradient[map==i,,j+1] <- (((y-eta[,j])/std)^2 -1)
            }
        }
    
        # convert from log to density
        yprob[map>0,] <- exp(yprob[map>0,])
        if (gradient) attr(yprob, "gradient") <- c(yprob) * ygrad
        yprob
    }

    # rfun will be called many times; set the defaults so that they don't
    #  have to be passed through the maximizer
    temp <- formals(rfun)
    temp$npeak <- npeak
    temp$map <- stateinfo$index
    formals(rfun) <- temp

    check <- function(y)
        if (!is.numeric(y)) stop("a gaussian marker must be numeric")
    list(name="gaussian", rfun=rfun, pname=pname, subset=subset,
         checkfun = check)
}

makedistlabels <- function(stateinfo, parms) {
    nlev <- length(stateinfo$levels)
    nparm <-length(parms)
    # I don't need to say "state(dementia)", "dementia" will do
    if (stateinfo$sname == "state") 
        rbind(rep(stateinfo$levels, each=nparm), rep(parms, nlev))
    else {
        temp <- paste0(stateinfo$sname, '(', stateinfo$levels, ')')
        rbind(rep(temp, each=nparm), rep(parms, nlev))
    }
}

# logistic, a bit fatter tails
logistic <- function(stateinfo, levels, param, ...) {
    npeak <- length(stateinfo$levels)
    pname <- makedistlabels(stateinfo, c("mean", "std"))
    if (missing(param)) subset= 1:ncol(pname)
    else {
        index <- match(param, c("mean", "std"))
        if (any(is.na(index))) 
            stop("unrecognized logistic parameter: ", param[is.na(index)])
        if (all(sort(index) == 1:2)) subset<- 1:ncol(pname)
        else if (index==1) subset <- 1:npeak * 2L - 1L
        else subset <- 2L* 1:npeak
    }
    
    rfun <- function(y, eta, gradient=FALSE,  npeak, map) {
        # y a vector of m values, m= number of measurements of this biomarker
        # eta = matrix of values, m rows by (2*npeak) columns of mean, log(std),
        #    mean, log(std), etc.
        # map= a map of peak to state
        # value: nstate rows by m columns, row j= f(y| state=j)
        #
        nstate <- length(map)
        mcount <- table(map[map!=0])
        yprob <- matrix(0., nstate, length(y))
        if (gradient) ygrad <- array(0, dim=c(nstate, length(y), ncol(eta)))
   
        for (i in 1:npeak) {
            j <- i*2L - 1L  # npeak (mean, log(std)) pairs
            scale <- exp(eta[,j])* sqrt(3)/pi # R logis has "scale" not "std"
            yprob[map==i,] <- rep(dlogis(y, eta[,j], scale, log=TRUE),
                                  each= mcount[i])
            if (gradient) {
                temp1 <- (y- eta[,j])/scale
                temp2 <- exp(temp1)/(1+ exp(temp1))
                gradient[map==i,,j] = -temp1*(1 + temp2)
                gradient[map==i,,j+1] = -(1+ temp1*(1+temp2))
            }
            # convert from log to density
            yprob[map>0,] <- exp(yprob[map>0,])
            if (gradient) attr(yprob, "gradient") <- c(yprob) * ygrad
            yprob
        }
    }
    # rfun will be called many times; set the defaults so that they don't
    #  have to be passed through the maximizer
    temp <- formals(rfun)
    temp$npeak <- npeak
    temp$map <- stateinfo$index
    formals(rfun) <- temp

    check <- function(y)
        if (!is.numeric(y)) stop("a gaussian marker must be numeric")
    list(name="logistic", rfun=rfun, pname=pname, subset=subset, checkfun=check)
}

# beta distribution
beta <- function(stateinfo, markerlevel, param) {
    npeak <- length(stateinfo$levels)
    pname <- makedistlabels(stateinfo, c("shape1", "shape2"))
    if (missing(param)) subset= 1:ncol(pname)
    else {
        index <- match(param, c("shape1", "shape2"))
        if (any(is.na(index))) 
            stop("unrecognized beta parameter: ", param[is.na(index)])
        if (all(sort(index) == 1:2)) subset<- 1:ncol(pname)
        else if (index==1) subset <- 1:npeak * 2L - 1L
        else subset <- 2L* 1:npeak
    }
    
    rfun <- function(y, eta, gradient=FALSE,  npeak, map) {
        # y a vector of m values, m= number of measurements of this biomarker
        # eta = matrix of values, m rows by (2*npeak) columns of log(shape1),
        #    log(shape2) for first peak, then second, ...
        # map= a map of peak to state
        # value: nstate rows by m columns, row j= f(y| state=j) 
        #
        nstate <- length(map)
        mcount <- table(map[map!=0])
        yprob <- matrix(0., nstate, length(y))
        if (gradient) ygrad <- array(0, dim=c(nstate, length(y), ncol(eta)))
   
        for (i in 1:npeak) {
            j <- i*2L - 1L  # npeak (log(shape1), log(shape2)) pairs
            a <- exp(eta[,j])
            b <- exp(eta[,j+1])
            f <- dbeta(y, a, b, log=FALSE)
            yprob[map==i,] <- rep(f, each= mcount[i])
            if (gradient) {
                #see the derivation in the code vignette
                dga <- psi(a + b) - psi(a)
                dgb <- psi(a + b) - psi(b)
                g <- gamma(a+b)/(gamma(a)* gamma(b))
                dha <- (a-1)*y^(a-2)* (1-y)^(b-1)
                dhb <- -(y^(a-1) * (b-1)*(1-y)^(b-2))
                ygrad[map==i,,j] = a*(dga* f + g*dha)
                ygrad[map==i,,j+1] = b*(dgb*f + g*dhb)
                attr(yprob, "gradient") <- ygrad
            }

            yprob
        }
    }
    # rfun will be called many times; set the defaults so that they don't
    #  have to be passed through the maximizer
    temp <- formals(rfun)
    temp$npeak <- npeak
    temp$map <- stateinfo$index
    formals(rfun) <- temp

    check <- function(y) {
        if (any(y<0 | y>1)) stop("invalid marker value for beta distribution")
    }
    list(name="beta", rfun=rfun, pname=pname, subset=subset, checkfun= check)
}

# multivariate logit, first category is the reference. If there are k
#  groups, eta will have k-1 columns
# This is used by other functions
mlogit <- function(eta, gradient=FALSE) {
    m <- ncol(eta)
    denom <- 1 + rowSums(exp(eta))
    pi <- cbind(1, exp(eta)) /denom

    if (gradient) {
        if (nrow(eta) ==1) {  # only one subject
            pi2 <- drop(pi)
            #temp <- diag(pi) - outer(pi, pi) # next line is a touch faster
            temp <- diag(pi2) -  rep(pi2, m+1) * rep(pi2, each= m+1)
            dmat <- temp[,-1,drop=FALSE] 
            dim(dmat) <- c(1, dim(dmat))  # make it a "row" per subject
        }
        else {
            dmat <- array(0., dim=c(nrow(eta), m+1, m))
            dmat[,1,] <- -pi[,-1]/denom
            for (j in 1:m) { # for each column of eta
                dmat[,j+1,] <- -pi[,j+1]* pi[,-1]
                dmat[,j+1, j] <- pi[,j+1]*(1-pi[,j+1])
            }
        }
        attr(pi, "gradient") <- dmat
    }
    pi
}
# 
#  The rfun for discrete has the usual y and eta arguments, along with a
# matrix etamap. The response y is a factor, each level maps to a column
# of etamap, rows of etamap correspond to states, and elements of etamap index
# to the linear predictors.  The details of this corresondence are controled
# by a pattern argument.  If there are no linear predictors and instead
# an init matrix, that is used directly: a row per state, col per response.
#  This function sets this all up and creates rfun
#
discrete <- function(stateinfo, mlevel, init, pattern, static) {
    ngroup <- length(stateinfo$levels)  #number of calls to hmmlogit, later
    nstate <- length(stateinfo$index)
    if (length(mlevel) ==0)
        stop("a discrete marker must be a factor")
    nlevel <- length(mlevel)
    
    if (!missing(pattern)) {
        # check for a valid pattern matrix
        if (!is.matrix(pattern))
            stop("pattern argument must be a matrix")
        if (ncol(pattern) != nlevel)
            stop("pattern matrix's columns must match marker variable")
        if (!is.null(colnames(pattern))) {
            k <- match(colnames(pattern), mlevel)
            if (any(is.na(k)) || any(duplicated(k)))
                stop("pattern matrix's columns must match marker variable")
            pattern <- pattern[,k]
        }  # cols are now in the order of the marker variable (factor)
        if (nrow(pattern) != ngroup)
            stop("pattern should have one row per ", stateinfo$sname, " level")
        if (!is.null(rownames(pattern))) {
            k <- match(rownames(pattern), stateinfo$levels)
            if (any(is.na(k)) || any(duplicated(k)))
                stop("rows of pattern should match levels of ", stateinfo$sname)
        }
        if (any(pattern) != floor(pattern) | any(pattern < 0))
            stop("a pattern matrix for discrete must be integers >=0")
        n.prob <- apply(pattern, 1, function(x) length(unique(x[x>0])))
        # each row has to sum to 1, so there fewer linear predictors than probs

        # Create the emap matrix. There will be one vector of probabilities,
        #  which sum to 1, for postive value of statefig$index.
        # emap: ngroup rows and nlevel columns (same as pattern)
        #  if emap[i,j]= k, then eta[,k] is used for index 1 and mlevel j. 
        # Multiple emap elements can be associated with the same linear 
        #  predictor k. 
        emap <- matrix(0L, ngroup, nlevel)
        # any coefficients shared across rows? (other than 0)
        # first, make the values be 0,1,2,.. with no gaps
        # The user will often assign numbers by row, use unique rather than
        #  sort(unique( to keep things in the same order
        emap[,] <- match(pattern, unique(c(0L, pattern))) -1L
        maxp <- max(emap) # total number of unique phat values
        nz <- (emap >0)
        pcount <- table(row(emap)[nz], emap[nz])
        across.row <- (colSums(pcount>0) > 1)
        if (any(across.row)) emap[,across.row] <- emap[,across.row] + maxp
        # by default, choose the smallest non-shared index in each row as the
        #  reference value, but avoid those that are shared across rows
        for (i in 1:nrow(emap)) {
            j <- min(emap[i, nz[i,]])
            reference <- (emap[i,] ==j)  # reference cell(s) for this row
            emap[i,reference] <- -1    # special code for reference cells
        }
        uval <- unique(emap[emap>0]) # not counting reference cells
        emap[emap>0] <- match(emap[emap>0], uval) # re-number
        dimnames(emap) <- list(paste0(stateinfo$sname, stateinfo$levels),
                               marker=mlevel)
        # why match rather than emap>0?  I don't want duplicates
        pname <- outer(rownames(emap), colnames(emap),sep='_')[match(uval, emap)]
    } else if (!static) {
        # Assume the fully parameterized missclassification matrix
        uval <- seq.int(1, ngroup *(nlevel-1))
        temp <- matrix(uval, nrow=ngroup, byrow=TRUE)
        emap <- cbind(-1, temp) # first state is reference group
        dimnames(emap) <- list(paste0(stateinfo$sname, stateinfo$levels),
                               marker=mlevel)
        pname <- outer(rownames(emap), colnames(emap),sep='_')[match(uval, emap)]
    }  

    if (!missing(init)) {
        # check for a valid init matrix
        if (!is.matrix(init))
            stop("init must be a missclassification matrix")
        if (ncol(init) != nlevel)
            stop("init matrix's columns must match marker variable")
        if (!is.null(colnames(init))) {
            k <- match(colnames(init), mlevel)
            if (any(is.na(k)) || any(duplicated(k)))
                stop("init matrix's columns must match marker variable")
            init <- init[,k]
        }  # cols are now in the order of the marker variable (factor)
        if (!is.null(rownames(init))) {
            k <- match(rownames(init), stateinfo$levels)
            if (any(is.na(k)) || any(duplicated(k)))
                stop("rows of init should match levels of", stateinfo$sname)
        }
        if (!all(rowSums(init) ==1))
            stop("row sums of init matrix must be 1")
        if (nrow(init) != nstate) {
            # assume missing rows are states for which the marker is irrelevant
            i2 <- matrix(0, nstate, ncol(init))
            i2[k,] <- init
            init <- i2
        }
    }

    if (static) { # no parameters
        if (missing(init)) stop("init argument is needed discrete")
        # expand rows to one per state
        missmat <- matrix(0, nrow=nstate, ncol= nlevel)
        j <- stateinfo$index
        missmat[j>0,] <- init[j,]
        colnames(missmat) <- mlevel
        rfun <- function(y, missclass= missmat, deriv=FALSE) missclass[,y]
        return(list(name="discrete", rfun=rfun, pname=NULL))
    }
        
    # Remainder is the more common non-static case
    if (!missing(init)) 
        warning("init option not yet available for pattern matrix")
    # There is not necessarily a set of eta values that will exactly
    #  produce a user's desired initial probabilities, so implementation
    #  will require a non-linear maximization.  The result would
    #  be passed back and used as the initial intercept parameters.
          
    rfun <- function(y, eta, gradient=FALSE, emap=emap, index=stateinfo$index) {
        # y = a factor, , m= number of measurements of this biomarker
        # eta = matrix of values: m rows, one column per linear predictor
        # emap = map from lp to errors, -1= reference cell
        ny <- length(y)
        if (is.factor(y)) y <- as.integer(y)  # might already have been converted
        if (!is.matrix(eta)) eta <- matrix(eta, ncol=1) #single linear predictor
        nstate <- length(index)
        phat <- matrix(0., nrow= nstate, ncol=ny)
        if (any(index==0)) phat[index==0,] <- 0  #marker uninformative for state
        if (gradient) gmat <- array(0., dim=c(nstate, ny, ncol(eta)))

        # Do one row of  emap at a time
        # 
        for (i in 1:nrow(emap)) {
            etemp <- emat[i,]
            ecol <- sort(unque(etemp[etemp>0]))
            mtemp <- mlogit(eta[,ecol]) # result has col 1= ref, then others
            indx <- c(which(etemp== -1), 1L+ match(etemp, ecol, nomatch= -1))
            yindx <- indx[y]
            if (gradient) {
                for (j in which(index==i)) {
                    phat[j, yindx>0] <- mtemp[yindx]
                    gmat[j, yindx>0,ecol]<- attr(mtemp, "gradient")[yindx,]
                }
            }
        }
        if (gradient) attr(phat, "gradient") <- gmat
        phat
    }

    check <- function(y) {
        if (!is.factor(y))
            stop("marker for a discrete distribution must be categorical")
    }
            
    list(name="discrete", rfun=rfun, pname=pname, checkfun= check)
}
        
