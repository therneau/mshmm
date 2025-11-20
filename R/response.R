# Marker distribution functions
# hmm(...  marker= list(log(pib) ~ gaussian, 
#                             dx ~ multinomial(pattern=pmat)),..
# The parsemarker2 routine will create a call of
#           gaussian(stateinfo, ...) where 'gaussian' is the function below that
#              sets and returns a gaussian response function + other info
# likewise  multinomial(statepattern, nclass, pattern)
#
# hmm.dist is used by parsemarker to check for legal names
hmm.dist <- c("gaussian", "logistic", "beta", "multinomial", "noerror")

# The stateinfo has name and level information that is use to create labels
#  for the linear predictors, plus an index matching gaussian densities to
#  states, 
#  the markerlevel argument is used by categorical methods, pattern by 
#   the multinomial dist, and the param option for return information.
#
# return a list with 
#   rfun: the response function
#   pname: labels for the parameters
#   subset: which subset are referred to by the param argument
gaussian <- function(stateinfo, markerlevel, param) {
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

    list(name="gaussian", rfun=rfun, pname=pname, subset=subset)
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
logistic <- function(stateinfo, markerlevel, param) {
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

    list(name="logistic", rfun=rfun, pname=pname, subset=subset)
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

    checkfun <- function(y) {
        if (any(y<0 | y>1)) stop("invalid marker value for beta distribution")
    }
    list(name="beta", rfun=rfun, pname=pname, subset=subset, check= check)
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
# mulinomial distribution
#  Say the the marker has 5 levels and there were k states, for each
#  state there will 4 parameters to create the 5 probabilities of p=
#  (1, exp(eta2), exp(eta3), exp(eta4), exp(eta5)) / (1 + exp(eta2) + ...eta5))
# wlog eta1 is taken to be 0, leading to 4k linear predictors eta.
# The derivatives turn out to have the same form as a multinomial variance:
#  d p_i/d eta_i = p_i(1-p_i) and d p_i/ d eta_j = -p_ip_j; which is easy to
#  remember. See the code vignette for more explanation. 
# If the marker had 5 categories, the evaluation routine needs to report back
#  pr(observed marker value, given state), i.e., a single probability p for 
#  each state, for that observation (nstate by n.obs matrix), along with the 4
#  derivatives of that value wrt the 4 eta values (nstate by n.obs by 4).
#  The other p_i are just a tool to compute that set of 4 derivative values: 
#    there is no variance matrix.
# 
# For more complex cases the user can supply a pattern matrix with one row
#  per state and one col per value of the marker. A zero value in the matrix
#  indicates that that state/marker pairing will not occur; there is no need
#  to waste linear predictors for that eventuality. If there were m=5 marker
#  values but one of them can not occur when the true state is 'A', there will
#  be 3 eta vectors, not 4, for A.
# Other values in the matrix determine the mapping, i.e., the smallest 
#  non-zero value in a row identifies the reference category, the order of the
#  remaining matrix values determine the mapping of eta to state/prob. 
##  
multinomial <- function(stateinfo, nlevel, pattern) {
    nstate <- length(stateinfo$index)
    if (nlevel ==0) stop("marker must be a factor for multinomial distribution")
    ngroup <- max(stateinfo$index)  #number of predicted phat vectors

    if (missing(pattern)) {
        # the compute function gets two lists with one elment per state
        #  eindex = which columns of eta for this state
        #    our default is to use 1,2,..., k-1 for state 1, k, ... for 
        #    state 2, etc where k is the number of levels for the marker
        #  mindex = non-zero responses
        n.eta <- nstate * (nlevel-1)
        eindex <- split(1:n.eta, rep(1:nstate, each=nlevel-1))
        mindex <- lapply(1:nstate, function(x) 1:nlevel)
        nparm <- nstate
    }
    else {
        if (!is.matrix(pattern) || nrow(pattern) != ngroup ||
            ncol(pattern) != nlevel)
        stop("pattern must be a matrix with one row per group of states",
             " and one column per level of the marker")
        if (any(is.na(pattern))) stop("missing value in pattern matrix")

        
        nphat <- apply(pattern!=0, 1, sum)  # number of probabilities per row
        if (any(nphat ==0)) stop("pattern matrix has a zero row")
        n.eta  <- sum(nphat -1) # total number of linear predictors
    }
    p2 <- pattern # modify this into "standard" form
    ref <- apply(pattern, 1, function(x) min(which(x!=0)))
    p2[cbind(1:ngroup, ref)] <- -1
    index <- which(p2>0)
    p2[index] <- rank(p2[index], ties="first")
    # make pname
    
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

    checkfun <- function(y, nc=ncol(pattern)) {
        if (!is.factor(y) || length(levels(y))!= nc)
            stop("marker must be a factor with ", nc, 
                 " levels, to match the pattern matrix")
    }
            
    list(name="multinomial", rfun=rfun, pname=pname, subset=subset, 
         check= checkfun)
}


# the distribution for an state observed without error
noerror <- function(stateinfo, ...) {
    rfun <- function(y, nstate) {
        # eta should be 0 columns, gradient will be ignored
        temp <- diag(nstate)
        if (any(y < 1 | y>nstate | floor(y) !=y))
            stop("y must be an integer between 1 and number of states")
        temp[,y, drop = FALSE]
        }
    temp <- formals(rfun)
    temp$nstate <- length(stateinfo$index)
    formals(rfun) <- temp

    checkfun <- function(y, nstate) {
        if (any(y < 1 | y>nstate | floor(y) !=y))
            stop("y must be an integer between 1 and number of states")
    }
    temp <- formals(checkfun)
    temp$nstate <- length(stateinfo$index)
    formals(checkfun) <- temp
    list(name="noerror", rfun=rfun, pname=NULL, check= checkfun)
}
        
