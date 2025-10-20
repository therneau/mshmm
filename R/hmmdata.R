# Automatically generated from hmmcode
hmmdata <- function(formula, data, na.action=na.pass, id,
                    istate="istate", tstart="tstart") {
    # Start with standard model.frame stuff.  The right hand side of the
    #  formula will nearly always be "."    
    Call <- match.call()
    indx <- match(c("formula", "data", "id"), 
                  names(Call), nomatch=0)
    if (indx[1] ==0) stop("a formula argument is required")
    temp <- Call[c(1L, indx)]
    temp$na.action <- na.pass
    temp[[1L]] <- quote(stats::model.frame)
    mf <- eval.parent(temp)
    
    Y <- model.response(mf)
    if (!is.Surv(Y)) stop("left hand side of the formula must be a Surv call")
    if (ncol(Y) ==3) stop("data set already has a start-stop format")
    id <- model.extract(mf, "id")
    if (is.null(id)) stop("id variable not found")
    uid <- unique(id)  #unique id values, in the order in which they appear
    index <- order(match(id, uid), Y[,2])
    mf <- mf[index,]
    Y  <- Y[index,]
    id <- id[index]
    first <- match(uid, id)
    if (attr(Y, "type") == "mright") {
        states <- attr(Y, "states")
        status <- factor(Y[,2], 0:length(states), labels=c("censor", states))
    }
    else status <- Y[,2]

    last <- c(first[-1] -1, nrow(mf))
    newY <- data.frame(t1= Y[-last,1], t2 =Y[-first,1], 
                       status = status[-first], 
                       istate= status[match(id, id)][-first],
                       stringsAsFactors=FALSE)
    temp <- formula[[2]]  # find the old names
    if (class(temp[[2]])=="call" && temp[[2]][[1]] == as.name("Surv"))
        yname <- c(tstart, as.character(temp[[2]][[2]]), istate, 
                   as.character(temp[[2]][[3]]))
    else yname <- c(tstart, dimnames(Y)[[2]], istate)
    names(newY) <- yname

    toss <- c(1, which(names(mf) == "(id)"))
    newdata <- cbind(newY, mf[-last, -toss])  #new Y variables, then all the others
    row.names(newdata) <- NULL
    newdata
    }
