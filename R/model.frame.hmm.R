model.frame.hmm <- function (formula, ...) {
    dots <- list(...)
    nargs <- dots[match(c("data", "na.action", "subset"), 
            names(dots), 0)]

    if (length(nargs) || is.null(formula$model)) {
        fcall <- formula$call
        indx <- match(c("formula", "data", "weights", "subset", 
            "na.action"), names(fcall), nomatch = 0)
        if (indx[1] == 0) 
            stop("The hmm call is missing a formula!")
        temp <- fcall[c(1, indx)]
        temp[[1L]] <- quote(stats::model.frame) 
        temp$formula <- formula$terms 
        temp$xlev <- formula$xlevels
        if (length(nargs) > 0) 
            temp[names(nargs)] <- nargs
        if (is.null(environment(formula$terms))) 
            eval(temp, parent.frame())
        else eval(temp, environment(formula$terms), parent.frame())
     }
    else formula$model
}

model.matrix.hmm <- function(object, data,  ...) {
    if (missing(data) && !is.null(object[["x"]]))
        return(object[["x"]])

    Terms <- delete.response(object$terms)
    if (missing(data)) 
        mf <- stats::model.frame(object, ...)
    else {
        if (is.null(attr(data, "terms")))
            mf <- stats::model.frame(Terms, data, xlev=object$xlevels)
        else mf <- data  #assume we were given a model frame     
    } 
        
    model.matrix(Terms, mf, contrasts.arg= object$contrasts)
}   
