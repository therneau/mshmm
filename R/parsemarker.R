# Parse the list of markers for latent states.  This will be a list of
#  formulas, each left hand side a state:marker pair, a formula which is
#  usually ~1 but can have other covariates, and / followed by options.
# The code partially mimics parsecovar.
#
# Separate out the states, markers, and options.
parsemarker <- function(flist, statedata) {
    if (inherits(flist, "formula")) flist <- list(flist) # only one marker
    if (any(sapply(flist, function(x) !inherits(x, "formula"))))
        stop("an element of the marker list is not a formula")
    if (any(sapply(flist, length) != 3))
        stop("all formulas must have a left and right side")
    
    # split the formulas into a right hand and left hand side
    lhs <- lapply(flist, function(x) x[-3])   # keep the ~
    rhs <- lapply(flist, function(x) x[[3]])  # don't keep the ~
    
    rhs <- rightslash(rhs) # defer parsing the options until later 
    

    # deal with the left hand side of the formula
    # the next routine cuts at '+' signs
    pcut <- function(form) {
        if (length(form)==3) {
            if (form[[1]] == '+') 
                c(pcut(form[[2]]), pcut(form[[3]]))
            else if (form[[1]] == '~') pcut(form[[2]])
            else list(form)
        }
        else list(form)
    }

    # cut the LHS into parts
    lcut <- lapply(lhs, function(x) pcut(x[[2]]))
    env1 <- new.env(parent= parent.frame(2))
    env2 <- new.env(parent= env1)

    # We allow people to say "state('dementia', 'death'), so make
    # "state" a function that is a variant of c()
    if (missing(statedata)) {
        assign("state", function(...) list(stateid= "state", 
                                           values=c(...)), env1)
        assign("state", list(stateid="state"))
    }
    else {
        # if there is statedata, then every column becomes such a function
        for (i in statedata) {
            assign(i, eval(list(stateid=i)), env2)
            tfun <- eval(parse(text=paste0("function(...) list(stateid='"
                                           , i, "', values=c(...))")))
            assign(i, tfun, env1)
        }
    }
    lterm <- lapply(lcut, function(x) {
        lapply(x, function(z) {
            if (length(z)==1) {
                eval(z, envir= env2)
            }
            else if (length(z) ==3 && z[[1]]==':') {
                # The right hand side of a colon should be a single variable
                #   name
                if (!is.name(z[[3]])) 
                    stop("marker does not map to a unique variable name: ",
                         deparse(z))
                list(states = eval(z[[2]], envir=env2), 
                     marker = z[[3]])
            }
            else stop("invalid term: ", deparse(z))
        })
    })

    list(lhs= lterm, rhs=rightslash(rhs))
}
