link_weight <- function(classes, logliks, adjmatrix) {
    ## calculate link weight
    M <- nrow(adjmatrix)
    Dkl <- matrix(0, M, M)

    active_nodes <- which(apply(adjmatrix==1, 1, any))
    for (m in active_nodes) {
        classes_m <- classes == m
        neighs <- which(adjmatrix[m,]==1)

        ll.origin <- logliks[classes_m, m]
        ll.move <- logliks[classes_m, neighs, drop=F]

        Dkl[m, neighs] <- colMeans(ll.origin - ll.move)
        Dkl[m, m] <- mean(ll.origin)
    }

    weights <- 0.5 * Dkl + 0.5 * t(Dkl)

    return(weights)
}

link_cutting <- function(beta, weights, adjmatrix) {
    ## cut link
    adj.tmp <- adjmatrix
    M <- nrow(adjmatrix)
    active_nodes <- which(apply(adjmatrix==1, 1, any))

    ## threshold
    h <- max(-diag(weights))
    ## h <- ifelse(h < 0, 0.1, h) ## necessary?

    for (m in active_nodes) {
        Dml <- weights[m, ]
        cutp <- Dml > beta * h
        adj.tmp[m, which(cutp)] <- 0
        adj.tmp[which(cutp), m] <- 0
    }
    adjmatrix <- adj.tmp

    return(adjmatrix)
}

which_delete <- function(X, current_map, classes,
                         llconst, dtype, M, nsubc) {
    cmdl <- classif_mdl(X, current_map, classes, llconst, dtype)
    cmdl.prev <- cmdl
    selected_m <- 0
    selected_map <- NULL

    ## Determine start and end point for range of class search.
    ## Usually, start is 0 and end point is the number of clusters.
    ## However, if nsubc (Mk, the number of subclasses within class) is given,
    ## it is subclass detection task, so the algorithm should search the
    ## nearest subclass within the class which xi belongs to.
    ## start_m and end_m are arguments of the cpp function `find_nearest`
    ## to controll those search range.
    if (missing(nsubc)) {
        start_m <- rep(0, M)
        end_m <- rep(M-1, M)
    } else {
        K <- length(nsubc)
        start_m <- c(0, cumsum(nsubc[0:(K-1)]))
        end_m <- cumsum(nsubc) - 1

        names(start_m) <- 1:K
        start_m <- rep(start_m, nsubc)
        end_m <- rep(end_m, nsubc)
    }

    for (target_m in 1:M) {
        mrange <- c(start_m[target_m], end_m[target_m])
        if ((mrange[2] - mrange[1]) == 0) next
        candidate_map <- eval_without(target_m-1, X, mrange, current_map,
                                      classes-1, llconst, dtype)
        cmdl_c <- candidate_map$value
        if (cmdl > cmdl_c) {
            cmdl <- cmdl_c
            selected_m <- target_m
            selected_map <- candidate_map
        }
    }
    return(list(selected_m=selected_m, selected_map=selected_map,
                cmdl=cmdl, cmdl.prev=cmdl.prev))
}

accommodate_delete <- function(adjmatrix, selected_m) {
    neighs <- which(adjmatrix[selected_m, ]==1)
    for (i in neighs){
        adjmatrix[i, neighs[neighs != i]] <- 1
        adjmatrix[neighs[neighs != i], i] <- 1
    }
    adjmatrix <- adjmatrix[-selected_m, -selected_m]

    return(adjmatrix)
}

loglikelihood.mvn <- function(X, Mus, Sigmas) {
    if (ncol(X) != ncol(Mus))
        stop("ncol(X) and ncol(Mus) should be equal")
    if (ncol(X) != dim(Sigmas)[1L])
        stop("ncol(X) and dim(Mus)[1] should be equal")
    if (dim(Sigmas)[1L] != dim(Sigmas)[2L])
        stop("Each Sigmas should be a square matrix")

    M <- nrow(Mus)
    mu_list <- list()
    Sigma_list <- list()
    for (m in 1:M) {
        mu_list[[m]] <- Mus[m, ]
        Sigma_list[[m]] <- Sigmas[,,m]
    }

    params <- list(M=M, mu=mu_list, Sigma=Sigma_list)
    dtype <- match_dtype("mvnorm")
    llconst <- loglikelihood_const(X, dtype)
    logliks <- loglikelihood(X, params, llconst, dtype)

    return(logliks)
}

## =============================================================
##   functions for mdl calculation
## =============================================================

classif_mdl <- function(data, parameters, classes, llconst, dtype) {
    p <- ncol(data)
    M <- parameters$M
    if (dtype == 0)
        df <- dof.mvn(p, M)
    else if (dtype == 1)
        df <- dof.multinom(p, M)
    else if (dtype == 2)
        df <- dof.norms(p, M)
    else {
        stop("dtype error: such class not implemented")
    }

    n <- nrow(data)
    logliks_c <- classification_loglikelihood(data, parameters, classes-1, llconst, dtype)

    return(-sum(logliks_c) + 0.5 * df * log(n) + n*log(M))
}

dof.mvn <- function(p, M) {
    return(0.5 * M * p * (p + 3))
}

dof.multinom <- function(p, M) {
    return(M * (p-1))
}

dof.norms <- function(p, M) {
    return(2*M*p)
}


## =============================================================
##   functions for utils
## =============================================================

list_to_mvn.Mus <- function(params) {
    M <- params$M
    mu_list <- params$mu
    Mus <- c()
    for (m in 1:M)
        Mus <- rbind(Mus, mu_list[[m]])
    return(Mus)
}

list_to_mvn.Sigmas <- function(params) {
    M <- params$M
    p <- ncol(params$Sigma[[1]])
    Sigma_list <- params$Sigma
    Sigmas <- array(0, dim=c(p, p, M))
    for (m in 1:M)
        Sigmas[,,m] <- Sigma_list[[m]]
    return(Sigmas)
}
