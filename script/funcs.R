## =============================================================
##   functions for smlsom calculation
## =============================================================

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
    h <- ifelse(h < 0, 0.5, h) ## necessary?

    for (m in active_nodes) {
        Dml <- weights[m, ]
        cutp <- Dml > beta * h
        adj.tmp[m, which(cutp)] <- 0
        adj.tmp[which(cutp), m] <- 0
    }
    adjmatrix <- adj.tmp

    return(adjmatrix)
}

which_delete <- function(X, current_map, classes, dtype, Mhat) {
    cmdl <- classif_mdl(X, current_map, classes, dtype)
    cmdl.prev <- cmdl
    selected_m <- 0
    for (target_m in 1:Mhat) {
        candidate_map <- eval_without(target_m-1, X, current_map, classes-1, dtype)
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
    adjmatrix <- adjmatrix[-selected_m, ]
    adjmatrix <- adjmatrix[, -selected_m]

    return(adjmatrix)
}


## =============================================================
##   functions for mdl calculation
## =============================================================

classif_mdl <- function(data, parameters, classes, dtype) {
    p <- ncol(data)
    M <- parameters$M
    if (dtype == 0)
        df <- dof.mvn(p, M)
    else if (dtype == 1)
        df <- dof.multinom(p, M)
    else {
        stop("dtype error: such class not implemented")
    }

    n <- nrow(data)
    logliks_c <- classification_loglikelihood(data, parameters, classes-1, dtype)

    return(-sum(logliks_c) + 0.5 * df * log(n) + n*log(M))
}

dof.mvn <- function(p, M) {
    return(0.5 * M * p * (p + 3))
}

dof.multinom <- function(p, M) {
    return(M * (p-1))
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

plot_grid <- function(adjmatrix, grid, Mhat) {
    plot(grid$pts, pch="")
    text(grid$pts, labels=1:Mhat)
    for (m in 1:Mhat) {
        neigh <- which(adjmatrix[m, ]==1)
        if (length(neigh) == 0)
            next
        segments(x0=grid$pts[m, 1], y0=grid$pts[m, 2],
                 x1=grid$pts[neigh, 1], y1=grid$pts[neigh, 2])
    }
}

draw.nl <- function(code, neigh,
                    lwd=2, col="black",
                    spf.node=NULL, spf.col="blue") {
    in.draw.nl <- function(i,col) {
        tmp <- code[neigh[, i]==1,,drop=F]
        tmp.nr <- nrow(tmp)
        if (tmp.nr > 0)
            segments(x0=rep(code[i, 1], tmp.nr),
                     y0=rep(code[i, 2], tmp.nr),
                     x1=tmp[, 1],
                     y1=tmp[, 2],
                     lwd=lwd, col=col)
    }
    ng <- nrow(code)
    for (i in setdiff(1:ng, spf.node)) #i <- 4
        in.draw.nl(i, col)
    if (!is.null(spf.node))
        in.draw.nl(spf.node, spf.col)
}

dist_from_adj <- function(adjmatrix) {
    M <- nrow(adjmatrix)
    mat <- matrix(as.integer(adjmatrix), nrow=M)
    A <- mat
    ddmat <- ifelse(adjmatrix==1, 1, Inf)
    for (n in 2:M) {
        mat <- mat %*% A
        ddmat[(mat > 0) & (ddmat==Inf)] <- n
    }
    diag(ddmat) <- 0

    return(ddmat)
}
