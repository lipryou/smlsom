#' mlsom for multivariate normals
#'
#' This function provides mlsom clustering.
#'
#' @param data A numeric (n, p) matrix.
#' @param niter The number of iteration. Default is \code{nrow(data)}.
#' @param Mus A numeric (M, p) matrix. The rows are the mean vector of each cluster. Default is to cluster \code{data} randomly.
#' @param Sigmas A numeric (p, p, M) array. \code{Sigmas[,,m]} is the variance-covariance matrix of the mth cluster. Default is identity matrices.
#' @param nhbrdist Distance matrix between nodes. Note that "node" means a point in the SOM map space. Usually, a two-dimensional lattice space is used.
#' @param alpha Range of learning rates: \eqn{\alpha = (\alpha_1, \alpha_2)^t}. Note that \eqn{\alpha_1 > \alpha_2 > 0}. Monotonically decreasing from \eqn{\alpha_1} to \eqn{\alpha_2} for the number of iterations. Default is \code{alpha=c(0.05, 0.01)}. See \code{kohonen} package for details.
#' @param radii Range of neighbourhood radius: \eqn{r = (r_1, r_2)^t}. Note that \eqn{r_1 > r_2 > 0}. Monotonically decreasing from \eqn{r_1} to \eqn{r_2} for the number of iterations. Default is to start with a value that covers 2/3 of node distances. See \code{kohonen} package for detals.
#'
#' @return hogehoge
#'
#' @export
mlsom.mvn <- function(data, Mus, Sigmas, nhbrdist,
                      niter = nrow(data), alpha = c(0.05, 0.01),
                      radii = stats::quantile(nhbrdist,.67) * c(1, -1))
{
    dtype <- match_dtype("mvnorm")

    if (!is.numeric(data))
        stop("Argument data should be numeric")
    if (!(is.matrix(data) | is.table(data)))
        stop("Argument data should be matrix or table")

    if (!(is.matrix(nhbrdist)))
        stop("Argument data should be matrix")

    M <- nrow(nhbrdist)

    if (ncol(nhbrdist) != M)
        stop("Argument nhbrdist is not square")

    n <- nrow(data)
    p <- ncol(data)

    if (missing(Mus)) {
        cl <- sample(1:M, n, replace=T)
        Mus <- as.matrix(stats::aggregate(data, by=list(cl), mean)[,-1])
    } else {
        if (is.matrix(Mus))
            stop("Argument `Mus` should be matrix")
        if (nrow(Mus) != M)
            stop("The number of rows of `Mus` is not equal to M")
        if (ncol(Mus) != p)
            stop("The number of columns of `Mus` is not equal to p")
    }
    if (missing(Sigmas)) {
        Sigmas <- array(0, dim=c(p, p, M))
        for (m in 1:M)
            Sigmas[,,m] <- diag(1,p)
    } else {
        if (is.array(Sigmas))
            stop("Argument `Sigmas` should be array")
        if (all(dim(Sigmas) == c(p, p, M)))
            stop("`Sigmas` should be p x p x M array")
    }

    if (length(radii) == 1)
        radii <- sort(radii * c(1, -1), decreasing=TRUE)

    ## initial parameter list
    mu_list <- list()
    Sigma_list <- list()
    for (m in 1:M) {
        mu_list[[m]] <- Mus[m, ]
        Sigma_list[[m]] <- Sigmas[,,m]
    }
    params_init <- list(M=M, mu=mu_list, Sigma=Sigma_list)

    params <- mlsom(data, params_init, dtype, niter, nhbrdist, alpha, radii)

    ## get estimated parameters
    mu_list <- params$mu
    Sigma_list <- params$Sigma
    for (m in 1:M) {
        Mus[m, ] <- mu_list[[m]]
        Sigmas[,,m] <- Sigma_list[[m]]
    }

    logliks <- loglikelihood(data, params, dtype)

    classes <- apply(logliks, 1, which.max)

    structure(list(Mus=Mus, Sigmas=Sigmas, logliks=logliks, classes=classes,
                   nhbrdist=nhbrdist, alpha=alpha, radii=radii),
              class = "mlsom.mvn")
}

#' smlsom for multivariate normals
#'
#' This function provides mlsom clustering.
#'
#' @param data A numeric (n, p) matrix.
#' @param xdim hogehoge
#' @param ydim hogehoge
#' @param niter The number of iteration. Default is \code{nrow(data)}.
#' @param Mus A numeric (M, p) matrix. The rows are the mean vector of each cluster. Default is to cluster \code{data} randomly.
#' @param Sigmas A numeric (p, p, M) array. \code{Sigmas[,,m]} is the variance-covariance matrix of the mth cluster. Default is identity matrices.
#' @param alpha Range of learning rates: \eqn{\alpha = (\alpha_1, \alpha_2)^t}. Note that \eqn{\alpha_1 > \alpha_2 > 0}. Monotonically decreasing from \eqn{\alpha_1} to \eqn{\alpha_2} for the number of iterations. Default is \code{alpha=c(0.05, 0.01)}. See \code{kohonen} package for details.
#'
#' @return hogehoge
#'
#' @export
smlsom.mvn <- function(data, xdim, ydim, beta=5, niter = nrow(data),
                       Mus, Sigmas, alpha = c(0.05, 0.01),
                       radii = stats::quantile(nhbrdist,.67) * c(1, -1)) {
    dtype <- match_dtype("mvnorm")

    if (!is.numeric(data))
        stop("Argument data should be numeric")
    if (!(is.matrix(data) | is.table(data)))
        stop("Argument data should be matrix or table")

    if ((length(xdim)==1) & (length(ydim)==1)) {
        if (!(is.integer(xdim) | !(is.integer(ydim))))
            stop("Argument xdim and ydim are should be scalar integer")
    }

    M <- xdim * ydim
    n <- nrow(data)
    p <- ncol(data)

    if (missing(Mus)) {
        cl <- sample(1:M, n, replace=T)
        Mus <- as.matrix(stats::aggregate(data, by=list(cl), mean)[,-1])
    } else {
        if (is.matrix(Mus))
            stop("Argument `Mus` should be matrix")
        if (nrow(Mus) != M)
            stop("The number of rows of `Mus` is not equal to M")
        if (ncol(Mus) != p)
            stop("The number of columns of `Mus` is not equal to p")
    }
    if (missing(Sigmas)) {
        Sigmas <- array(0, dim=c(p, p, M))
        for (m in 1:M)
            Sigmas[,,m] <- diag(1,p)
    } else {
        if (is.array(Sigmas))
            stop("Argument `Sigmas` should be array")
        if (all(dim(Sigmas) == c(p, p, M)))
            stop("`Sigmas` should be p x p x M array")
    }

    ## create som map
    grid <- kohonen::somgrid(xdim, ydim, topo="h")
    adjmatrix <- as.matrix(dist(grid$pts)) <= 1.5
    diag(adjmatrix) = 0

    nhbrdist <- dist_from_adj(adjmatrix)
    if (length(radii) == 1)
        radii <- sort(radii * c(1, -1), decreasing=TRUE)

    ## create params
    mu_list <- list()
    Sigma_list <- list()
    for (m in 1:M) {
        mu_list[[m]] <- Mus[m, ]
        Sigma_list[[m]] <- Sigmas[,,m]
    }
    current_map <- list(M=M, mu=mu_list, Sigma=Sigma_list)

    llconst <- loglikelihood_const(data, dtype)

    while(M > 1) {
        result <- do_smlsom(data, current_map, adjmatrix, nhbrdist,
                            beta, niter, alpha, radii, llconst, dtype)
        current_map <- result$current_map
        adjmatrix_new <- result$adjmatrix

        if (all(dim(adjmatrix) == dim(adjmatrix_new)))
            break

        adjmatrix <- adjmatrix_new
        nhbrdist <- dist_from_adj(adjmatrix)
        M <- current_map$M
    }

    ## get estimated parameters
    Mus <- list_to_mvn.Mus(current_map)
    Sigmas <- list_to_mvn.Sigmas(current_map)

    logliks <- loglikelihood(data, current_map, llconst, dtype)
    classes <- apply(logliks, 1, which.max)

    cmdl <- classif_mdl(data, current_map, classes, llconst, dtype)

    structure(list(M=M, Mus=Mus, Sigmas=Sigmas, logliks=logliks, classes=classes,
                   adjmatrix=adjmatrix, beta=beta, alpha=alpha, radii=radii, cmdl=cmdl),
              class = "smlsom.mvn")
}

do_smlsom <- function(data, current_map, adjmatrix, nhbrdist,
                      beta, niter, alpha, radii, llconst, dtype) {
    M <- current_map$M

    ## mlsom
    current_map <- mlsom(data, current_map, dtype, niter, nhbrdist, alpha, radii)

    ## link cutting
    logliks <- loglikelihood(data, current_map, llconst, dtype)
    classes <- apply(logliks, 1, which.max)

    current_map <- onebatch(data, current_map, classes-1, dtype)

    weight <- link_weight(classes, logliks, adjmatrix)
    adjmatrix <- link_cutting(beta, weight, adjmatrix)

    ## node deletion
    selected_list <- which_delete(data, current_map, classes, llconst, dtype, M)
    selected_m <- selected_list$selected_m
    selected_map <- selected_list$selected_map

    ## map update.
    if (selected_m != 0) {
        current_map <- selected_map
        ## compensation for deleted nodes
        adjmatrix <- accommodate_delete(adjmatrix, selected_m)
    }

    return(list(current_map=current_map, adjmatrix=adjmatrix))
}
