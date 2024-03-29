#' mlsom for multivariate normals
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
#' @param radii Range of neighbourhood radius: \eqn{r = (r_1, r_2)^t}. Note that \eqn{r_1 > r_2 > 0}. Monotonically decreasing from \eqn{r_1} to \eqn{r_2} for the number of iterations. Default is to start with a value that covers 2/3 of node distances. See \code{kohonen} package for detals.
#' @param cov.type hogehoge
#'
#' @return hogehoge
#'
#' @export
mlsom.mvn <- function(data, xdim, ydim, Mus, Sigmas,
                      niter = nrow(data), alpha = c(0.05, 0.01),
                      radii = stats::quantile(nhbrdist,.67) * c(1, -1),
                      cov.type = c("full", "diag"))
{
    cov.type <- match.arg(cov.type)
    if (cov.type == "full")
        dtype <- match_dtype("mvnorm")
    else if (cov.type == "diag")
        dtype <- match_dtype("norms")

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
        if (!is.matrix(Mus))
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
        if (!is.array(Sigmas))
            stop("Argument `Sigmas` should be array")
        if (any(dim(Sigmas) != c(p, p, M)))
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
    for (m in 1:M)
        mu_list[[m]] <- Mus[m, ]
    if (cov.type == "full") {
        for (m in 1:M)
            Sigma_list[[m]] <- Sigmas[,,m]
    } else if (cov.type == "diag") {
        for (m in 1:M)
            Sigma_list[[m]] <- diag(Sigmas[,,m])
    }

    params_init <- list(M=M, mu=mu_list, Sigma=Sigma_list)
    current_map <- mlsom(data, params_init, dtype, niter, nhbrdist, alpha, radii)

    ## get estimated parameters
    Mus <- list_to_mvn.Mus(current_map)
    if (cov.type == "full")
        Sigmas <- list_to_mvn.Sigmas(current_map)
    else if (cov.type == "diag") {
        Sigmas <- array(0, dim=c(p, p, M))
        for (m in 1:M)
            Sigmas[,,m] <- diag(current_map$Sigma[[m]])
    }

    llconst <- loglikelihood_const(data, dtype)
    logliks <- loglikelihood(data, current_map, llconst, dtype)
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
#' @param beta hogehoge
#' @param niter The number of iteration. Default is \code{nrow(data)}.
#' @param Mus A numeric (M, p) matrix. The rows are the mean vector of each cluster. Default is to cluster \code{data} randomly.
#' @param Sigmas A numeric (p, p, M) array. \code{Sigmas[,,m]} is the variance-covariance matrix of the mth cluster. Default is identity matrices.
#' @param alpha Range of learning rates: \eqn{\alpha = (\alpha_1, \alpha_2)^t}. Note that \eqn{\alpha_1 > \alpha_2 > 0}. Monotonically decreasing from \eqn{\alpha_1} to \eqn{\alpha_2} for the number of iterations. Default is \code{alpha=c(0.05, 0.01)}. See \code{kohonen} package for details.
#' @param radii Range of neighbourhood radius: \eqn{r = (r_1, r_2)^t}. Note that \eqn{r_1 > r_2 > 0}. Monotonically decreasing from \eqn{r_1} to \eqn{r_2} for the number of iterations. Default is to start with a value that covers 2/3 of node distances. See \code{kohonen} package for detals.
#' @param cov.type hogehoge
#'
#' @return hogehoge
#'
#' @export
smlsom.mvn <-
    function(data, xdim, ydim, beta=5, niter = nrow(data),
             Mus, Sigmas, alpha = c(0.05, 0.01),
             radii = stats::quantile(nhbrdist,.67) * c(1, -1),
             cov.type = c("full", "diag"), verbose=F)
{
    cov.type <- match.arg(cov.type)
    if (cov.type == "full")
        dtype <- match_dtype("mvnorm")
    else if (cov.type == "diag")
        dtype <- match_dtype("norms")

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
        if (!is.matrix(Mus))
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
        if (!is.array(Sigmas))
            stop("Argument `Sigmas` should be array")
        if (any(dim(Sigmas) != c(p, p, M)))
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
    for (m in 1:M)
        mu_list[[m]] <- Mus[m, ]
    if (cov.type == "full") {
        for (m in 1:M)
            Sigma_list[[m]] <- Sigmas[,,m]
    } else if (cov.type == "diag") {
        for (m in 1:M)
            Sigma_list[[m]] <- diag(Sigmas[,,m])
    }

    current_map <- list(M=M, mu=mu_list, Sigma=Sigma_list)

    llconst <- loglikelihood_const(data, dtype)

    while(M > 1) {
        if (verbose)
            cat(sprintf("Current M = %d\n", M))

        result <- do_smlsom(data, current_map, adjmatrix, nhbrdist,
                            beta, niter, alpha, radii, llconst, dtype, verbose)
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
    if (cov.type == "full")
        Sigmas <- list_to_mvn.Sigmas(current_map)
    else if (cov.type == "diag") {
        Sigmas <- array(0, dim=c(p, p, M))
        for (m in 1:M)
            Sigmas[,,m] <- diag(current_map$Sigma[[m]])
    }

    logliks <- loglikelihood(data, current_map, llconst, dtype)
    classes <- apply(logliks, 1, which.max)

    cmdl <- classif_mdl(data, current_map, classes, llconst, dtype)

    structure(list(M=M, Mus=Mus, Sigmas=Sigmas, logliks=logliks, classes=classes,
                   adjmatrix=adjmatrix, beta=beta, alpha=alpha, radii=radii, cmdl=cmdl),
              class = "smlsom.mvn")
}

#' smlsom for classification problem using multivariate normals
#'
#' This function provides smlsom clustering.
#'
#' @param X A numeric (n, p) matrix.
#' @param y class labels. n integer vector.
#' @param K The number of classes.
#' @param Mk The number of subclasses within each class.
#' @param beta hogehoge
#' @param niter The number of iteration. Default is \code{nrow(data)}.
#' @param Mus A numeric (M, p) matrix. The rows are the mean vector of each cluster. Default is to cluster \code{data} randomly.
#' @param Sigmas A numeric (p, p, M) array. \code{Sigmas[,,m]} is the variance-covariance matrix of the mth cluster. Default is identity matrices.
#' @param alpha Range of learning rates: \eqn{\alpha = (\alpha_1, \alpha_2)^t}. Note that \eqn{\alpha_1 > \alpha_2 > 0}. Monotonically decreasing from \eqn{\alpha_1} to \eqn{\alpha_2} for the number of iterations. Default is \code{alpha=c(0.05, 0.01)}. See \code{kohonen} package for details.
#' @param radii Range of neighbourhood radius: \eqn{r = (r_1, r_2)^t}. Note that \eqn{r_1 > r_2 > 0}. Monotonically decreasing from \eqn{r_1} to \eqn{r_2} for the number of iterations. Default is to start with a value that covers 2/3 of node distances. See \code{kohonen} package for detals.
#' @param cov.type hogehoge
#'
#' @return hogehoge
#'
#' @export
smlsom_clf.mvn <- function(X, y, K, Mk = 5, beta = 5, niter = nrow(X),
                           Mus, Sigmas, alpha = c(0.05, 0.01),
                           radii, cov.type = c("full", "diag")) {
    cov.type <- match.arg(cov.type)
    if (cov.type == "full")
        dtype <- match_dtype("mvnorm")
    else if (cov.type == "diag")
        dtype <- match_dtype("norms")

    if (!is.numeric(X))
        stop("Argument X should be numeric")
    if (!(is.matrix(X) | is.table(X)))
        stop("Argument X should be matrix or table")

    n <- nrow(X)
    p <- ncol(X)

    if (!is.integer(y))
        stop("Argument y should be integer")
    if (!is.vector(y))
        stop("Argument y should be vector")

    class_size <- table(y)
    if (length(class_size) != K)
        stop("Class label y contains larger than K classes")

    M <- Mk * K
    subclasses <- rep(1:K, rep(Mk, K))
    nsubc <- table(subclasses)

    if (missing(Mus)) {
        Mus <- matrix(0, M, p)
        for (m in 1:M)
            Mus[m, ] <- X[sample(which(y==subclasses[m]), 1), ]
    } else {
        if (!is.matrix(Mus))
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
        if (!is.array(Sigmas))
            stop("Argument `Sigmas` should be array")
        if (any(dim(Sigmas) != c(p, p, M)))
            stop("`Sigmas` should be p x p x M array")
    }

    ## create som map
    adjmatrix <- matrix(0, M, M, dimnames=list(1:M, 1:M))
    for (k in 1:K) {
        grid <- kohonen::somgrid(xdim=nsubc[k], ydim=1)
        adjmatrix_sub <- as.matrix(dist(grid$pts)) <= 1.5
        diag(adjmatrix_sub) = 0
        adjmatrix[subclasses==k, subclasses==k] <- adjmatrix_sub
    }

    nhbrdist <- dist_from_adj(adjmatrix)
    if (missing(radii))
        radii <- stats::quantile(nhbrdist[nhbrdist!=Inf],.67) * c(1, -1)
    if (length(radii) == 1)
        radii <- sort(radii * c(1, -1), decreasing=TRUE)

    ## create params
    mu_list <- list()
    Sigma_list <- list()
    for (m in 1:M)
        mu_list[[m]] <- Mus[m, ]
    if (cov.type == "full") {
        for (m in 1:M)
            Sigma_list[[m]] <- Sigmas[,,m]
    } else if (cov.type == "diag") {
        for (m in 1:M)
            Sigma_list[[m]] <- diag(Sigmas[,,m])
    }

    current_map <- list(M=M, mu=mu_list, Sigma=Sigma_list)

    llconst <- loglikelihood_const(X, dtype)

    while(M > 1) {
        if (verbose)
            cat(sprintf("Current M = %d\n", M))

        cumnsubc <- cumsum_rshift(nsubc)

        result <- do_smlsom_clf(X, y, subclasses, nsubc, cumnsubc, current_map, adjmatrix,
                                nhbrdist, beta, niter, alpha, radii, llconst, dtype, verbose)
        current_map <- result$current_map
        adjmatrix_new <- result$adjmatrix
        nsubc <- result$nsubc
        subclasses <- result$subclasses

        if (all(dim(adjmatrix) == dim(adjmatrix_new)))
            break

        adjmatrix <- adjmatrix_new
        nhbrdist <- dist_from_adj(adjmatrix)
        M <- current_map$M
    }

    ## get estimated parameters
    cumnsubc <- cumsum_rshift(nsubc)
    Mus <- list_to_mvn.Mus(current_map)
    if (cov.type == "full")
        Sigmas <- list_to_mvn.Sigmas(current_map)
    else if (cov.type == "diag") {
        Sigmas <- array(0, dim=c(p, p, M))
        for (m in 1:M)
            Sigmas[,,m] <- diag(current_map$Sigma[[m]])
    }

    logliks <- loglikelihood(X, current_map, llconst, dtype)
    classes <- classifsubc_within_class(X, y-1, nsubc, cumnsubc, current_map, dtype)

    cmdl <- classif_mdl(X, current_map, classes, llconst, dtype)

    ## calcluate prior probs.
    cl_prior <- as.vector(table(y) / n)
    sc_prob <- table(classes) / n
    sc_prior <- as.vector(sc_prob / cl_prior[subclasses])

    structure(list(M=M, Mus=Mus, Sigmas=Sigmas, logliks=logliks, classes=classes,
                   adjmatrix=adjmatrix, beta=beta, alpha=alpha, radii=radii, cmdl=cmdl,
                   subclasses=subclasses, nsubc=nsubc,
                   class_prior=cl_prior, subclass_prior=sc_prior),
              class = "smlsom_clf.mvn")
}


do_smlsom <-
    function(data, current_map, adjmatrix, nhbrdist,
             beta, niter, alpha, radii, llconst, dtype, verbose)
{
    M <- current_map$M

    ## mlsom
    current_map <- mlsom(data, current_map, dtype,
                         niter, nhbrdist, alpha, radii)

    ## link cutting
    logliks <- loglikelihood(data, current_map, llconst, dtype)
    classes <- apply(logliks, 1, which.max)

    if(verbose) {
        cat(" Distribution of cluster sizes\n")
        print(summary(as.vector(table(classes))))
    }

    current_map <- onebatch(data, current_map, classes-1, dtype)

    weight <- link_weight(classes, logliks, adjmatrix)
    adjmatrix <- link_cutting(beta, weight, adjmatrix)

    ## node deletion
    selected_list <- which_delete(data, current_map, classes, llconst, dtype, M)
    selected_m <- selected_list$selected_m
    selected_map <- selected_list$selected_map

    if (verbose)
        cat(sprintf(" current_mdl = %.3f previous mdl = %.3f\n",
                selected_list$cmdl, selected_list$cmdl.prev))

    ## map update.
    if (selected_m != 0) {
        current_map <- selected_map
        ## compensation for deleted nodes
        adjmatrix <- accommodate_delete(adjmatrix, selected_m)
    }

    return(list(current_map=current_map, adjmatrix=adjmatrix))
}

do_smlsom_clf <-
    function(X, Y, subclasses, nsubc, cumnsubc,
             current_map, adjmatrix, nhbrdist,
             beta, niter, alpha, radii, llconst, dtype, verbose)
{
    M <- current_map$M

    ## mlsom
    current_map <- mlsom_clf(X, Y-1, nsubc, cumnsubc, current_map, dtype, niter, nhbrdist, alpha, radii)

    ## link cutting
    logliks <- loglikelihood(X, current_map, llconst, dtype)
    classes <- classifsubc_within_class(X, Y-1, nsubc, cumnsubc, current_map, dtype)

    if(verbose) {
        cat(" Distribution of cluster sizes\n")
        print(summary(as.vector(table(classes))))
    }

    weight <- link_weight(classes, logliks, adjmatrix)
    adjmatrix <- link_cutting(beta, weight, adjmatrix)

    ## node deletion
    selected_list <- which_delete(X, current_map, classes, llconst, dtype, M, nsubc)
    selected_m <- selected_list$selected_m
    selected_map <- selected_list$selected_map

    if (verbose)
        cat(sprintf(" current_mdl = %.3f previous mdl = %.3f\n",
                    selected_list$cmdl, selected_list$cmdl.prev))

    ## map update.
    if (selected_m != 0) {
        M <- M - 1

        dec_class <- subclasses[selected_m]
        nsubc[dec_class] <- nsubc[dec_class] - 1

        subclasses <- subclasses[-selected_m]

        current_map <- selected_map

        adjmatrix <- accommodate_delete(adjmatrix, selected_m)
    }

    return(list(current_map=current_map, adjmatrix=adjmatrix,
                nsubc=nsubc, subclasses=subclasses))
}

#' prediction function for smlsom_clf.mvn
#'
#' This function provides smlsom clustering.
#'
#' @param object smlsom_clf.mvn object.
#' @param newdata a numeric (n, p) matrix.
#' @param type type of return. \code{class} gives a vector of length n whose elements are the classes with the largest posterior probability. \code{probability} provides the posterior probabilities matrix.
#'
#' @return hogehoge
#'
#' @export
predict.smlsom_clf.mvn <-
    function(object, newdata, type=c("class", "probability"))
{
    type <- match.arg(type)

    n <- nrow(newdata)
    K <- length(object$class_prior)

    sc <- object$subclasses
    scores <- matrix(0, K, n)
    logliks <- loglikelihood.mvn(newdata, object$Mus, object$Sigmas)
    for (k in 1:K)
        scores[k, ] <- colSums(exp(t(logliks[, sc==k])) * object$subclass_prior[sc==k])

    scores <- t(log(scores) + log(object$class_prior))

    if (type == "class")
        pred <- apply(scores, 1, which.max)
    else if (type == "probability")
        pred <- exp(scores) / rowSums(exp(scores))

    return(pred)
}
