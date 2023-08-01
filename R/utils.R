match_dtype <- function(model_type=c("mvnorm", "multinomial", "norms")) {
    model_type <- match.arg(model_type)

    if (model_type == "mvnorm")
        dtype <- 0
    else if (model_type == "multinomial")
        dtype <- 1
    else if (model_type == "norms")
        dtype <- 2

    return(dtype)
}

LinearInit <- function(data, xdim=8, ydim=6) {
    pcm <- prcomp(data)
    pc <- pcm$rotation[, 1:2]
    sd <- pcm$sdev[1:2]
    mn <- apply(data, 2, mean)
    ans <- matrix(NA, xdim * ydim, dim(data)[2])
    colnames(ans) <- colnames(data)
    if (xdim >= ydim) {
        xtick <- sd[1] * pc[, 1]
        ytick <- sd[2] * pc[, 2]
    }
    else {
        xtick <- sd[2] * pc[, 2]
        ytick <- sd[1] * pc[, 1]
    }
    if (xdim == 1)
        xis <- rep(0, xdim)
    else xis <- seq(-2, 2, length = xdim)
    if (ydim == 1)
        yis <- rep(0, ydim)
    else yis <- seq(-2, 2, length = ydim)
    for (i in 1:(xdim * ydim)) {
        xi <- (i - 1)%%xdim + 1
        yi <- (i - 1)%/%xdim + 1
        ans[i, ] <- mn + xis[xi] * xtick + yis[yi] * ytick
    }
    ans
}

dist_from_adj <- function(adjmatrix) {
    if (length(adjmatrix) == 1)
        return(0)
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

cumsum_rshift <- function(x, left=0) {
    if (length(x) < 2)
        stop("the length of x should be longer than 1")

    x.cum <- cumsum(x)
    x.cum[2:length(x)] <- x.cum[1:(length(x)-1)]
    x.cum[1] <- left

    return(x.cum)
}
