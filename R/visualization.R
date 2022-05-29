plot.ellipses <- function(Mus, Sigs, labs, level=.8, col="gray") {
    if (length(col) == 1) {
        col <- rep(col,nrow(Mus))
        names(col) <- labs
    }
    for (m in labs) {
        points(ellipse::ellipse(Sigs[,,m],
                                centre=Mus[m,],
                                level=level),
               type="l",col=col[m])
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
