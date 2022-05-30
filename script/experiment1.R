devtools::load_all("..")

library(MixSim)
library(mvtnorm)

source("~/Document/Rfuncs/SMLSOM.R")

### data gen
n <- 3000
M <- 6
p <- 2

dtype <- match_dtype("mvnorm")

beta <- 5
niter <- 3e+3
alpha <- c(0.05, 0.01)
radii <- c(2, -2)

### smlsom
ndat <- 20
nrep <- 100
performs <- array(0, dim=c(nrep, 6, ndat))

for (i in 1:ndat) {
    seed1 <- i * 100
    set.seed(seed1)
    Q <- MixSim(BarOmega=0.05, K=M, p=p, int=c(-3, 3), sph=TRUE)
    A <- simdataset(n, Q$Pi, Q$Mu, Q$S)
    X <- A$X

    cat("dataset", i, ": ", date(), "\n", sep="")
    for (tt in 1:nrep) {
        seed2 <- i*100 + tt

        set.seed(seed2)
        tm1 <- system.time(ret1 <- smlsom::smlsom.mvn(X, 3, 3, beta, niter, radii=radii))
        set.seed(seed2)
        tm2 <- system.time(ret2 <- smlsom.mvn(X, 3, 3, rlen=1, beta=beta, radius=radii,
                                              init.type="random", silent=1))

        performs[tt, 1, i] <- tm1["elapsed"]
        performs[tt, 2, i] <- RandIndex(ret1$classes, A$id)$AR
        performs[tt, 3, i] <- ret1$M

        performs[tt, 4, i] <- tm2["elapsed"]
        performs[tt, 5, i] <- RandIndex(as.integer(ret2$unit.classif), A$id)$AR
        performs[tt, 6, i] <- ret2$nd
    }
    save(list="performs", file="test.Rdata")
    cat("...dataset",i," done: ", date(), "\n", sep="")
}

