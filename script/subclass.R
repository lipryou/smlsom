devtools::load_all("..")

library(MixSim)
library(mvtnorm)

n <- 1000
M <- 8
K <- 3
p <- 2

set.seed(123)
Q <- MixSim(MaxOmega=0.5, K=M, p=p, sph=T, hom=T, int=c(-3, 3))
A <- simdataset(n, Q$Pi, Q$Mu, Q$S)
X <- A$X

plot(X, col=A$id+1, type="n")
points(Q$Mu, pch=as.character(1:M))
plot.ellipses(Q$Mu, Q$S, 1:M)

plot(X, col=A$id+1)
points(Q$Mu, pch=16)
plot.ellipses(Q$Mu, Q$S, 1:M)

dd <- as.dist(1-Q$OmegaMap)
dd.hc <- hclust(dd, method="ward.D2")
plot(dd.hc)

classif <- cutree(dd.hc, K)

Y <- rep(classif, table(A$id))
plot(X, col=Y+1)
points(Q$Mu, pch=16)
plot.ellipses(Q$Mu, Q$S, 1:M)

nsubc <- table(classif)
cumnsubc <- cumsum(nsubc)
cumnsubc[2:K] <- cumnsubc[1:(K-1)]
cumnsubc[1] <- 0

dtype <- match_dtype("mvnorm")

beta <- 5
niter <- 1e+4
alpha <- c(0.05, 0.01)
radii <- c(2, -2)

Mhat <- 8
subclasses <- c(1, 1, 1, 1, 2, 2, 3, 3)
mu_list <- list()
Sigma_list <- list()
for (m in 1:Mhat) {
    mu_list[[m]] <- X[sample(which(Y==subclasses[m]), 1), ]
    Sigma_list[[m]] <- diag(p)
}

current_map <- list(M=Mhat, mu=mu_list, Sigma=Sigma_list)
adjmatrix <- matrix(0, M, M, dimnames=list(1:M, 1:M))
for (k in 1:K) {
    grid <- kohonen::somgrid(xdim=1, ydim=nsubc[k])
    adjmatrix_sub <- as.matrix(dist(grid$pts)) <= 1.5
    diag(adjmatrix_sub) = 0
    adjmatrix[subclasses==k, subclasses==k] <- adjmatrix_sub
}

## plot
Mus <- list_to_mvn.Mus(current_map)
Sigmas <- list_to_mvn.Sigmas(current_map)

plot(X, col=Y+1)
points(Mus, pch=16)
draw.nl(Mus, adjmatrix)

llconst <- loglikelihood_const(X, dtype)

## 1. do mlsom
nhbrdist <- dist_from_adj(adjmatrix)
current_map <- mlsom_clf(X, Y-1, nsubc, cumnsubc, current_map, dtype, niter, nhbrdist, alpha, radii)

## plot
Mus <- list_to_mvn.Mus(current_map)
Sigmas <- list_to_mvn.Sigmas(current_map)

plot(X, col=Y+1)
points(Mus, pch=16)
points(Q$Mu, pch=17)
plot.ellipses(Mus, Sigmas, 1:M)
draw.nl(Mus, adjmatrix)

## 2. link
logliks <- loglikelihood(X, current_map, llconst, dtype)
classes <- classifsubc_within_class(X, Y-1, nsubc, cumnsubc, current_map, dtype)

weight <- link_weight(classes, logliks, adjmatrix)
adjmatrix <- link_cutting(beta, weight, adjmatrix)

## plot
plot(X, col=Y+1)
points(Mus, pch=16)
points(Q$Mu, pch=17)
plot.ellipses(Mus, Sigmas, 1:M)
draw.nl(Mus, adjmatrix)

## 3. node deletion
