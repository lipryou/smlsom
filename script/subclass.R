devtools::load_all("..")

library(MixSim)
library(mvtnorm)

n <- 1000
MM <- 8
K <- 3
p <- 2

set.seed(123)
Q <- MixSim(MaxOmega=0.5, K=MM, p=p, sph=T, hom=T, int=c(-3, 3))
A <- simdataset(n, Q$Pi, Q$Mu, Q$S)
X <- A$X

plot(X, col=A$id+1, type="n")
points(Q$Mu, pch=as.character(1:MM))
plot.ellipses(Q$Mu, Q$S, 1:MM)

plot(X, col=A$id+1)
points(Q$Mu, pch=16)
plot.ellipses(Q$Mu, Q$S, 1:MM)

dd <- as.dist(1-Q$OmegaMap)
dd.hc <- hclust(dd, method="ward.D2")
plot(dd.hc)

classif <- cutree(dd.hc, K)

Y <- rep(classif, table(A$id))
plot(X, col=Y+1)
points(Q$Mu, pch=16)
plot.ellipses(Q$Mu, Q$S, 1:MM)
dtype <- match_dtype("mvnorm")

beta <- 5
niter <- 1e+4
alpha <- c(0.05, 0.01)

Mk <- 5
subclasses <- rep(1:K, rep(Mk, K))
M <- length(subclasses)

nsubc <- table(subclasses)
cumnsubc <- cumsum_rshift(nsubc)

mu_list <- list()
Sigma_list <- list()
set.seed(456)
for (m in 1:M) {
    mu_list[[m]] <- X[sample(which(Y==subclasses[m]), 1), ]
    Sigma_list[[m]] <- diag(p)
}

current_map <- list(M=M, mu=mu_list, Sigma=Sigma_list)
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

par(mar=c(1,1,1,1))
plot(X, col=Y+1, xlab="", ylab="", xaxt="n", yaxt="n")
points(Mus, pch=16)
draw.nl(Mus, adjmatrix)
points(Q$Mu, pch=17, col="red")

##########

llconst <- loglikelihood_const(X, dtype)
nhbrdist <- dist_from_adj(adjmatrix)
radii <- c(1, -1) * quantile(nhbrdist[nhbrdist!=Inf], .67)

## 1. do mlsom
nhbrdist <- dist_from_adj(adjmatrix)
current_map <- mlsom_clf(X, Y-1, nsubc, cumnsubc, current_map, dtype, niter, nhbrdist, alpha, radii)

## plot
Mus <- list_to_mvn.Mus(current_map)
Sigmas <- list_to_mvn.Sigmas(current_map)

plot(X, col=Y+1)
points(Mus, pch=16)
plot.ellipses(Mus, Sigmas, 1:M)
draw.nl(Mus, adjmatrix)
points(Q$Mu, pch=17, col="red")

## 2. link
logliks <- loglikelihood(X, current_map, llconst, dtype)
classes <- classifsubc_within_class(X, Y-1, nsubc, cumnsubc, current_map, dtype)

weight <- link_weight(classes, logliks, adjmatrix)
adjmatrix <- link_cutting(beta, weight, adjmatrix)

## plot
png(file="smlsomDAproc01.png", width=6, height=6, units="in", res=300)

par(mar=c(1,1,1,1))
plot(X, col=Y+1, xlab="", ylab="", xaxt="n", yaxt="n")
points(Mus, pch=16)
draw.nl(Mus, adjmatrix)
points(Q$Mu, pch=17, col="red")

dev.off()

## 3. node deletion

#### 3-1. current mdl
selected_list <- which_delete(X, current_map, classes, llconst, dtype, M, nsubc)
selected_m <- selected_list$selected_m
selected_map <- selected_list$selected_map

#### 3-2. mdl without each node
if (selected_m != 0) {
    M <- M - 1

    dec_class <- subclasses[selected_m]
    nsubc[dec_class] <- nsubc[dec_class] - 1
    cumnsubc <- cumsum_rshift(nsubc)

    subclasses <- subclasses[-selected_m]

    current_map <- selected_map

    adjmatrix <- accommodate_delete(adjmatrix, selected_m)
}

## plot
Mus <- list_to_mvn.Mus(current_map)
Sigmas <- list_to_mvn.Sigmas(current_map)

plot(X, col=Y+1)
points(Mus, pch=16)
points(Q$Mu, pch=17)
plot.ellipses(Mus, Sigmas, 1:M)
draw.nl(Mus, adjmatrix)

### ===========================================
###  test implemented smlsom_clf.mvn
### ===========================================

n <- 3000
MM <- 6
K <- 3
p <- 2

set.seed(122)
Q <- MixSim(BarOmega=0.05, K=MM, p=p, sph=F, hom=T, int=c(-3, 3))
A <- simdataset(n, Q$Pi, Q$Mu, Q$S)
X <- A$X

dd <- as.dist(1-Q$OmegaMap)
dd.hc <- hclust(dd, method="complete")
classif <- cutree(dd.hc, K)
Y <- rep(classif, table(A$id))

## plot
plot(X, col=Y+1)
points(Q$Mu, pch=16)
plot.ellipses(Q$Mu, Q$S, 1:MM)


beta <- 5
niter <- 1e+4
alpha <- c(0.05, 0.01)
Mk <- 5

system.time(clf.mvn <- smlsom_clf.mvn(X, Y, K=K, Mk=Mk, beta=beta,
                                      niter=niter, alpha=alpha, radii=c(2, -2)))

## plot
M <- clf.mvn$M
Mus <- clf.mvn$Mus
Sigmas <- clf.mvn$Sigmas
adjmatrix <- clf.mvn$adjmatrix

plot(X, col=Y+1)
points(Mus, pch=16)
points(Q$Mu, pch=17)
plot.ellipses(Mus, Sigmas, 1:M)
draw.nl(Mus, adjmatrix)

## prediction
sc <- clf.mvn$subclasses
cl_prior <- as.vector(table(Y) / n)
sc_prob <- table(clf.mvn$classes) / n
sc_prior <- as.vector(scratio / class_prob[sc])

### for train
scores <- matrix(0, K, n)
logliks <- clf.mvn$logliks
for (k in 1:K)
    scores[k, ] <- colSums(exp(t(logliks[, sc==k])) * sc_prior[sc==k])
scores <- t(log(scores) + log(cl_prior))
pred <- apply(scores, 1, which.max)

### for test
n.test <- 1000
B <- simdataset(n.test, Q$Pi, Q$Mu, Q$S)
Xte <- B$X
Yte <- classif[B$id]

## plot check
plot(X, col=Y+1)
plot(Xte, col=Yte+1)

### predict

pred <- predict.smlsom_clf.mvn(clf.mvn, Xte, type="c")
