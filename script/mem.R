devtools::load_all("..")

library(MixSim)
library(mvtnorm)

## ======================================================================
##   loglikelihood:: multivariate normal
## ======================================================================

n <- 1000
M <- 2
p <- 10
Q <- MixSim(BarOmega=0.1, K=M, p=p)
A <- simdataset(n, Q$Pi, Q$Mu, Q$S)

mu_list <- list()
Sigma_list <- list()
for (m in 1:M) {
    mu_list[[m]] <- Q$Mu[m, ]
    Sigma_list[[m]] <- Q$S[,,m]
}

dtype <- match_dtype("mvnorm")
params <- list(M=M, mu=mu_list, Sigma=Sigma_list)
llconst <- loglikelihood_const(A$X, dtype)
logliks1 <- loglikelihood(A$X, params, llconst, dtype)

logliks2 <- sapply(1:M, function(m) dmvnorm(A$X, Q$Mu[m, ], Q$S[,,m], log=T))

norm(logliks1 - logliks2)

## ======================================================================
##   loglikelihood:: multinomial
## ======================================================================

M <- 2

theta1 <- c(0.7, 0.3)
theta2 <- c(0.3, 0.7)

n1 <- rpois(10, 10) + 1
x1 <- t(sapply(n1, function(n) rmultinom(1, n, theta1)))
n2 <- rpois(10, 10) + 1
x2 <- t(sapply(n1, function(n) rmultinom(1, n, theta2)))

n <- c(n1, n2)
X <- rbind(x1, x2)

dtype <- match_dtype("multinom")
params <- list(M=M, theta=list(theta1, theta2))
llconst <- loglikelihood_const(X, dtype)
logliks1 <- loglikelihood(X, params, llconst, dtype)

logliks2 <- sapply(1:length(n), function(i) dmultinom(X[i, ], prob=theta1, log=T))
logliks2 <- cbind(logliks2, sapply(1:length(n), function(i) dmultinom(X[i, ], prob=theta2, log=T)))

norm(logliks1 - logliks2)

## ======================================================================
##   mlsom.mvn
## ======================================================================

M <- 4
Q <- MixSim(BarOmega=0.05, K=M, p=2)
A <- simdataset(1000, Q$Pi, Q$Mu, Q$S)

mu_list <- list()
Sigma_list <- list()
for (m in 1:M) {
    mu_list[[m]] <- Q$Mu[m, ]
    Sigma_list[[m]] <- Q$S[,,m]
}

dtype <- match_dtype("mvnorm")
params <- list(M=M, mu=mu_list, Sigma=Sigma_list)

grid <- kohonen::somgrid(xdim=2, ydim=2, topo="h")
adjmatrix <- as.matrix(dist(grid$pts)) <= 1.5
diag(adjmatrix) = 0
nhbrdist <- dist_from_adj(adjmatrix)

mvn <- mlsom.mvn(A$X, niter=1e+4, nhbrdist=nhbrdist, radii=2)

plot(A$X, col=A$id+1)
points(mvn$Mus, pch=16)
plot.ellipses(mvn$Mus, mvn$Sigmas, 1:M)

grid <- kohonen::somgrid(xdim=3, ydim=3, topo="h")
adjmatrix <- as.matrix(dist(grid$pts)) <= 1.5
diag(adjmatrix) = 0

nhbrdist <- dist_from_adj(adjmatrix)

plot_grid(adjmatrix, grid, 9)

## ===========================================
##   method test
## ==========================================
n <- 1000
M <- 4
p <- 2
Q <- MixSim(BarOmega=0.05, K=M, p=p)
A <- simdataset(n, Q$Pi, Q$Mu, Q$S)

mu_list <- list()
Sigma_list <- list()
for (m in 1:M) {
    mu_list[[m]] <- Q$Mu[m,]
    Sigma_list[[m]] <- Q$S[,,m]
}

dtype <- match_dtype("mvnorm")
params <- list(M=M, mu=mu_list, Sigma=Sigma_list)

llconst <- loglikelihood_const(A$X, dtype)
logliks <- loglikelihood(A$X, params, llconst, dtype)

classes <- apply(logliks, 1, which.max)

target_m <- 2

## method test1
params2 <- method_test1(params, target_m-1, p)

## method test2
classes2 <- method_test2(A$X, params, classes-1, target_m-1) + 1
### check
logliks_tmp <- loglikelihood(A$X, params2, llconst, dtype)
classes_tmp <- classes
classes_tmp[classes > target_m] <- classes_tmp[classes > target_m] - 1
classes_tmp[classes==target_m] <- apply(logliks_tmp[classes==target_m, ], 1, which.max)
table(classes_tmp, classes2)

## method test3
params_batch <- method_test3(A$X, params, classes-1)
### check
mus <- by(A$X, classes, colMeans)
Sigmas <- by(A$X, classes, function(x) (nrow(x)-1)/ nrow(x) * cov(x))

## eval_without
mrange <- c(0, params$M-1)
candidate_map <- eval_without(target_m-1, A$X, mrange, params, classes-1, llconst, dtype)
classes_c <- candidate_map$classes + 1
cmdl <- classif_mdl(A$X, candidate_map, classes_c, dtype)
cmdl - candidate_map$value

## =========================================================
##    draft smlsom for mvn
## =========================================================

### data gen
n <- 1000
M <- 4
p <- 2
Q <- MixSim(BarOmega=0.05, K=M, p=p, int=c(-3, 3))
A <- simdataset(n, Q$Pi, Q$Mu, Q$S)
X <- A$X

plot(X, col=A$id+1)

### select model
dtype <- match_dtype("mvnorm")

### mlsom
beta <- 5
niter <- 1e+4
alpha <- c(0.05, 0.01)
radii <- c(2, -2)

Mhat <- 6
mu_list <- list()
Sigma_list <- list()
for (m in 1:Mhat) {
    mu_list[[m]] <- rnorm(p)
    Sigma_list[[m]] <- diag(p)
}

current_map <- list(M=Mhat, mu=mu_list, Sigma=Sigma_list)
grid <- kohonen::somgrid(xdim=3, ydim=2, topo="h")
adjmatrix <- as.matrix(dist(grid$pts)) <= 1.5
diag(adjmatrix) = 0

llconst <- loglikelihood_const(X, dtype)

### 1. do mlsom
nhbrdist <- dist_from_adj(adjmatrix)
current_map <- mlsom(X, current_map, dtype, niter, nhbrdist, alpha, radii)

#### plot
Mus <- list_to_mvn.Mus(current_map)
Sigmas <- list_to_mvn.Sigmas(current_map)

plot(X, col=A$id+1)
points(Mus, pch=16)
plot.ellipses(Mus, Sigmas, 1:Mhat)
draw.nl(Mus, adjmatrix)

### 2. link cutting
logliks <- loglikelihood(X, current_map, llconst, dtype)
classes <- apply(logliks, 1, which.max)

weight <- link_weight(classes, logliks, adjmatrix)
adjmatrix <- link_cutting(beta, weight, adjmatrix)

#### plot
plot(X, col=A$id+1)
points(Mus, pch=16)
plot.ellipses(Mus, Sigmas, 1:Mhat)
draw.nl(Mus, adjmatrix)

### 3. node deletion

#### 3-1. current mdl
selected_list <- which_delete(X, current_map, classes, llconst, dtype, Mhat)
selected_m <- selected_list$selected_m
selected_map <- selected_list$selected_map

#### 3-2. mdl without each node
if (selected_m != 0) {
    Mhat <- Mhat - 1
    current_map <- selected_map
    adjmatrix <- accommodate_delete(adjmatrix, selected_m)
}

#### plot
Mus <- list_to_mvn.Mus(current_map)
Sigmas <- list_to_mvn.Sigmas(current_map)

plot(X, col=A$id+1)
points(Mus, pch=16)
plot.ellipses(Mus, Sigmas, 1:Mhat)
draw.nl(Mus, adjmatrix)

### goto 1.

## =========================================================
##    smlsom function for mvn
## =========================================================

### experiment1.R
load("test.Rdata")

#### plot
Mhat <- ret1$M
Mus <- ret1$Mus
Sigmas <- ret1$Sigmas
adjmatrix <- ret1$adjmatrix

plot(X, col=A$id+1)
points(Mus, pch=16)
plot.ellipses(Mus, Sigmas, 1:Mhat)
draw.nl(Mus, adjmatrix)

Mus2 <- ret2$mu1[ret2$lives, ]
mu2 <- ret2$mu2[ret2$lives]
Sigmas2 <- array(0, dim=c(p,p,ret2$nd))
for (m in 1:ret2$nd)
    Sigmas2[,,m] <- mu2[[m]] - Mus2[m, ] %*% t(Mus2[m, ])

plot(X, col=A$id+1)
points(Mus2, pch=16)
plot.ellipses(Mus2, Sigmas2, 1:Mhat)
