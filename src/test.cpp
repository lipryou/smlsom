#include "model.h"
#include "model_io.h"
#include "clustering.h"
#include "likelihood.h"
#include "method.h"

// [[Rcpp::export]]
NumericMatrix model_test1(NumericMatrix A) {
  /*R--------------------------------
    A <- matrix(c(3, 2, 2, 5), ncol=2)
    solve(chol(A))
    ---------------------------------
  */

  return choldc_inv(A);
}

// [[Rcpp::export]]
double model_test2(NumericVector x, NumericVector mu, NumericMatrix Sigma) {
  /*R--------------------------------
    library(mvtnorm)

    mu <- c(1, 2)
    S <- matrix(c(3, 1, 1, 5), ncol=2)

    x <- c(0, 0)
    dmvnorm(x, mu, S, log=T) + log(2*pi)
    ---------------------------------
  */

  gaussian g (mu, Sigma);

  return g.loglikelihood(x);
}

// [[Rcpp::export]]
double model_test3(NumericVector x,
                    NumericVector mu1, NumericMatrix Sigma1,
                    NumericVector mu2, NumericMatrix Sigma2) {
  /*R--------------------------------
    library(mvtnorm)

    mu <- c(1, 2)
    S <- matrix(c(3, 1, 1, 5), ncol=2)

    mu2 <- c(0.5, 0.5)
    S2 <- matrix(c(2, .5, .5, 3), ncol=2)

    x <- c(0, 0)
    dmvnorm(x, mu1, S1, log=T) + log(2*pi)
    dmvnorm(x, mu2, S2, log=T) + log(2*pi)
    ---------------------------------
  */

  gaussian g (mu1, Sigma1);

  List params = List::create(Named("mu") = mu2, _["Sigma"] = Sigma2);

  g.update(params);

  return g.loglikelihood(x);
}

// [[Rcpp::export]]
List model_test4(double alpha, NumericVector x,
                      NumericVector mu, NumericMatrix Sigma) {
  /*R--------------------------------
    mu <- c(1, 2)
    S <- matrix(c(3, 1, 1, 5), ncol=2)

    x <- c(0, 0)
    alpha <- 0.05

    tmp <- x - mu
    mu + alpha * tmp
    S + alpha * ((1-alpha) * tmp %*% t(tmp) - S)
    ---------------------------------
  */

  gaussian g (mu, Sigma);

  g.mom_by_sample(x, alpha);

  return List::create(g.mu, g.Sigma);
}

// [[Rcpp::export]]
void model_io_test1(List parameters) {
  /*R--------------------------------
    M <- 2
    mu1 <- c(0, 0)
    mu2 <- c(1, 1)
    S1 <- matrix(c(3, 1, 1, 5), ncol=2)
    S2 <- matrix(c(3, 1, 1, 5), ncol=2)

    params <- list(M=M, mu=list(mu1, mu2), Sigma=list(S1, S2))

    model_test1(params)
    ---------------------------------
  */

  int M = parameters["M"];

  model** mlist = model_import_gaussian(parameters);

  NumericVector x(2);

  for (int k = 0; k < M; k++)
    Rprintf("%d: %.5f\n",k, mlist[k]->loglikelihood(x));
}

// [[Rcpp::export]]
void model_io_test2(List parameters) {
  /*R--------------------------------
    M <- 2
    theta1 <- c(0.1, 0.9)
    theta2 <- c(0.9, 0.1)

    params <- list(M=M, theta=list(theta1, theta2))

    model_test2(params)
    ---------------------------------
  */
  int M = parameters["M"];

  model** mlist = model_import_multinomial(parameters);

  NumericVector x = {2, 3};

  for (int k = 0; k < M; k++)
    Rprintf("%d: %.5f\n",k, mlist[k]->loglikelihood(x));
}

// [[Rcpp::export]]
List model_io_test3(List parameters, int dtype) {
  /*R--------------------------------
    M <- 2
    mu1 <- c(0, 0)
    mu2 <- c(1, 1)
    S1 <- matrix(c(3, 1, 1, 5), ncol=2)
    S2 <- matrix(c(1, 0, 0, 1), ncol=2)

    params_i <- list(M=M, mu=list(mu1, mu2), Sigma=list(S1, S2))

    dtype <- match_dtype("mvnorm")
    params_o <- model_io_test3(params_i, dtype)

    theta1 <- c(0.9, 0.1)
    theta2 <- c(0.1, 0.9)

    params_i <- list(M=M, theta=list(theta1, theta2))

    dtype <- match_dtype("multinom")
    params_o <- model_io_test3(params_i, dtype)
    ---------------------------------
  */
  int M = parameters["M"];
  model_class mclass = match_model_class(dtype);
  model** mlist = model_import(parameters, mclass);

  return model_export(mlist, M, mclass);
}

// [[Rcpp::export]]
List model_io_test4(List parameters, int p) {
  /*R--------------------------------
    M <- 2
    mu1 <- c(0, 0)
    mu2 <- c(1, 1)
    S1 <- matrix(c(3, 1, 1, 5), ncol=2)
    S2 <- matrix(c(1, 0, 0, 1), ncol=2)

    params_i <- list(M=M, mu=list(mu1, mu2), Sigma=list(S1, S2))

    params_j <- model_io_test4(params_i, 2)
    ---------------------------------
  */

  int M = parameters["M"];

  model** mlist = model_import(parameters, model_class::mvnorm);
  model** mlist_c = model_dummy(model_class::mvnorm, M, p);

  // model copy
  for (int m = 0; m < M; m++) {
    List params = mlist[m]->get_parameters();
    mlist_c[m]->update(params);
  }

  return model_export(mlist_c, M, model_class::mvnorm);
}


// [[Rcpp::export]]
IntegerVector clustering_test1(NumericMatrix data, List parameters) {
  /*R--------------------------------
    library(MixSim)

    M <- 5
    Q <- MixSim(BarOmega=0.1, K=M, p=10)
    A <- simdataset(1000, Q$Pi, Q$Mu, Q$S)

    mu_list <- list()
    Sigma_list <- list()
    for (m in 1:M) {
      mu_list[[m]] <- Q$Mu[m, ]
      Sigma_list[[m]] <- Q$S[,,m]
    }

    params <- list(M=M, mu=mu_list, Sigma=Sigma_list)
    classes <- mlsom_test1(A$X, params)

    logliks <- sapply(1:M, function(m) dmvnorm(A$X, Q$Mu[m, ], Q$S[,,m], log=T))
    cls <- apply(logliks, 1, which.max)

    table(classes, cls)
    ---------------------------------
  */
  int n = data.nrow();
  int M = parameters["M"];

  NumericVector xi;

  IntegerVector classes(n);

  model** mlist = model_import_gaussian(parameters);

  for (int i = 0; i < n; i++) {
    xi = data(i, _);

    classes[i] = find_nearest(xi, mlist, 0, M);
  }

  return classes;
}

// [[Rcpp::export]]
double clustering_test2(int t, int niter, NumericVector radii) {
  /*R--------------------------------------
    radii <- c(2, -2)
    niter <- 1000
    rs <- sapply(1:niter, function(t) clustering_test2(t, niter, radii))

    plot(1:niter, rs, type="l")
  */
  return radius(t, niter, radii);
}

// [[Rcpp::export]]
double clustering_test3(int t, int niter, NumericVector alphas) {
  /*R--------------------------------------
    alpha <- c(0.05, 0.01)
    niter <- 1000
    alphas <- sapply(1:niter, function(t) clustering_test3(t, niter, alpha))

    plot(1:niter, alphas, type="l")
  */
  return learning_rate(t, niter, alphas);
}

// [[Rcpp::export]]
void list_test(NumericMatrix X, IntegerVector classes, int M) {
  int n = X.nrow();
  int p = X.ncol();

  std::vector<std::list<int>*> Sm(M);
  for (int m = 0; m < M; m++)
    Sm[m] = new std::list<int>();

  for (int i = 0; i < n; i++)
    Sm[classes[i]]->push_back(i);

  for (int m = 0; m < M; m++)
    Rprintf("%dth length: %d\n", m, Sm[m]->size());

  List list;
  for (int m = 0; m < M; m++) {
    int im = 0;
    NumericMatrix Xm(Sm[m]->size(), p);
    Rprintf("nrow = %d\n", Xm.nrow());
    for (auto itr = Sm[m]->begin(); itr != Sm[m]->end(); ++itr) {
      Xm(im, _) = X(*itr, _);
      im++;
    }
    for (int j = 0; j < p; j++)
      Rprintf(" %.3f", (double)sum(Xm(_, j)) );
    Rprintf("\n");
  }
}

// [[Rcpp::export]]
List method_test1(List parameters, int target_m, int p) {
  int M = parameters["M"];
  int Mhat = M-1;

  model** mlist = model_import(parameters, model_class::mvnorm);
  model** mlist_c = model_drop_copy(mlist, model_class::mvnorm, M, target_m, p);

  return model_export(mlist_c, Mhat, model_class::mvnorm);
}

// [[Rcpp::export]]
IntegerVector method_test2(NumericMatrix data, List parameters,
                           IntegerVector classes, int target_m) {
  int M = parameters["M"];
  int Mhat = M-1;
  int p = data.ncol();

  model** mlist = model_import(parameters, model_class::mvnorm);
  model** mlist_c = model_drop_copy(mlist, model_class::mvnorm, M, target_m, p);

  return redistribute(data, mlist_c, classes, target_m, 0, Mhat);
}

// [[Rcpp::export]]
List method_test3(NumericMatrix data, List parameters, IntegerVector classes) {
  int M = parameters["M"];

  model** mlist = model_import(parameters, model_class::mvnorm);

  model_batch(data, mlist, classes, M);

  return model_export(mlist, M, model_class::mvnorm);
}
