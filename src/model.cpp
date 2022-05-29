#include "model.h"

/*
   m : mean vector
   S : covariance matrix
     -> S = LL^t
   log f (x | p) = - p/2 log(2 * pi) + log(|S|^-1/2) - 1/2 (x-mu)^t S^-1 (x-mu)
                 = const + log(|L^-1|) - 1/2 || L^-1 (x-mu) ||^2
*/
double gaussian::loglikelihood(NumericVector x) {
  double tmp;
  double exp = 0.0, logdet = 0.0;
  int p = x.length();

  for (int i = 0; i < p; i++) {
    tmp = 0.0;
    for (int j = 0; j <= i; j++)
      tmp += invL(i, j) * (x[j] - mu[j]); //L^-1 (x-mu)
    exp += tmp * tmp;
    logdet += log(invL(i, i));
  }

  return (logdet - 0.5 * exp);
}

void gaussian::update(List parameters) {
  NumericVector m = parameters["mu"];
  NumericMatrix s = parameters["Sigma"];

  mu = m;
  Sigma = s;
  invL = choldc_inv(s);
}


void gaussian::mom_by_sample(NumericVector x, double alpha) {
  int p = x.length();
  double delta;

  NumericVector tmp;
  NumericVector new_mu;
  NumericMatrix new_Sigma (p, p);

  tmp = x - mu;

  new_mu = mu + alpha * tmp;

  for (int i = 0; i < p; i++) {
    for (int j = 0; j < p; j++) {
      delta = (1 - alpha) * tmp[i] * tmp[j] - Sigma(i, j);
      new_Sigma(i, j) = Sigma(i, j) + alpha * delta;
    }
  }

  List parameters = List::create(Named("mu") = new_mu, _["Sigma"] = new_Sigma);
  update(parameters);
}

void gaussian::batch(NumericMatrix X) {
  int n = X.nrow();
  int p = X.ncol();

  NumericVector new_mu(p);
  NumericMatrix new_Sigma(p, p);

  for (int j = 0; j < p; j++)
    new_mu[j] = sum(X(_, j)) / n;

  for (int j = 0; j < p; j++) {
    for (int l = 0; l < p; l++)
      new_Sigma(j, l) = 1.0 / n * sum((X(_, j) - new_mu[j]) * (X(_, l) - new_mu[l]));
  }

  List parameters = List::create(Named("mu") = new_mu, _["Sigma"] = new_Sigma);
  update(parameters);
}

List gaussian::get_parameters() {
  return List::create(Named("mu") = mu, _["Sigma"] = Sigma);
}

/*
   p : multinomial distribution parameter
   log f (x | p) = lgamma(sum_j(xj) + 1) - sum_j(lgamma(xj+1)) + sum_j (xj log(pj))
                 = const + sum_j (xj log(pj))
*/
double multinomial::loglikelihood(NumericVector x) {
  double loglik = 0.0;

  for (int j = 0; j < x.length(); j++)
    loglik += x[j] * log(theta[j]);

  return loglik;
}

void multinomial::update(List parameters) {
  NumericVector p = parameters["theta"];
  theta = p;
}

void multinomial::mom_by_sample(NumericVector x, double alpha) {
  int p = x.length();
  double delta;

  double sum = 0.0;
  NumericVector new_theta (p);

  for (int j = 0; j < p; j++) {
    sum += x[j];
  }
  if (sum == 0.0)
    return;
  for(int j = 0; j < p; j++) {
    delta = (x[j]) / sum - theta[j];
    new_theta[j] = theta[j] + alpha * delta;
  }

  List parameters = List::create(Named("theta") = new_theta);
  update(parameters);
}

void multinomial::batch(NumericMatrix X) {
  int p = X.ncol();
  NumericVector col_total(p);
  NumericVector new_theta(p);

  for (int j = 0; j < p; j++)
    col_total[j] = sum(X(_, j));

  double total = sum(col_total);
  for (int j = 0; j < p; j++)
    new_theta[j] = col_total[j] / total;

  List parameters = List::create(Named("theta") = new_theta);
  update(parameters);
}

List multinomial::get_parameters() {
  return List::create(Named("theta") = theta);
}

NumericMatrix choldc_inv(NumericMatrix A) {
  /* A should be n x n matrix */

  double sum;
  int i, j, k;
  int n = A.nrow();

  NumericMatrix L(n, n);
  NumericVector diagL(n);

  // cholesky decomposition
  for (i = 0; i < n; i++) {
    for (j = 0; j < n; j++) {
      sum = A(i, j);
      for (k = i-1; k >= 0; k--)
        sum -= L(i, k) * L(j, k);
      if (i == j) {
        if (sum <= 0) {
          stop("'A' must be a positive definite.");
        }
        diagL[i] = sqrt(sum);
      } else {
        L(j, i) = sum / diagL[i];
      }
    }
  }

  // inverse L
  for (i = 0; i < n; i++) {
    L(i, i) = 1.0 / diagL[i];
    for (j = i+1; j < n; j++) {
      sum = 0.0;
      for (k = i; k < j; k++)
        sum -= L(j, k) * L(k, i);
      L(j, i) = sum / diagL[j];
    }
  }

  return L;
}
