#include "model.h"
#include "model_io.h"
#include "method.h"
#include "clustering.h"

// [[Rcpp::export]]
List mlsom(NumericMatrix data, List parameters, int dtype, int niter,
           NumericMatrix nhbrdist, NumericVector alphas, NumericVector radii) {
  int n = data.nrow();
  int M = parameters["M"];

  NumericVector xi, nhbrdist_m;
  int i, nearest;
  double r, alpha;

  model_class mclass = match_model_class(dtype);

  model** mlist = model_import(parameters, mclass);

  RANDIN;

  for (int t = 0; t < niter; t++) {
    i = (int)(n * UNIF);

    xi = data(i, _);

    nearest = find_nearest(xi, mlist, 0, M);

    r = radius(t, niter, radii);
    alpha = learning_rate(t, niter, alphas);

    nhbrdist_m = nhbrdist(nearest, _);

    for (int m = 0; m < M; m++) {
      if (nhbrdist_m[m] > r)
        continue;
      mlist[m]->mom_by_sample(xi, alpha);
    }
  }

  RANDOUT;

  return model_export(mlist, M, mclass);
}

// [[Rcpp::export]]
List mlsom_clf(NumericMatrix X, IntegerVector y,
               IntegerVector nsubc, IntegerVector cumnsubc, List parameters,
               int dtype, int niter, NumericMatrix nhbrdist, NumericVector alphas, NumericVector radii) {
  int n = X.nrow();
  int M = parameters["M"];

  NumericVector xi, nhbrdist_m;
  int i, nearest, yi;
  double r, alpha;

  model_class mclass = match_model_class(dtype);

  model** mlist = model_import(parameters, mclass);

  RANDIN;

  for (int t = 0; t < niter; t++) {
    i = (int)(n * UNIF);

    xi = X(i, _);
    yi = y[i];

    nearest = find_nearest(xi, mlist, cumnsubc[yi], nsubc[yi]+cumnsubc[yi]);

    r = radius(t, niter, radii);
    alpha = learning_rate(t, niter, alphas);

    nhbrdist_m = nhbrdist(nearest, _);

    for (int m = 0; m < M; m++) {
      if (nhbrdist_m[m] > r)
        continue;
      mlist[m]->mom_by_sample(xi, alpha);
    }
  }

  RANDOUT;

  return model_export(mlist, M, mclass);
}

// [[Rcpp::export]]
List onebatch(NumericMatrix data, List parameters, IntegerVector classes, int dtype) {
  int M = parameters["M"];

  model_class mclass = match_model_class(dtype);

  model** mlist = model_import(parameters, mclass);

  model_batch(data, mlist, classes, M);

  return model_export(mlist, M, mclass);
}

// [[Rcpp::export]]
IntegerVector classifsubc_within_class(NumericMatrix X, IntegerVector y,
                                       IntegerVector nsubc, IntegerVector cumnsubc, List parameters,
                                       int dtype) {
  int n = X.nrow();
  int yi;
  IntegerVector classes(n);

  model_class mclass = match_model_class(dtype);

  model** mlist = model_import(parameters, mclass);

  for (int i = 0; i < n; i++) {
    yi = y[i];
    classes[i] = find_nearest(X(i, _), mlist, cumnsubc[yi], nsubc[yi]+cumnsubc[yi]);
  }

  return classes + 1;
}

int find_nearest(NumericVector xi, model** mlist, int start_m, int end_m) {
  int nind = 0;
  int nearest = -1;
  double dm = DOUBLE_XMAX;
  double dist;

  // find 'winner' among  nodes
  for (int m = start_m; m < end_m; m++) {
    dist = (-1) * mlist[m]->loglikelihood(xi);
    if (dist <= dm * (1 + EPS)) {
      if (dist < dm * (1 - EPS)) {
        nind = 0;
        nearest = m;
      } else {
        if(++nind * UNIF < 1.0) nearest = m;
      }
      dm = dist;
    }
  }
  if (nearest < 0)
    stop("No nearest neighbour found...");

  return nearest;
}

double radius(int t, int niter, NumericVector radii) {
  double threshold = radii[0] - (radii[0] - radii[1]) * (double)t/(double)niter;
  if (threshold < 1.0) threshold = 0.5;

  return threshold;
}

double learning_rate(int t, int niter, NumericVector alphas) {
  double alpha = alphas[0] - (alphas[0] - alphas[1]) * (double)t/(double)niter;

  return alpha;
}
