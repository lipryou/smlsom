#include "model.h"
#include "model_io.h"
#include "likelihood.h"

// [[Rcpp::export]]
NumericMatrix loglikelihood(NumericMatrix data, List parameters, NumericVector llconst, int dtype) {
  int n = data.nrow();
  int M = parameters["M"];
  NumericMatrix lliks(n, M);

  model_class mclass = match_model_class(dtype);

  model** mlist = model_import(parameters, mclass);

  //NumericVector llconst = loglikelihood_const(data, mclass);

  for (int i = 0; i < n; i++) {
    for (int m = 0; m < M; m++){
      lliks(i, m) = mlist[m]->loglikelihood(data(i, _));
      lliks(i, m) += llconst[i];
    }
  }

  return lliks;
}

// [[Rcpp::export]]
NumericVector classification_loglikelihood(NumericMatrix data, List parameters,
                                           IntegerVector classes, NumericVector llconst, int dtype) {
  int n = data.nrow();
  NumericVector clliks(n);

  model_class mclass = match_model_class(dtype);

  model** mlist = model_import(parameters, mclass);

  //NumericVector llconst = loglikelihood_const(data, mclass);

  for (int i = 0; i < n; i++) {
    clliks[i] = mlist[classes[i]]->loglikelihood(data(i, _));
    clliks[i] += llconst[i];
  }

  return clliks;
}

// [[Rcpp::export]]
NumericVector loglikelihood_const(NumericMatrix data, int dtype) {
  int n = data.nrow();
  int p = data.ncol();
  NumericVector llconst(n);

  model_class mclass = match_model_class(dtype);

  switch (mclass) {

  case model_class::mvnorm:
    for (int i = 0; i < n; i++)
      llconst[i] = (-1) * p / 2.0 * std::log(2*M_PI);
    return llconst;

  case model_class::multinom:
    {
      double constr = 0.0;

      for (int i = 0; i < n; i++) {
        constr = 0.0;
        for (int j = 0; j < p; j++) constr += std::lgamma(data(i, j) + 1);
        llconst[i] = std::lgamma(sum(data(i, _)) + 1) - constr;
      }
    }

    return llconst;

  default:
    stop("model_class error: procedure for such class not exist");
  }
}
