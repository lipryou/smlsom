#ifndef INCLUDE_likelihood_h_
#define INCLUDE_likelihood_h_

#include "base.h"

NumericMatrix loglikelihood(NumericMatrix, List, NumericVector, int);
NumericVector classification_loglikelihood(NumericMatrix, List, IntegerVector, NumericVector, int);
NumericVector loglikelihood_const(NumericMatrix, int);

#endif //INCLUDE_likelihood_h_
