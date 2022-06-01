#ifndef INCLUDE_base_h_
#define INCLUDE_base_h_

#include <cmath>
#include <float.h>
#include <Rcpp.h>

using namespace Rcpp;

#define RANDIN  GetRNGstate()
#define RANDOUT PutRNGstate()
#define UNIF unif_rand() //一様乱数

#define EPS 1e-4

enum class model_class
{
  mvnorm,
  multinom,
  norms
};

model_class match_model_class(int);

#endif //INCLUDE_base_h_
