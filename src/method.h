#ifndef INCLUDE_method_h_
#define INCLUDE_method_h_

#include "base.h"

model** model_drop_copy(model**, model_class, int, int, int);
IntegerVector redistribute(NumericMatrix, model**, IntegerVector, int, int, int);
void model_batch(NumericMatrix, model**, IntegerVector, int);
double cmdl(NumericMatrix, model**, IntegerVector, NumericVector, int);

#endif //INCLUDE_method_h_
