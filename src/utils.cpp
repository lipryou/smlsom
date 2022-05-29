#include "base.h"

model_class match_model_class(int dtype) {
  switch (dtype) {
  case 0:
    return model_class::mvnorm;
  case 1:
    return model_class::multinom;
  default:
    stop("Such `dtype` is not exist.");
  }
}
