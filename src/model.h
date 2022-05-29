#ifndef INCLUDE_model_h_
#define INCLUDE_model_h_

#include "base.h"

NumericMatrix choldc_inv(NumericMatrix);

class model
{
public:
  double df;

  model() {};

  virtual double loglikelihood(NumericVector) = 0;

  virtual void update(List) = 0;

  virtual void mom_by_sample(NumericVector, double) = 0;

  virtual void batch(NumericMatrix) = 0;

  virtual List get_parameters(void) = 0;

};

class gaussian : public model
{
public:
  NumericVector mu;
  NumericMatrix Sigma;
  NumericMatrix invL;

  gaussian(NumericVector m, NumericMatrix s) {
    int p = m.length();

    mu = m;
    Sigma = s;
    invL = choldc_inv(s);

    df = p + 0.5 * p * (p + 1);
  }

  double loglikelihood(NumericVector);

  void update(List);

  void mom_by_sample(NumericVector, double);

  void batch(NumericMatrix);

  List get_parameters(void);
};

class multinomial : public model
{
public:
  NumericVector theta;

  multinomial(NumericVector p) {
    theta = p;
    df = p.length() - 1;
  }

  double loglikelihood(NumericVector);

  void update(List);

  void mom_by_sample(NumericVector, double);

  void batch(NumericMatrix);

  List get_parameters(void);
};

#endif //INCLUDE_model_h_
