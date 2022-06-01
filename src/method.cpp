#include "model.h"
#include "model_io.h"
#include "clustering.h"
#include "likelihood.h"
#include "method.h"

// [[Rcpp::export]]
List eval_without(int target_m, NumericMatrix data,
                  IntegerVector search_mrange, List parameters,
                  IntegerVector classes, NumericVector llconst, int dtype) {
  int M = parameters["M"];
  int Mhat = M-1;
  int p = data.ncol();
  int start_m = search_mrange[0];
  int end_m = search_mrange[1];

  model_class mclass = match_model_class(dtype);

  model** mlist = model_import(parameters, mclass);

  // drop
  model** mlist_c = model_drop_copy(mlist, mclass, M, target_m, p);

  // redistribute
  IntegerVector classes_c = redistribute(data, mlist_c, classes, target_m, start_m, end_m);

  // update
  model_batch(data, mlist_c, classes_c, Mhat);

  // evaluate
  double value = cmdl(data, mlist_c, classes_c, llconst, Mhat);

  // prepare return
  List candidate = model_export(mlist_c, Mhat, mclass);
  candidate.push_back(classes_c, "classes");
  candidate.push_back(value, "value");

  return candidate;
}

double cmdl(NumericMatrix data, model** mlist,
            IntegerVector classes, NumericVector llconst, int M) {
  double value = 0.0;
  double df = 0.0;
  int n = data.nrow();

  for (int i = 0; i < n; i++)
    value += mlist[classes[i]]->loglikelihood(data(i, _));
  value += sum(llconst);

  for (int m = 0; m < M; m++)
    df += mlist[m]->df;

  return (-value + 0.5*df*log(n) + n*log(M));
}

model** model_drop_copy(model** mlist, model_class mclass,
                        int M, int target_m, int p) {
  int Mhat = M-1;
  model** mlist_c = model_dummy(mclass, Mhat, p);

  // model copy
  for (int m = 0; m < M; m++) {
    if (m == target_m)
      continue;
    List params = mlist[m]->get_parameters();
    if (m < target_m)
      mlist_c[m]->update(params);
    if (m > target_m)
      mlist_c[m-1]->update(params);
  }

  return mlist_c;
}

IntegerVector redistribute(NumericMatrix data, model** mlist,
                           IntegerVector classes, int target_m, int start_m, int end_m) {
  int n = classes.length();
  IntegerVector classes_c(n);

  for (int i = 0; i < n; i++) {
    if (classes[i] < target_m)
      classes_c[i] = classes[i];
    if (classes[i] == target_m)
      classes_c[i] = find_nearest(data(i, _), mlist, start_m, end_m);
    if (classes[i] > target_m)
      classes_c[i] = classes[i] - 1;
  }

  return classes_c;
}

void model_batch(NumericMatrix data, model** mlist, IntegerVector classes, int M) {
  int n = data.nrow();
  int p = data.ncol();

  std::vector<std::list<int>*> Sm(M);

  for (int m = 0; m < M; m++)
    Sm[m] = new std::list<int>();

  for (int i = 0; i < n; i++)
    Sm[classes[i]]->push_back(i);

  for (int m = 0; m < M; m++) {
    int im = 0;
    NumericMatrix Xm(Sm[m]->size(), p);
    for (auto itr = Sm[m]->begin(); itr != Sm[m]->end(); ++itr) {
      Xm(im, _) = data(*itr, _);
      im++;
    }
    mlist[m]->batch(Xm);
  }

  for (int m = 0; m < M; m++)
    Sm[m]->clear();
}
