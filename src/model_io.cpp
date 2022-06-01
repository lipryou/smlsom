#include "model.h"
#include "model_io.h"

model** model_import(List parameters, model_class mclass) {
  switch (mclass) {

  case model_class::mvnorm:
    return model_import_gaussian(parameters);

  case model_class::multinom:
    return model_import_multinomial(parameters);

  case model_class::norms:
    return model_import_norms(parameters);

  default:
    stop("model_class error: procedure for such class not exist");
  }
}

List model_export(model** mlist, int M, model_class mclass) {
  switch (mclass) {

  case model_class::mvnorm:
    return model_export_gaussian(mlist, M);

  case model_class::multinom:
    return model_export_multinomial(mlist, M);

  case model_class::norms:
    return model_export_norms(mlist, M);

  default:
    stop("model_class error: procedure for such class not exist");
  }
}

model** model_dummy(model_class mclass, int M, int p) {
  switch (mclass) {

  case model_class::mvnorm:
    return model_dummy_gaussian(M, p);

  case model_class::multinom:
    return model_dummy_multinomial(M, p);

  case model_class::norms:
    return model_dummy_norms(M, p);

  default:
    stop("model_class error: procedure for such class not exist");
  }
}

// model import

model** model_import_gaussian(List parameters) {
  int M = parameters["M"];
  List mu_list = parameters["mu"];
  List Sigma_list = parameters["Sigma"];

  model** mlist = new model* [M];

  for (int m = 0; m < M; m++) {
    NumericVector mu = mu_list[m];
    NumericMatrix Sigma = Sigma_list[m];

    mlist[m] = new gaussian(mu, Sigma);
  }

  return mlist;
}

model** model_import_multinomial(List parameters) {
  int M = parameters["M"];
  List theta_list = parameters["theta"];

  model** mlist = new model* [M];

  for (int m = 0; m < M; m++) {
    NumericVector theta = theta_list[m];

    mlist[m] = new multinomial(theta);
  }

  return mlist;
}

model** model_import_norms(List parameters) {
  int M = parameters["M"];
  List mu_list = parameters["mu"];
  List sigma_list = parameters["Sigma"];

  model** mlist = new model* [M];

  for (int m = 0; m < M; m++) {
    NumericVector mu = mu_list[m];
    NumericVector sigma = sigma_list[m];

    mlist[m] = new norms(mu, sigma);
  }

  return mlist;
}


// model export

List model_export_gaussian(model** mlist, int M) {
  List tmp_list;
  List mu_list = List::create();
  List Sigma_list = List::create();

  for (int m = 0; m < M; m++) {
    tmp_list = mlist[m]->get_parameters();
    mu_list.push_back(tmp_list["mu"]);
    Sigma_list.push_back(tmp_list["Sigma"]);
  }

  return List::create(Named("M")=M, _["mu"]=mu_list, _["Sigma"]=Sigma_list);
}

List model_export_multinomial(model** mlist, int M) {
  List tmp_list;
  List theta_list = List::create();

  for (int m = 0; m < M; m++) {
    tmp_list = mlist[m]->get_parameters();
    theta_list.push_back(tmp_list["theta"]);
  }

  return List::create(Named("M")=M, _["theta"]=theta_list);
}

List model_export_norms(model** mlist, int M) {
  List tmp_list;
  List mu_list = List::create();
  List sigma_list = List::create();

  for (int m = 0; m < M; m++) {
    tmp_list = mlist[m]->get_parameters();
    mu_list.push_back(tmp_list["mu"]);
    sigma_list.push_back(tmp_list["Sigma"]);
  }

  return List::create(Named("M")=M, _["mu"]=mu_list, _["Sigma"]=sigma_list);
}

// model dummy

model** model_dummy_gaussian(int M, int p) {
  List mu_list = List::create();
  List Sigma_list = List::create();

  for (int m = 0; m < M; m++) {
    NumericVector mu(p);
    NumericMatrix sigma(p, p);

    sigma.fill_diag(1);

    mu_list.push_back(mu);
    Sigma_list.push_back(sigma);
  }

  List params = List::create(Named("M")=M, _["mu"]=mu_list, _["Sigma"]=Sigma_list);

  return model_import_gaussian(params);
}

model** model_dummy_multinomial(int M, int p) {
  List theta_list = List::create();

  for (int m = 0; m < M; m++) {
    NumericVector theta(p);
    theta_list.push_back(theta);
  }

  List params = List::create(Named("M")=M, _["theta"]=theta_list);
  return model_import_multinomial(params);
}

model** model_dummy_norms(int M, int p) {
  List mu_list = List::create();
  List sigma_list = List::create();

  for (int m = 0; m < M; m++) {
    NumericVector mu(p);
    NumericVector sigma(p);

    for (int j=0; j < p; j++) sigma[j] = 1;

    mu_list.push_back(mu);
    sigma_list.push_back(sigma);
  }

  List params = List::create(Named("M")=M, _["mu"]=mu_list, _["Sigma"]=sigma_list);

  return model_import_norms(params);
}
