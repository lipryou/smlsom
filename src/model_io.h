#ifndef INCLUDE_model_io_h_
#define INCLUDE_model_io_h_

#include "base.h"

model** model_import(List, model_class);
List model_export(model**, int, model_class);
model** model_dummy(model_class, int, int);

model** model_import_gaussian(List);
model** model_import_multinomial(List);
model** model_import_norms(List);

List model_export_gaussian(model**, int);
List model_export_multinomial(model**, int);
List model_export_norms(model**, int);

model** model_dummy_gaussian(int, int);
model** model_dummy_multinomial(int, int);
model** model_dummy_norms(int, int);

#endif //INCLUDE_model_io_h_
