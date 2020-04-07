#ifndef RUNMODEL_H
#define RUNMODEL_H

#include "define.h"
#include "io.h"
#include "solver.h"
#include "structures.h"

void runModel(model_data *data, const model_parameters *const params,
              const time_series *const time_series,
              const model_functions *const functions);

#endif /* !RUNMODEL_H */
