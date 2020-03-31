#ifndef RUNMODEL_H
#define RUNMODEL_H

#include "io.h"
#include "define.h"
#include "solver.h"
#include "structures.h"

void runModel(model_data *data,const model_parameters * const params,
              const time_series * const time_series);

#endif  /* !RUNMODEL_H */
