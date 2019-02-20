#ifndef RUNMODEL_H
#define RUNMODEL_H

#include "io.h"
#include "define.h"
#include "solver.h"
#include "structures.h"

void runModel(double* tnew, double* surfaceTemp,double* iceThickness,double* acc,double* acc2,
              double* melt,double* freeze, double** temperature, double** density,
              const double* const surfaceTempLoad,const double* const iceThicknessLoad,
              const double* const accLoad,
              const double mw, const double QG, const double tCor, const double tCor2, const double pCor,
              const double deltaH, const double len,const double flat);

#endif  /* !RUNMODEL_H */
