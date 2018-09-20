#ifndef STRUCTURES_H
#define STRUCTURES_H

#include <stdbool.h>
#include <stdlib.h>

typedef struct _time_series {
  double* surfaceTempLoad;
  double* iceThicknessLoad;
  double* accLoad;
  double* age;
  double* borehole_temp;
} time_series;

#endif  /* !STRUCTURES_H */

bool initTimeSeries(time_series *ts);
