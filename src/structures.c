#include "structures.h"
#include "define.h"

bool initTimeSeries(time_series *ts)
{
  ts->surfaceTempLoad=calloc(T,sizeof(double));
  if(!ts->surfaceTempLoad)
  {
    return false;
  }
  ts->iceThicknessLoad=calloc(T,sizeof(double));
  if(!ts->iceThicknessLoad)
  {
    return false;
  }
  ts->accLoad=calloc(T,sizeof(double));
  if(!ts->accLoad)
  {
    return false;
  }
  ts->age=calloc(Z,sizeof(double));
  if(!ts->age)
  {
    return false;
  }
  ts->borehole_temp=calloc(Z,sizeof(double));
  if(!ts->borehole_temp)
  {
    return false;
  }

  return true;
}
