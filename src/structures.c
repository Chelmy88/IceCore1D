#include "structures.h"
#include "io.h"

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

void deleteTimeSeries(time_series *ts)
{
  free(ts->surfaceTempLoad);
  free(ts->iceThicknessLoad);
  free(ts->accLoad);
  free(ts->age);
  free(ts->borehole_temp);
	printf("[I] Time series correctly deleted\n");
}

bool initModelParameters(model_parameters *params, char* fileName)
{
  params->Z1=-1;
  params->T1=-1;
  params->S1=-1;
  params->OUTPUT_PATH;
  params->SAVE_TYPE1=ST_UNSET;
  params->SCHEME1=SC_UNSET;
  params->rhoSnow1=-1;
  params->THERMAL1=TH_UNSET;
  params->FIRN1=FI_UNSET;
  params->RHO1=RHO_UNSET;
  params->VERTICAL1=VP_UNSET;
  params->INTERNAL_ENERGY1=IE_UNSET;
  params->MELTING1=ME_UNSET;
  return readInitFile(params, fileName);
}

char *SAVE_TYPE_STR[3] = {"ST_MATRIX", "ST_VECTOR", "ST_UNSET"};
char *SCHEME_STR[3] = {"SC_CN", "SC_EXPL", "SC_UNSET"};
char *THERMAL_STR[3] = {"TH_CP", "TH_GO", "TH_UNSET"};
char *FIRN_STR[4] = {"FI_SC", "FI_CP", "FI_FI", "FI_UNSET"};
char *RHO_TYPE_STR[3] = {"RHO_FIRN", "RHO_CONST", "RHO_UNSET"};
char *VERTICAL_PROFILE_STR[3] = {"VP_FI", "VP_PA", "VP_UNSET"};
char *INTERNAL_ENERGY_STR[3] = {"IE_ON", "IE_OFF", "IE_UNSET"};
char *MELTING_STR[4] = {"ME_FREE_MELT", "ME_FREEZING_NO_ICE", "ME_FREEZING", "ME_UNSET"};

void printModelParameters(model_parameters *params)
{
  printf("\t--Model Parameters--\n");
  printf("\tZ\t\t:\t%d\n", params->Z1);
  printf("\tT\t\t:\t%d\n", params->T1);
  printf("\tS\t\t:\t%d\n", params->S1);
  printf("\tOUTPUT_PATH\t:\t%s\n", params->OUTPUT_PATH);
  printf("\tSAVE_TYPE\t:\t%s\n", SAVE_TYPE_STR[params->SAVE_TYPE1]);
  printf("\tSCHEME\t\t:\t%s\n", SCHEME_STR[params->SCHEME1]);
  printf("\trhoSnow\t\t:\t%d\n", params->rhoSnow1);
  printf("\tTHERMAL\t\t:\t%s\n", THERMAL_STR[params->THERMAL1]);
  printf("\tFIRN1\t\t:\t%s\n", FIRN_STR[params->FIRN1]);
  printf("\tRHO1\t\t:\t%s\n", RHO_TYPE_STR[params->RHO1]);
  printf("\tVERTICAL1\t:\t%s\n", VERTICAL_PROFILE_STR[params->VERTICAL1]);
  printf("\tINTERNAL_ENERGY1:\t%s\n", INTERNAL_ENERGY_STR[params->INTERNAL_ENERGY1]);
  printf("\tMELTING1\t:\t%s\n", MELTING_STR[params->MELTING1]);
}


void deleteModelParameters(model_parameters *params)
{
  free(params->OUTPUT_PATH);
  printf("[I] Parameters correctly deleted\n");

}
