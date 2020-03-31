#include "structures.h"
#include "io.h"
#include "physicalFunctions.h"

bool initTimeSeries(time_series *ts, model_parameters *params)
{
  ts->surfaceTempLoad=calloc(params->T1,sizeof(real));
  if(!ts->surfaceTempLoad)
  {
    return false;
  }
  if(!readTable(ts->surfaceTempLoad,params->TEMPERATURE_FILE))
  {
    return false;
  }

  ts->iceThicknessLoad=calloc(params->T1,sizeof(real));
  if(!ts->iceThicknessLoad)
  {
    return false;
  }
  if(!readTable(ts->iceThicknessLoad,params->ICE_THICKNESS_FILE))
  {
    return false;
  }

  ts->accLoad=calloc(params->T1,sizeof(real));
  if(!ts->accLoad)
  {
    return false;
  }
  if(!readTable(ts->accLoad,params->ACCUMULATION_FILE))
  {
    return false;
  }

  ts->age=calloc(params->Z1,sizeof(real));
  if(!ts->age)
  {
    return false;
  }
  if(!readTable(ts->age,params->AGE_FILE))
  {
    return false;
  }

  ts->borehole_temp=calloc(params->Z1,sizeof(real));
  if(!ts->borehole_temp)
  {
    return false;
  }
  if(!readTable(ts->borehole_temp,params->BOREHOLE_TEMPERATURE_FILE))
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
  params->SAVE_TYPE1=ST_UNSET;
  params->SCHEME1=SC_UNSET;
  params->rhoSnow1=-1;
  params->THERMAL1=TH_UNSET;
  params->FIRN1=FI_UNSET;
  params->RHO1=RHO_UNSET;
  params->VERTICAL_PROFILE1=VP_UNSET;
  params->INTERNAL_ENERGY1=IE_UNSET;
  params->MELTING1=ME_UNSET;
  params->OUTPUT_PATH=NULL;
  params->values.mw=NULL;
  params->values.mw_n=0;
  params->values.QG=NULL;
  params->values.QG_n=0;
  params->values.TCor=NULL;
  params->values.TCor_n=0;
  params->values.TCor2=NULL;
  params->values.TCor2_n=0;
  params->values.PCor=NULL;
  params->values.PCor_n=0;
  params->values.deltaH=NULL;
  params->values.deltaH_n=0;
  params->values.len=NULL;
  params->values.len_n=0;
  params->values.flat=NULL;
  params->values.flat_n=0;
  params->values.tot=0;

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
  printf("\tTEMPERATURE_FILE\t:\t%s\n", params->TEMPERATURE_FILE);
  printf("\tACCUMULATION_FILE\t:\t%s\n", params->ACCUMULATION_FILE);
  printf("\tICE_THICKNESS_FILE\t:\t%s\n", params->ICE_THICKNESS_FILE);
  printf("\tAGE_FILE\t:\t%s\n", params->AGE_FILE);
  printf("\tBOREHOLE_TEMPERATURE_FILE\t:\t%s\n", params->BOREHOLE_TEMPERATURE_FILE);
  printf("\tOUTPUT_PATH\t:\t%s\n", params->OUTPUT_PATH);
  printf("\tSAVE_TYPE\t:\t%s\n", SAVE_TYPE_STR[params->SAVE_TYPE1]);
  printf("\tSCHEME\t\t:\t%s\n", SCHEME_STR[params->SCHEME1]);
  printf("\tRHOSNOW\t\t:\t%d\n", params->rhoSnow1);
  printf("\tTHERMAL\t\t:\t%s\n", THERMAL_STR[params->THERMAL1]);
  printf("\tFIRN\t\t:\t%s\n", FIRN_STR[params->FIRN1]);
  printf("\tRHO\t\t:\t%s\n", RHO_TYPE_STR[params->RHO1]);
  printf("\tVERTICAL_PROFILE:\t%s\n", VERTICAL_PROFILE_STR[params->VERTICAL_PROFILE1]);
  printf("\tINTERNAL_ENERGY\t:\t%s\n", INTERNAL_ENERGY_STR[params->INTERNAL_ENERGY1]);
  printf("\tMELTING\t\t:\t%s\n", MELTING_STR[params->MELTING1]);
  printf("\tMW\t\t:\t%.2f", params->values.mw[0]);
  for(size_t i=1; i<params->values.mw_n;++i)
  {
    printf(" %.2f", params->values.mw[i]);
  }
  printf("\n");
  printf("\tQG\t\t:\t%.2f", params->values.QG[0]);
  for(size_t i=1; i<params->values.QG_n;++i)
  {
    printf(" %.2f", params->values.QG[i]);
  }
  printf("\n");
  printf("\tTCOR\t\t:\t%.2f", params->values.TCor[0]);
  for(size_t i=1; i<params->values.TCor_n;++i)
  {
    printf(" %.2f", params->values.TCor[i]);
  }
  printf("\n");
  printf("\tTCOR2\t\t:\t%.2f", params->values.TCor2[0]);
  for(size_t i=1; i<params->values.TCor2_n;++i)
  {
    printf(" %.2f", params->values.TCor2[i]);
  }
  printf("\n");
  printf("\tPCOR\t\t:\t%.2f", params->values.PCor[0]);
  for(size_t i=1; i<params->values.PCor_n;++i)
  {
    printf(" %.2f", params->values.PCor[i]);
  }
  printf("\n");
  printf("\tDELTAH\t\t:\t%.2f", params->values.deltaH[0]);
  for(size_t i=1; i<params->values.deltaH_n;++i)
  {
    printf(" %.2f", params->values.deltaH[i]);
  }
  printf("\n");
  printf("\tLEN\t\t:\t%.2f", params->values.len[0]);
  for(size_t i=1; i<params->values.len_n;++i)
  {
    printf(" %.2f", params->values.len[i]);
  }
  printf("\n");
  printf("\tFLAT\t\t:\t%.2f", params->values.flat[0]);
  for(size_t i=1; i<params->values.flat_n;++i)
  {
    printf(" %.2f", params->values.flat[i]);
  }
  printf("\n");


}


void deleteModelParameters(model_parameters *params)
{
  free(params->OUTPUT_PATH);
  free(params->values.mw);
  free(params->values.QG);
  free(params->values.TCor);
  free(params->values.TCor2);
  free(params->values.PCor);
  free(params->values.deltaH);
  free(params->values.len);
  free(params->values.flat);

  printf("[I] Parameters correctly deleted\n");

}


bool initModelData(model_data *data,const model_parameters * const params,size_t mwL,
                   size_t QGL,size_t TcorL, size_t TcorL2,size_t PcorL,size_t deltaHL,
                   size_t lenL,size_t flatL)
{
  data->temperature = calloc( params->Z1, sizeof(real*));
  if(!data->temperature)
  {
    return false;
  }
  for (int li = 0; li < params->Z1; li++)
  {
    data->temperature[li] = calloc(params->T1, sizeof(real));
    if(!data->temperature[li])
    {
      return false;
    }
  }

  data->density = calloc( params->Z1, sizeof(real*));
  if(!data->density)
  {
    return false;
  }
  for (int li = 0; li < params->Z1; li++)
  {
    data->density[li] = calloc(params->T1, sizeof(real));
    if(!data->density[li])
    {
      return false;
    }
  }

  data->surfaceTemp = calloc( params->T1, sizeof(real*));
  if(!data->surfaceTemp)
  {
    return false;
  }

  data->iceThickness = calloc( params->T1, sizeof(real*));
  if(!data->iceThickness)
  {
    return false;
  }

  data->acc = calloc( params->T1, sizeof(real*));
  if(!data->acc)
  {
    return false;
  }

  data->acc2 = calloc( params->T1, sizeof(real*));
  if(!data->acc2)
  {
    return false;
  }

  data->melt = calloc( params->T1, sizeof(real*));
  if(!data->melt)
  {
    return false;
  }

  data->freeze = calloc( params->T1, sizeof(real*));
  if(!data->freeze)
  {
    return false;
  }

  data->tnew = calloc( params->Z1, sizeof(real*));
  if(!data->tnew)
  {
    return false;
  }

  data->mw=params->values.mw[mwL];
  data->QG=params->values.QG[QGL];
  data->tCor=params->values.TCor[TcorL];
  data->tCor2=params->values.TCor2[TcorL2];
  data->pCor=params->values.PCor[PcorL];
  data->deltaH=params->values.deltaH[deltaHL];
  data->len=params->values.len[lenL];
  data->flat=params->values.flat[flatL];

  return true;
}

void deleteModelData(model_data *data,const model_parameters * const params)
{
  for (int li = 0; li < params->Z1; li++)
  {
    double* currentIntPtr = data->temperature[li];
    free(currentIntPtr);
  }
  for (int li = 0; li < params->Z1; li++)
  {
    double* currentIntPtr = data->density[li];
    free(currentIntPtr);
  }
  free(data->temperature);
  free(data->density);
  free(data->tnew);
  free(data->surfaceTemp);
  free(data->iceThickness);
  free(data->acc);
  free(data->acc2);
  free(data->melt);
  free(data->freeze);
}

bool initModelFunctions(model_functions *functions,const model_parameters * const params)
{
  if(params->RHO1==RHO_FIRN)
  {
    functions->setRho=&setRho_FIRN;
  }
  else if(params->RHO1==RHO_CONST)
  {
    functions->setRho=&setRho_CONST;
  }
  else
  {
    return false;
  }
  return true;
}

void deleteModelFunctions(model_functions *functions)
{
  free(functions->setRho);
  // free(functions->setHeatVar);
  // free(functions->computeMelt);
  // free(functions->wDef);
  // free(functions->setABW);
  // free(functions->setSe);
  // free(functions->getDwdz);
  // free(functions->getA);
  // free(functions->getDudz);

}
