#include "structures.h"
#include "io.h"
#include "physicalFunctions.h"

bool initTimeSeries(time_series *ts, model_parameters *params) {
  ts->surfaceTempLoad = calloc(params->T, sizeof(real));
  if (!ts->surfaceTempLoad) {
    return false;
  }
  if (!readTable(ts->surfaceTempLoad, params->TEMPERATURE_FILE)) {
    return false;
  }

  ts->iceThicknessLoad = calloc(params->T, sizeof(real));
  if (!ts->iceThicknessLoad) {
    return false;
  }
  if (!readTable(ts->iceThicknessLoad, params->ICE_THICKNESS_FILE)) {
    return false;
  }

  ts->accLoad = calloc(params->T, sizeof(real));
  if (!ts->accLoad) {
    return false;
  }
  if (!readTable(ts->accLoad, params->ACCUMULATION_FILE)) {
    return false;
  }

  ts->age = calloc(params->Z, sizeof(real));
  if (!ts->age) {
    return false;
  }
  if (!readTable(ts->age, params->AGE_FILE)) {
    return false;
  }

  ts->borehole_temp = calloc(params->Z, sizeof(real));
  if (!ts->borehole_temp) {
    return false;
  }
  if (!readTable(ts->borehole_temp, params->BOREHOLE_TEMPERATURE_FILE)) {
    return false;
  }
  return true;
}

void deleteTimeSeries(time_series *ts) {
  free(ts->surfaceTempLoad);
  free(ts->iceThicknessLoad);
  free(ts->accLoad);
  free(ts->age);
  free(ts->borehole_temp);
  printf("[I] Time series correctly deleted\n");
}

bool initModelParameters(model_parameters *params, char *fileName) {
  params->strings = (char ***)malloc(DATA_ENUM_SIZE * sizeof(char **));

  for (size_t i = 0; i < DATA_ENUM_SIZE; i++) {
    params->strings[i] = (char **)malloc(5 * sizeof(char *));
    for (size_t j = 0; j < 5; j++) {
      params->strings[i][j] = (char *)malloc(50 * sizeof(char));
    }
  }

  params->strings[SAVE_TYPE][ST_MATRIX] = "MATRIX";
  params->strings[SAVE_TYPE][ST_VECTOR] = "VECTOR";

  params->strings[SCHEME][SC_CN] = "CN";
  params->strings[SCHEME][SC_EXPL] = "EXPL";

  params->strings[THERMAL_ICE][TI_CP] = "CP";
  params->strings[THERMAL_ICE][TI_GO] = "GO";

  params->strings[THERMAL_FIRN][TF_SC] = "SC";
  params->strings[THERMAL_FIRN][TF_CP] = "CP";
  params->strings[THERMAL_FIRN][TF_CP_LIN] = "CP_LIN";
  params->strings[THERMAL_FIRN][TF_SC_LIN] = "SC_LIN";

  params->strings[HEAT_CAPACITY][CP_CP] = "CP";
  params->strings[HEAT_CAPACITY][CP_AL] = "AL";

  params->strings[RHO_FIRN][RHO_FIRN] = "HL";
  params->strings[RHO_FIRN][RF_CONST] = "CONST";

  params->strings[VERTICAL_PROFILE][VP_FI] = "FI";
  params->strings[VERTICAL_PROFILE][VP_PA] = "PA";

  params->strings[INTERNAL_ENERGY][IE_ON] = "ON";
  params->strings[INTERNAL_ENERGY][IE_OFF] = "OFF";

  params->strings[MELTING][SC_EXPL] = "FREE_MELT";
  params->strings[MELTING][SC_EXPL] = "FREEZING_NO_ICE";
  params->strings[MELTING][SC_EXPL] = "FREEZING";

  params->Z = -1;
  params->T = -1;
  params->S = -1;
  params->SAVE_TYPE = ST_UNSET;
  params->SCHEME = SC_UNSET;
  params->RHO_SNOW = -1;
  params->THERMAL_ICE = TI_UNSET;
  params->THERMAL_FIRN = TF_UNSET;
  params->RHO_FIRN = RF_UNSET;
  params->VERTICAL_PROFILE = VP_UNSET;
  params->INTERNAL_ENERGY = IE_UNSET;
  params->MELTING = ME_UNSET;
  params->OUTPUT_PATH = NULL;
  params->values.mw = NULL;
  params->values.mw_n = 0;
  params->values.QG = NULL;
  params->values.QG_n = 0;
  params->values.TCor = NULL;
  params->values.TCor_n = 0;
  params->values.TCor2 = NULL;
  params->values.TCor2_n = 0;
  params->values.PCor = NULL;
  params->values.PCor_n = 0;
  params->values.deltaH = NULL;
  params->values.deltaH_n = 0;
  params->values.len = NULL;
  params->values.len_n = 0;
  params->values.flat = NULL;
  params->values.flat_n = 0;
  params->values.tot = 0;

  return readInitFile(params, fileName);
}

void printModelParameters(model_parameters *params) {
  printf("\t--Model Parameters--\n");
  printf("\tZ\t\t:\t%d\n", params->Z);
  printf("\tT\t\t:\t%d\n", params->T);
  printf("\tS\t\t:\t%d\n", params->S);
  printf("\tTEMPERATURE_FILE\t:\t%s\n", params->TEMPERATURE_FILE);
  printf("\tACCUMULATION_FILE\t:\t%s\n", params->ACCUMULATION_FILE);
  printf("\tICE_THICKNESS_FILE\t:\t%s\n", params->ICE_THICKNESS_FILE);
  printf("\tAGE_FILE\t:\t%s\n", params->AGE_FILE);
  printf("\tBOREHOLE_TEMPERATURE_FILE\t:\t%s\n",
         params->BOREHOLE_TEMPERATURE_FILE);
  printf("\tOUTPUT_PATH\t:\t%s\n", params->OUTPUT_PATH);
  printf("\tSAVE_TYPE\t:\t%s\n", params->strings[SAVE_TYPE][params->SAVE_TYPE]);
  printf("\tSCHEME\t\t:\t%s\n", params->strings[SCHEME][params->SCHEME]);
  printf("\tRHOSNOW\t\t:\t%d\n", params->RHO_SNOW);
  printf("\tTHERMAL\t\t:\t%s\n",
         params->strings[THERMAL_ICE][params->THERMAL_ICE]);
  printf("\tFIRN\t\t:\t%s\n",
         params->strings[THERMAL_FIRN][params->THERMAL_FIRN]);
  printf("\tRHO\t\t:\t%s\n", params->strings[RHO_FIRN][params->RHO_FIRN]);
  printf("\tVERTICAL_PROFILE:\t%s\n",
         params->strings[VERTICAL_PROFILE][params->VERTICAL_PROFILE]);
  printf("\tINTERNAL_ENERGY\t:\t%s\n",
         params->strings[INTERNAL_ENERGY][params->INTERNAL_ENERGY]);
  printf("\tMELTING\t\t:\t%s\n", params->strings[MELTING][params->MELTING]);
  printf("\tMW\t\t:\t%.2f", params->values.mw[0]);
  for (size_t i = 1; i < params->values.mw_n; ++i) {
    printf(" %.2f", params->values.mw[i]);
  }
  printf("\n");
  printf("\tQG\t\t:\t%.2f", params->values.QG[0]);
  for (size_t i = 1; i < params->values.QG_n; ++i) {
    printf(" %.2f", params->values.QG[i]);
  }
  printf("\n");
  printf("\tTCOR\t\t:\t%.2f", params->values.TCor[0]);
  for (size_t i = 1; i < params->values.TCor_n; ++i) {
    printf(" %.2f", params->values.TCor[i]);
  }
  printf("\n");
  printf("\tTCOR2\t\t:\t%.2f", params->values.TCor2[0]);
  for (size_t i = 1; i < params->values.TCor2_n; ++i) {
    printf(" %.2f", params->values.TCor2[i]);
  }
  printf("\n");
  printf("\tPCOR\t\t:\t%.2f", params->values.PCor[0]);
  for (size_t i = 1; i < params->values.PCor_n; ++i) {
    printf(" %.2f", params->values.PCor[i]);
  }
  printf("\n");
  printf("\tDELTAH\t\t:\t%.2f", params->values.deltaH[0]);
  for (size_t i = 1; i < params->values.deltaH_n; ++i) {
    printf(" %.2f", params->values.deltaH[i]);
  }
  printf("\n");
  printf("\tLEN\t\t:\t%.2f", params->values.len[0]);
  for (size_t i = 1; i < params->values.len_n; ++i) {
    printf(" %.2f", params->values.len[i]);
  }
  printf("\n");
  printf("\tFLAT\t\t:\t%.2f", params->values.flat[0]);
  for (size_t i = 1; i < params->values.flat_n; ++i) {
    printf(" %.2f", params->values.flat[i]);
  }
  printf("\n");
}

void deleteModelParameters(model_parameters *params) {
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

bool initModelData(model_data *data, const model_parameters *const params,
                   size_t mwL, size_t QGL, size_t TcorL, size_t TcorL2,
                   size_t PcorL, size_t deltaHL, size_t lenL, size_t flatL) {
  data->temperature = calloc(params->Z, sizeof(real *));
  if (!data->temperature) {
    return false;
  }
  for (int li = 0; li < params->Z; li++) {
    data->temperature[li] = calloc(params->T, sizeof(real));
    if (!data->temperature[li]) {
      return false;
    }
  }

  data->density = calloc(params->Z, sizeof(real *));
  if (!data->density) {
    return false;
  }
  for (int li = 0; li < params->Z; li++) {
    data->density[li] = calloc(params->T, sizeof(real));
    if (!data->density[li]) {
      return false;
    }
  }

  data->surfaceTemp = calloc(params->T, sizeof(real *));
  if (!data->surfaceTemp) {
    return false;
  }

  data->iceThickness = calloc(params->T, sizeof(real *));
  if (!data->iceThickness) {
    return false;
  }

  data->acc = calloc(params->T, sizeof(real *));
  if (!data->acc) {
    return false;
  }

  data->acc2 = calloc(params->T, sizeof(real *));
  if (!data->acc2) {
    return false;
  }

  data->melt = calloc(params->T, sizeof(real *));
  if (!data->melt) {
    return false;
  }

  data->freeze = calloc(params->T, sizeof(real *));
  if (!data->freeze) {
    return false;
  }

  data->tnew = calloc(params->Z, sizeof(real *));
  if (!data->tnew) {
    return false;
  }

  data->mw = params->values.mw[mwL];
  data->QG = params->values.QG[QGL];
  data->tCor = params->values.TCor[TcorL];
  data->tCor2 = params->values.TCor2[TcorL2];
  data->pCor = params->values.PCor[PcorL];
  data->deltaH = params->values.deltaH[deltaHL];
  data->len = params->values.len[lenL];
  data->flat = params->values.flat[flatL];

  return true;
}

void deleteModelData(model_data *data, const model_parameters *const params) {
  for (int li = 0; li < params->Z; li++) {
    double *currentIntPtr = data->temperature[li];
    free(currentIntPtr);
  }
  for (int li = 0; li < params->Z; li++) {
    double *currentIntPtr = data->density[li];
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

bool initModelFunctions(model_functions *functions,
                        const model_parameters *const params) {

  if (params->RHO_FIRN == RF_HL) {
    functions->setRho = &setRho_HL;
  } else if (params->RHO_FIRN == RF_CONST) {
    functions->setRho = &setRho_CONST;
  } else {
    functions->setRho = NULL;
    return false;
  }

  if (params->THERMAL_ICE == TI_CP) {
    functions->setThermalIce = &setThermalIce_CP;
  } else if (params->THERMAL_ICE == TI_GO) {
    functions->setThermalIce = &setThermalIce_GO;
  } else {
    functions->setThermalIce = NULL;
    return false;
  }

  if (params->THERMAL_FIRN == TF_CP) {
    functions->setThermalFirn = &setThermalFirn_CP;
  } else if (params->THERMAL_FIRN == TF_SC) {
    functions->setThermalFirn = &setThermalFirn_SC;
  } else if (params->THERMAL_FIRN == TF_CP_LIN) {
    functions->setThermalFirn = &setThermalFirn_CP_LIN;
  } else if (params->THERMAL_FIRN == TF_SC_LIN) {
    functions->setThermalFirn = &setThermalFirn_SC_LIN;
  } else {
    functions->setThermalFirn = NULL;
    return false;
  }
  if (params->HEAT_CAPACITY == CP_CP) {
    functions->setHeatCapacity = &setHeatCapacity_CP;
  } else if (params->HEAT_CAPACITY == CP_AL) {
    functions->setHeatCapacity = &setHeatCapacity_AL;
  } else {
    functions->setHeatCapacity = NULL;
    return false;
  }

  if (params->MELTING == ME_FREE_MELT) {
    functions->computeMelt = &computeMelt_FREE_MELT;
  } else if (params->MELTING == ME_FREEZING_NO_ICE) {
    functions->computeMelt = &computeMelt_FREEZING_NO_ICE;
  } else if (params->MELTING == ME_FREEZING) {
    functions->computeMelt = &computeMelt_FREEZING;
  } else {
    functions->computeMelt = NULL;
    return false;
  }

  if (params->VERTICAL_PROFILE == VP_FI) {
    functions->wDef = &wDef_FI;
  } else if (params->VERTICAL_PROFILE == VP_PA) {
    functions->wDef = &wDef_PA;
  } else {
    functions->computeMelt = NULL;
    return false;
  }

  if (params->INTERNAL_ENERGY == IE_ON) {
    functions->setSe = &setSe_ON;
  } else if (params->INTERNAL_ENERGY == IE_OFF) {
    functions->setSe = &setSe_OFF;
  } else {
    functions->computeMelt = NULL;
    return false;
  }

  return true;
}
