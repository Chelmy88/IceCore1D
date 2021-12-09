// Written by Adrien Michel
// adrien.michel@epfl.ch
// For the purpose of a Master thesis at the
// Climate and Environmental Group
// And
// Oeschger Center for Climate Change Research
// University of Bern
// February 2016
// Modified during PhD at EPFL, ENAC, CRYOS
// February 2019

// License: distribted under GNU GPL V3 license

// This code is developed to run on Linux machine, but should run on Windows or
// Mac OS
// The only parameters defined in main.c are the free parameters of the model,
// the correction to the boundary condition time series, and the initial
// temperature profile for spin up For any other modification see the main.h
// file

// Include the header file main.h, containing all the key functions and
// parameters.

#include "define.h"
#include "io.h"
#include "runModel.h"
#include "solver.h"
#include "structures.h"
#include <math.h>
#include <time.h>
#include <unistd.h>

bool mainLoop(model_parameters *params, time_series *ts,
              model_functions *functions);
double **computeAge(const model_data *const data,
                    const model_functions *const functions, size_t ageVerRes,
                    size_t ageHorRes, size_t ageCor);
void saveData(const model_data *const data,
              const model_parameters *const params, double **age,
              size_t ageVerRes, size_t ageHorRes);
int main() {
  // Set internal timer
  double begin = omp_get_wtime();

  // Read the model parameters
  model_parameters params;
  if (!initModelParameters(&params, "init.txt")) {
    deleteModelParameters(&params);
    printf("[E] Error reading the ini file. Exiting now\n");
    exit(EXIT_FAILURE);
  }

  // Creating the output directories.
  if (!createOutputDirs(&params)) {
    deleteModelParameters(&params);
    printf("[E] Error creating output directories. Exiting now\n");
    exit(EXIT_FAILURE);
  }

  // Tables to load data
  time_series ts;
  if (!initTimeSeries(&ts, &params)) {
    deleteTimeSeries(&ts);
    deleteModelParameters(&params);
    printf("[E] Error reading the time series files. Exiting now\n");
    exit(EXIT_FAILURE);
  }
  // Tables to load data
  model_functions functions;
  if (!initModelFunctions(&functions, &params)) {
    deleteTimeSeries(&ts);
    deleteModelParameters(&params);
    printf("[E] Error Initializing model functions. Exiting now\n");
    exit(EXIT_FAILURE);
  }

  char path[120];
  sprintf(path, "%s/init.txt", params.OUTPUT_PATH);
  copyFile("init.txt", path);

  size_t exit_status;
  if (mainLoop(&params, &ts, &functions)) {
    printf("[I] Simulation done in %f seconds\n", (double)(omp_get_wtime() - begin));
    exit_status = 0;
  } else {
    exit_status = EXIT_FAILURE;
  }

  deleteModelParameters(&params);
  deleteTimeSeries(&ts);

  return exit_status;
}

bool mainLoop(model_parameters *params, time_series *ts,
              model_functions *functions) {

  bool succes = true;
  int tot =
      params->values
          .tot; // Number of run, the size of the summary table should be bigger
  // Loops for the values of the free parameter
  int count = 0;
// Initialize the core splitting, should be placed in front of a loop having if
// possible a number of elements corresponding to a multiple of the core numbers
#pragma omp parallel for collapse(8) shared(succes)
  for (size_t mwL = 0; mwL < params->values.mw_n; mwL++) {
    for (size_t QGL = 0; QGL < params->values.QG_n; QGL++) {
      for (size_t TcorL = 0; TcorL < params->values.TCor_n; TcorL++) {
        for (size_t TcorL2 = 0; TcorL2 < params->values.TCor2_n; TcorL2++) {
          for (size_t PcorL = 0; PcorL < params->values.PCor_n; PcorL++) {
            for (size_t deltaHL = 0; deltaHL < params->values.deltaH_n;
                 deltaHL++) {
              for (size_t lenL = 0; lenL < params->values.len_n; lenL++) {
                for (size_t flatL = 0; flatL < params->values.flat_n; flatL++) {
                  if (!succes) {
                    continue;
                  }
                  model_data data;

                  if (!initModelData(&data, params, mwL, QGL, TcorL, TcorL2,
                                     PcorL, deltaHL, lenL, flatL)) {
                    deleteModelData(&data, params);
                    printf("[E] Error Initializing mode data. Exiting when all "
                           "threads are done\n");
                    succes = false;
                    continue;
                  }

                  runModel(&data, params, ts, functions);

                  size_t ageVerRes = 0;
                  size_t ageHorRes = 0;
                  size_t ageCor = 0;
                  if (params->SAVE_TYPE == ST_MATRIX) {
                    if (params->T == 10001) {
                      ageVerRes = (size_t)(params->Z / 5);
                      ageHorRes = (size_t)(params->T / 10);
                      ageCor = 0;
                    } else if (params->T == 40001) {
                      ageVerRes = (size_t)(params->Z / 5);
                      ageHorRes = (size_t)(params->T / 20);
                      ageCor = 20000;
                    }
                  } else if (params->SAVE_TYPE == ST_VECTOR) {
                    ageVerRes = (size_t)(params->Z / 5);
                    ageHorRes = 1;
                    ageCor = params->T - 11;
                  }

                  double **age = computeAge(&data, functions, ageVerRes,
                                            ageHorRes, ageCor);
                  saveData(&data, params, age, ageVerRes, ageHorRes);

#pragma omp atomic
                  count++;
#pragma omp flush(count)
                  printf("[I] Finished simulation: %d/%d\n", count, tot);
                  fflush(stdout);
                  deleteModelData(&data, params);
                }
              }
            }
          }
        }
      }
    }
  }
  return succes;
}

double **computeAge(const model_data *const data,
                    const model_functions *const functions, size_t ageVerRes,
                    size_t ageHorRes, size_t ageCor) {
  // // POST PROCESSING//
  // //Compute the difference with the age profile for 213 points and compute
  // the mean of the difference
  printf("[i] Computing age ...\n");
  fflush(stdout);
  double age_timer = omp_get_wtime();
  double **ageRel;

  ageRel = malloc((size_t)ageVerRes * sizeof(*ageRel));
  for (size_t i = 0; i < ageVerRes; i++) {
    ageRel[i] = malloc((ageHorRes + 1) * sizeof(ageRel));
  }
  for (size_t a = ageHorRes; a > 0; a--) {
    for (size_t h = 1; h < ageVerRes; h++) {
      float height = h * 5;
      int age = a * 10 + ageCor;
      while (height < data->iceThickness[age] && age > 0) {
        height += ((data->acc[age] * SEC_YEAR * 100 - data->melt[age] * 100 -
                    (data->iceThickness[age] - data->iceThickness[age - 1])) *
                       functions->wDef((double)height, data->iceThickness[age],
                                       data->mw) + data->melt[age] * 100);
        age--;
      }
      ageRel[h - 1][0] = h * 5;
      ageRel[h - 1][a] = (a * 10 + ageCor - age) * 100;
    }
  }
  printf("[i] Age profile computed in: %f secondes\n", (double)(omp_get_wtime() - age_timer));

  return (ageRel);
}

void saveData(const model_data *const data,
              const model_parameters *const params, double **age,
              size_t ageVerRes, size_t ageHorRes) {
  // //Compute the difference with the borehole  temperature profile below 600 m
  // deep (because upper part of the measurements are affected by seasonality)
  // double tempDiff=0;
  // double tnew2[Z]= {0};
  // for(size_t li=0; li<=(size_t)iceThickness[T-1]-7; li++)
  // {
  //   tnew2[li]=temperature[li][T-1];
  //   if(li<((size_t)iceThickness[T-1]-600))
  //   {
  //       tempDiff+=fabs(tnew2[li]-ts->borehole_temp[li]);
  //   }
  // }
  // tempDiff/=(iceThickness[T-1]);
  //
  // // Generate a file name with the free parameters

  char fileName[400] = "";
  char path[400] = "";
  sprintf(path,
          "%s/"
          "m_%.3f_Q_%.2f_Pcor_%.0f_Tcor_%.1f_Tcor2_%.1f_dH_%.0f_len_%.0f_flat_%"
          ".0f_%s_Rho_Snow_%d_Thermal_Ice_%s_Thermal_Firn_%s_Heat_Capacity_%s_"
          "Rho_Firn_%s_"
          "Internal_Energy_%s_Scheme_%s",
          params->OUTPUT_PATH, data->mw, data->QG * 1000, data->pCor,
          data->tCor, data->tCor2, data->deltaH, data->len, data->flat, "EDC",
          params->RHO_SNOW, params->strings[THERMAL_ICE][params->THERMAL_ICE],
          params->strings[THERMAL_FIRN][params->THERMAL_FIRN],
          params->strings[HEAT_CAPACITY][params->HEAT_CAPACITY],
          params->strings[RHO_FIRN][params->RHO_FIRN],
          params->strings[INTERNAL_ENERGY][params->INTERNAL_ENERGY],
          params->strings[SCHEME][params->SCHEME]);
  // Save the temperature profile, the melt rate and the age scale
  sprintf(fileName, "%s.dat", "melt_rate");
  saveTable(data->melt, fileName, path, params->T);
  sprintf(fileName, "%s.dat", "frozen_ice");
  saveTable(data->freeze, fileName, path, params->T);
  if (params->SAVE_TYPE == ST_MATRIX) {
    sprintf(fileName, "%s.dat", "age_matrix");
    save2DTable(age, fileName, path, ageVerRes, ageHorRes + 1, 1, 1, 0);
    sprintf(fileName, "%s.dat", "temperature_matrix");
    save2DTable(data->temperature, fileName, path, params->Z, params->T, 1, 1,
                0);
    sprintf(fileName, "%s.dat", "density_matrix");
    save2DTable(data->density, fileName, path, params->Z, params->T, 1, 1, 0);
  } else if (params->SAVE_TYPE == ST_VECTOR) {
    sprintf(fileName, "%s.dat", "age_profile");
    save2DTable(age, fileName, path, ageVerRes, ageHorRes + 1, 1, 1, 0);
    sprintf(fileName, "%s.dat", "temp_profile");
    saveTable(data->tnew, fileName, path, params->Z);
    real density[params->Z];
    memset(density, 0, params->Z * sizeof(real));
    for (int i = 0; i < params->Z; i++) {
      density[i] = data->density[i][params->T - 1];
    }
    sprintf(fileName, "%s.dat", "density_profile");
    saveTable(density, fileName, path, params->Z);
    real ice_density[params->Z];
    memset(ice_density, 0, params->Z * sizeof(real));
    for (int i = 0; i < params->Z; i++) {
      ice_density[i] = data->ice_density[i][params->T - 1];
    }
    sprintf(fileName, "%s.dat", "pure_ice_density_profile");
    saveTable(ice_density, fileName, path, params->Z);
  }
}
