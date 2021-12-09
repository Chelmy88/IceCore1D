#include "runModel.h"
#include <time.h>

void runModel(model_data *data, const model_parameters *const params,
              const time_series *const ts,
              const model_functions *const functions) {
  size_t li = 0;
  const size_t Z = params->Z;
  const size_t T = params->T;

  real *spin_up_temp, *spin_up_temp2, *temperatureBorder, *dens, *ice_density;

  spin_up_temp = calloc(params->Z, sizeof(real));
  spin_up_temp2 = calloc(params->Z, sizeof(real));
  temperatureBorder = calloc(params->Z, sizeof(real));
  dens = calloc(params->Z, sizeof(real));
  ice_density = calloc(params->Z, sizeof(real));

  // SET BC VALUES //
  // Loop over the time steps to implement the correction on the time series
  for (li = 0; li < T; li++) {
    data->surfaceTemp[li] = ts->surfaceTempLoad[li];
    data->iceThickness[li] = ts->iceThicknessLoad[li];
    data->acc[li] = ts->accLoad[li] * 3600 * 24 * 365 / SEC_YEAR;
    // surfaceTemp[li]=surfaceTemp[li]+(iceThicknessLoad[T-1]-iceThicknessLoad[li])/100;
  }
  data->iceThickness[T] = data->iceThickness[T - 1];

  real t0 = data->surfaceTemp[T - 1];
  real tLGM = data->surfaceTemp[T - 257];

  for (li = 0; li < T; li++) {
    data->surfaceTemp[li] =
        (data->surfaceTemp[li] - t0) * (tLGM - t0 + data->tCor2) / (tLGM - t0) +
        t0;
    if (T == 40001) {
      if (li > 39980 && li < 40001) {
        data->surfaceTemp[li] = data->surfaceTemp[li] - data->tCor;
      }
    } else if (T == 10001) {
      if (li > 9980 && li < 10001) {
        data->surfaceTemp[li] = data->surfaceTemp[li] - data->tCor;
      }
    }
    data->acc[li] += data->acc[li] * data->pCor / 100.;
    data->acc2[li] = data->acc[li] * SEC_YEAR;
  }
  // set the initial temperature profile
  real Tmelt = 273.16 - 9.8 * 7.42 * 1E-8 * 921 * data->iceThickness[0];
  real Tmelt2 =
      273.16 - 9.8 * 7.42 * 1E-8 * 921 * (data->iceThickness[0] - data->deltaH);
  real Tsurf = data->surfaceTemp[0];

  // SPIN UP MODEL//

	double spin_up_timer = omp_get_wtime();


  // Spin up the second profile used for the valley effect
  if (data->deltaH != 0) {
    for (li = 0; li <= (size_t)(data->iceThickness[0] - data->deltaH); li++) {
      spin_up_temp2[li] =
          Tmelt2 + (Tsurf - Tmelt2) *
                       pow(li / (data->iceThickness[0] - data->deltaH), 1);
    }

    spin_up(functions, params, spin_up_temp2,
            data->iceThickness[0] - data->deltaH, data->surfaceTemp[0],
            data->acc[0], data->QG, data->mw, temperatureBorder, data->deltaH,
            1, data->len, data->flat, data->melt, dens);
    for (li = 0; li <= (size_t)(data->iceThickness[0] - data->deltaH); li++) {
      temperatureBorder[li] = spin_up_temp2[li];
    }
  }

  // Spin up the main profile first without the valley effect and a second time
  // if the valley effect is enabled
  for (li = 0; li <= (size_t)data->iceThickness[0]; li++) {
    spin_up_temp[li] =
        Tmelt + (Tsurf - Tmelt) * pow(li / data->iceThickness[0], 1);
  }
  spin_up(functions, params, spin_up_temp, data->iceThickness[0],
          data->surfaceTemp[0], data->acc[0], data->QG, data->mw,
          temperatureBorder, data->deltaH, 1, data->len, data->flat, data->melt,
          dens);

  if (data->deltaH != 0) {
    spin_up(functions, params, spin_up_temp, data->iceThickness[0],
            data->surfaceTemp[0], data->acc[0], data->QG, data->mw,
            temperatureBorder, data->deltaH, 0, data->len, data->flat,
            data->melt, dens);
  }

	printf("[I] Spin up done in %f seconds\n", (double)(omp_get_wtime() - spin_up_timer));


  for (li = 0; li < Z; li++) {
    data->temperature[li][0] = spin_up_temp[li];
    data->density[li][0] = dens[li];
    data->tnew[li] = spin_up_temp[li];
  }

  // RUN MODEL//

  // Loop over all the time steps
  size_t time = 1;
  real time_for_loop = 0;
  for (time = 1; time < T; time++) {
    real begin2 = omp_get_wtime();
    for (li = 0; li < Z; li++) {
      ice_density[li] = 0.;
      dens[li] = 0.;
    }

    for (li = 0; li <= (size_t)data->iceThickness[time - 1]; li++) {
      data->tnew[li] = data->temperature[li][time - 1];
    }

    tempScale(data->tnew, data->iceThickness[time - 1],
              data->iceThickness[time], data->surfaceTemp[time - 1], Z);

    if (data->deltaH != 0) {
      tempScale(temperatureBorder, data->iceThickness[time - 1] - data->deltaH,
                data->iceThickness[time] - data->deltaH,
                data->surfaceTemp[time - 1], Z);
    }

    if (data->deltaH != 0) {

      t_solve(functions, params, temperatureBorder, time,
              data->iceThickness[time] - data->deltaH,
              data->iceThickness[time + 1] - data->deltaH,
              data->surfaceTemp[time], data->acc[time], data->melt, data->QG,
              data->mw, temperatureBorder, data->deltaH, 1, data->len,
              data->flat, data->freeze, dens, ice_density);
    }

    t_solve(functions, params, data->tnew, time, data->iceThickness[time],
            data->iceThickness[time + 1], data->surfaceTemp[time],
            data->acc[time], data->melt, data->QG, data->mw, temperatureBorder,
            data->deltaH, 0, data->len, data->flat, data->freeze, dens,
            ice_density);

    for (li = 0; li < Z; li++) {
      data->temperature[li][time] = data->tnew[li];
      data->density[li][time] = dens[li];
      data->ice_density[li][time] = ice_density[li];
    }
    time_for_loop += (real)(omp_get_wtime() - begin2); // Store the loop time
  }

  free(spin_up_temp);
  free(spin_up_temp2);
  free(temperatureBorder);
  free(dens);
  free(ice_density);

  // melt[0]=melt[1];
  // Print the time for the run and the individual loop mean time
  printf("[I] Integration ok in: %f secondes \n", time_for_loop);
  printf("[I] Mean run time: %f miliseconds\n", time_for_loop / (T - 1.) * 1000.);
}
