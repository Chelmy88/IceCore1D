#ifndef STRUCTURES_H
#define STRUCTURES_H

#include <stdbool.h>
#include <stdlib.h>
#include "define.h"

// Structure to store time series informations
typedef struct _time_series {
  double* surfaceTempLoad;
  double* iceThicknessLoad;
  double* accLoad;
  double* age;
  double* borehole_temp;
} time_series;

enum SAVE_TYPE_ENUM{ST_MATRIX, ST_VECTOR, ST_UNSET};
enum SCHEME_ENUM {SC_CN, SC_EXPL, SC_UNSET};
enum THERMAL_ENUM {TH_CP, TH_GO, TH_UNSET};
enum FIRN_ENUM {FI_SC, FI_CP, FI_FI, FI_UNSET};
enum RHO_TYPE_ENUM {RHO_FIRN, RHO_CONST, RHO_UNSET};
enum VERTICAL_PROFILE_ENUME {VP_FI, VP_PA, VP_UNSET};
enum INTERNAL_ENERGY_ENUM {IE_ON, IE_OFF, IE_UNSET};
enum MELTING_ENUM {ME_FREE_MELT, ME_FREEZING_NO_ICE, ME_FREEZING, ME_UNSET};

bool initTimeSeries(time_series *ts);

void deleteTimeSeries(time_series *ts);


// Structure to store time series informations
typedef struct _model_parameters {
  int Z1; //3400//height of the table
  int T1; //10001//10001 //width of the table
  int S1; //1500 //Length of the spin up in hYr
  enum SAVE_TYPE_ENUM SAVE_TYPE1; //"VECTOR"//MATRIX or VECTOR to save age and temperature
  enum SCHEME_ENUM SCHEME1; //"CN"//Scheme used, values can be CN or expl
  int rhoSnow1; //350 //Value of the snow density used in the computation of the density profile
  enum THERMAL_ENUM THERMAL1;// "CP"//Model used for themal parameters, can be CP or GO
  enum FIRN_ENUM FIRN1;// "SC" // Correction for the firn thermal conductivity, can be CP,SC or FI
  enum RHO_TYPE_ENUM RHO1;// "FIRN" // Set the density profile to realistic (FIRN) or constant (CONST)
  enum VERTICAL_PROFILE_ENUME VERTICAL1;// "FI" // Set the flux shape function to FI or PA
  enum INTERNAL_ENERGY_ENUM INTERNAL_ENERGY1;// "OFF" //Decide wether internal energy should be included or not
  enum MELTING_ENUM MELTING1;// "FREE_MELT" //Basal malting-refeezing handeling : FREE_MELT->no basal refreezing, temperature decreases if there is no melt, FREEZING_NO_ICE -> some refreezing is possible, but the ice dissapear (i.e. bottom temp is always tmelt, no other difference), FREEZING -> some water is allowed to refreez, when melting comes back, first this ice is melted before real melting occures (refreezing and melting of frozen ice have no inflence on vertical velocity).
  char* OUTPUT_PATH; //"output"//Directory to store data file into
}model_parameters;

bool initModelParameters(model_parameters *params, char* fileName);

void printModelParameters(model_parameters *params);

void deleteModelParameters(model_parameters *params);


#endif  /* !STRUCTURES_H */
