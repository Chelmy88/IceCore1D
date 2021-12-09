#ifndef STRUCTURES_H
#define STRUCTURES_H

#include "define.h"
#include <stdbool.h>
#include <stdlib.h>

typedef struct _model_values {
  real *mw;
  size_t mw_n;
  real *QG;
  size_t QG_n;
  real *TCor;
  size_t TCor_n;
  real *TCor2;
  size_t TCor2_n;
  real *PCor;
  size_t PCor_n;
  real *deltaH;
  size_t deltaH_n;
  real *len;
  size_t len_n;
  real *flat;
  size_t flat_n;
  size_t tot;
} model_values;

enum SAVE_TYPE_ENUM { ST_MATRIX, ST_VECTOR, ST_UNSET };
enum SCHEME_ENUM { SC_CN, SC_EXPL, SC_UNSET };
enum THERMAL_ICE_ENUM { TI_CP, TI_GO, TI_UNSET };
enum THERMAL_FIRN_ENUM {
  TF_SC,
  TF_CP,
  TF_CP_LIN,
  TF_SC_LIN,
  TF_CP_AL,
  TF_SC_AL,
  TF_SC_ST,
  TF_CP_ST,
  TF_CP_WE_ADD,
  TF_CP_WE_LIN,
  TF_SC_WE_ADD,
  TF_SC_WE_LIN,
  TF_NULL,
  TF_UNSET
};
enum RHO_FIRN_ENUM { RF_HL, RF_CONST, RF_UNSET };
enum HEAT_CAPACITY_ENUM { CP_CP, CP_CP_AL, CP_UNSET };
enum VERTICAL_PROFILE_ENUME { VP_FI, VP_PA, VP_UNSET };
enum INTERNAL_ENERGY_ENUM { IE_ON, IE_OFF, IE_UNSET };
enum MELTING_ENUM { ME_FREE_MELT, ME_FREEZING_NO_ICE, ME_FREEZING, ME_UNSET };

enum DATA_ENUM {
  SAVE_TYPE,
  SCHEME,
  THERMAL_ICE,
  THERMAL_FIRN,
  RHO_FIRN,
  HEAT_CAPACITY,
  VERTICAL_PROFILE,
  INTERNAL_ENERGY,
  MELTING,
  DATA_ENUM_SIZE
};

// Structure to store time series informations
typedef struct _model_parameters {
  int Z; // 3400//height of the table
  int T; // 10001//10001 //width of the table
  int S; // 1500 //Length of the spin up in hYr
  enum SAVE_TYPE_ENUM
      SAVE_TYPE; //"VECTOR"//MATRIX or VECTOR to save age and temperature
  enum SCHEME_ENUM SCHEME; //"CN"//Scheme used, values can be CN or expl
  int RHO_SNOW; // 350 //Value of the snow density used in the computation of
                // the density profile
  enum THERMAL_ICE_ENUM
      THERMAL_ICE; // "CP"//Model used for themal parameters, can be CP or GO
  enum THERMAL_FIRN_ENUM
      THERMAL_FIRN; // "SC" // Correction for the firn thermal
  // conductivity, can be CP,SC or FI
  enum HEAT_CAPACITY_ENUM
      HEAT_CAPACITY;           // "SC" // Correction for the firn thermal
                               // conductivity, can be CP,SC or FI
  enum RHO_FIRN_ENUM RHO_FIRN; // "FIRN" // Set the density profile to realistic
                               // (FIRN) or constant (CONST)
  enum VERTICAL_PROFILE_ENUME
      VERTICAL_PROFILE; // "FI" // Set the flux shape function to FI or PA
  enum INTERNAL_ENERGY_ENUM INTERNAL_ENERGY; // "OFF" //Decide wether internal
                                             // energy should be included or not
  enum MELTING_ENUM
      MELTING; // "FREE_MELT" //Basal malting-refeezing handeling :
               // FREE_MELT->no basal refreezing, temperature decreases if
               // there is no melt, FREEZING_NO_ICE -> some refreezing is
               // possible, but the ice dissapear (i.e. bottom temp is always
               // tmelt, no other difference), FREEZING -> some water is
               // allowed to refreez, when melting comes back, first this ice
               // is melted before real melting occures (refreezing and melting
               // of frozen ice have no inflence on vertical velocity).
  char *TEMPERATURE_FILE;          //"output"//Directory to store data file into
  char *ACCUMULATION_FILE;         //"output"//Directory to store data file into
  char *ICE_THICKNESS_FILE;        //"output"//Directory to store data file into
  char *AGE_FILE;                  //"output"//Directory to store data file into
  char *BOREHOLE_TEMPERATURE_FILE; //"output"//Directory to store data file into
  char *OUTPUT_PATH;               //"output"//Directory to store data file into
  char ***strings;
  model_values values;
} model_parameters;

bool initModelParameters(model_parameters *params, char *fileName);

void printModelParameters(model_parameters *params);

void deleteModelParameters(model_parameters *params);

// Structure to store time series informations
typedef struct _time_series {
  real *surfaceTempLoad;
  real *iceThicknessLoad;
  real *accLoad;
  real *age;
  real *borehole_temp;
} time_series;

bool initTimeSeries(time_series *ts, model_parameters *params);

void deleteTimeSeries(time_series *ts);

// Structure to store model data while running
typedef struct _model_data {
  real **temperature;
  real **density;
  real **ice_density;
  real *tnew;
  real *surfaceTemp;
  real *iceThickness;
  real *acc;
  real *acc2;
  real *melt;
  real *freeze;
  real mw;
  real QG;
  real tCor;
  real tCor2;
  real pCor;
  real deltaH;
  real len;
  real flat;
} model_data;

bool initModelData(model_data *data, const model_parameters *const params,
                   size_t mwL, size_t QGL, size_t TcorL, size_t TcorL2,
                   size_t PcorL, size_t deltaHL, size_t lenL, size_t flatL);

void deleteModelData(model_data *data, const model_parameters *const params);

typedef struct _model_functions {
  void (*setRho)(double rhoSnowConst, double *rho, double *rhoIce, double *temp,
                 int thickness, double acc);
  // Compute the density profile

  // Compute the values of the K and c thermal variables, called by spin_up()
  // and t_solve()
  void (*setThermalIce)(double *K, double *temperature, int thickness);
  void (*setThermalFirn)(double *K, double *rho, double *rhoIce, double *cp,
                         int thickness, double *temperature);
  void (*setHeatCapacity)(double *cp, double *rho, double *rhoIce,
                          double *temperature, int thickness);

  void (*computeMelt)(double diff, double tmelt, double *m, double *tground,
                      double *rho, double L, double K0, double cp0,
                      double told1, double told0, double delz, double QG,
                      double *f);
  // Compute the melt rate, called by spin_up() and t_solve()

  double (*wDef)(double z, double thickness, double mw);
  // Compute the flux shape function values, called by spin_up() and t_solve()

  void (*setSe)(double *se, double *rho, double *w, double *cp, int thickness,
                double *told, double delt);
  // Compute the internal energy production and the lateral heat flux (valley
  // effect)

} model_functions;

bool initModelFunctions(model_functions *functions,
                        const model_parameters *const params);

#endif /* !STRUCTURES_H */
