#ifndef STRUCTURES_H
#define STRUCTURES_H

#include <stdbool.h>
#include <stdlib.h>
#include "define.h"

typedef struct _model_values {
  real* mw;
  size_t mw_n;
  real* QG;
  size_t QG_n;
  real* TCor;
  size_t TCor_n;
  real* TCor2;
  size_t TCor2_n;
  real* PCor;
  size_t PCor_n;
  real* deltaH;
  size_t deltaH_n;
  real* len;
  size_t len_n;
  real* flat;
  size_t flat_n;
  size_t tot;
}model_values;

enum SAVE_TYPE_ENUM{ST_MATRIX, ST_VECTOR, ST_UNSET};
enum SCHEME_ENUM {SC_CN, SC_EXPL, SC_UNSET};
enum THERMAL_ENUM {TH_CP, TH_GO, TH_UNSET};
enum FIRN_ENUM {FI_SC, FI_CP, FI_FI, FI_UNSET};
enum RHO_TYPE_ENUM {RHO_FIRN, RHO_CONST, RHO_UNSET};
enum VERTICAL_PROFILE_ENUME {VP_FI, VP_PA, VP_UNSET};
enum INTERNAL_ENERGY_ENUM {IE_ON, IE_OFF, IE_UNSET};
enum MELTING_ENUM {ME_FREE_MELT, ME_FREEZING_NO_ICE, ME_FREEZING, ME_UNSET};
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
  enum VERTICAL_PROFILE_ENUME VERTICAL_PROFILE1;// "FI" // Set the flux shape function to FI or PA
  enum INTERNAL_ENERGY_ENUM INTERNAL_ENERGY1;// "OFF" //Decide wether internal energy should be included or not
  enum MELTING_ENUM MELTING1;// "FREE_MELT" //Basal malting-refeezing handeling : FREE_MELT->no basal refreezing, temperature decreases if there is no melt, FREEZING_NO_ICE -> some refreezing is possible, but the ice dissapear (i.e. bottom temp is always tmelt, no other difference), FREEZING -> some water is allowed to refreez, when melting comes back, first this ice is melted before real melting occures (refreezing and melting of frozen ice have no inflence on vertical velocity).
  char* TEMPERATURE_FILE; //"output"//Directory to store data file into
  char* ACCUMULATION_FILE; //"output"//Directory to store data file into
  char* ICE_THICKNESS_FILE; //"output"//Directory to store data file into
  char* AGE_FILE; //"output"//Directory to store data file into
  char* BOREHOLE_TEMPERATURE_FILE; //"output"//Directory to store data file into
  char* OUTPUT_PATH; //"output"//Directory to store data file into
  model_values values;
}model_parameters;

bool initModelParameters(model_parameters *params, char* fileName);

void printModelParameters(model_parameters *params);

void deleteModelParameters(model_parameters *params);


// Structure to store time series informations
typedef struct _time_series {
  real* surfaceTempLoad;
  real* iceThicknessLoad;
  real* accLoad;
  real* age;
  real* borehole_temp;
} time_series;


bool initTimeSeries(time_series *ts, model_parameters *params);

void deleteTimeSeries(time_series *ts);

// Structure to store model data while running
typedef struct _model_data {
  real** temperature;
  real** density;
  real* tnew;
  real* surfaceTemp;
  real* iceThickness;
  real* acc;
  real* acc2;
  real* melt;
  real* freeze;
  real mw;
  real QG;
  real tCor;
  real tCor2;
  real pCor;
  real deltaH;
  real len;
  real flat;
} model_data;

bool initModelData(model_data *data,const model_parameters * const params,size_t mwL,
                   size_t QGL,size_t TcorL, size_t TcorL2,size_t PcorL,size_t deltaHL,size_t lenL,size_t flatL);

void deleteModelData(model_data *data,const model_parameters * const params);


typedef struct _model_functions {
  void (*setRho)(double* rho, double *rhoIce, double* temp, int thickness,double acc);
  //Compute the density profile

  void (*setHeatVar)(double *K,double *cp,double *told,int thickness,double *rho, double* rhoIce);
  //Compute the values of the K and c thermal variables, called by spin_up() and t_solve()

  void (*computeMelt)(double* m,double* tground,double* rho,double L,double K0,double cp0, double told1,double told0,double thickness,double delz,double QG, double* f);
  //Compute the melt rate, called by spin_up() and t_solve()

  double (*wDef)(double z, double thickness,double mw);
  //Compute the flux shape function values, called by spin_up() and t_solve()

  void (*setABW)(double* a,double* b,double* w,double* cp,double* K,double* rho,double delt,double delz,double acc,double m,double dhdt,double* w_def,int thickness, double* rhoIce);
  //Compute the vertical velocity and the a,b (explicit scheme) or alpha,beta(CN scheme) values, called by spin_up() and t_solve()

  void (*setSe)(double *se,double *rho,double *w, double *cp, double *K,double delt, int thickness, double* told, double deltaH,double dhdt,double * tborder,int border,double len, double flat);
  //Compute the internal energy production and the lateral heat flux (valley effect)

  double (*getDwdz)(double*w,int z,int thickness);
  //Compute the vertical derivative of the vertical velocity profile, called by setSe()

  double (*getA)(double t);
  //Compute the creep factor A values from piecewise linear approximation, called by setSe()

  double (*getDudz)(double zh);
  //Compute the vertical derivative of horizontal velocity profile from piecewise linear approximation, called by setSe()

} model_functions;

bool initModelFunctions(model_functions *functions,const model_parameters * const params);

void deleteModelFunctions(model_functions *functions);

#endif  /* !STRUCTURES_H */
