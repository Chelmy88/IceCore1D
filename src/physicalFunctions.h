#ifndef PHYSICALFUNCTIONS_H
#define PHYSICALFUNCTIONS_H
// Written by Adrien Michel
// adrien.michel@no-log.org
// For the purpose of a Master thesis at the
// Climate and Environmental Group
// And
// Oeschger Center for Climate Change Research
// University of Bern
// February 2016

// This code is developed to run on Linux machine, but should run on Windows or Mac OS

//*******LIBRARY**********
// Open MP should be installed (http://openmp.org/wp/)
// -lm and -fopenmp flags must be used for compilation

#include <unistd.h>
#include <math.h>
#include <string.h>
#include "define.h"
#include "structures.h"


//*******Definition the  main variables defined in main.h***********


//*************File management functions*************

// char* name= --> Used to store the name of the created file
// char fileName[120] --> Used to combine the relative path and the file name of the created file




//*******FUNCTIONS PROTOTYPE***********

//*************Computational functions*************


//Compute the density profile
void setRho_FIRN(double* rho, double *rhoIce, double* temp, int thickness,double acc);
void setRho_CONST(double* rho, double *rhoIce, double* temp, int thickness,double acc);

//Compute the values of the K and c thermal variables, called by spin_up() and t_solve()
void setHeatVar(const model_functions * const functions, double *K,double *cp,double *told,int thickness,double *rho, double* rhoIce);

void setThermalIce_CP(double *K,double *temperature,int thickness);

void setThermalIce_GO(double *K,double *temperature,int thickness);

void setThermalFirn_CP(double *K,double *rho, double* rhoIce,int thickness);

void setThermalFirn_SC(double *K,double *rho, double* rhoIce,int thickness);

void setHeatCapacity(double *cp,double *temperature,int thickness);






void computeMelt(double* m,double* tground,double* rho,double L,double K0,double cp0, double told1,double told0,double thickness,double delz,double QG, double* f);
//Compute the melt rate, called by spin_up() and t_solve()

double wDef(double z, double thickness,double mw);
//Compute the flux shape function values, called by spin_up() and t_solve()

void setABW(double* a,double* b,double* w,double* cp,double* K,double* rho,double delt,double delz,double acc,double m,double dhdt,double* w_def,int thickness, double* rhoIce);
//Compute the vertical velocity and the a,b (explicit scheme) or alpha,beta(CN scheme) values, called by spin_up() and t_solve()

void setSe(double *se,double *rho,double *w, double *cp, double *K,double delt, int thickness, double* told, double deltaH,double dhdt,double * tborder,int border,double len, double flat);
//Compute the internal energy production and the lateral heat flux (valley effect)

double getDwdz(double*w,int z,int thickness);
//Compute the vertical derivative of the vertical velocity profile, called by setSe()

double getA(double t);
//Compute the creep factor A values from piecewise linear approximation, called by setSe()

double getDudz(double zh);
//Compute the vertical derivative of horizontal velocity profile from piecewise linear approximation, called by setSe()

#endif  /* !PHYSICALFUNCTIONS_H */
