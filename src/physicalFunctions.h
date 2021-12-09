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

// This code is developed to run on Linux machine, but should run on Windows or
// Mac OS

//*******LIBRARY**********
// Open MP should be installed (http://openmp.org/wp/)
// -lm and -fopenmp flags must be used for compilation

#include "define.h"
#include "structures.h"
#include <math.h>
#include <string.h>
#include <unistd.h>

//*******Definition the  main variables defined in main.h***********

//*************File management functions*************

// char* name= --> Used to store the name of the created file
// char fileName[120] --> Used to combine the relative path and the file name of
// the created file

//*******FUNCTIONS PROTOTYPE***********

//*************Computational functions*************

// Compute the density profile
void setRho_HL(const double rhoSnowConst, double *rho, double *rhoIce, const double *temp,
               const int thickness, double acc);
void setRho_CONST(const double rhoSnowConst, double *rho, double *rhoIce,
                  const double *temp, const int thickness, double acc);

// Compute the values of the K and c thermal variables, called by spin_up() and
// t_solve()
void setHeatVar(const model_functions *const functions, double *K, double *cp,
                double *told, int thickness, double *rho, double *rhoIce);

void setHeatCapacity_CP(double *cp, double *rho, double *rhoIce,
                        double *temperature, int thickness);

void setHeatCapacity_CP_AL(double *cp, double *rho, double *rhoIce,
                           double *temperature, int thickness);

void setThermalIce_CP(double *K, double *temperature, int thickness);

void setThermalIce_GO(double *K, double *temperature, int thickness);

void setThermalFirn_CP(double *K, double *rho, double *rhoIce, double *cp,
                       int thickness, double *temperature);

void setThermalFirn_SC(double *K, double *rho, double *rhoIce, double *cp,
                       int thickness, double *temperature);

void setThermalFirn_AL(double *K, double *rho, double *rhoIce, double *cp,
                       int thickness, double *temperature);

void setThermalFirn_CP_LIN(double *K, double *rho, double *rhoIce, double *cp,
                           int thickness, double *temperature);

void setThermalFirn_SC_LIN(double *K, double *rho, double *rhoIce, double *cp,
                           int thickness, double *temperature);

void setThermalFirn_CP_AL(double *K, double *rho, double *rhoIce, double *cp,
                          int thickness, double *temperature);

void setThermalFirn_SC_AL(double *K, double *rho, double *rhoIce, double *cp,
                          int thickness, double *temperature);

void setThermalFirn_ST(double *K, double *rho, double *rhoIce, double *cp,
                       int thickness, double *temperature);

void setThermalFirn_SC_ST(double *K, double *rho, double *rhoIce, double *cp,
                          int thickness, double *temperature);

void setThermalFirn_CP_ST(double *K, double *rho, double *rhoIce, double *cp,
                          int thickness, double *temperature);

void setThermalFirn_CP_WE_ADD(double *K, double *rho, double *rhoIce,
                              double *cp, int thickness, double *temperature);

void setThermalFirn_CP_WE_LIN(double *K, double *rho, double *rhoIce,
                              double *cp, int thickness, double *temperature);

void setThermalFirn_SC_WE_ADD(double *K, double *rho, double *rhoIce,
                              double *cp, int thickness, double *temperature);

void setThermalFirn_SC_WE_LIN(double *K, double *rho, double *rhoIce,
                              double *cp, int thickness, double *temperature);

void computeMelt(const model_functions *const functions, double *m,  double prev_melt,
                 double *tground, double *rho, double L, double K0, double cp0,
                 double told1, double told0, double thickness, double delz,
                 double QG, double *f);

void computeMelt_FREE_MELT(double diff, double tmelt, double *m, double prev_melt,
                           double *tground, double *rho, double L, double K0,
                           double cp0, double told1, double told0, double delz,
                           double QG, double *f);
void computeMelt_FREEZING_NO_ICE(double diff, double tmelt, double *m, double prev_melt,
                                 double *tground, double *rho, double L,
                                 double K0, double cp0, double told1,
                                 double told0, double delz, double QG,
                                 double *f);
void computeMelt_FREEZING(double diff, double tmelt, double *m, double prev_melt,
                          double *tground, double *rho, double L,
                          double K0, double cp0, double told1,
                          double told0, double delz, double QG,
                          double *f);

// Compute the flux shape function values, called by spin_up() and t_solve()
double wDef_FI(double z, double thickness, double mw);
double wDef_PA(double z, double thickness, double mw);

void setABW(double *a, double *b, double *w, double *cp, double *K, double *rho,
            double delt, double delz, double acc, double m, double dhdt,
            double *w_def, int thickness, double *rhoIce);
// Compute the vertical velocity and the a,b (explicit scheme) or alpha,beta(CN
// scheme) values, called by spin_up() and t_solve()

void setInternal(const model_functions *const functions, double *se,
                 double *rho, double *w, double *cp, double *K, double delt,
                 int thickness, double *told, double deltaH, double *tborder,
                 int border, double len, double flat);
// Compute the internal energy production and the lateral heat flux (valley
// effect)
void setSe_ON(double *se, double *rho, double *w, double *cp, int thickness,
              double *told, double delt);
void setSe_OFF(double *se, double *rho, double *w, double *cp, int thickness,
               double *told, double delt);

double getDwdz(double *w, int z, int thickness);
// Compute the vertical derivative of the vertical velocity profile, called by
// setSe()

double getA(double t);
// Compute the creep factor A values from piecewise linear approximation, called
// by setSe()

double getDudz(double zh);
// Compute the vertical derivative of horizontal velocity profile from piecewise
// linear approximation, called by setSe()

#endif /* !PHYSICALFUNCTIONS_H */
