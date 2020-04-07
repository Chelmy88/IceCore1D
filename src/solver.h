#ifndef SOLVER_H
#define SOLVER_H

#include "define.h"
#include "physicalFunctions.h"
#include "structures.h"
#include <math.h>
#include <omp.h>
#include <stdio.h>
#include <time.h>
#include <unistd.h>

//*************Computational function*************

// int thickness --> Ice thickness obtained from main.c, transformed to an
// integer double tnew[Z] --> Table used to store the new temperature computed
// double told[Z] --> Table used to store the temperature at the begining of the
// time step obtained from main.c int i,li --> Variables used in loops double L
// =333500 --> Latent heat of ice in J/kg double rho[Z] --> Table used to store
// the computed density profile double rhoIce[Z] --> Table used to compute the
// pure ice density profile (for actual temperature and pressure) double
// a[Z],b[Z] --> double a2[Z],b2[Z] --> double m --> Melt rate double tground -->
// Ground temperature double K[Z]  --> Ice thermal conductivity double cp[Z] -->
// Ice specific heat capacity double w[Z] --> Table used to store the velocity
// profile double w_def[Z] --> Table used to store the flux shape function values
// double delt=31556926.*100 --> Time step (100kyr)
// double delz=1 -->  Height Step
// double tmelt --> Melt temperature computed with the bottom pressure
// double dhdt  --> Thickness time derivative
// double se[Z]  --> Internal heat (valley effect + internal heat production)
// power density

// double rhoIceConst=917  -->  Pure ice density for a first approximation of
// the density profile double rhoSnowConst  --> Snow density (value defined in
// header) double R=8.3144  --> Gaz constant double k0  --> Value computed in the
// H-L density model double k1  --> Value computed in the H-L density model
// double z55  --> Value computed in the H-L density model
// double z0[Z]  --> Values computed in the H-L density model

// double c1[Z]  --> Sub-diagonal matrix element computed in the explicit scheme
// double c2[Z]  --> Diagonal matrix element computed in the explicit scheme
// double c3[Z]  --> Sub-diagonal matrix element computed in the explicit scheme

// double l[Z]  --> Sub-diagonal matrix element computed in the C-N scheme
// double d[Z]  --> Diagonal matrix element computed in the C-N scheme
// double r[Z]  --> Sub-diagonal matrix element computed in the C-N scheme
// double b[Z]  --> Vector to be multiplied with he inverse matrix in the C-N
// scheme (explicit part)

void spin_up(const model_functions *const functions,
             const model_parameters *const params, double *temperature,
             double thick, double tsurf, double acc, double QG, double mw,
             double *tborder, double deltaH, int border, double len,
             double flat, double *melt, double *density);
// Perform the spin_up of the model for the given time using the CN scheme

void t_solve(const model_functions *const functions,
             const model_parameters *const params, double *temperature,
             int time, double thickness, double thicknessFuture, double tsurf,
             double acc, double *melt, double QG, double mw, double *tborder,
             double deltaH, int border, double len, double flat, double *freeze,
             double *density);
// Getting the 1D array temperature a t-1, return the temperature at T.
// This function calls various function to compute all the needed parameters
// Finally, this function calls the defined algorithm to compute the temperature

void integrate_CN(double *tint, double *told, double *alpha, double *beta,
                  double *alpha1, double *beta1, double tground, double tsurf,
                  int thickness, double *se); //, double* se);
// Compute the temperature using the CN scheme, called by spin_up() and
// t_solve()

void integrate_expl(double *told, double *a, double *b, double tground,
                    double tsurf, double tsurf_old, int thickness, double *se);
// Compute the temperature using the explicit scheme, called by  t_solve()

void tempScale(double *told, double thickness, double thicknessFuture,
               double tsurf);
// Scale the temperature profile to the thickness value of the next step

#endif /* !SOLVER_H */
