#include "physicalFunctions.h"
#include <stdio.h>
//*************DEFINTION OF THE FUNCTIONS*************

static const real K_eff_WE[11] = {0.32, 0.44, 0.51, 0.64, 0.76, 0.78,
                                  0.77, 0.71, 0.65, 0.68, 0.71};

static const real K_eff_WE_delta[11] = {0.056, 0.023, 0.071, 0.154,
                                        0.238, 0.250, 0.196, 0.112,
                                        0.027, 0.014, 0.000};
//*************Computational functions*************

void setRho_HL(double rhoSnowConst, double *rho, double *rhoIce, double *temp,
               int thickness, double acc) {
  const double rhoIceConst = 917;
  const double R = 8.3144;
  const double g = 9.81;

  const double k0 = 11 * exp(-10160 / (R * temp[thickness]));
  const double k1 = 575 * exp(-21400 / (R * temp[thickness]));
  acc = acc * 31556926.;
  const real z55 = 1 / (rhoIceConst / 1000 * k0) *
                   (log(0.55 / (rhoIceConst / 1000 - 0.55)) -
                    log(rhoSnowConst / (rhoIceConst - rhoSnowConst)));
  real z0[thickness + 1];
  real pressure[thickness + 2];
  memset(pressure, 0, (thickness + 2) * sizeof(real));
  // First pass without pressure
  for (int li = thickness; li >= 0; li--) {
    rhoIce[li] = 916.5 - 0.14438 * (temp[li] - 271.16) -
                 0.00015175 * (temp[li] - 273.16) * (temp[li] - 273.16);
    if ((thickness - li) < z55) {
      z0[li] = exp(rhoIce[li] / 1000 * k0 * (thickness - li)) * rhoSnowConst /
               (rhoIce[li] - rhoSnowConst);
      rho[li] = rhoIce[li] * z0[li] / (1 + z0[li]);
    } else if (li > 2000) {
      z0[li] =
          exp(rhoIce[li] / 1000 * k1 * (thickness - li - z55) / sqrt(acc)) *
          0.55 / (rhoIce[li] / 1000 - 0.55);
      rho[li] = rhoIce[li] * z0[li] / (1 + z0[li]);
    } else {
      rho[li] = rhoIce[li];
    }
    if (rho[li] > rhoIce[li]) {
      rho[li] = rhoIce[li];
    }
    pressure[li] += g * rhoIce[li] + pressure[li + 1];
  }

  // second pass with pressure
  for (int li = 0; li <= thickness; li++) {
    rhoIce[li] += 1.1e-7 * pressure[li];
    if ((thickness - li) < z55) {
      z0[li] = exp(rhoIce[li] / 1000 * k0 * (thickness - li)) * rhoSnowConst /
               (rhoIce[li] - rhoSnowConst);
      rho[li] = rhoIce[li] * z0[li] / (1 + z0[li]);
    } else if (li > 2000) {
      z0[li] =
          exp(rhoIce[li] / 1000 * k1 * (thickness - li - z55) / sqrt(acc)) *
          0.55 / (rhoIce[li] / 1000 - 0.55);
      rho[li] = rhoIce[li] * z0[li] / (1 + z0[li]);
    } else {
      rho[li] = rhoIce[li];
    }
    if (rho[li] > rhoIce[li]) {
      rho[li] = rhoIce[li];
    }
  }
}

void setRho_CONST(double rhoSnowConst, double *rho, double *rhoIce,
                  double *temp, int thickness, double acc) {
  UNUSED(rhoSnowConst);
  UNUSED(rhoIce);
  UNUSED(temp);
  UNUSED(acc);
  for (int li = 0; li <= thickness; li++) {
    rho[li] = 921;
  }
}

void setHeatVar(const model_functions *const functions, double *K, double *cp,
                double *temperature, int thickness, double *rho,
                double *rhoIce) {


  functions->setHeatCapacity(cp, rho, rhoIce, temperature, thickness);

  functions->setThermalIce(K, temperature, thickness);

  if(functions->setThermalFirn){
    functions->setThermalFirn(K, rho, rhoIce, cp, thickness, temperature);
  }
}

void setThermalIce_CP(double *K, double *temperature, int thickness) {
  for (int li = 0; li <= thickness; li++) {
    K[li] = 9.828 * exp(-0.0057 * temperature[li]);
  }
}

void setThermalIce_GO(double *K, double *temperature, int thickness) {
  for (int li = 0; li <= thickness; li++) {
    K[li] = 2.22 * (1 - 0.0067 * (temperature[li] - 273.15));
  }
}

/**
 * @name THERMAL CONDUCTIVITY OF ICE
 * @brief Based on master thesis of Tobias Hipp, who used relationships by Ling
 * & Yhang (2005).
 * @version 11.03
 * @param Temperature Temperature (K)
 * @return Thermal conductivity of ice
 */
// double conductivity_ice(const double &Temperature) {
//   const double ki = 0.4685 + 488.19 / Temperature;
//   return ki;
// }

void setThermalFirn_CP(double *K, double *rho, double *rhoIce, double *cp,
                       int thickness, double *temperature) {
  UNUSED(cp);
  UNUSED(temperature);
  for (int li = 0; li <= thickness; li++) {
    K[li] = 2. * K[li] * rho[li] / (3 * rhoIce[li] - rho[li]);
  }
}

void setThermalFirn_SC(double *K, double *rho, double *rhoIce, double *cp,
                       int thickness, double *temperature) {
  UNUSED(cp);
  UNUSED(temperature);
  for (int li = 0; li <= thickness; li++) {
    K[li] = K[li] * pow((rho[li] / rhoIce[li]), 2 - 0.5 * rho[li] / rhoIce[li]);
  }
}

void setThermalFirn_AL(double *K, double *rho, double *rhoIce, double *cp,
                       int thickness, double *temperature) {
  UNUSED(rhoIce);
  for (int li = 0; li <= thickness; ++li) {
    K[li] = rho[li] * cp[li] * (1 - 0.00882 * (temperature[li] + 30 - 273.15)) *
            (-1.229e-14 * rho[li] * rho[li] * rho[li] +
             2.1312 * 1e-11 * rho[li] * rho[li] - 9.4e-9 * rho[li] + 1.779e-6);
  }
}

void setThermalFirn_CP_LIN(double *K, double *rho, double *rhoIce, double *cp,
                           int thickness, double *temperature) {
  real cp_top = cp[thickness];
  setThermalFirn_CP(K, rho, rhoIce, cp, thickness, temperature);
  real k0 = 25 * rho[thickness] * cp_top / (365.25 * 24 * 3600);
  real k100 = K[thickness - 100];
  for (int li = thickness; li >= thickness - 100; --li) {
    K[li] = k0 + (thickness - li) / 100. * (k100 - k0);
  }
}

void setThermalFirn_SC_LIN(double *K, double *rho, double *rhoIce, double *cp,
                           int thickness, double *temperature) {
  real cp_top = cp[thickness];
  setThermalFirn_SC(K, rho, rhoIce, cp, thickness, temperature);
  real k0 = 25 * rho[thickness] * cp_top / (365.25 * 24 * 3600);
  real k100 = K[thickness - 100];
  for (int li = thickness; li >= thickness - 100; --li) {
    K[li] = k0 + (thickness - li) / 100. * (k100 - k0);
  }
}

void setThermalFirn_CP_AL(double *K, double *rho, double *rhoIce, double *cp,
                          int thickness, double *temperature) {
  double K1[thickness + 1];
  setThermalFirn_AL(K1, rho, rhoIce, cp, thickness, temperature);
  setThermalFirn_CP(K, rho, rhoIce, cp, thickness, temperature);
  for (int li = thickness; li >= thickness - 300; --li) {
    K[li] = (K1[li] - K[li]) *
                ((rho[li] - rhoIce[li]) / (rho[thickness] - rhoIce[li])) +
            K[li];
  }
}

void setThermalFirn_SC_AL(double *K, double *rho, double *rhoIce, double *cp,
                          int thickness, double *temperature) {
  double K1[thickness + 1];
  setThermalFirn_AL(K1, rho, rhoIce, cp, thickness, temperature);
  setThermalFirn_SC(K, rho, rhoIce, cp, thickness, temperature);
  for (int li = thickness; li >= thickness - 300; --li) {
    K[li] = (K1[li] - K[li]) *
                ((rho[li] - rhoIce[li]) / (rho[thickness] - rhoIce[li])) +
            K[li];
  }
}

void setThermalFirn_ST(double *K, double *rho, double *rhoIce, double *cp,
                       int thickness, double *temperature) {
  UNUSED(rhoIce);
  UNUSED(cp);
  UNUSED(temperature);
  for (int li = 0; li <= thickness; ++li) {
    K[li] = 0.138 - 1.01e-3 * rho[li] / 1000 + 3.233e-6 * rho[li] * rho[li];
  }
}

void setThermalFirn_CP_ST(double *K, double *rho, double *rhoIce, double *cp,
                          int thickness, double *temperature) {
  double K1[thickness + 1];
  setThermalFirn_ST(K1, rho, rhoIce, cp, thickness, temperature);
  setThermalFirn_CP(K, rho, rhoIce, cp, thickness, temperature);
  int li600 = thickness;
  while (rho[li600] < 600 && li600 > 0) {
    K[li600] = K1[li600];
    --li600;
  }
  for (int li = li600; li >= 0; --li) {
    K[li] = (K1[li] - K[li]) *
                ((rho[li] - rhoIce[li]) / (rho[li600] - rhoIce[li])) +
            K[li];
  }
}

void setThermalFirn_SC_ST(double *K, double *rho, double *rhoIce, double *cp,
                          int thickness, double *temperature) {

  double K1[thickness + 1];
  setThermalFirn_ST(K1, rho, rhoIce, cp, thickness, temperature);
  setThermalFirn_SC(K, rho, rhoIce, cp, thickness, temperature);
  int li600 = thickness;
  while (rho[li600] < 600 && li600 > 0) {
    K[li600] = K1[li600];
    --li600;
  }
  for (int li = li600; li >= 0; --li) {
    K[li] = (K1[li] - K[li]) *
                ((rho[li] - rhoIce[li]) / (rho[li600] - rhoIce[li])) +
            K[li];
  }
}

void setThermalFirn_CP_WE_ADD(double *K, double *rho, double *rhoIce,
                              double *cp, int thickness, double *temperature) {
  setThermalFirn_CP(K, rho, rhoIce, cp, thickness, temperature);
  for (int li = 0; li <= 10; ++li) {
    K[thickness - li] += K_eff_WE_delta[li];
  }
}

void setThermalFirn_CP_WE_LIN(double *K, double *rho, double *rhoIce,
                              double *cp, int thickness, double *temperature) {
  setThermalFirn_CP(K, rho, rhoIce, cp, thickness, temperature);

  for (int li = thickness - 10; li >= thickness - 300; --li) {
    real factor = ((rho[li] - rhoIce[li]) / (rho[thickness - 10] - rhoIce[li]));
    K[li] = (K_eff_WE[10] - K[li]) * factor * factor * factor * factor + K[li];
  }
  for (int li = 0; li <= 10; ++li) {
    K[thickness - li] = K_eff_WE[li];
  }
}

void setThermalFirn_SC_WE_ADD(double *K, double *rho, double *rhoIce,
                              double *cp, int thickness, double *temperature) {
  setThermalFirn_SC(K, rho, rhoIce, cp, thickness, temperature);
  for (int li = 0; li <= 10; ++li) {
    K[thickness - li] += K_eff_WE_delta[li];
  }
}

void setThermalFirn_SC_WE_LIN(double *K, double *rho, double *rhoIce,
                              double *cp, int thickness, double *temperature) {
  setThermalFirn_SC(K, rho, rhoIce, cp, thickness, temperature);
  for (int li = thickness - 10; li >= thickness - 300; --li) {
    real factor = ((rho[li] - rhoIce[li]) / (rho[thickness - 10] - rhoIce[li]));
    K[li] = (K_eff_WE[10] - K[li]) * factor * factor * factor * factor + K[li];
  }
  for (int li = 0; li <= 10; ++li) {
    K[thickness - li] = K_eff_WE[li];
  }
}

void setHeatCapacity_CP(double *cp, double *rho, double *rhoIce,
                        double *temperature, int thickness) {
  UNUSED(rho);
  UNUSED(rhoIce);
  for (int li = 0; li <= thickness; li++) {
    cp[li] = 152.5 + 7.122 * temperature[li];
  }
}

void setHeatCapacity_CP_AL(double *cp, double *rho, double *rhoIce,
                           double *temperature, int thickness) {
  for (int li = 0; li <= thickness; li++) {
    real cp_al = -27.796 + 7.7752 * temperature[li];
    real cp_sc = 152.5 + 7.122 * temperature[li];
    cp[li] = (cp_al - cp_sc) *
                 ((rho[li] - rhoIce[li]) / (rho[thickness] - rhoIce[li])) +
             cp_sc;
  }
}

void computeMelt(const model_functions *const functions, double *m,
                 double *tground, double *rho, double L, double K0, double cp0,
                 double told1, double told0, double thick, double delz,
                 double QG, double *f) {
  int li = 0;
  double tmelt = 0;
  double pressure = 0;
  // Computation of the pressure and the melting point
  for (li = 0; li <= (int)thick; li++) {
    pressure += rho[li];
  }
  pressure += (thick - (int)thick) * rho[(int)thick];
  tmelt = 273.16 - 7.2 * pow(10, -8) * (pressure)*9.8;
  double diff = QG + K0 * (told1 - tmelt) / delz;

  functions->computeMelt(diff, tmelt, m, tground, rho, L, K0, cp0, told1, told0,
                         delz, QG, f);
}

void computeMelt_FREE_MELT(double diff, double tmelt, double *m,
                           double *tground, double *rho, double L, double K0,
                           double cp0, double told1, double told0, double delz,
                           double QG, double *f) {
  if (diff > 0) // If enough energy is available to melt ice
  {
    *m = 1 /
         (rho[0] * (L - cp0 * (told0 - tmelt)) + cp0 * (tmelt - told1) / 2) *
         (-rho[0] * cp0 * (tmelt - told0) / (2. * 31556926. * 100.) + diff);
    *tground = tmelt;
    *f = 0;
  } else // If not enough energy is available, bottom temperature is decreased
  {
    *m = 0;
    *tground = QG * delz / K0 + told1;
    *f = 0;
  }
}

void computeMelt_FREEZING_NO_ICE(double diff, double tmelt, double *m,
                                 double *tground, double *rho, double L,
                                 double K0, double cp0, double told1,
                                 double told0, double delz, double QG,
                                 double *f) {
  UNUSED(K0);
  UNUSED(delz);
  UNUSED(QG);

  if (diff > 0) // If enough energy is available to melt ice
  {
    *m = 1 /
         (rho[0] * (L - cp0 * (told0 - tmelt)) + cp0 * (tmelt - told1) / 2) *
         (-rho[0] * cp0 * (tmelt - told0) / (2. * 31556926. * 100.) + diff);
    *tground = tmelt;
    *f = 0;
  } else // If not enough energy is available, temperature stays at t_melt but
         // no ice is added
  {
    *m = 0;
    *tground = tmelt;
    *f = 0;
  }
}

void computeMelt_FREEZING(double diff, double tmelt, double *m, double *tground,
                          double *rho, double L, double K0, double cp0,
                          double told1, double told0, double delz, double QG,
                          double *f) {
  UNUSED(K0);
  UNUSED(delz);
  UNUSED(QG);

  if (diff > 0) // If enough energy is available to melt ice
  {
    *m = 1 /
         (rho[0] * (L - cp0 * (told0 - tmelt)) + cp0 * (tmelt - told1) / 2) *
         (-rho[0] * cp0 * (tmelt - told0) / (2. * 31556926. * 100.) + diff);
    *tground = tmelt;
    if (*m >= *f / 31556926.) {
      *m -= *f / 31556926;
      *f = 0;
    } else {
      *f -= *m * 31556926;
      *m = 0;
    }
  } else // If not enough energy is available, temperature stays at t_melt and
         // ice is added
  {
    *m = 0;
    *tground = tmelt;
    *f = diff / (rho[0] * L) * 31556926;
  }
}

double wDef_FI(double z, double thickness, double mw) {
  const double zeta = z / thickness;
  if (mw == 0.5) {
    return (zeta) * sqrt(zeta);
  } else {
    return pow(((double)z / thickness), (1 + mw));
  }
}

double wDef_PA(double z, double thickness, double mw) {
  const double p = mw;
  const double zeta = z / thickness;
  return (1 - (p + 2) / (p + 1) * (1 - zeta) +
            1 / (p + 1) * pow(1 - zeta, p + 2));
}
// Compute the matrix element a and b used in explicit and CN scheme and the
// velocity profile
void setABW(double *a, double *b, double *w, double *cp, double *K, double *rho,
            double delt, double delz, double acc, double m, double dhdt,
            double *w_def, int thickness, double *rhoIce) {
  int li = 0;

  for (li = 0; li <= thickness; li++) {
    b[li] = delt * K[li] / (rho[li] * cp[li] * delz * delz);
    w[li] = -rhoIce[li] / rho[li] * (acc - m - dhdt) * w_def[li] - m;
    a[li] = delt / (delz * 2) *
            (1 / (rho[li] * cp[li] * 2 * delz) * (K[li + 1] - K[li - 1]) - w[li]);
  }
}

void setInternal(const model_functions *const functions, double *se,
                 double *rho, double *w, double *cp, double *K, double delt,
                 int thickness, double *told, double dH, double *tborder,
                 int border, double len, double flat) {
  int li = 0;
  const int deltaH = (int)dH;

  // Internal energy production
  functions->setSe(se, rho, w, cp, thickness, told, delt);

  // Valley effect
  if (deltaH > 0 && border == 0) {
    for (li = 0; li <= thickness - deltaH; li++) {
      se[li + deltaH] += 4 * K[li + deltaH] * 2 *
                         (tborder[li] - told[li + deltaH]) * delt /
                         (rho[li + deltaH] * cp[li + deltaH] * len * len);
    }
    for (li = 1; li < deltaH; li++) {
      double l = ((double)li * (double)li / ((double)deltaH * (double)deltaH));
      se[li] += K[li] * 2 * 4 * (told[0] + li * 6.50E-4 - told[li]) /
                ((flat + l * (len - flat)) * (flat + l * (len - flat)) *
                 rho[li] * cp[li]) *
                delt;
    }
  }
}

void setSe_ON(double *se, double *rho, double *w, double *cp, int thickness,
              double *told, double delt) {
  double P = 0;
  int li = 0;
  for (li = thickness - 1; li >= 1; li--) {
    P += rho[li - 1] * 9.81;
    se[li] = 0;
    double cr =
        cbrt(getDwdz(w, li, thickness) + getDudz((double)li / thickness));
    se[li] = 2 * cr * cr * getA(told[li]) * delt / (rho[li] * cp[li]) +
             (w[li] * P / rho[li] * (rho[li + 1] - rho[li - 1]) / 2) * delt /
                 (rho[li] * cp[li]);
  }
}
void setSe_OFF(double *se, double *rho, double *w, double *cp, int thickness,
               double *told, double delt) {
  UNUSED(rho);
  UNUSED(w);
  UNUSED(cp);
  UNUSED(told);
  UNUSED(delt);
  int li = 0;
  for (li = thickness - 1; li >= 1; li--) {
    se[li] = 0;
  }
}

// Linear piecewise approximation of the horizontal velocity vertical derivative
// profile squared
double getDudz(double zh) {
  const double us = 2.2594e-19;
  double dudz = 0;
  if (zh >= 0 && zh <= 0.05) {
    dudz = 1.942488E-06 - 1.952373E-05 * zh;
  } else if (zh <= 0.1) {
    dudz = 1.357489E-06 - 8.359790E-06 * zh;
  } else if (zh <= 0.2) {
    dudz = 8.170050E-07 - 2.934782E-06 * zh;
  } else if (zh <= 0.3) {
    dudz = 4.579192E-07 - 1.076580E-06 * zh;
  } else if (zh <= 0.5) {
    dudz = 2.656468E-07 - 4.343890E-07 * zh;
  } else if (zh <= 0.75) {
    dudz = 1.327851E-07 - 1.664704E-07 * zh;
  } else if (zh <= 1) {
    dudz = 4.236023E-08 - 4.340220E-08 * zh;
  }
  return (us * dudz);
}

// Square of the vertical derivative of the vertical velocity profile
double getDwdz(double *w, int z, int thickness) {
  double dwdz;
  if (z == 0) {
    dwdz = w[1] - w[0];
  }
  if (z == thickness) {
    dwdz = w[thickness] - w[thickness - 1];
  } else {
    dwdz = (w[z + 1] - w[z - 1]) / 2;
  }
  return dwdz * dwdz;
}

// Linear piecewise approximation of A power -1/3 value profile
double getA(double t) {
  double A = 0;
  if (t <= 210) {
    A = 31519576971 - 141986135 * t;
  } else if (t <= 220) {
    A = 19796924504 - 86163980 * t;
  } else if (t <= 230) {
    A = 8301327112 - 33868831 * t;
  } else if (t <= 240) {
    A = 4956807495 - 19262412 * t;
  } else if (t <= 255) {
    A = 2755424785 - 10096091 * t;
  } else if (t <= 275) {
    A = 1802280733 - 6324074 * t;
  }
  return A;
}
