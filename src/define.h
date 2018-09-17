#ifndef DEFINE_H
#define DEFINE_H
//*******PRE PROCESSOR VARIABLES**********
// The preprocessor variables allow to change parameters of the model. If non accepted values are entered a fatal error may occur, values are not verified.
#define Z 3400//height of the table
#define T 10001//10001 //width of the table
#define S 1500 //Length of the spin up in hYr
#define DIR_PATH "output"//Directory to store data file into
#define SAVE_TYPE "VECTOR"//MATRIX or VECTOR to save age and temperature
#define TYPE "CN"//Scheme used, values can be CN or expl
#define rhoSnow 350 //Value of the snow density used in the computation of the density profile
#define THERMAL "CP"//Model used for themal parameters, can be CP or GO
#define FIRN "SC" // Correction for the firn thermal conductivity, can be CP,SC or FI
#define RHO "FIRN" // Set the density profile to realistic (FIRN) or constant (CONST)
#define VERTICAL "FI" // Set the flux shape function to FI or PA
#define INTERNAL_ENERGY "OFF" //Decide wether internal energy should be included or not
#define MELTING "FREE_MELT" //Basal malting-refeezing handeling : FREE_MELT->no basal refreezing, temperature decreases if there is no melt, FREEZING_NO_ICE -> some refreezing is possible, but the ice dissapear (i.e. bottom temp is always tmelt, no other difference), FREEZING -> some water is allowed to refreez, when melting comes back, first this ice is melted before real melting occures (refreezing and melting of frozen ice have no inflence on vertical velocity).

#endif  /* !DEFINE_H */
