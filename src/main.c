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

// This code is developed to run on Linux machine, but should run on Windows or Mac OS
//The only parameters defined in main.c are the free parameters of the model, the correction to the boundary condition time series,
//and the initial temperature profile for spin up
//For any other modification see the main.h file

//Include the header file main.h, containing all the key functions and parameters.

#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>
#include <math.h>
#include <time.h>
#include "io.h"
#include "define.h"
#include "solver.h"
#include "structures.h"
#include "runModel.h"

int main()
{

    //Set internal timer
  double begin = omp_get_wtime();

	//Read the model parameters
  model_parameters params;
  if(!initModelParameters(&params, "init.txt"))
  {
    deleteModelParameters(&params);
    printf("[E] Error reading the ini file. Exiting now\n");
    exit(EXIT_FAILURE);
  }

  //Tables to load data
  time_series ts;
  initTimeSeries(&ts,&params);

  double ageGRIP[Z],tGRIP[Z]= {0};
  int i=0;

  //DEfine a table to store a summary of the various runs performed
  double** summary;
  summary = malloc( 15000  * sizeof(*summary));
  for (i = 0; i < 15000; i++)
  {
      summary[i] = malloc(11*sizeof(summary));
  }
  //Load borehole temperature and age profile for comparison
  readTable(ageGRIP,"time_series/EDC_age_forC.dat");
  readTable(tGRIP,"time_series/EDC_temp_forC.dat");

  //The following array contain the value of the different free parameters
  double mwArr[] = {0.5,0.6}; //Form factor
  int mwN=sizeof(mwArr) / sizeof(mwArr[0]);
 // int mwL=0;
  double QGArr[] = {0.054}; //Ground heat flux
  int QGN=sizeof(QGArr) / sizeof(QGArr[0]);
  double TcorArr[] = {1}; //Correction for temperature between 6000 and 2000 BP in K
  int TcorN=sizeof(TcorArr) / sizeof(TcorArr[0]);
  double TcorArr2[] = {2}; //Correction for cold periods temperature in K
  int TcorN2=sizeof(TcorArr2) / sizeof(TcorArr2[0]);
  double PcorArr[] = {10}; //Correction of accumulation time series in %
  int PcorN=sizeof(PcorArr) / sizeof(PcorArr[0]);
  double deltaHArr[] = {100}; //Depth of the valley
  int deltaHN=sizeof(deltaHArr) / sizeof(deltaHArr[0]);
  double lenArr[] = {5000}; //Width of the valley
  int lenN=sizeof(lenArr) / sizeof(lenArr[0]);
  double flatArr[] = {500}; //Flat arrea at the bottom of the valley
  int flatN=sizeof(flatArr) / sizeof(flatArr[0]);

  int tot=mwN*QGN*TcorN*PcorN*deltaHN*lenN*flatN*TcorN2; //Number of run, the size of the summary table should be bigger

  // Loops for the values of the free parameter
  int count=0;
  //Initialize the core splitting, should be placed in front of a loop having if possible a number of elements corresponding to a multiple of the core numbers
  #pragma omp parallel for collapse(8)
  for (int mwL=0; mwL<mwN; mwL++)
  {
  //int QGL=0;
  for (int QGL=0; QGL<QGN; QGL++)
  {
  // int TcorL=0;
  for (int TcorL=0; TcorL<TcorN; TcorL++)
  {
  //int TcorL2=0;
  for (int TcorL2=0; TcorL2<TcorN2; TcorL2++)
  {
  //int PcorL=0;
  for (int PcorL=0; PcorL<PcorN; PcorL++)
  {
  //int deltaHL=0;
  for (int deltaHL=0; deltaHL<deltaHN; deltaHL++)
  {
  //int lenL=0;
  for (int lenL=0; lenL<lenN; lenL++)
  {
  //int flatL=0;
  for (int flatL=0; flatL<flatN; flatL++)
  {
  const double mw=mwArr[mwL];
	const double QG=QGArr[QGL];
  const double tCor=TcorArr[TcorL];
  const double tCor2=TcorArr2[TcorL2];
  const double pCor=PcorArr[PcorL];
  const double deltaH=deltaHArr[deltaHL];
  const double len=lenArr[lenL];
  const double flat=flatArr[flatL];

  printf("%f %f %f",ts.surfaceTempLoad[0],ts.iceThicknessLoad[0],ts.accLoad[0]);

  printf("\n I'm running with : %f %f %f %f %f %f %f %f \n",mw, QG, tCor, tCor2, pCor, deltaH, len, flat);

  int li,co=0;

  // ALLOCATE //
  double** temperature;
  temperature = malloc( Z  * sizeof(*temperature));
  for (li = 0; li < Z; li++)
  {
    temperature[li] = malloc(T*sizeof(temperature));
  }
  double** density;
  density = malloc( Z  * sizeof(*density));
  for (li = 0; li < Z; li++)
  {
    density[li] = malloc(T*sizeof(density));
  }
  for (li = 0; li < Z; li++)
  {
    for (co=0; co<T; co++)
    {
      density[li][co] = 0;
    }
  }

  double surfaceTemp[T],iceThickness[T],acc[T],acc2[T],melt[T],freeze[T]= {0};
  double tnew[Z]={0};

  runModel(tnew, surfaceTemp, iceThickness, acc, acc2, melt, freeze, temperature,
           density, ts.surfaceTempLoad, ts.iceThicknessLoad, ts.accLoad,
           mw, QG, tCor, tCor2, pCor, deltaH, len, flat);
  // POST PROCESSING//
  //Compute the difference with the age profile for 213 points and compute the mean of the difference
 	printf("\nComputing age ... ");
 	fflush(stdout);
 	double begin3=omp_get_wtime();
 	double** ageRel;
 	int ageVerRes = 0;
  int ageHorRes = 0;
  int ageCor= 0;
  if(strcmp(SAVE_TYPE,"MATRIX")==0){
    if (T==10001){
      ageVerRes = (int) (Z/5);
      ageHorRes = (int) (T/10);
      ageCor=0;
    }
    else if (T==40001){
      ageVerRes = (int) (Z/5);
      ageHorRes = (int) (T/20);
      ageCor=20000;
    }
  }
  else if (strcmp(SAVE_TYPE,"VECTOR")==0){
    ageVerRes = (int) (Z/5);
    ageHorRes = 1;
    ageCor=T-11;
  }

  ageRel = malloc( (int) ageVerRes * sizeof(*ageRel));
  for (i = 0; i < ageVerRes; i++)
  {
    ageRel[i] = malloc((ageHorRes+1)*sizeof(ageRel));
  }
  double ageDiff=0;
  int a=0;
   //   #pragma omp parallel for
  for(a=ageHorRes; a>0; a--){
    int h=0;
    for(h=1; h<ageVerRes; h++)
    {
      float height=h*5;
      int age=a*10+ageCor;
      while (height<iceThickness[age] && age>0)
      {
        height+=((acc[age]*31556926*100-melt[age]*100-(iceThickness[age]-iceThickness[age-1]))*wDef((double)height,iceThickness[age],mw)+melt[age]*100);
        age--;
      }
      ageRel[h-1][0]=h*5;
      ageRel[h-1][a]=(a*10+ageCor-age)*100;
    }
  }

  printf("Age ok in in: %f secondes\n",(double)(omp_get_wtime() - begin3));
  //Compute the difference with the borehole  temperature profile below 600 m deep (because upper part of the measurements are affected by seasonality)
  double tempDiff=0;
  double tnew2[Z]= {0};
  for(li=0; li<=(int)iceThickness[T-1]-7; li++)
  {
    tnew2[li]=temperature[li][T-1];
    if(li<((int)iceThickness[T-1]-600))
    {
        tempDiff+=fabs(tnew2[li]-tGRIP[li]);
    }
  }
  tempDiff/=(iceThickness[T-1]);

  // Generate a file name with the free parameters
  char fileName[180]="";
  char path[180]="";
  sprintf(path, "%s/m_%.3f_Q_%.2f_Pcor_%.0f_Tcor_%.1f_Tcor2_%.1f_dH_%.0f_len_%.0f_flat_%.0f_%s_Thermal_%s_Firn_%s_Internal_Energy_%s_Scheme_%s",DIR_PATH,mw,QG*1000,pCor,tCor,tCor2,deltaH,len,flat,"EDC",THERMAL,FIRN,INTERNAL_ENERGY,TYPE);

  //Check if the man export directory and the subdirectory are already existing and create them if missing
  struct stat st = {0};
  if (stat(DIR_PATH, &st) == -1) {
  mkdir(DIR_PATH, 0700);
  }
  if (stat(path, &st) == -1) {
    mkdir(path, 0700);
  }
  // Save the temperature profile, the melt rate and the age scale
  //sprintf(fileName, "%s.dat","temp_profile");
  //saveTable(tnew,fileName,path,Z);
  sprintf(fileName, "%s.dat","melt_rate");
  saveTable(melt,fileName,path,T);
  //sprintf(fileName, "%s.dat","frozen_ice");
  //saveTable(freeze,fileName,path,T);
  if(strcmp(SAVE_TYPE,"MATRIX")==0){
    sprintf(fileName, "%s.dat","age_matrix");
    save2DTable(ageRel,fileName,path,ageVerRes,ageHorRes+1,1,1,0);
    sprintf(fileName, "%s.dat","temp_matrix");
    if (T==10001){
    	save2DTable(temperature,fileName,path,Z,T,1,10,0);
    }
    else if(T==40001){
     	save2DTable(temperature,fileName,path,Z,T,1,10,20000);
    }
  }
  else if(strcmp(SAVE_TYPE,"VECTOR")==0){
    sprintf(fileName, "%s.dat","age_profile");
    save2DTable(ageRel,fileName,path,ageVerRes,ageHorRes+1,1,1,0);
    sprintf(fileName, "%s.dat","temperature_today");
    saveTable(tnew,fileName,path,Z);
  }
  sprintf(fileName, "%s.dat","density_matrix_top");
  if (T==10001){
      save2DTable_top(density,fileName,path,iceThickness,200,T);
  }
  else if(T==40001){
    save2DTable_top(density,fileName,path,iceThickness,200,T);
  }
  sprintf(fileName, "%s.dat","temperature_matrix_top");
  if (T==10001){
    save2DTable_top(temperature,fileName,path,iceThickness,200,T);
  }
  else if(T==40001){
    save2DTable_top(temperature,fileName,path,iceThickness,200,T);
  }
  sprintf(fileName, "%s.dat","surface_temperature");
  saveTable(surfaceTemp,fileName,path,T);
  sprintf(fileName, "%s.dat","accumulation_rate");
  saveTable(acc2,fileName,path,T);
  sprintf(fileName, "%s.dat","ice_thickness");
  saveTable(iceThickness,fileName,path,T);
  //Create a summary table for all the runs
  summary[count][0]=mw;
  summary[count][1]=QG;
  summary[count][2]=pCor;
  summary[count][3]=tCor;
  summary[count][4]=tCor2;
  summary[count][5]=deltaH;
  summary[count][6]=len;
  summary[count][7]=flat;
  summary[count][8]=tempDiff;
  summary[count][9]=ageDiff;
  summary[count][10]=tnew2[0]-273.15;
  #pragma omp atomic
  count++;
  #pragma omp flush (count)
  printf("END LOOP : %d/%d\n\n",count,tot);
      //Uncomment the line below to save the whole temperature matrix
     // saveTemp(temperature);
  for (li = 0; li < Z; li++)
  {
    double* currentIntPtr = temperature[li];
    free(currentIntPtr);
  }
  fflush(stdout);
}}}}}}}}

    //Save the summary table
    //save2DTable(summary,"Summary",path,15000,11,1,1);
  printf("Run OK -- in %f seconds\n",(double)(omp_get_wtime() - begin));

  return 0;
}
