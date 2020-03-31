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

bool mainLoop(model_parameters* params,time_series* ts, model_functions* functions, double**  summary);

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
  if(!initTimeSeries(&ts,&params))
  {
    deleteTimeSeries(&ts);
    deleteModelParameters(&params);
    printf("[E] Error reading the time series files. Exiting now\n");
    exit(EXIT_FAILURE);
  }
  //Tables to load data
  model_functions functions;
  if(!initModelFunctions(&functions,&params))
  {
    deleteModelFunctions(&functions);
    deleteTimeSeries(&ts);
    deleteModelParameters(&params);
    printf("[E] Error Initializing model functions. Exiting now\n");
    exit(EXIT_FAILURE);
  }


  //Define a table to store a summary of the various runs performed
  double** summary;
  summary = malloc( 15000  * sizeof(*summary));
  for (size_t i = 0; i < 15000; i++)
  {
      summary[i] = malloc(11*sizeof(summary));
  }

  size_t exit_status;
  if(mainLoop(&params, &ts, &functions, summary))
  {
    printf("Run OK -- in %f seconds\n",(double)(omp_get_wtime() - begin));
    exit_status=0;
  }
  else{
    exit_status=EXIT_FAILURE;
  }

  deleteModelParameters(&params);
  deleteTimeSeries(&ts);
  deleteModelFunctions(&functions);

  return exit_status;
}

bool mainLoop(model_parameters* params, time_series* ts, model_functions* functions, double**  summary){

  bool succes=true;
  int tot=params->values.tot; //Number of run, the size of the summary table should be bigger
  // Loops for the values of the free parameter
  int count=0;
  //Initialize the core splitting, should be placed in front of a loop having if possible a number of elements corresponding to a multiple of the core numbers
  #pragma omp parallel for collapse(8) shared(succes)
  for (size_t mwL=0; mwL<params->values.mw_n; mwL++)
  {
  for (size_t QGL=0; QGL<params->values.QG_n; QGL++)
  {
  for (size_t TcorL=0; TcorL<params->values.TCor_n; TcorL++)
  {
  for (size_t TcorL2=0; TcorL2<params->values.TCor2_n; TcorL2++)
  {
  for (size_t PcorL=0; PcorL<params->values.PCor_n; PcorL++)
  {
  for (size_t deltaHL=0; deltaHL<params->values.deltaH_n; deltaHL++)
  {
  for (size_t lenL=0; lenL<params->values.len_n; lenL++)
  {
  for (size_t flatL=0; flatL<params->values.flat_n; flatL++)
  {
  if(!succes)
  {
    continue;
  }
  model_data data;

  if(!initModelData(&data,params,mwL,QGL,TcorL,TcorL2,PcorL,deltaHL,lenL,flatL))
  {
    deleteModelData(&data,params);
    printf("[E] Error Initializing mode data. Exiting when all threads are done\n");
    succes = false;
    continue;
  }

  //printf("\n I'm running with : %f %f %f %f %f %f %f %f \n",mw, QG, tCor, tCor2, pCor, deltaH, len, flat);

  runModel(&data,params,ts,functions);

  // // POST PROCESSING//
  // //Compute the difference with the age profile for 213 points and compute the mean of the difference
 	printf("\nComputing age ... ");
 	fflush(stdout);
 	double begin3=omp_get_wtime();
 	double** ageRel;
 	size_t ageVerRes = 0;
  size_t ageHorRes = 0;
  size_t ageCor= 0;
  if(strcmp(SAVE_TYPE,"MATRIX")==0){
    if (T==10001){
      ageVerRes = (size_t) (Z/5);
      ageHorRes = (size_t) (T/10);
      ageCor=0;
    }
    else if (T==40001){
      ageVerRes = (size_t) (Z/5);
      ageHorRes = (size_t) (T/20);
      ageCor=20000;
    }
  }
  else if (strcmp(SAVE_TYPE,"VECTOR")==0){
    ageVerRes = (size_t) (Z/5);
    ageHorRes = 1;
    ageCor=T-11;
  }

   ageRel = malloc( (size_t) ageVerRes * sizeof(*ageRel));
   for (size_t i = 0; i < ageVerRes; i++)
   {
     ageRel[i] = malloc((ageHorRes+1)*sizeof(ageRel));
  }
  double ageDiff=0;   //   #pragma omp parallel for
  for(size_t a=ageHorRes; a>0; a--){
    for(size_t h=1; h<ageVerRes; h++)
    {
      float height=h*5;
      int age=a*10+ageCor;
      while (height<data.iceThickness[age] && age>0)
      {
        height+=((data.acc[age]*31556926*100-data.melt[age]*100-
                 (data.iceThickness[age]-data.iceThickness[age-1]))*wDef((double)height,
                 data.iceThickness[age],data.mw)+data.melt[age]*100);
        age--;
      }
      ageRel[h-1][0]=h*5;
      ageRel[h-1][a]=(a*10+ageCor-age)*100;
    }
  }

  printf("Age ok in in: %f secondes\n",(double)(omp_get_wtime() - begin3));
  // //Compute the difference with the borehole  temperature profile below 600 m deep (because upper part of the measurements are affected by seasonality)
  // double tempDiff=0;
  // double tnew2[Z]= {0};
  // for(size_t li=0; li<=(size_t)iceThickness[T-1]-7; li++)
  // {
  //   tnew2[li]=temperature[li][T-1];
  //   if(li<((size_t)iceThickness[T-1]-600))
  //   {
  //       tempDiff+=fabs(tnew2[li]-ts->borehole_temp[li]);
  //   }
  // }
  // tempDiff/=(iceThickness[T-1]);
  //
  // // Generate a file name with the free parameters
   char fileName[180]="";
   char path[180]="";
   sprintf(path, "%s/m_%.3f_Q_%.2f_Pcor_%.0f_Tcor_%.1f_Tcor2_%.1f_dH_%.0f_len_%.0f_flat_%.0f_%s_Thermal_%s_Firn_%s_Internal_Energy_%s_Scheme_%s",
           DIR_PATH,data.mw,data.QG*1000,data.pCor,data.tCor,data.tCor2,data.deltaH,
           data.len,data.flat,"EDC",THERMAL,FIRN,INTERNAL_ENERGY,TYPE);
  //
  //Check if the man export directory and the subdirectory are already existing and create them if missing
  struct stat st = {0};
  if (stat(DIR_PATH, &st) == -1) {
  mkdir(DIR_PATH, 0700);
  }
  if (stat(path, &st) == -1) {
    mkdir(path, 0700);
  }
  // Save the temperature profile, the melt rate and the age scale
  sprintf(fileName, "%s.dat","temp_profile");
  saveTable(data.tnew,fileName,path,Z);
  sprintf(fileName, "%s.dat","melt_rate");
  saveTable(data.melt,fileName,path,T);
  sprintf(fileName, "%s.dat","frozen_ice");
  saveTable(data.freeze,fileName,path,T);
  if(strcmp(SAVE_TYPE,"MATRIX")==0){
    // sprintf(fileName, "%s.dat","age_matrix");
    // save2DTable(ageRel,fileName,path,ageVerRes,ageHorRes+1,1,1,0);
    // sprintf(fileName, "%s.dat","temp_matrix");
    // if (T==10001){
    // 	save2DTable(temperature,fileName,path,Z,T,1,10,0);
    // }
    // else if(T==40001){
    //  	save2DTable(temperature,fileName,path,Z,T,1,10,20000);
    // }
  }
  else if(strcmp(SAVE_TYPE,"VECTOR")==0){
    sprintf(fileName, "%s.dat","age_profile");
    save2DTable(ageRel,fileName,path,ageVerRes,ageHorRes+1,1,1,0);
    sprintf(fileName, "%s.dat","temperature_today");
    saveTable(data.tnew,fileName,path,Z);
  }
  // // sprintf(fileName, "%s.dat","density_matrix_top");
  // // if (T==10001){
  // //     save2DTable_top(density,fileName,path,iceThickness,200,T);
  // // }
  // // else if(T==40001){
  // //   save2DTable_top(density,fileName,path,iceThickness,200,T);
  // // }
  // // sprintf(fileName, "%s.dat","temperature_matrix_top");
  // // if (T==10001){
  // //   save2DTable_top(temperature,fileName,path,iceThickness,200,T);
  // // }
  // // else if(T==40001){
  // //   save2DTable_top(temperature,fileName,path,iceThickness,200,T);
  // // }
  // sprintf(fileName, "%s.dat","surface_temperature");
  // saveTable(surfaceTemp,fileName,path,T);
  // sprintf(fileName, "%s.dat","accumulation_rate");
  // saveTable(acc2,fileName,path,T);
  // sprintf(fileName, "%s.dat","ice_thickness");
  // saveTable(iceThickness,fileName,path,T);
  // //Create a summary table for all the runs
  // summary[count][0]=mw;
  // summary[count][1]=QG;
  // summary[count][2]=pCor;
  // summary[count][3]=tCor;
  // summary[count][4]=tCor2;
  // summary[count][5]=deltaH;
  // summary[count][6]=len;
  // summary[count][7]=flat;
  // summary[count][8]=tempDiff;
  // summary[count][9]=ageDiff;
  // summary[count][10]=tnew2[0]-273.15;
  #pragma omp atomic
  count++;
  #pragma omp flush (count)
  printf("END LOOP : %d/%d\n\n",count,tot);
      //Uncomment the line below to save the whole temperature matrix
     // saveTemp(temperature);
  fflush(stdout);
  deleteModelData(&data,params);


}}}}}}}}
  for (size_t i = 0; i < 15000; i++)
  {
      double* currentIntPtr = summary[i];
      free(currentIntPtr);
  }
  char path[120];
  sprintf(path,"%s/init.txt",params->OUTPUT_PATH);
  copyFile("init.txt",path);
  return succes;

}
