#include "runModel.h"
#include <time.h>

void runModel(double* tnew, double* surfaceTemp,double* iceThickness,double* acc,double* acc2,
              double* melt,double* freeze, double** temperature, double** density,
              const double* const surfaceTempLoad,const double* const iceThicknessLoad,
              const double* const accLoad,
              const double mw, const double QG, const double tCor, const double tCor2, const double pCor,
              const double deltaH, const double len,const double flat)
{
  size_t li=0;

  double spin_up_temp[Z],spin_up_temp2[Z],temperatureBorder[Z]= {0};
  //Dynamically allow the memory for the temperature matrix

  double dens[Z]= {0};
  // SET BC VALUES //
  //Loop over the time steps to implement the correction on the time series
  for(li=0; li<T; li++)
  {
    surfaceTemp[li]=surfaceTempLoad[li];
    iceThickness[li]=iceThicknessLoad[li];
    acc[li]=accLoad[li]*3600*24*365/31556926;
    //surfaceTemp[li]=surfaceTemp[li]+(iceThicknessLoad[T-1]-iceThicknessLoad[li])/100;
  }
  double t0=surfaceTemp[T-1];
  double tLGM=surfaceTemp[T-257];

  for(li=0; li<T; li++){
    surfaceTemp[li]=(surfaceTemp[li]-t0)*(tLGM-t0+tCor2)/(tLGM-t0)+t0;
    if(T==40001){
      if (li>39980 && li<40001)
        {
        surfaceTemp[li]=surfaceTemp[li]-tCor;
      }
    }
    else if(T==10001){
      if (li>9980 && li<10001)
      {
        surfaceTemp[li]=surfaceTemp[li]-tCor;
      }
    }
    acc[li]+=acc[li]*pCor/100.;
    acc2[li]=acc[li]*31556926;
  }
  // set the initial temperature profile
  double Tmelt=273.16-9.8*7.42*1E-8*921*iceThickness[0];
  double Tmelt2=273.16-9.8*7.42*1E-8*921*(iceThickness[0]-deltaH);
  double Tsurf=surfaceTemp[0];

  // SPIN UP MODEL//

  //Spin up the second profile used for the valley effect
  if(deltaH!=0)
  {
    for(li=0; li <=(size_t)(iceThickness[0]-deltaH); li++)
    {
      spin_up_temp2[li]=Tmelt2+(Tsurf-Tmelt2)*pow(li/(iceThickness[0]-deltaH),1);
    }
    spin_up(spin_up_temp2,iceThickness[0]-deltaH,surfaceTemp[0],acc[0],QG,mw,temperatureBorder,deltaH,1,len,flat,melt,dens);
    for(li=0; li <=(size_t)(iceThickness[0]-deltaH); li++)
    {
      temperatureBorder[li]=spin_up_temp2[li];
    }
  }
  //Spin up the main profile first without the valley effect and a second time if the valley effect is enabled
  for(li=0; li <=(size_t)iceThickness[0]; li++)
  {
    spin_up_temp[li]=Tmelt+(Tsurf-Tmelt)*pow(li/iceThickness[0],1);
  }
  spin_up(spin_up_temp,iceThickness[0],surfaceTemp[0],acc[0],QG,mw,temperatureBorder,deltaH,1,len,flat,melt, dens);

  if(deltaH!=0){
    spin_up(spin_up_temp,iceThickness[0],surfaceTemp[0],acc[0],QG,mw,temperatureBorder,deltaH,0,len,flat,melt, dens);
  }

  for(li=0; li<Z; li++)
  {
    temperature[li][0]=spin_up_temp[li];
    density[li][0]=dens[li];
    tnew[li]=spin_up_temp[li];
  }

  // RUN MODEL//

  //Loop over all the time steps
  int time=1;
  float time_for_loop=0;
  for (time=1; time<T; time++)
  {
    double begin2=omp_get_wtime();
    for(li=0; li <Z; li++)
    {
      tnew[li]=0;
      dens[li]=0;
    }
    for(li=0; li <=(size_t)iceThickness[time-1]; li++)
    {
      tnew[li]=temperature[li][time-1];
    }
    if(deltaH!=0)
    {
      t_solve(temperatureBorder,time, iceThickness[time-1]-deltaH, iceThickness[time]-deltaH, surfaceTemp[time],acc[time],melt,QG,mw,temperatureBorder,deltaH,1,len,flat,freeze,dens);
    }
    t_solve(tnew,time, iceThickness[time-1], iceThickness[time], surfaceTemp[time],acc[time],melt,QG,mw,temperatureBorder,deltaH,0,len,flat,freeze,dens);
    tempScale(tnew,  iceThickness[time-1], iceThickness[time], surfaceTemp[time]);
    if(deltaH!=0)
    {
      tempScale(temperatureBorder,iceThickness[time-1]-deltaH, iceThickness[time]-deltaH, surfaceTemp[time]);
    }
    for(li=0; li <=(size_t)iceThickness[time]; li++)
    {
      temperature[li][time]=tnew[li];
      density[li][time]=dens[li];
    }
    time_for_loop+=(double)(omp_get_wtime() - begin2); //Store the loop time
  }
//melt[0]=melt[1];
//Print the time for the run and the individual loop mean time
  printf("\nIntegration ok in: %f secondes  ",time_for_loop);
  printf("Mean run time: %f miliseconds\n",time_for_loop/(T-1.)*1000.);

}
