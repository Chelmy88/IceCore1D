// Written by Adrien Michel
// adrien.michel@no-log.org
// For the purpose of a Master thesis at the
// Climate and Environmental Group
// And
// Oeschger Center for Climate Change Research
// University of Bern
// February 2016

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

int main()
{

    //Set internal timer
    double begin = omp_get_wtime();

		//Read the model parameters
    model_parameters params;
    if(!initModelParameters(&params, "init.txt"))
		{
			deleteModelParameters(&params);
			exit(EXIT_FAILURE);
		}

    //Tables to load data
    time_series ts;
    initTimeSeries(&ts);
    double surfaceTempLoad[T],iceThicknessLoad[T],accLoad[T]= {0};
    double ageGRIP[Z],tGRIP[Z]= {0};
    int i=0;

    //DEfine a table to store a summary of the various runs performed
    double** summary;
    summary = malloc( 15000  * sizeof(*summary));
    for (i = 0; i < 15000; i++)
    {
        summary[i] = malloc(11*sizeof(summary));
    }

    //Load boundary condition LR04EDC time series. Data file exist for 1Myr of 4Myr
    if(T==10001){
	    readTable(surfaceTempLoad,"time_series/LR04-EDC_temp_1Myr.dat");
	    readTable(iceThicknessLoad,"time_series/LR04-EDC_thickness_1Myr.dat");
	    readTable(accLoad,"time_series/LR04-EDC_acc_1Myr.dat");
      //readTable(ts.surfaceTempLoad,"../time_series/LR04-EDC_temp_1Myr.dat");
      //readTable(ts.iceThicknessLoad,"../time_series/LR04-EDC_thickness_1Myr.dat");
      //readTable(ts.accLoad,"../time_series/LR04-EDC_acc_1Myr.dat");
    }
    else if(T==40001){
    	readTable(surfaceTempLoad,"time_series/LR04-EDC_temp_4Myr.dat");
	    readTable(iceThicknessLoad,"time_series/LR04-EDC_thickness_4Myr.dat");
	    readTable(accLoad,"time_series/LR04-EDC_acc_4Myr.dat");
      //readTable(ts.surfaceTempLoad,"../time_series/LR04-EDC_temp_1Myr.dat");
      //readTable(ts.iceThicknessLoad,"../time_series/LR04-EDC_thickness_1Myr.dat");
      //readTable(ts.accLoad,"../time_series/LR04-EDC_acc_1Myr.dat");
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
		double mw=mwArr[mwL];
		double QG=QGArr[QGL];
		double tCor=TcorArr[TcorL];
		double tCor2=TcorArr2[TcorL2];
		double pCor=PcorArr[PcorL];
		double deltaH=deltaHArr[deltaHL];
		double len=lenArr[lenL];
    	double flat=flatArr[flatL];

    	printf("\n I'm running with : %f %f %f %f %f %f %f %f \n",mw, QG, tCor, tCor2, pCor, deltaH, len, flat);

        double surfaceTemp[T],iceThickness[T],acc[T],acc2[T],melt[T],freeze[T]= {0};
        double spin_up_temp[Z],spin_up_temp2[Z],temperatureBorder[Z],tnew[Z]= {0};

        //Dynamically allow the memory for the temperature matrix
        double** temperature;
        int li,co=0;

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

        double dens[Z]= {0};

        //Loop over the time steps to implement the correction on the time series
        for(li=0; li<T; li++)
        {
            surfaceTemp[li]=surfaceTempLoad[li];
            iceThickness[li]=iceThicknessLoad[li];
            acc[li]=accLoad[li]*3600*24*365/31556926;
          //  surfaceTemp[li]=surfaceTemp[li]+(iceThicknessLoad[T-1]-iceThicknessLoad[li])/100;
        }
        double t0=surfaceTemp[T-1];
		double tLGM=surfaceTemp[T-257];

	//	double acc0=acc[T-1];
		//double accLGM=acc[T-257];

        for(li=0; li<T; li++){
			surfaceTemp[li]=(surfaceTemp[li]-t0)*(tLGM-t0+tCor2)/(tLGM-t0)+t0;
			//acc[li]=(1.5*pow(2,(surfaceTemp[li]-273.15)/10)-0.005)/31556926;

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

       //    	acc[li]=(acc[li]-acc0)*(accLGM-acc0+accLGM*(pCor/100))/(accLGM-acc0)+acc0;
            acc[li]+=acc[li]*pCor/100.;
            acc2[li]=acc[li]*31556926;

//            if((li>38900 && li<39850) ||( li< 38650 && li > 38100) || ( li< 37800 && li > 37650) || ( li< 37500 && li > 36800) || ( li< 36600 && li > 36100)|| ( li< 35700 && li > 35200)|| ( li< 34700 && li > 34400)|| ( li< 33700 && li > 33300)|| ( li< 32600 && li > 32400)|| ( li< 32050 && li > 31900) || ( li< 31300 && li > 30700) )
  //          {
            //    surfaceTemp[li]+=tCor2;
    //        }
      //      else
        //    {
        }


        // set the initial temperature profile
        double Tmelt=273.16-9.8*7.42*1E-8*921*iceThickness[0];
        double Tmelt2=273.16-9.8*7.42*1E-8*921*(iceThickness[0]-deltaH);
        double Tsurf=surfaceTemp[0];

        //Spin up the second profile used for the valley effect
        if(deltaH!=0)
        {
            for(li=0; li <=(int)(iceThickness[0]-deltaH); li++)
            {
                spin_up_temp2[li]=Tmelt2+(Tsurf-Tmelt2)*pow(li/(iceThickness[0]-deltaH),1);
            }
            spin_up(spin_up_temp2,iceThickness[0]-deltaH,surfaceTemp[0],acc[0],QG,mw,temperatureBorder,deltaH,1,len,flat,melt,dens);
            for(li=0; li <=(int)(iceThickness[0]-deltaH); li++)
            {
                temperatureBorder[li]=spin_up_temp2[li];
            }

        }
        //Spin up the main profile first without the valley effect and a second time if the valley effect is enabled
        for(li=0; li <=(int)iceThickness[0]; li++)
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
            for(li=0; li <=(int)iceThickness[time-1]; li++)
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

            for(li=0; li <=(int)iceThickness[time]; li++)
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
