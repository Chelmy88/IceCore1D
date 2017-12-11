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
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>
#include <math.h>
#include <time.h>
#include <omp.h>

//*******PREÊPROCESSORÊVARIABLES**********
// The preprocessor variables allow to change parameters of the model. If non accepted values are entered a fatal error may occur, values are not verified.
#define Z 3400//height of the table
#define T 10001//10001 //width of the table
#define S 1500 //Length of the spin up in hYr
#define DIR_PATH "For_daniel"//Directory to store data file into
#define SAVE_TYPE "VECTOR"//MATRIX or VECTOR to save age and temperature
#define TYPE "CN"//Scheme used, values can be CN or expl
#define rhoSnow 350 //Value of the snow density used in the computation of the density profile
#define THERMAL "CP"//Model used for themal parameters, can be CP or GO
#define FIRN "SC" // Correction for the firn thermal conductivity, can be CP,SC or FI
#define RHO "FIRN" // Set the density profile to realistic (FIRN) or constant (CONST)
#define VERTICAL "FI" // Set the flux shape function to FI or PA
#define INTERNAL_ENERGY "OFF" //Decide wether internal energy should be included or not
#define MELTING "FREE_MELT" //Basal malting-refeezing handeling : FREE_MELT->no basal refreezing, temperature decreases if there is no melt, FREEZING_NO_ICE -> some refreezing is possible, but the ice dissapear (i.e. bottom temp is always tmelt, no other difference), FREEZING -> some water is allowed to refreez, when melting comes back, first this ice is melted before real melting occures (refreezing and melting of frozen ice have no inflence on vertical velocity).


//*******Definition the  main variables defined in main.h***********


//*************File management functions*************

// char* name= --> Used to store the name of the created file
// char fileName[120] --> Used to combine the relative path and the file name of the created file


//*************Computational function*************

//int thickness --> Ice thickness obtained from main.c, transformed to an integer
//double tnew[Z] --> Table used to store the new temperature computed
//double told[Z] --> Table used to store the temperature at the begining of the time step obtained from main.c
//int i,li --> Variables used in loops
//double L =333500 --> Latent heat of ice in J/kg
//double rho[Z] --> Table used to store the computed density profile
//double rhoIce[Z] --> Table used to compute the pure ice density profile (for actual temperature and pressure)
//double a[Z],b[Z] -->
//double a2[Z],b2[Z] -->
//double m --> Melt rate
//double tground --> Ground temperature
//double K[Z]  --> Ice thermal conductivity
//double cp[Z] --> Ice specific heat capacity
//double w[Z] --> Table used to store the velocity profile
//double w_def[Z] --> Table used to store the flux shape function values
//double delt=31556926.*100 --> Time step (100kyr)
//double delz=1 -->  Height Step
//double tmelt --> Melt temperature computed with the bottom pressure
//double dhdt  --> Thickness time derivative
//double se[Z]  --> Internal heat (valley effect + internal heat production) power density

//double rhoIceConst=917  -->  Pure ice density for a first approximation of the density profile
//double rhoSnowConst  --> Snow density (value defined in header)
//double R=8.3144  --> Gaz constant
//double k0  --> Value computed in the H-L density model
//double k1  --> Value computed in the H-L density model
//double z55  --> Value computed in the H-L density model
//double z0[Z]  --> Values computed in the H-L density model

//double c1[Z]  --> Sub-diagonal matrix element computed in the explicit scheme
//double c2[Z]  --> Diagonal matrix element computed in the explicit scheme
//double c3[Z]  --> Sub-diagonal matrix element computed in the explicit scheme

//double l[Z]  --> Sub-diagonal matrix element computed in the C-N scheme
//double d[Z]  --> Diagonal matrix element computed in the C-N scheme
//double r[Z]  --> Sub-diagonal matrix element computed in the C-N scheme
//double b[Z]  --> Vector to be multiplied with he inverse matrix in the C-N scheme (explicit part)



//*******FUNCTIONSÊPROTOTYPE***********

//*************File management functions*************

void readTable(double* table,char* fileName);
// Read the indicated data file and store it to the given table. Size is not controlled, to avoid error the table should be large enough

void saveTable(double *table, char *name, char* path, int tabSize);
// Save general 1D table containing doubles.

void save2DTable(double **table, char *name, char* path, int nRow, int nCol, int skipR, int skipC, int startC);
// Save general 2D table containing doubles, skipR and skipC allow to save respectively only every "skipR" lines and every "skipC" column, startC gives the startig column to be saved

void save2DTable_top(double **table, char *name, char* path, double* thickness, int nrow, int ncol);
// TO DO



//*************Computational functions*************

void spin_up(double *temperature, double thick, double tsurf,double acc, double QG, double mw,double* tborder, double deltaH,int border,double len, double flat, double* melt, double *density);
//Perform the spin_up of the model for the given time using the CN scheme

void t_solve(double *temperature, int time, double thickness, double thicknessFuture, double tsurf,double acc,double* melt, double QG, double mw, double* tborder,double deltaH, int border,double len, double flat, double* freeze, double *density);
//Getting the 1D array temperature a t-1, return the temperature at T.
//This function calls various function to compute all the needed parameters
//Finally, this function calls the defined algorithm to compute the temperature

void setRho(double* rho, double *rhoIce, double* temp, int thickness,double acc);
//Compute the density profile

void setHeatVar(double *K,double *cp,double *told,int thickness,double *rho, double* rhoIce);
//Compute the values of the K and c thermal variables, called by spin_up() and t_solve()

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

void integrate_CN(double* tint, double* told,double* alpha,double* beta,double* alpha1,double* beta1,double tground,double tsurf,int thickness,int step,double* se);//, double* se);
//Compute the temperature using the CN scheme, called by spin_up() and t_solve()

void integrate_expl(double* told, double* a, double* b, double tground, double tsurf, double tsurf_old, int thickness, double* se);
//Compute the temperature using the explicit scheme, called by  t_solve()

void tempScale(double* told, double thickness,double thicknessFuture,double tsurf);
//Scale the temperature profile to the thickness value of the next step



//*************DEFINTIONÊOFÊTHEÊFUNCTIONS*************

//*************File management functions*************

void readTable(double* table,char* fileName)
{
    // Read a table from the file called "filename" and store it into the double table called "table". Table should be 1D. "filename" can contain directory path.
    FILE *fp;
    int li=0;
    double a=0;
    if((fp=fopen(fileName, "r"))==NULL)
    {
        printf("Cannot open file: %s\n",fileName);
        exit(0);
    }
    else
    {
        printf("File: %s opened in reading mode\n",fileName);
        while(fscanf(fp,"%lf",&a)==1)
        {
            table[li]=a;
            li++;
        }
        if(fclose(fp)==0)
        {
            printf("File imported successfully (%d data) and closed: %s \n\n",li,fileName);
        }
        else
        {
            printf("Not able to close: %s \n\n",fileName);
        }
    }
}


void saveTable(double *table,char *name, char* path,int tabSize)
{
    // Read the 1D double table called "table" and store it a file called "filename".  Table should be tab delimited. "filename" can contain directory path. The default path is a folder called "export". The table size should be passed as parameter.
    FILE *fp;
    int li=0;
    char full_path[180]="";
    sprintf(full_path,"%s/%s",path,name);
    if((fp=fopen(full_path, "w+"))==NULL)
    {
        printf("Cannot open file: %s\n",full_path);
    }
    else
    {
        printf("File oppened:%s\n...writing...\n",full_path);
        for (li=0; li<tabSize; li++)
        {
            fprintf(fp,"%f \n",table[li]);
        }
        fclose(fp);
        printf("File closed: %s \n\n",full_path);
    }
}

void save2DTable(double **table,char *name, char* path,int nRow,int nCol, int skipR, int skipC,int startC)
{
    // Read the 2D double table called "table" and store it a tab delimited file called "filename". The default path is a folder called "export". The table dimensions should be passed as parameter.
    FILE *fp;
    int li=0;
    int co=0;
    char full_path[180]="";
    sprintf(full_path,"%s/%s",path,name);
    if((fp=fopen(full_path, "w+"))==NULL)
    {
        printf("Cannot open file: %s\n",full_path);
    }
    else
    {
        printf("File oppened:%s\n...writing...\n",full_path);
        for (li=0; li<nRow; li+=skipR)
        {
            for (co=startC; co<nCol; co+=skipC)
            {
                if(co!=startC)
                {
                    fprintf(fp,"\t");
                }
                if(table[li][co]>0)
                {
                    fprintf(fp,"%f",table[li][co]);
                }
                else
                {
                    fprintf(fp,"NaN");
                }
            }
            fprintf(fp,"\n");
        }
        fclose(fp);
        printf("File closed: %s \n\n",full_path);
    }
}

//save2DTable_top(density,fileName,path,iceThickness,100,T);

void save2DTable_top(double **table,char *name, char* path,double* thickness, int nRow,int nCol)
{
    // Read the 2D double table called "table" and store it a tab delimited file called "filename". The default path is a folder called "export". The table dimensions should be passed as parameter.
    FILE *fp;
    int li=0;
    int co=0;
    int startC=0;
    char full_path[180]="";
    sprintf(full_path,"%s/%s",path,name);
    if((fp=fopen(full_path, "w+"))==NULL)
    {
        printf("Cannot open file: %s\n",full_path);
    }
    else
    {
        printf("File oppened:%s\n...writing...\n",full_path);
        for (li=0; li<nRow; li+=1)
        {
            for (co=startC; co<nCol; co+=1)
            {
				int real_li=0;
                if(co!=startC)
                {
                    fprintf(fp,"\t");
                    real_li = (int)(thickness[co-1])-li;

                }
                else{
					 real_li = (int)(thickness[co])-li;
				}
                if(table[real_li][co]>0)
                {
                    fprintf(fp,"%f",table[real_li][co]);
                }
                else
                {
                    fprintf(fp,"NaN");
                }
            }
            fprintf(fp,"\n");
        }
        fclose(fp);
        printf("File closed: %s \n\n",full_path);
    }
}



//*************Computational functions*************

void spin_up(double *told, double thick, double tsurf,double acc,double QG,double mw,double* tborder, double deltaH,int border,double len, double flat,double* melt, double *density)
{
    // Perform a spin up for the time indicated in header (time in hyr). The spin up is done with a 2-passes implicit scheme.
    double begin3=omp_get_wtime();
    int thickness=(int)thick;
    int i,li=0;
    double L=333500;
    double rho[Z],rhoIce[Z]= {[0 ... Z-1] = 921};
    double a[Z],b[Z],a2[Z],b2[Z]= {0};
    double m,f,tground=0;
    double K[Z],cp[Z], w[Z], w_def[Z]= {0};
    double delt=31556926.*100.;
    double delz=1.;
    double dhdt=0;
    double se[Z]= {0};
    for(li=0; li <=thickness; li++)
    {
        w_def[li]=wDef ((double) li, (double) thickness,mw);
    }
    for(i=0; i<S; i++)
    {
        double tint[Z],rho_first[Z],se_first[Z]= {0};
        setRho(rho,rhoIce, told, thickness, acc);
        setHeatVar(K, cp, told, thickness,rho,rhoIce);
        computeMelt(&m,&tground,rho,L,K[1],cp[0],told[1],told[0],thick,delz,QG,&f);
        told[0]=tground;
        setABW(a,b,w,cp,K,rho,delt,delz,acc,m,dhdt,w_def,thickness,rhoIce);
        setSe(se,rho,w,cp,K,delt,thickness,told,deltaH,dhdt,tborder, border,len,flat);
        integrate_CN(tint,told, a, b, a, b, tground, tsurf, thickness, 1,se);
        for(li=0; li <=thickness; li++)
        {
            rho_first[li]=rho[li];
        }
        // Second pass
        setRho(rho,rhoIce, tint, thickness, acc);
        setHeatVar(K, cp, tint, thickness,rho,rhoIce);
        computeMelt(&m,&tground,rho_first,L,K[1],cp[0],tint[1],tint[0],thick,delz,QG,&f);
        tint[0]=tground;
        setABW(a2,b2,w,cp,K,rho,delt,delz,acc,m,dhdt,w_def,thickness,rhoIce);
        setSe(se,rho,w,cp,K,delt,thickness,told,deltaH,dhdt,tborder, border,len,flat);
        for(li=0; li <=thickness; li++)
            {
                if(se[li]>0)
                {
                    se[li]=(se[li]+se_first[li])/2;
                }
            }
        integrate_CN(tint,told, a, b, a2, b2, tground,tsurf, thickness, 1,se);
        for(li=0; li <=thickness; li++)
        {
			density[li] = rho[li];
            told[li]=tint[li];
        }
    }
    melt[0]=m;
    printf("Spin up OK -- in %f seconds\n",(double)(omp_get_wtime() - begin3));
}


void t_solve(double *temperature, int time, double thick, double thickFuture, double tsurf, double acc, double* melt, double QG, double mw,double* tborder,double deltaH,int border,double len, double flat,double* freeze, double *density)
{

    int thickness=(int)thick;
    double told[Z]= {0};
    int i,li=0;
    double L=333500;
    double rho[Z],rhoIce[Z]= {[0 ... Z-1] = 921};
    double a[Z],b[Z],a2[Z],b2[Z]= {0};
    double m,f,tground=0;
    double K[Z],cp[Z], w[Z], w_def[Z]= {0};
    double delt=31556926.*100.;
    double delz=1.;
    double dhdt=(thickFuture-thick)/delt;
    double se[Z]= {0};
    
    f=freeze[time-1];

    for(li=0; li <=thickness; li++)
    {
        told[li]=temperature[li];
        w_def[li]=wDef ((double) li, (double) thick,mw);
    }
    double tsurf_old=told[thickness];

    if(strcmp(TYPE,"CN")==0)//C-N scheme
    {
        double rep=1; //Define the number of passes-1 in the C-N scheme
        double tint[Z],rho_first[Z],rho_mean[Z],se_first[Z]= {0};
        double cp0,K1=0;

        setRho(rho,rhoIce, told, thickness, acc);
        setHeatVar(K, cp, told, thickness,rho,rhoIce);
        computeMelt(&m,&tground,rho,L,K[1],cp[0],told[1],told[0],thick,delz,QG,&f);
        setABW(a,b,w,cp,K,rho,delt,delz,acc,m,dhdt,w_def,thickness,rhoIce);
        setSe(se,rho,w,cp,K,delt,thickness,told,deltaH,dhdt,tborder, border,len,flat);
        integrate_CN(tint,told, a, b, a, b, tground, tsurf, thickness, 1,se);
        for(li=0; li <=thickness; li++)
        {
            rho_first[li]=rho[li];
        }
        cp0=cp[0];
        K1=K[1];
        melt[time]=m*31556926.;
        freeze[time]=f;
        f=freeze[time];

        for (i=0; i<(int)rep; i++)
        {
            setRho(rho,rhoIce, tint, thickness, acc);
            setHeatVar(K, cp, tint, thickness,rho,rhoIce);
            for(li=0; li <=thickness; li++)
            {
                rho_mean[li]=(rho[li]+rho_first[li])/2;
            }
            computeMelt(&m,&tground,rho_mean,L,(K[1]+K1)/2,(cp[0]+cp0)/2,(tint[1]+told[1])/2,(tint[0]+told[0])/2,thick,delz,QG,&f);
            tint[0]=tground;
            setABW(a2,b2,w,cp,K,rho,delt,delz,acc,m,dhdt,w_def,thickness,rhoIce);
            setSe(se,rho,w,cp,K,delt,thickness,tint,deltaH,dhdt,tborder, border,len,flat);
            for(li=0; li <=thickness; li++)
            {
                if(se[li]>0)
                {
                    se[li]=(se[li]+se_first[li])/2;
                }
            }
            integrate_CN(tint,told, a, b, a2, b2, tground, tsurf, thickness, 1,se);
        }
        for(li=0; li <=thickness; li++)
        {
            told[li]=tint[li];
        }
        melt[time]+=m*31556926;
        melt[time]/=2;
        freeze[time]+=f;
        freeze[time]/=2;
    }
    else if(strcmp(TYPE,"EXPL")==0) // Explicit scheme
    {
        setRho(rho,rhoIce, told, thickness, acc);
        setHeatVar(K, cp, told, thickness,rho,rhoIce);
        computeMelt(&m,&tground,rho,L,K[1],cp[0],told[1],told[0],thick,delz,QG,&f);
        setABW(a,b,w,cp,K,rho,delt,delz,acc,m,dhdt,w_def,thickness,rhoIce);
        setSe(se,rho,w,cp,K,delt,thickness,told,deltaH,dhdt,tborder, border,len,flat);
        integrate_expl(told, a,b, tground, tsurf,  tsurf_old, thickness, se);
        melt[time]+=m*31556926;
        freeze[time]=f;
    }
    for(li=0; li <=thickness; li++)
    {
		density[li] = rho[li];
        temperature[li]=told[li];
    }
}

void setRho(double* rho, double *rhoIce,double* temp, int thickness,double acc)
{
    int li=0;
    double rhoIceConst=917;
    double rhoSnowConst=350;
    double R=8.3144;
    double k0=11*exp(-10160/(R*temp[thickness]));
    double k1=575*exp(-21400/(R*temp[thickness]));
    acc=acc*31556926.;
    double z55 = 1/(rhoIceConst/1000*k0)*(log(0.55/(rhoIceConst/1000-0.55))-log(rhoSnowConst/(rhoIceConst-rhoSnowConst)));
    double z0[Z]= {0};
    for (li=0; li<=thickness; li++)
    {
        if(strcmp(RHO,"FIRN")==0){
            rhoIce[li]=916.5-0.14438*(temp[li]-271.16)-0.00015175*(temp[li]-273.16)*(temp[li]-273.16);
            if((thickness-li)<z55)
            {
                z0[li]=exp(rhoIce[li]/1000*k0*(thickness-li))*rhoSnowConst/(rhoIce[li]-rhoSnowConst);
                rho[li]=rhoIce[li]*z0[li]/(1+z0[li]);
            }
            else if (li>2000)
            {
                z0[li]=exp(rhoIce[li]/1000*k1*(thickness-li-z55)/sqrt(acc))*0.55/(rhoIce[li]/1000-0.55);
                rho[li]=rhoIce[li]*z0[li]/(1+z0[li]);
            }
            else
            {
                rho[li]=rhoIce[li];
            }
            if(rho[li]>rhoIce[li])
            {
                rho[li]=rhoIce[li];
            }
        }
        else if(strcmp(RHO,"CONST")==0){
            rho[li]=921;
        }
    }
}

void setHeatVar(double *K,double *cp,double *told,int thickness, double *rho, double *rhoIce)
{
    int li=0;
    for(li=0; li<=thickness; li++)
    {
    	//if(strcmp(THERMAL,"CP")==0){
          //  K[li]=9.828*exp(-0.0057*told[li]);
           // cp[li]=152.5 + 7.122*270.4;
        //}
        //else{
          //  K[li]=9.828*exp(-0.0057*told[li]);
           // if(strcmp(THERMAL,"SC")==0){
             //   K[li]=K[li]*pow((rho[li]/rhoIce[li]),2-0.5*rho[li]/rhoIce[li]);
            //}
            //else if (strcmp(THERMAL,"FI")==0){
              //  K[li]=2.*K[li]*rho[li]/(3*rhoIce[li]-rho[li]);
                //K[li]=1.1*K[li];
            //}
            //cp[li]=152.5 + 7.122*told[li];
       // }
        
    	// Ice thermal conductivity
        if(strcmp(THERMAL,"CP")==0){
            K[li]=9.828*exp(-0.0057*told[li]);
            cp[li]=152.5 + 7.122*told[li];
        }
        else if(strcmp(THERMAL,"GO")==0){
            K[li]=2.22*(1-0.0067*(told[li]-273.15));
            cp[li]=152.5 + 7.122*told[li];
        }
      
      	// Firn thermal conductivity
        if(strcmp(FIRN,"SC")==0){
            K[li]=K[li]*pow((rho[li]/rhoIce[li]),2-0.5*rho[li]/rhoIce[li]);
        }
        else if(strcmp(FIRN,"CP")==0){
            K[li]=2.*K[li]*rho[li]/(3*rhoIce[li]-rho[li]);
        }
        
    }
}

void computeMelt(double* m,double* tground,double* rho,double L,double K0,double cp0, double told1,double told0,double thick,double delz,double QG,double* f)
{
    int li=0;
    double tmelt=0;
    double pressure=0;
    //Computation of the pressure and the melting point
    for(li=0; li<=(int)thick; li++)
    {
        pressure+=rho[li];
    }
    pressure+=(thick-(int)thick)*rho[(int)thick];
    tmelt=273.16-7.2*pow(10,-8)*(pressure)*9.8;
    double diff=QG+K0*(told1-tmelt)/delz;

    if(diff>0) //If enough energy is available to melt ice
    {
    	if (strcmp(MELTING,"FREE_MELT")==0){
        	*m= 1/(rho[0]*(L-cp0*(told0-tmelt))+cp0*(tmelt-told1)/2)* (-rho[0]*cp0*(tmelt-told0)/(2.*31556926.*100.) +diff);
	        *tground=tmelt;
	        *f=0;
    	}
    	if (strcmp(MELTING,"FREEZING_NO_ICE")==0){
        	*m= 1/(rho[0]*(L-cp0*(told0-tmelt))+cp0*(tmelt-told1)/2)* (-rho[0]*cp0*(tmelt-told0)/(2.*31556926.*100.) +diff);
	        *tground=tmelt;
	        *f=0;
    	}
    	if (strcmp(MELTING,"FREEZING")==0){
        	*m= 1/(rho[0]*(L-cp0*(told0-tmelt))+cp0*(tmelt-told1)/2)* (-rho[0]*cp0*(tmelt-told0)/(2.*31556926.*100.) +diff);
        	if(*m >= *f/31556926.){
        		*m-=*f/31556926;
        		*f=0;	
        	}
        	else{
        		*f-=*m*31556926;
        		*m=0;
        	}
	        *tground=tmelt;
    	}
    }
    else if(diff<=0) //If not enough energy is available, bottom temperature is decreased
    {
   		if (strcmp(MELTING,"FREE_MELT")==0){
   			*m= 0;
           	*tground=QG*delz/K0+told1;
           	*f=0;
   		}
   		else if (strcmp(MELTING,"FREEZING_NO_ICE")==0){
   			*m= 0;
           	*tground=tmelt;
           	*f=0;
   		}
   		if (strcmp(MELTING,"FREEZING")==0){
   			*m= 0;
           	*tground=tmelt;
   			*f= diff/(rho[0]*L)*31556926;
   		}
            
    }
}

double wDef (double z, double thickness,double mw)
{
    if (strcmp(VERTICAL,"FI")==0){
    	//	if (mw==0.5){
    	//		return ((double)z/ thickness)*sqrt((double)z/ thickness);
    	//		}
    	//	else{
            	return pow(((double)z/ thickness),(1+mw));
    	//	}
    }
    else if (strcmp(VERTICAL,"PA")==0){
        double p=mw;
        return ((z/ thickness)*0+(1.-0)*(1-(p+2)/(p+1)*(1-(z/ thickness))+1/(p+1)*pow(1-(z/ thickness),p+2)));
    }

}

//Compute the matrix element a and b used in explicit and CN scheme and the velocity profile
void setABW(double* a,double* b,double* w,double* cp,double* K,double* rho,double delt,double delz,double acc,double m,double dhdt,double* w_def,int thickness,double* rhoIce)
{
	clock_t t0, t1;
	t0=clock();
    int li=0;
    
    for(li=0; li<=thickness; li++)
    { 
        b[li]=delt*K[li]/(rho[li]*cp[li]*delz*delz);  
        w[li]=-rhoIce[li]/rho[li]*(acc-m-dhdt)*w_def[li]-m;
        a[li]=delt/(delz*2)*(1/(rho[li]*cp[li]*2*delz)*(K[li+1]-K[li-1])-w[li]); 
    }
}

void setSe(double *se,double *rho,double *w, double *cp, double *K, double delt,int thickness, double* told, double dH,double dhdt,double * tborder, int border,double len, double flat)
{
    double P=0;
    int li=0;
    int deltaH=(int)dH;
    //Internal energy production
    if (strcmp(INTERNAL_ENERGY,"ON")==0){
	    for(li=thickness-1; li>=1; li--)
	    {
	        P+=rho[li-1]*9.81;
	        se[li]=0;
	        double cr=cbrt(getDwdz(w,li,thickness)+getDudz((double)li/thickness));
	        se[li]=2*cr*cr*getA(told[li])*delt/(rho[li]*cp[li])+(w[li]*P/rho[li]*(rho[li+1]-rho[li-1])/2)*delt/(rho[li]*cp[li]);
	    }
    }
    else if (strcmp(INTERNAL_ENERGY,"OFF")==0){
    	for(li=thickness-1; li>=1; li--)
	    {
	        se[li]=0;
	    }
    }
    //Valley effect
    if(deltaH>0 && border==0)
    {
        for(li=0; li<=thickness-deltaH; li++)
        {
            se[li+deltaH]+=4*K[li+deltaH]*2*(tborder[li]-told[li+deltaH])*delt/(rho[li+deltaH]*cp[li+deltaH]*len*len);
        }
        for(li=1; li<deltaH; li++)
        {
            double l=((double)li*(double)li/((double)deltaH* (double)deltaH));
            se[li]+=K[li]*2*4*(told[0]+li*6.50E-4-told[li])/((flat+l*(len-flat))*(flat+l*(len-flat))*rho[li]*cp[li])*delt;
        }
    }
}

//Linear piecewise approximation of the horizontal velocity vertical derivative profile squared
double getDudz(double zh)
{
    double us=2.2594e-19;
    double dudz=0;
    if(zh>=0 && zh<=0.05)
    {
        dudz=1.942488E-06 - 1.952373E-05*zh;
    }
    else if(zh<=0.1)
    {
        dudz=1.357489E-06 - 8.359790E-06*zh;
    }
    else if(zh<=0.2)
    {
        dudz=8.170050E-07 - 2.934782E-06*zh;
    }
    else if(zh<=0.3)
    {
        dudz=4.579192E-07 - 1.076580E-06 *zh;
    }
    else if(zh<=0.5)
    {
        dudz=2.656468E-07 - 4.343890E-07*zh;
    }
    else if(zh<=0.75)
    {
        dudz=1.327851E-07 - 1.664704E-07*zh;
    }
    else if(zh<=1)
    {
        dudz=4.236023E-08 - 4.340220E-08*zh;
    }
    return (us*dudz);
}

//Square of the vertical derivative of the vertical velocity profile
double getDwdz(double*w,int z,int thickness)
{
    double dwdz;
    if(z==0)
    {
        dwdz=w[1]-w[0];
    }
    if(z==thickness)
    {
        dwdz=w[thickness]-w[thickness-1];
    }
    else
    {
        dwdz=(w[z+1]-w[z-1])/2;
    }
    return dwdz*dwdz;
}

//Linear piecewise approximation of A power -1/3 value profile
double getA(double t)
{
    double A=0;
    if(t<=210)
    {
        A=31519576971 - 141986135 *t;
    }
    else if(t<=220)
    {
        A=19796924504 - 86163980 *t;
    }
    else if(t<=230)
    {
        A=8301327112 - 33868831*t;
    }
    else if(t<=240)
    {
        A=4956807495 - 19262412*t;
    }
    else if(t<=255)
    {
        A=2755424785 - 10096091*t;
    }
    else if(t<=275)
    {
        A=1802280733 - 6324074*t;
    }
    return A;
}

//Explicit integrations scheme
void integrate_expl(double* told, double* a, double* b, double tground, double tsurf, double tsurf_old, int thickness, double* se)
{
    int li,loop=0;
    double c1[Z]= {0};
    double c2[Z]= {0};
    double c3[Z]= {0};
    for(li=0; li <=thickness; li++)
    {
        c1[li]=(b[li]-a[li])/(365*100);
        c2[li]=-2*b[li]/(365*100)+1;
        c3[li]=(b[li]+a[li])/(365*100);
    }
    double tnew[Z]= {0};
    //internal loop for a daily time step
    for(loop=0; loop<100*365; loop++)
    {
        for(li=1; li<thickness; li++)
        {
            tnew[li]=told[li+1]*c3[li]+told[li]*c2[li]+told[li-1]*c1[li]+se[li]/(365.*100.);
        }
        //Set boundary conditions for the next loop
        tnew[0]=tground;
        tnew[thickness]=tsurf_old+(tsurf-tsurf_old)*(loop+1)/(365.0*100.0);
        for(li=0; li <=thickness; li++)
        {
            told[li]=tnew[li];
        }
    }
}

//CN integrations scheme
void integrate_CN(double* tint, double* told,double* alpha,double* beta,double* alpha1,double* beta1,double tground,double tsurf,int thickness,int step,double* se)
{

    double l[Z]= {0};
    double d[Z]= {0};
    double r[Z]= {0};
    double b[Z]= {0};
    int li,i=0;
    double fact=0.7;
    for(li=1; li <thickness; li++)
    {
        l[li]=fact*(-beta1[li]+alpha1[li]);
        d[li]=fact*2*beta1[li]+1;
        r[li]=fact*(-beta1[li]-alpha1[li]);
        b[li]=told[li-1]*(1-fact)*(beta[li]-alpha[li])+told[li]*(1-(1-fact)*2*beta[li])+told[li+1]*(1-fact)*(beta[li]+alpha[li])+se[li];
        se[li]=se[li];
    }
    b[1]=told[0]*(1-fact)*(beta[1]-alpha[1])+told[1]*(1-(1-fact)*2*beta[1])+told[2]*(1-fact)*(beta[1]+alpha[1])-tground*fact*(-beta1[1]+alpha1[1])+se[1];
    b[thickness-1]=told[thickness-2]*(1-fact)*(beta[thickness-1]-alpha[thickness-1])+told[thickness-1]*(1-(1-fact)*2*beta[thickness-1])+told[thickness]*(1-fact)*(beta[thickness-1]+alpha[thickness-1])-tsurf*fact*(-beta1[thickness-1]-alpha1[thickness-1])+se[thickness-1];
    double dp[Z],bp[Z],x[Z]= {0};
    dp[1]=d[1];
    bp[1]=b[1];
    for(i=1; i<thickness; i++)
    {
        dp[i+1]=d[i+1]-l[i+1]/dp[i]*r[i];
        bp[i+1]=b[i+1]-l[i+1]/dp[i]*bp[i];
    }
    x[thickness-1]=bp[thickness-1]/dp[thickness-1];
    for (i=thickness-2; i>0; i--)
    {
        x[i]=(bp[i]-r[i]*x[i+1])/dp[i];
    }
    tint[0]=tground,
            tint[thickness]=tsurf;
    for (i=1; i<thickness; i++)
    {
        tint[i]=x[i];
    }
    
    
}

//Scale the temperature profile to the next thickness value
void tempScale(double* told, double thick,double thickFuture,double tsurf)
{
    double temperature[Z]= {0};
    int thickness=(int)thick;
    int thicknessFuture=(int) thickFuture;
    int li=0;
    double deltaThick=thicknessFuture-thickness;
    if (deltaThick>0)//If next thickness is bigger, ass some layers at surface temperature
    {
        for(li=0; li <=thickness; li++)
        {
            temperature[li]=told[li];
        }
        for(li=1; li <=deltaThick; li++)
        {
            temperature[li+thickness]=tsurf;
        }
    }
    else if (deltaThick<0)//If the next thickness is smaller, linearly scale the temperature profile
    {
        for(li=0; li <thicknessFuture; li++)
        {
            int oldLi=li*thickness/thicknessFuture;
            int oldLiF=floor(oldLi);
            temperature[li]= told[oldLiF]+(told[oldLiF+1]-told[oldLiF])*(oldLi-oldLiF);
        }
        temperature[thicknessFuture]=told[thickness];
    }
    else
    {
        for(li=0; li <=thickness; li++)
        {
            temperature[li]=told[li];
        }
    }
    for(li=0; li <Z; li++)
    {
        told[li]=temperature[li];
    }
    
}

