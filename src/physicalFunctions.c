#include "physicalFunctions.h"

//*************DEFINTION OF THE FUNCTIONS*************

//*************Computational functions*************


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
