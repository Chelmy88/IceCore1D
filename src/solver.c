#include "solver.h"
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


void spin_up(const model_functions * const functions, double *told, double thick,
            double tsurf,double acc,double QG,double mw,double* tborder,
            double deltaH,int border,double len, double flat,
            double* melt, double *density)
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
        w_def[li]=functions->wDef ((double) li, (double) thickness,mw);
    }
    for(i=0; i<S; i++)
    {
        double tint[Z],rho_first[Z],se_first[Z]= {0};
        functions->setRho(rho,rhoIce, told, thickness, acc);
        setHeatVar(functions, K, cp, told, thickness,rho,rhoIce);
        computeMelt(functions, &m,&tground,rho,L,K[1],cp[0],told[1],told[0],thick,delz,QG,&f);
        told[0]=tground;
        setABW(a,b,w,cp,K,rho,delt,delz,acc,m,dhdt,w_def,thickness,rhoIce);
        setInternal(functions,se,rho,w,cp,K,delt,thickness,told,deltaH,dhdt,tborder, border,len,flat);
        integrate_CN(tint,told, a, b, a, b, tground, tsurf, thickness, 1,se);
        for(li=0; li <=thickness; li++)
        {
            rho_first[li]=rho[li];
        }
        // Second pass
        functions->setRho(rho,rhoIce, tint, thickness, acc);
        setHeatVar(functions,K, cp, tint, thickness,rho,rhoIce);
        computeMelt(functions,&m,&tground,rho_first,L,K[1],cp[0],tint[1],tint[0],thick,delz,QG,&f);
        tint[0]=tground;
        setABW(a2,b2,w,cp,K,rho,delt,delz,acc,m,dhdt,w_def,thickness,rhoIce);
        setInternal(functions,se,rho,w,cp,K,delt,thickness,told,deltaH,dhdt,tborder, border,len,flat);
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


void t_solve(const model_functions * const functions, double *temperature, int time, double thick, double thickFuture, double tsurf, double acc, double* melt, double QG, double mw,double* tborder,double deltaH,int border,double len, double flat,double* freeze, double *density)
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
        w_def[li]=functions->wDef ((double) li, (double) thick,mw);
    }
    double tsurf_old=told[thickness];

    if(strcmp(TYPE,"CN")==0)//C-N scheme
    {
        double rep=1; //Define the number of passes-1 in the C-N scheme
        double tint[Z],rho_first[Z],rho_mean[Z],se_first[Z]= {0};
        double cp0,K1=0;

        functions->setRho(rho,rhoIce, told, thickness, acc);
        setHeatVar(functions,K, cp, told, thickness,rho,rhoIce);
        computeMelt(functions, &m,&tground,rho,L,K[1],cp[0],told[1],told[0],thick,delz,QG,&f);
        setABW(a,b,w,cp,K,rho,delt,delz,acc,m,dhdt,w_def,thickness,rhoIce);
        setInternal(functions,se,rho,w,cp,K,delt,thickness,told,deltaH,dhdt,tborder, border,len,flat);
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
            functions->setRho(rho,rhoIce, tint, thickness, acc);
            setHeatVar(functions,K, cp, tint, thickness,rho,rhoIce);
            for(li=0; li <=thickness; li++)
            {
                rho_mean[li]=(rho[li]+rho_first[li])/2;
            }
            computeMelt(functions,&m,&tground,rho_mean,L,(K[1]+K1)/2,(cp[0]+cp0)/2,(tint[1]+told[1])/2,(tint[0]+told[0])/2,thick,delz,QG,&f);
            tint[0]=tground;
            setABW(a2,b2,w,cp,K,rho,delt,delz,acc,m,dhdt,w_def,thickness,rhoIce);
            setInternal(functions,se,rho,w,cp,K,delt,thickness,tint,deltaH,dhdt,tborder, border,len,flat);
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
        functions->setRho(rho,rhoIce, told, thickness, acc);
        setHeatVar(functions,K, cp, told, thickness,rho,rhoIce);
        computeMelt(functions, &m,&tground,rho,L,K[1],cp[0],told[1],told[0],thick,delz,QG,&f);
        setABW(a,b,w,cp,K,rho,delt,delz,acc,m,dhdt,w_def,thickness,rhoIce);
        setInternal(functions,se,rho,w,cp,K,delt,thickness,told,deltaH,dhdt,tborder, border,len,flat);
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
