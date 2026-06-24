#include <stdlib.h>
#include <math.h>
#include "Plasma.h"
#include "IDLinterface.h"
#include "ExtMath.h"
#include "DF.h"
#include "GS.h"
#include "FF.h"
#include "Messages.h"

void ComputeMW(int *Lparms, double *Rparms, double *Sparms, double *Bparms, double *Flux)
{
 int Nf=max(Lparms[i_Nf], 0);
 double *f=(Nf>0) ? (double*)malloc(sizeof(double)*Nf) : 0;

 int Ncomp=max(Lparms[i_Ncomp], 0);
 int Nnodes=max(Lparms[i_Nnodes], 16);
 int FF_on=(Lparms[i_FFflag]==0);

 double Sang=Rparms[i_S]/(sqr(AU)*sfu);

 double f0=Rparms[i_f0];
 if (f0>0)
 {
  double df=pow(10.0, Rparms[i_df]);
  f[0]=f0;
  for (int i=1; i<Nf; i++) f[i]=f[i-1]*df;
 }
 else for (int i=0; i<Nf; i++) f[i]=Flux[D2(OutSize, 0, i)]*1e9;

 double dz=Sparms[i_dz];

 double nth=Sparms[i_nth];
 double ne_total=nth;
 for (int i=0; i<Ncomp; i++) if (Bparms[D2(BpSize, i_AddN, i)]) ne_total+=Bparms[D2(BpSize, i_nb, i)];
 double fpe=qe*sqrt(ne_total/(M_PI*me));

 double B=Sparms[i_B];
 double fBe=qe*B/(2*M_PI*me*c);

 double T0=Sparms[i_T0];
 double theta=Sparms[i_theta]*M_PI/180;
 double ct=cos(theta);
 double st=sin(theta);
 if (st<0)
 {
  st=-st;
  ct=-ct;
 }

 DF **DF_arr=(Ncomp>0) ? (DF**)malloc(sizeof(DF*)*Ncomp) : 0;

 for (int i=0; i<Ncomp; i++) DF_arr[i]=new DF(Bparms+BpSize*i);

 for (int i=0; i<Nf; i++) 
 {
  Flux[D2(OutSize, 0, i)]=f[i]/1e9;

  for (int sigma=-1; sigma<=1; sigma+=2)
  {
   double I=dNaN;

   double N, T, L, Zh;
   int em_valid=FindPlasmaDispersion(f[i], fpe, fBe, st, ct, sigma, &N, &T, &L, &Zh);
   if (em_valid)
   {
    double J=0, K=0;

    for (int j=0; j<Ncomp; j++)
    {
     double Jloc, Kloc;
     FindGSjk(f[i], B, st, ct, N, T, L, DF_arr[j], Nnodes, &Jloc, &Kloc);
     J+=Jloc;
     K+=Kloc;
    }

    if (FF_on)
    {
     double Jloc, Kloc;
     FindFFjk(f[i], nth, T0, N, Zh, &Jloc, &Kloc);
     J+=Jloc;
     K+=Kloc;
    }

    double Rtau=-K*dz; //-tau = -kappa * dz
    double expRtau=(Rtau<700) ? exp(Rtau) : 0.0;
    I=(Rtau>700) ? 0.0 : ((K==0.0) ? J*dz : J/K*((1.0-expRtau) ? 1.0-expRtau : -Rtau));
   }
  
   Flux[D2(OutSize, (sigma+3)/2, i)]=I*Sang;
  }
 }

 for (int i=0; i<Ncomp; i++) delete DF_arr[i];

 if (DF_arr) free(DF_arr);
 if (f) free(f);
}