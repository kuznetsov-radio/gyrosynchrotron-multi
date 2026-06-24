#include <math.h>
#include <float.h>
#include "ExtMath.h"

int FindPlasmaDispersion(double f, double fpe, double fBe, double st, double ct, int sigma, double *N, double *T, double *L, double *Zh)
{
 int res=1;

 double u=sqr(fBe/f);
 double v=sqr(fpe/f);

 double D2=sqr(u*sqr(st))+4.0*u*sqr((1.0-v)*ct);
 if (D2<0) res=0;
 else
 {
  double D=sqrt(D2);

  double N2=1.0-2.0*v*(1.0-v)/(2.0*(1.0-v)-u*sqr(st)+sigma*D);
  if (N2<0) res=0;
  else
  {
   *N=sqrt(N2);
   *T=2.0*sqrt(u)*(1.0-v)*ct/(u*sqr(st)-sigma*D);
   *L=(v*sqrt(u)*st+(*T)*u*v*st*ct)/(1.0-u-v+u*v*sqr(ct));
   *Zh=u ? 2.0*(u*sqr(st)+2.0*sqr(1.0-v)-sigma*sqr(u*sqr(st))/D)/sqr(2.0*(1.0-v)-u*sqr(st)+sigma*D) : 1.0;

   res=finite(*N) && finite(*T) && finite(*L) && finite(*Zh);
  }
 }

 return res;
}