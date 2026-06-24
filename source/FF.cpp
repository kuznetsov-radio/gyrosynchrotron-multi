#include <math.h>
#include "ExtMath.h"
#include "Plasma.h"

void FindFFjk(double f, double nth, double T0, double N, double Zh, double *J, double *K)
{
 double Z=(T0>3.5e4) ? 1.146 : 1.0;
 double lnL=(T0<2e5) ? 18.2+1.5*log(T0)-log(f) : 24.573+log(T0/f);

 *K=9.786e-3*sqr(nth)*lnL/(N*sqr(f)*T0*sqrt(T0))*Zh*Z;
 *J=(*K)*sqr(N*f/c)*kB*T0; 
}