#include <math.h>
#include <float.h>
#include "ExtMath.h"

double SecantRoot(Integrable1F *F, double x1, double x2, double acc)
{
 double fl, f, dx, swap, xl, rts;

 fl=F->F(x1);
 f=F->F(x2);
 if (fabs(fl)<fabs(f)) 
 { 
  rts=x1; 
  xl=x2;
  swap=fl;
  fl=f;
  f=swap;
 } 
 else 
 {
  xl=x1;
  rts=x2;
 }

 int j=0;
 do
 { 
  dx=(xl-rts)*f/(f-fl); 
  xl=rts;
  fl=f;
  rts+=dx;
  f=F->F(rts);
  j++;
 }
 while (fabs(dx)>acc && f!=0.0 && j<SecantMAXIT);

 return (j<SecantMAXIT) ? rts : dNaN;
}

double BrentRoot(Integrable1F *F, double x1, double x2, double acc)
{
 int iter;
 double a=x1, b=x2, c=x2, d, e, min1, min2;
 double fa=F->F(a), fb=F->F(b), fc, p, q, r, s, tol1, xm;

 if (!finite(fa) || !finite(fb)) return dNaN;

 if (fa*fb>0.0) return dNaN;
 else
 {
  fc=fb;
  for (iter=1; iter<=BrentMAXIT; iter++) 
  {
   if ((fb>0.0 && fc>0.0) || (fb<0.0 && fc<0.0)) 
   {
    c=a; 
    fc=fa;
    e=d=b-a;
   }
   if (fabs(fc)<fabs(fb)) 
   {
    a=b;
    b=c;
    c=a;
    fa=fb;
    fb=fc;
    fc=fa;
   }
   tol1=0.5*acc; 
   xm=0.5*(c-b);
   if (fabs(xm)<=tol1 || fb==0.0) return b;
   if (fabs(e)>=tol1 && fabs(fa)>fabs(fb)) 
   {
    s=fb/fa; 
    if (a==c) 
	{
     p=2.0*xm*s;
     q=1.0-s;
    } 
	else 
	{
     q=fa/fc;
     r=fb/fc;
     p=s*(2.0*xm*q*(q-r)-(b-a)*(r-1.0));
     q=(q-1.0)*(r-1.0)*(s-1.0);
    }
    if (p>0.0) q=-q; 
    p=fabs(p);
    min1=3.0*xm*q-fabs(tol1*q);
    min2=fabs(e*q);
    if (2.0*p<(min1<min2 ? min1 : min2)) 
	{
     e=d; 
     d=p/q; 
    } 
	else 
	{
     d=xm; 
     e=d;
    }
   } 
   else 
   { 
    d=xm;
    e=d;
   }
   a=b; 
   fa=fb;
   if (fabs(d)>tol1) b+=d;
   else b+=(fabs(tol1)*sign(xm));
   fb=F->F(b); 
   if (!finite(fb)) return dNaN;
  }
  return dNaN; 
 }
}

void trapzd2(Integrable2F *F, double a, double b, int N, double *I1, double *I2)
{
 double s1=0, s2=0;
 double x=a;
 double dx=(b-a)/N;
 for (int i=0; i<=N; i++)
 {
  double y1, y2;
  F->F(x, &y1, &y2);
  if ((i==0) || (i==N))
  {
   y1/=2;
   y2/=2;
  }
  s1+=y1;
  s2+=y2;
  x+=dx;
 }
 *I1=s1*dx;
 *I2=s2*dx;
}

void trapzdLog2(Integrable2F *F, double a, double b, int N, double *I1, double *I2)
{
 double s1=0, s2=0;
 double t=log(a);
 double dt=log(b/a)/N;
 for (int i=0; i<=N; i++)
 {
  double y1, y2;
  double x=exp(t);
  F->F(x, &y1, &y2);
  y1*=x;
  y2*=x;
  if ((i==0) || (i==N))
  {
   y1/=2;
   y2/=2;
  }
  s1+=y1;
  s2+=y2;
  t+=dt;
 }
 *I1=s1*dt;
 *I2=s2*dt;
}

double trapzdQ(Integrable1F *F, double a, double b, int n, double *s)
{
 double x, tnm, sum, del;
 int it, j;
 if (n==1)
 {
  return ((*s)=0.5*(b-a)*(F->F(a)+F->F(b)));
 }
 else
 {
  for (it=1, j=1; j<n-1; j++) it<<=1;
  tnm=it;
  del=(b-a)/tnm;
  x=a+0.5*del;
  for (sum=0.0, j=1; j<=it; j++, x+=del) sum+=F->F(x);
  *s=0.5*(*s+(b-a)*sum/tnm);
  return *s;
 }
}

void polint(double *xa, double *ya, int n, double x, double *y, double *dy)
{
 int i, m, ns=1;
 double den, dif, dift, ho, hp, w;
 double c[10], d[10]; //fixed-size arrays; n must be <=10
 dif=fabs(x-xa[1]);
 for (i=1; i<=n; i++)
 {
  if ((dift=fabs(x-xa[i]))<dif)
  {
   ns=i;
   dif=dift;
  }
  c[i]=ya[i];
  d[i]=ya[i];
 }
 *y=ya[ns--];
 for (m=1; m<n; m++)
 {
  for (i=1; i<=n-m; i++)
  {
   ho=xa[i]-x;
   hp=xa[i+m]-x;
   w=c[i+1]-d[i];
   den=ho-hp; 
   den=w/den;
   d[i]=hp*den;
   c[i]=ho*den;
  }
  *y+=(*dy=(2*ns<(n-m) ? c[ns+1] : d[ns--]));
 }
}

#define JMAX 20
#define JMAXP (JMAX+1)
#define romK 6

double qromb(Integrable1F *F, double a, double b, double acc, int *err)
{
 if (err) *err=0;               
 double ts, ss, dss;
 double s[JMAXP], h[JMAXP+1];
 int j;
 h[1]=1.0;
 for (j=1; j<=JMAX; j++)
 {
  s[j]=trapzdQ(F, a, b, j, &ts);
  if (j>=romK)
  {
   polint(&h[j-romK], &s[j-romK], romK, 0.0, &ss, &dss);
   if (fabs(dss)<=acc*fabs(ss) || !finite(ss)) return ss;
  }
  h[j+1]=0.25*h[j];
 }
 if (err) *err=1; 
 return ss;
}

class Integrable1Flog : public Integrable1F
{
 public:
 Integrable1F *oldF=0;
 double F(double t);
};

double Integrable1Flog :: F(double t)
{
 double x=exp(t);
 return oldF->F(x)*x;
}

double qrombLog(Integrable1F *F, double a, double b, double acc, int *err)
{
 Integrable1Flog *ifl=new Integrable1Flog;
 ifl->oldF=F;
 return qromb(ifl, log(a), log(b), acc, err);
 delete ifl;
}

double ErfC(double x)
{
 double t, z, ans;
 z=fabs(x);
 t=1.0/(1.0+0.5*z);
 ans=t*exp(-z*z-1.26551223+t*(1.00002368+t*(0.37409196+t*(0.09678418+
     t*(-0.18628806+t*(0.27886807+t*(-1.13520398+t*(1.48851587+
     t*(-0.82215223+t*0.17087277)))))))));
 return (x>=0.0) ? ans : 2.0-ans;
}

double Erf(double x)
{
 return 1.0-ErfC(x);
}