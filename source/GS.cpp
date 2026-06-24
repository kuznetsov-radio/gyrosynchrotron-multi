#include <math.h>
#include <float.h>
#include "ExtMath.h"
#include "DF.h"
#include "Plasma.h"

class MuSolveFunction : public Integrable1F
{
 public:
 DF *df=0;
 double f, fB, N, T, L, st, ct;
 double E, G, beta;

 int Qcorr;
 double mu_arr[2], lnQ1_arr[2];
 int m=0;

 double F(double mu);
};

double MuSolveFunction :: F(double mu)
{
 double g1=df->g1(mu);

 double sa=sqrt(1.0-mu*mu);
 double x=N*beta*st*sa/(1.0-N*beta*mu*ct);
 double xx=sqrt(1.0-x*x);
 double lnZ=xx+log(x/(1.0+xx));
 double h1mu=2.0*G*f/fB*((N*beta*ct-mu)/(1.0-mu*mu)*xx-N*beta*ct*lnZ)+g1;

 if (Qcorr)
 {
  double s=f/fB*G*(1.0-N*beta*mu*ct);
  double a1=xx*xx*xx+0.503297/s;
  double b1=xx*xx*xx+1.193000/s;
  double b2=1.0-1.0/5/pow(s, 2.0/3);
  double a=pow(a1, 1.0/6);
  double b=pow(b1, 1.0/6)*b2;
  double chi=3.0*x*xx*(mu-N*beta*ct)/sqr(1.0-N*beta*mu*ct)*N*beta*st/sa;
  double eta=f/fB*G*N*beta*ct/(s*s);
  double a1a=(chi+0.503297*eta)/(6.0*a1);
  double b1b=(chi+1.193000*eta)/(6.0*b1)-2.0/15*eta*pow(s, 1.0/3)/b2;
  double lnQ1=2.0*(df->eps*a*b*(N*beta*ct-(a1a+b1b)*(1.0-N*beta*mu*ct))-T*N*beta)/
	          (T*(ct-N*beta*mu)+L*st-df->eps*a*b*(1-N*beta*mu*ct))-2.0*a1a+N*beta*ct/(1.0-N*beta*mu*ct);

  mu_arr[m]=mu;
  lnQ1_arr[m]=lnQ1;
  m^=1;

  h1mu+=lnQ1;
 }

 return h1mu;
}

class GSIntegrandApprox : public Integrable2F
{
 public:
 DF *df=0;
 double f, fB, N, T, L, st, ct;
 double E, G, beta, p;

 MuSolveFunction *msf=0;

 void F(double _E, double *IntJ, double *IntK);
 void Find_mu0(double *mu0, double *lnQ2);
 GSIntegrandApprox();
 ~GSIntegrandApprox();
};

GSIntegrandApprox :: GSIntegrandApprox()
{
 msf=new MuSolveFunction;
}

GSIntegrandApprox :: ~GSIntegrandApprox()
{
 delete msf;
}

void GSIntegrandApprox :: Find_mu0(double *mu0, double *lnQ2)
{
 msf->df=df;
 msf->f=f;
 msf->fB=fB; 
 msf->N=N; 
 msf->T=T;
 msf->L=L;
 msf->st=st;
 msf->ct=ct;
 msf->E=E;
 msf->G=G;
 msf->beta=beta;

 *lnQ2=0;

 msf->Qcorr=0;
 double mu1=BrentRoot(msf, -1.0+1e-5, 1.0-1e-5, df->acc_mu);

 if (finite(mu1))
 {
  msf->Qcorr=1;

  double mu2=SecantRoot(msf, mu1, mu1-1e-4*sign(mu1), df->acc_mu);
  if (finite(mu2) && fabs(mu2)<1.0) mu1=mu2;
  else
  {
   mu2=BrentRoot(msf, -1.0+1e-5, mu1, df->acc_mu);
   if (finite(mu2)) mu1=mu2;
   else mu1=BrentRoot(msf, mu1, 1.0-1e-5, df->acc_mu);
  }
 }

 *mu0=mu1;

 if (finite(mu1) && fabs(mu1)<1.0) *lnQ2=(msf->lnQ1_arr[1]-msf->lnQ1_arr[0])/(msf->mu_arr[1]-msf->mu_arr[0]);
}

void GSIntegrandApprox :: F(double _E, double *IntJ, double *IntK)
{
 *IntJ=*IntK=0;
 E=_E;

 if (E>0.0)
 {
  G=E/(df->m*c*c)+1.0;
  beta=sqrt(sqr(G)-1.0)/G;
  p=beta*G*df->m*c;

  double mu0, lnQ2;
  Find_mu0(&mu0, &lnQ2);

  if (finite(mu0))
  {
   double f0, df0_dE, df0_dmu, g1, g2;
   df->F(E, mu0, &f0, &df0_dE, &df0_dmu, &g1, &g2);

   double s=f/fB*G*(1.0-N*beta*mu0*ct);
   double sa=sqrt(1.0-mu0*mu0);
   double x=N*beta*st*sa/(1.0-N*beta*mu0*ct);
   double xx=sqrt(1.0-x*x);
   double Z=x*exp(xx)/(1.0+xx);
   double Z2s=pow(Z, 2.0*s);
   double a=pow(xx*xx*xx+0.503297/s, 1.0/6);
   double b=pow(xx*xx*xx+1.193000/s, 1.0/6)*(1.0-1.0/5/pow(s, 2.0/3));
   double Q=sqr(T*(ct-N*beta*mu0)+L*st-df->eps*a*b*(1-N*beta*mu0*ct))/(a*a*(1.0-N*beta*mu0*ct));

   double h2mu=-2.0*G*f/fB*(1.0+mu0*mu0-4.0*N*beta*mu0*ct+N*N*beta*beta*(mu0*mu0+ct*ct-st*st))/
	           (sqr(1.0-mu0*mu0)*(1.0-N*beta*mu0*ct)*xx)+lnQ2+g2-g1*g1;
   double LpFactor=(h2mu<0) ? sqrt(-2.0*M_PI/h2mu) : 0;

   double R=df0_dE-(1.0+beta*beta)/(c*p*beta)*f0+(N*beta*ct-mu0)/(c*p*beta)*df0_dmu;

   *IntJ=Z2s*Q*f0*LpFactor;
   *IntK=Z2s*Q*R*LpFactor;
  }
 }
}

void FindGSjk(double f, double B, double st, double ct, double N, double T, double L, DF *df, int Nnodes, double *J, double *K)
{
 GSIntegrandApprox *gsi=new GSIntegrandApprox();
 gsi->df=df;
 gsi->f=f;
 gsi->fB=df->q*B/(2*M_PI*df->m*c);
 gsi->N=N;
 gsi->T=T;
 gsi->L=L;
 gsi->st=st;
 gsi->ct=ct;

 double IntJ=0, IntK=0;
 for (int i=0; i<df->N_intervals; i++)
 {
  double IntJloc, IntKloc;
  if (df->logScale[i]) trapzdLog2(gsi, df->Ex[i], df->Ex[i+1], Nnodes, &IntJloc, &IntKloc);
  else trapzd2(gsi, df->Ex[i], df->Ex[i+1], Nnodes, &IntJloc, &IntKloc);

  IntJ+=IntJloc;
  IntK+=IntKloc;
 }

 *J=2.0*M_PI*sqr(df->q)/c*f/(N*(1.0+T*T)*st*st)*IntJ;
 *K=-2.0*M_PI*sqr(df->q)*c/(N*N*N*f*(1.0+T*T)*st*st)*IntK;

 delete gsi;
}