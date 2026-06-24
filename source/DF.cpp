#include <stdlib.h>
#include <math.h>
#include <float.h>
#include "DF.h"
#include "Plasma.h"
#include "IDLinterface.h"
#include "ExtMath.h"
#include "Messages.h"

DF_energy :: ~DF_energy()
{
 if (Ex) free(Ex);
 if (logScale) free(logScale);
}

//------------------------------------------------------------------

class PLWdf : public DF_energy // PLW=3
{
 double A, delta;
 public:
 void FE(double E, double *F, double *dF_dE);
 PLWdf(double *Bparms);
};

void PLWdf :: FE(double E, double *F, double *dF_dE)
{
 *F=A*pow(E, -delta);
 *dF_dE=(-delta/E)*(*F);
}

PLWdf :: PLWdf(double *Bparms)
{
 double E1=Bparms[i_Emin]*eV*1e6;
 double E2=Bparms[i_Emax]*eV*1e6;
 double nb=Bparms[i_nb];

 delta=Bparms[i_delta1];
 A=nb/(2.0*M_PI)*(delta-1.0)/(pow(E1, 1.0-delta)-pow(E2, 1.0-delta));

 N_intervals=1;
 Ex=(double*)malloc(sizeof(double)*(N_intervals+1));
 Ex[0]=E1;
 Ex[1]=E2;
 logScale=(int*)malloc(sizeof(int)*N_intervals);
 logScale[0]=1;

 valid=finite(A) && (A>0) && (E2>E1);
}

//--------------------------------------------------------------------

class DPLdf : public DF_energy // DPL=4
{
 double E_break, A1, A2, delta1, delta2;
 public:
 void FE(double E, double *F, double *dF_dE);
 DPLdf(double *Bparms);
};

void DPLdf :: FE(double E, double *F, double *dF_dE)
{
 if (E<E_break)
 {
  *F=A1*pow(E, -delta1);
  *dF_dE=(-delta1/E)*(*F);
 }
 else
 {
  *F=A2*pow(E, -delta2);
  *dF_dE=(-delta2/E)*(*F);
 }
}

DPLdf :: DPLdf(double *Bparms)
{
 double E1=Bparms[i_Emin]*eV*1e6;
 double E2=Bparms[i_Emax]*eV*1e6;
 double nb=Bparms[i_nb];

 E_break=Bparms[i_Ebreak]*eV*1e6;
 delta1=Bparms[i_delta1];
 delta2=Bparms[i_delta2];

 A1=nb/(2.0*M_PI)/((pow(E1, 1.0-delta1)-pow(E_break, 1.0-delta1))/(delta1-1.0)+
	pow(E_break, delta2-delta1)*(pow(E_break, 1.0-delta2)-pow(E2, 1.0-delta2))/(delta2-1.0));
 A2=A1*pow(E_break, delta2-delta1);

 N_intervals=2;
 Ex=(double*)malloc(sizeof(double)*(N_intervals+1));
 Ex[0]=E1;
 Ex[1]=E_break;
 Ex[2]=E2;
 logScale=(int*)malloc(sizeof(int)*N_intervals);
 logScale[0]=1;
 logScale[1]=1;

 valid=finite(A1) && (A1>0) && finite(A2) && (A2>0) && (E2>E_break) && (E_break>E1);
}

//----------------------------------------------------------------

class PLFdf : public DF_energy // PLF=13
{
 double A, delta, m;
 public:
 void FE(double E, double *F, double *dF_dE);
 PLFdf(double *Bparms);
};

void PLFdf :: FE(double E, double *F, double *dF_dE)
{
 double G=E/(m*c*c)+1.0;
 double beta=sqrt(sqr(G)-1.0)/G;

 *F=A*pow(E, -delta)/beta;
 *dF_dE=-(*F)*(delta/E+1.0/(beta*beta*G*G*G*m*c*c));
}

class PLFIntegrand : public Integrable1F
{
 public:
 double delta, m;
 double F(double E);
};

double PLFIntegrand :: F(double E)
{
 double G=E/(m*c*c)+1.0;
 double beta=sqrt(sqr(G)-1.0)/G;
 return pow(E, -delta)/beta;
}

PLFdf :: PLFdf(double *Bparms)
{
 double E1=Bparms[i_Emin]*eV*1e6;
 double E2=Bparms[i_Emax]*eV*1e6;
 double nb=Bparms[i_nb];

 delta=Bparms[i_delta1];
 m=Bparms[i_m]*me;

 PLFIntegrand *pi=new PLFIntegrand;
 pi->delta=delta;
 pi->m=m;
 A=nb/(2.0*M_PI)/qrombLog(pi, E1, E2, 1e-6, 0);
 delete pi;

 N_intervals=1;
 Ex=(double*)malloc(sizeof(double)*(N_intervals+1));
 Ex[0]=E1;
 Ex[1]=E2;
 logScale=(int*)malloc(sizeof(int)*N_intervals);
 logScale[0]=1;

 valid=finite(A) && (A>0) && (E2>E1);
}

//----------------------------------------------------------------

class DPFdf : public DF_energy // DPF=14
{
 double E_break, A1, A2, delta1, delta2, m;
 public:
 void FE(double E, double *F, double *dF_dE);
 DPFdf(double *Bparms);
};

void DPFdf :: FE(double E, double *F, double *dF_dE)
{
 double G=E/(m*c*c)+1.0;
 double beta=sqrt(sqr(G)-1.0)/G;

 if (E<E_break)
 {
  *F=A1*pow(E, -delta1)/beta;
  *dF_dE=-(*F)*(delta1/E+1.0/(beta*beta*G*G*G*m*c*c));
 }
 else
 {
  *F=A2*pow(E, -delta2)/beta;
  *dF_dE=-(*F)*(delta2/E+1.0/(beta*beta*G*G*G*m*c*c));
 }
}

DPFdf :: DPFdf(double *Bparms)
{
 double E1=Bparms[i_Emin]*eV*1e6;
 double E2=Bparms[i_Emax]*eV*1e6;
 double nb=Bparms[i_nb];

 E_break=Bparms[i_Ebreak]*eV*1e6;
 delta1=Bparms[i_delta1];
 delta2=Bparms[i_delta2];
 m=Bparms[i_m]*me;

 PLFIntegrand *pi=new PLFIntegrand;
 pi->m=m;
 pi->delta=delta1;
 double I1=qrombLog(pi, E1, E_break, 1e-6, 0);
 pi->delta=delta2;
 double I2=qrombLog(pi, E_break, E2, 1e-6, 0);
 A1=nb/(2.0*M_PI)/(I1+I2*pow(E_break, delta2-delta1));
 A2=A1*pow(E_break, delta2-delta1);
 delete pi;

 N_intervals=2;
 Ex=(double*)malloc(sizeof(double)*(N_intervals+1));
 Ex[0]=E1;
 Ex[1]=E_break;
 Ex[2]=E2;
 logScale=(int*)malloc(sizeof(int)*N_intervals);
 logScale[0]=1;
 logScale[1]=1;

 valid=finite(A1) && (A1>0) && finite(A2) && (A2>0) && (E2>E_break) && (E_break>E1);
}

//----------------------------------------------------------------

class ISOdf : public DF_angle // ISO=0, ISO1=1
{
 public:
 void Fmu(double mu, double *F, double *dF_dmu, double *g1, double *g2);
 double g1(double mu);
 ISOdf();
};

void ISOdf :: Fmu(double mu, double *F, double *dF_dmu, double *g1, double *g2)
{
 *F=0.5;
 *dF_dmu=0.0;
 *g1=0.0;
 *g2=0.0;
}

double ISOdf :: g1(double mu)
{
 return 0;
}

ISOdf :: ISOdf()
{
 acc_mu=1e-3;
 valid=1;
}

//----------------------------------------------------------------

class GAUdf : public DF_angle // GAU=3
{
 double A, mu_c, dmu;
 public:
 void Fmu(double mu, double *F, double *dF_dmu, double *g1, double *g2);
 double g1(double mu);
 GAUdf(double *Bparms);
};

void GAUdf :: Fmu(double mu, double *F, double *dF_dmu, double *g1, double *g2)
{
 double amu=fabs(mu);

 if (amu<mu_c)
 {
  *F=A;
  *dF_dmu=0.0;
  *g1=0.0;
  *g2=0.0;
 }
 else
 {
  *F=A*exp(-sqr((amu-mu_c)/dmu));
  *g1=-2.0*(amu-mu_c)/sqr(dmu)*sign(mu);
  *dF_dmu=(*F)*(*g1);
  *g2=4.0*sqr((amu-mu_c)/sqr(dmu))-2.0/sqr(dmu);
 }
}

double GAUdf :: g1(double mu)
{
 double amu=fabs(mu);
 return (amu<mu_c) ? 0.0 : -2.0*(amu-mu_c)/sqr(dmu)*sign(mu);
}

GAUdf :: GAUdf(double *Bparms)
{
 mu_c=cos(Bparms[i_alpha0]*M_PI/180);
 dmu=Bparms[i_dmu];

 A=0.5/(mu_c+dmu*sqrt(M_PI)/2*Erf((1.0-mu_c)/dmu));

 acc_mu=min(dmu*dmu/30, 1e-3);
 valid=finite(A) && (A>0.0);
}

//------------------------------------------------------------------

class GABdf : public DF_angle // GAB=4
{
 double A, mu0, dmu;
 public:
 void Fmu(double mu, double *F, double *dF_dmu, double *g1, double *g2);
 double g1(double mu);
 GABdf(double *Bparms);
};

void GABdf :: Fmu(double mu, double *F, double *dF_dmu, double *g1, double *g2)
{
 *F=A*exp(-sqr((mu-mu0)/dmu));
 *g1=-2.0*(mu-mu0)/sqr(dmu);
 *dF_dmu=(*F)*(*g1);
 *g2=4.0*sqr((mu-mu0)/sqr(dmu))-2.0/sqr(dmu);
}

double GABdf :: g1(double mu)
{
 return -2.0*(mu-mu0)/sqr(dmu);
}

GABdf :: GABdf(double *Bparms)
{
 mu0=cos(Bparms[i_alpha0]*M_PI/180);
 dmu=Bparms[i_dmu];

 A=2.0/(sqrt(M_PI)*dmu)/(Erf((1.0-mu0)/dmu)+Erf((1.0+mu0)/dmu));

 acc_mu=min(dmu*dmu/30, 1e-3);
 valid=finite(A) && (A>0.0);
}

//----------------------------------------------------------------------

void DF :: F(double E, double mu, double *F, double *dF_dE, double *dF_dmu, double *g1, double *g2)
{
 double f1, df1_dE;
 F1->FE(E, &f1, &df1_dE);

 double f2, df2_dmu;
 F2->Fmu(mu, &f2, &df2_dmu, g1, g2);

 *F=f1*f2;
 *dF_dE=df1_dE*f2;
 *dF_dmu=f1*df2_dmu;
}

double DF :: g1(double mu)
{
 return F2->g1(mu);
}

DF :: DF(double *Bparms)
{
 q=fabs(Bparms[i_q]*qe);
 m=Bparms[i_m]*me;
 eps=-sign(Bparms[i_q]);

 int E_id=(int)Bparms[i_Eid];
 switch(E_id)
 {
  case PLW: F1=new PLWdf(Bparms);
	        break;
  case DPL: F1=new DPLdf(Bparms);
	        break;
  case PLF: F1=new PLFdf(Bparms);
	        break;
  case DPF: F1=new DPFdf(Bparms);
	        break;
  default: F1=new PLWdf(Bparms);
 }

 N_intervals=F1->N_intervals;
 Ex=F1->Ex;
 logScale=F1->logScale;

 int mu_id=(int)Bparms[i_muid];
 switch(mu_id)
 {
  case ISO:
  case ISO1: F2=new ISOdf();
	         break;
  case GAU: F2=new GAUdf(Bparms);
	        break;
  case GAB: F2=new GABdf(Bparms);
	        break;
  default: F2=new ISOdf();
 }

 acc_mu=F2->acc_mu;

 valid=F1->valid && F2->valid;
}

DF :: ~DF()
{
 if (F1) delete F1;
 if (F2) delete F2;
}