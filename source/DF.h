#pragma once

#define THM 2  //relativistic thermal distribution
#define PLW 3  //power-law distribution on energy
#define DPL 4  //double power-law distribution on energy
#define PLF 13 //power-law distribution of the electron flux on energy
#define DPF 14 //double power-law distribution of the electron flux on energy

#define ISO 0  //isotropic distribution
#define ISO1 1 //same as isotropic
#define GAU 3  //gaussian (on pitch-angle cosine) distribution
#define GAB 4  //directed gaussian beam

class DF_energy
{
 public:
 int N_intervals=0;
 double *Ex=0;
 int *logScale=0;
 int valid;

 virtual void FE(double E, double *F, double *dF_dE)=0;
 ~DF_energy();
};

class DF_angle
{
 public:
 double acc_mu;
 int valid;

 virtual void Fmu(double mu, double *F, double *dF_dmu, double *g1, double *g2)=0;
 virtual double g1(double mu)=0;
 ~DF_angle() {};
};

class DF
{
 DF_energy *F1=0;
 DF_angle *F2=0;

 public:
 double q, m; 
 int eps, valid;

 int N_intervals;
 double *Ex;
 int *logScale;
 double acc_mu;

 void F(double E, double mu, double *F, double *dF_dE, double *dF_dmu, double *g1, double *g2);
 double g1(double mu);
 DF(double *Bparms);
 ~DF();
};