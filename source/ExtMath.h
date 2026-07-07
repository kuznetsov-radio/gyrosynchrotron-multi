#pragma once

inline double sqr(double x)
{
 return x*x;
}

inline int min(int a, int b)
{
 return (a<b) ? a : b;
}

inline int max(int a, int b)
{
 return (a>b) ? a : b;
}

inline double max(double a, double b)
{
 return (a>b) ? a : b;
}

inline double min(double a, double b)
{
 return (a<b) ? a : b;
}

inline int sign(double x)
{
 return (x<0.0) ? -1 : 1;
}

#ifndef LINUX
#define finite _finite
#else
#define finite isfinite
#endif

#define dNaN (double(NAN))

class Integrable1F
{
 public:
 virtual double F(double x)=0;
};

#define SecantMAXIT 20
#define BrentMAXIT 100

double SecantRoot(Integrable1F *F, double x1, double x2, double acc);
double BrentRoot(Integrable1F *F, double x1, double x2, double acc);
double qrombLog(Integrable1F *F, double a, double b, double acc, int *err);

class Integrable2F
{
 public:
 virtual void F(double x, double *y1, double *y2)=0;
};

void trapzd2(Integrable2F *F, double a, double b, int N, double *I1, double *I2);
void trapzdLog2(Integrable2F *F, double a, double b, int N, double *I1, double *I2);

double Erf(double);
double ExpBesselK(int, double);