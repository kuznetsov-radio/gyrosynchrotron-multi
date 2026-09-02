#include "IDLinterface.h"

#ifndef LINUX
extern "C" __declspec(dllexport) int pyGET_MW(int *Lparms, double *Rparms, double *Sparms, double *Bparms, double *Flux)
#else
extern "C" double pyGET_MW(int *Lparms, double *Rparms, double *Sparms, double *Bparms, double *Flux)
#endif
{
 void *ARGV[5];
 ARGV[0]=(void*)Lparms;
 ARGV[1]=(void*)Rparms;
 ARGV[2]=(void*)Sparms;
 ARGV[3]=(void*)Bparms;
 ARGV[4]=(void*)Flux;

 return GET_MW(5, ARGV);
}