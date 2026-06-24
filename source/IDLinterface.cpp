#include "MWmain.h"

#ifndef LINUX
extern "C" __declspec(dllexport) int GET_MW(int argc, void **argv)
#else
extern "C" double GET_MW(int argc, void **argv)
#endif
{
 int res=0;

 if (argc<5) res=-1;
 else
 {
  int *Lparms=(int*)argv[0];
  double *Rparms=(double*)argv[1];
  double *Sparms=(double*)argv[2];
  double *Bparms=(double*)argv[3];
  double *Flux=(double*)argv[4];

  ComputeMW(Lparms, Rparms, Sparms, Bparms, Flux);
 }

 return res;
}