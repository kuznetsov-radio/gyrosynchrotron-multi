#pragma once

#define me 9.1093837015e-28
#define c 2.99792458e10
#define qe (1.602176634e-19*c/10)
#define eV 1.602176634e-12
#define kB 1.380649e-16
#define AU 1.495978707e13 
#define sfu 1e-19

int FindPlasmaDispersion(double f, double fpe, double fBe, double st, double ct, int sigma, double *N, double *T, double *L, double *Zh);