#pragma once

#define acc_mu_default 1e-3

void FindGSjk(double f, double B, double st, double ct, double N, double T, double L, DF *df, int Nnodes, double *J, double *K);