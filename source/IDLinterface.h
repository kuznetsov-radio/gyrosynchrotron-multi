#pragma once

#define D2(s1, i1, i2)                 ((i1)+(i2)*(s1))
#define D3(s1, s2, i1, i2, i3)         ((i1)+((i2)+(i3)*(s2))*(s1))

#define i_Nf 0
#define i_Ncomp 1
#define i_Nnodes 2
#define i_FFflag 3
#define LpSize 4

#define i_S 0
#define i_f0 1
#define i_df 2
#define RpSize 3

#define i_dz 0
#define i_nth 1
#define i_T0 2
#define i_B 3
#define i_theta 4
#define SpSize 5

#define i_q 0
#define i_m 1
#define i_nb 2
#define i_Eid 3
#define i_Emin 4
#define i_Emax 5
#define i_Ebreak 6
#define i_delta1 7
#define i_delta2 8
#define i_muid 9
#define i_alpha0 10
#define i_dmu 11
#define i_AddN 13
#define BpSize 14

#define OutSize 3