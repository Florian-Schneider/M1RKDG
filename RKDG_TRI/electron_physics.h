#ifndef _electron_physics_h
#define _electron_physics_h
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
//#include<stdlib.h>
//#include<stdio.h>
void stopping_power_electrons(double*,int,int,int,double*);
void trans_coef_electrons(double*,int,int,int,int,int,double*);
double* data_sp_elec(int,int);
double* data_trans_coef_elec(int,int,int);
#endif

