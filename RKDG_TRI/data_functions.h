#ifndef _data_functions_h_
#define _data_functions_h_
#include<math.h>
#include<stdlib.h>
#include<stdio.h>
#include"data_structure.h"
#include"hyperbolic_solver.h"
#define CHAR_EPSILON 1.0e-05
void u_zero_val(double,double,parameters*,double*);
void flux_function(double*,double,double,double,double,parameters*,int,double*);
void bdry_function_val(double,double,double,double,double,int,int,element*,node*,edge*,double*,parameters*,double*);
void advection_v(double,double*,parameters*,int,double*);
void transf_to_char_var(double*,double*,double*,parameters*,double*);
void transform_to_char_var(double*,double*,double*,parameters*,double*,int);
void transf_from_char_var(double*,double*,double*,parameters*,double*);
void ensure_realizability_ortho_proj(double,double,double*);
void ensure_realizability_vert_proj(double,double*);
void ensure_positivity_sol(double,int,element*,node*,int,double*);
void boundary_beam(double,double,double,double,double,double,double,double,double,double,double,double*);
double sign(double);
double estimate_largest_eval(double*,double*,parameters*);
void realizability_hack(double, double*);
#endif
