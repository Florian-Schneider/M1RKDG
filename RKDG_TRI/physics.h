#ifndef _physics_h_
#define _physics_h_
#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include"data_structure.h"
#include"dg_monomial_basis.h"
#include"data_functions.h"
#include"mesh_utils.h"
#include"gauss_pt.h"
#include"hyperbolic_solver.h"

void assemble_source(int,element*,node*,parameters*,double*,int,double,double*);
void scattering(int,double*,element*,node*,parameters*,int,double,double*,double*);
void source_function(double,double*,parameters*,double*);
void init_sigma_vectors(int,double*,element*,node*,parameters*,double,double,double*,double*);
double eval_sigma_s(double,double,double,parameters*);
double eval_sigma_t(double,double,double,parameters*);
#endif
