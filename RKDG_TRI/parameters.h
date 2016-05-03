#ifndef _parameters_h_
#define _parameters_h_

#include"data_structure.h"
#include"mesh_utils.h"
#include"gauss_pt.h"
#include"gl_weight.h"
#include"hyperbolic_solver.h"
#include<fstream>
#include<iostream>
void init_parameters(parameters*);
void init_parameters(parameters*,std::ifstream*);
void init_density(element*,node*,edge*,parameters*,int);
void init_gpts(element*,node*,edge*,int,int);
void init_M1_trafo_matrices(parameters*);
#endif
