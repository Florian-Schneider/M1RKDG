#ifndef _dg_error_h_
#define _dg_error_h_
#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include"data_structure.h"
#include"gauss_pt.h"
#include"dg_monomial_basis.h"
#include"mesh_utils.h"
#include"data_functions.h"
#include"hyperbolic_solver.h"
void l1_dg_error(int,element*,node*,parameters*,double*,double,int,int,double*,double*);
void check_dg_sol_realizability(element*,int,node*,parameters*,double*,double,double*,double*);
void true_solution_val(double,double,double,double,double,int,element*,node*,double*,parameters*,double*);
void check_symmetry(element*,int,node*,edge*,double*,int,int,double,double*);
int find_symmetric_element(element*,int,node*,int);
int find_symmetric_edge(edge*,element*,node*,int,int,int);
void symmetry_correction(element*,node*,double*,int,double,parameters*);
#endif
