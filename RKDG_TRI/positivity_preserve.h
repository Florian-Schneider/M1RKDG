#ifndef _positivity_preserve_h_
#define _positivity_preserve_h_
#ifdef _OPENMP
	#include <omp.h>
#else
	#define omp_get_thread_num() 0
	#define omp_get_num_threads() 1
#endif
#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include"data_structure.h"
#include"dg_monomial_basis.h"
#include"data_functions.h"
#include"mesh_utils.h"
#include"gauss_pt.h"
#include"gl_weight.h"
#include"output.h"
#include"physics.h"
#include"compute_reference_mesh_error.h"
#include<assert.h>
#include"parameters.h"
#include"hyperbolic_solver.h"
void pos_assemble_source(int,element*,node*,parameters*,double*,int,double,double*);
void pos_scattering(int,double*,element*,node*,parameters*,int,double,double*,double*);
void init_limiting_gpts(int,int,element*,edge*,node*);
void pos_clocal_edge_loc(int,int,int,element*,node*,edge*,parameters*,double*,int,double,double*);
void pos_clocal_bdry_edge_loc(int,int,int,element*,node*,edge*,parameters*,double*,int,double,double*);
void pos_flocal_vec1(int,double,int,element*,node*,edge*,parameters*,double*,int,double*);
void pos_project_true_sol_to_dg_space(int,node*,element*,edge*,parameters*,int,double*,double*);
void ensure_realizability(int,node*,element*,edge*,double*,double*,double*,double,double*);
int is_realizable(double*);
void compute_solution_average(int,node*,element*,edge*,double*,double*);
void check_realizability_on_gpts(int,double*,double*,node*,element*,edge*,double,double*);
#endif

