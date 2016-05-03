#ifndef _shu_osher_limiter_h_
#define _shu_osher_limiter_h_

#include"data_structure.h"
#include"mesh_utils.h"
#include"gauss_pt.h"
#include"gl_weight.h"
#include"parameters.h"
#include"hyperbolic_solver.h"
#include"mesh_utils.h"
#include"refine.h"
#include"dg_monomial_basis.h"
#include"data_functions.h"
void min_mod_limiter(int,element*,node*,edge*,parameters*,int,int,double*,double*);
void project_solution_to_linear(int,element*,node*,parameters*,double*,double*);
void limiting_linear_basis_functions(double,double,double*);
void quad_sol_linear_phi(int,element*,node*,parameters*,double*,int,double*);
void linear_Aloc_mat(int,element*,node*,parameters*,double*);
void limiting_linear_cell_average(int,element*,node*,double*,double*);
void dg_quadratic_basis_cell_average(int,element*,node*,double*,double*);
void get_barycenter(int,element*,node*,double*);
void get_approx_sol_linear_basis(int,element*,node*,double*,double,double,double*);
void compute_edge_alphas(double*,double*,double*,double*,double*,int,double*,int*);
void get_edge_mid_point(edge*,int,node*,double*);
double min_mod_bar(double,double,parameters*);
double min_mod(double,double);
double min_fun(double,double);
double max_fun(double,double);
double limiting_linear_function_val(double,double,double,double,double,double);
void project_linear_to_dg_quad(int,element*,node*,parameters*,double,double,double,double,double*);
void linear_sol_phi(int,element*,node*,parameters*,double,double,double,double,double*);
void Aloc_mat_linear_mon(int,element*,node*,parameters*,double*);
void project_linear_to_dg_mon(int,element*,node*,parameters*,double,double,double,double,double*);
#endif
