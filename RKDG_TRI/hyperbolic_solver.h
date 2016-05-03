#ifndef _hyperbolic_solver_h_
#define _hyperbolic_solver_h_
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
#include"positivity_preserve.h"
#include"dg_error.h"
void rk_dg_time_step(element*,node*,edge*,FILE*,double,int,int,int,int,double,parameters*,double*);
void discrete_pde_rhs(int,int,node*,element*,edge*,parameters*,double*,int,double,double*);
void Alocal_mat(int,element*,node*,parameters*,double*);
void init_initial_conditions(int,node*,element*,edge*,parameters*,int,double*);
void get_approx_solution(double*,int,element*,node*,double,double,double*);
void get_approx_solution_lin(double*,int,element*,node*,double,double,double*);
void get_approx_solution_one_mom(double*,int,element*,node*,double,double,double*);
void flocal_vec1(int,double,int,element*,node*,parameters*,double*,int,double*);
void calc_normal_face(int,node*,edge*,int,double*);
double calc_length_face(edge*,int,node*);
void get_end_points(int,int,node*,element*,edge*,double*);
void get_phys_coords_edge(int,node*,element*,edge*,double,double*);
void clocal_edge_loc(int,int,int,element*,node*,edge*,parameters*,double*,int,double,double*);
void clocal_bdry_edge_loc(int,int,int,element*,node*,edge*,parameters*,double*,int,double,double*);
void lax_friedrichs_flux(edge*,int,double,double*,double*,double*,double,double,parameters*,int,double,double,double,double,double,double*);
void copy_matrix(double*,int,double*);
void shu_osher_limit_scalar(element*,node*,edge*,int,int,double,double,double*,double*);
void upwind_flux(int,int,double*,element*,node*,parameters*,double*,int,double,double*,double*,double*,double*);
void find_barycenter(int,element*,node*,double[2]);
void compute_alphas(double[2],double[2],double[2],double[2],double[2]);
double minmod(double,double,double);
double modified_minmod(double,double,double,double,double,double);
double cell_average_one_mom(double*,int,element*,node*);
void cell_average(double*,int,element*,node*,double*);
void cell_average_lin(double*,int,element*,node*,double*);
void qrdcmp(double**,int,double*,double*,int*);
void qrsolv(double**,int,double[],double[], double[]);
void rsolv(double**,int,double[],double[]);
void assemble_solve_qr(double*,double*,int,double*,double*,double*,double**,int*);
void assemble_qr(double*,int,double*,double*,double**,int*);
void shu_osher_limit(int,element*,node*,edge*,parameters*,int,int,double*,double*);
void normalize_vector(int,double*);
void fractional_step(int,node*,element*,edge*,parameters*,double,double,double*,double*,double,double*,double,double,int);
void sum_vectors(int,double*,double*,double,double*,double,int);
void project_true_sol_to_dg_space(int,node*,element*,edge*,parameters*,int,double*,double*);
void project_quadr_sol_vec_to_lin(int,double*,double*);
void quad_sol_linear_phi(int,element*,node*,parameters*,double*,int,double*);
void limiting_linear_basis_functions(double,double,double*);
void quad_sol_linear_phi(int,element*,node*,parameters*,double*,int,double*);
void linear_Aloc_mat(int,element*,node*,parameters*,double*);
void project_quad_to_linear(int,element*,node*,parameters*,double*,double*);
void project_linear_to_dg_mon(int,element*,node*,parameters*,double,double,double,double,double*);
void linear_sol_phi(int,element*,node*,parameters*,double,double,double,double,double*);
void Aloc_mat_linear_mon(int,element*,node*,parameters*,double*);
double limiting_linear_function_val(double,double,double,double,double,double);
void limiting_linear_basis_functions(double,double,double*);
void dg_sol_linear_phi(int,element*,node*,parameters*,double*,int,double*);
void project_solution_to_linear(int,element*,node*,parameters*,double*,double*);
#endif
