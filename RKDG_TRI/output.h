#ifndef _output_h_
#define _output_h_
#include<stdlib.h>
#include<stdio.h>
#include<string.h>
#include<fstream>
#include"data_structure.h"
#include"mesh_utils.h"
#include"hyperbolic_solver.h"
void store_solution(element*,node*,double*,int,int,int,char*,parameters*);
void save_solution_coefficients(double*,element*,node*,parameters*,int,int,double,char*);
void output_sol_at_vec(double*,node*,element*,int,double*,double*,int);
void ctrl_output_line(double*,node*,element*,parameters*,int);
void save_error_coeffs(double*,double*,parameters*,int,int); 
void save_cell_solution_vtp(double*,element*,node*,parameters*,int,int,char*);
void save_point_solution_vtp(double*,element*,node*,parameters*,int,int,char*);
void save_point_solution_error_vtp(double*,element*,node*,parameters*,int,int,double,char*);
void init_point_collector(parameters*);
void init_cell_collector(parameters*);
void write_point_collector(parameters*,double,char*);
void write_cell_collector(parameters*,double,char*);
void end_point_collector(parameters*);
void end_cell_collector(parameters*);
void save_ref_error(double*,double*,node*,node*,element*,element*,voxel*,int,int,int,int,parameters*);
void plot_realizability(double*,element*,node*,parameters*,int,int,double,char*);
#endif
