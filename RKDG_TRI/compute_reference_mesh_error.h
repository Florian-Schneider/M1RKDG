#ifndef _compute_reference_mesh_error_h_
#define _compute_reference_mesh_error_h_
#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<string.h>
#include"data_structure.h"
#include"gauss_pt.h"
#include"mesh_utils.h"
#include"data_functions.h"
#include"hyperbolic_solver.h"
void compute_relative_error(double*,double*,node*,node*,element*,element*,voxel*,int,int,int,double*,double*);
void compute_cell_error(double*,double*,int,node*,node*,element*,element*,voxel*,int,int,double*);
void load_sol_coeffs(double*,element*,node*,FILE*,int,int);
void create_voxel(double,element*,node*,double,double,double,double,voxel*,int*);
int check_triangle_in_voxel(int,int,voxel*,node*,element*);
void add_elements_to_voxel(int,int,voxel*,element*,node*,int*);
int check_point_in_voxel(double*,voxel*,int);
int check_point_in_bounding_box(double*,double*);
int check_point_in_element(int,element*,node*,double*);
double det(double*,double*);
void get_bounding_box(int,element*,node*,double*);
double max(double,double);
double min(double,double);
#endif
