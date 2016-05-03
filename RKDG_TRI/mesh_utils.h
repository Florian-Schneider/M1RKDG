#ifndef _mesh_utils_h_
#define _mesh_utils_h_
#include<stdlib.h>
#include<stdio.h>
#include<string.h>
#include<fstream>
#include<iostream>
#include<math.h>
#include<assert.h>
#include"data_structure.h"
#include"dg_error.h"
void init_mesh_data_structure(element*,node*,edge*,int,int,int*,int*);
void init_zero_d(double*,int);
void init_zero_m(double**,int,int);
void init_zero_int(int*,int);
void new_edge(new_edge_check*,int,int,int*);
void map_to_reference_element(element*,node*,int,double*,double,double);
double get_Element_data(double*,element*,node*,int);
void map_to_physical_element(element*,node*,int,double*,double,double);
void load_mesh_data_structure(FILE*,node*,element*,int,int*,int*,parameters*);
void make_boundary_periodic(edge*,node*,element*,int,parameters*);
int search_periodic_edge(int,edge*,node*,element*,int);
void label_reflective_boundary(int,edge*,node*,parameters*);
void check_symmetry_edges(element*,node*,edge*,int);
void check_symmetry_of_gpts(element*,node*,edge*,int,int);
#endif
