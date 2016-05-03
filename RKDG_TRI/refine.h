#ifndef _refine_h_
#define _refine_h_
#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#include"data_structure.h"
#include"mesh_utils.h"
void refine_mesh(element*,node*,edge*,int*,int*);
void reindex_mesh(element*,node*,int,int*);
void get_mid_point(edge*,int,node*,double*);
void set_neighbour(edge*,int,int);
#endif
