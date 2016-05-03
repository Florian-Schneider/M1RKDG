#ifndef _dg_monomial_basis_h_
#define _dg_monomial_basis_h_
#include<stdlib.h>
#include<stdio.h>
#include"data_structure.h"
#include"mesh_utils.h"

void init_monomial_basis(int,element*,node*,int,double,double,double*);
void init_monomial_basis_limit(int,element*,node*,int,double,double,double*);
void init_monomial_basis_deriv(int,element*,node*,int,double,double,double*,double*);
#endif
