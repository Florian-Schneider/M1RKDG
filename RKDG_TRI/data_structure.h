#ifndef _data_structure_h_
#define _data_structure_h_

#include <fstream>
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_cblas.h>

using namespace std;



#define TRANS_TO_CHAR 1
#define TRANS_FROM_CHAR 2

#define NUM_EDGES 3
#define NUM_NODES 3

#define INTERIOR 0
#define EXTERIOR 1

#define VERTICAL 1
#define HORIZONTAL 2

#define BOTTOM 1
#define TOP 2
#define LEFT 3
#define RIGHT 4

#define TO_BE_INDEXED 50

#define NUM_ELMTS_MAX 1000000
#define NUM_NODES_MAX 5000000
#define NUM_EDGES_MAX 1500000
#define VOXEL_MAX_NUM_ELTS 100000
#define NUM_VOXEL_MAX 10000

#define linear_NLOC 3
#define NLOC 6
#define Nlocmax 21

#define ngpts 7
#define nlgpts 3

#define nfmax 1000
#define nphimax 1000

#define PI 3.1415926535897932384626433832795028841971693993

#define T 1.0
#define EPSILON 1.0e-12
#define EPSILON_R 1.0e-06

#define EPSILON_S 1.0e-12
#define SYSTEM_DIM 3

#define grad_phi_mat(i,j,ldmat) grad_phi_mat[((j)-1) * (ldmat) + (i)-1]
#define inv_BE_T(i,j,ldmat) inv_BE_T[((j)-1) *(ldmat)+(i)-1]
#define Aloc_matrix(i,j,ldmat) Aloc_matrix[((j)-1) *(ldmat)+(i)-1]
#define temp_Aloc(i,j,ldmat) temp_Aloc[((j)-1) *(ldmat)+(i)-1]
#define linear_Aloc_matrix(i,j,ldmat) linear_Aloc_matrix[((j)-1) *(ldmat)+(i)-1]
//element data stucture
typedef struct{
	int vertex[NUM_NODES+1];  
	int domain;
	int degree;
	int edge[4];
	double permeability[5];
	int elmt_index;
	int refined;
	int parent;
	int children[5];
	double det;
	double el_density[ngpts+1];
	double el_gpts_x[ngpts+1];
	double el_gpts_y[ngpts+1];
	double el_gpts_ref_x[ngpts+1];
	double el_gpts_ref_y[ngpts+1];
	double el_gpts_w[ngpts+1];
	double el_gpts_basis[Nlocmax*(ngpts+1)];
    int limited;
    int flag1;
    int flag2;
    //interior points
    double interior_gpts_x[9];
    double interior_gpts_y[9];
    double interior_gpts_ref_x[9];
    double interior_gpts_ref_y[9];
    double interior_gpts_w[9];
    double interior_phi_vals[9*NLOC];
    //volume gpts
    double vol_edge_gpts_x[9];
    double vol_edge_gpts_y[9];
    double vol_edge_gpts_ref_x[9];
    double vol_edge_gpts_ref_y[9];
    double vol_edge_w[9];
    double vol_edge_phi_vals[9*NLOC];
    int symmetric_element;
}element;
typedef struct{
	int vertex[3];
	int neighbour[3];
	int child[3];
	int reftype;
	int edge_type;
	int domain;
	int midpoint;
	int edge_index;
	int refined;
	double length;
	int computed;
	int slope;
	int reflective;
	int periodic_found;
	int periodic_nedge;
	int boundary_side; 
	double ed_density[nlgpts+1];
	double E1_ed_gpts_x[nlgpts+1]; 	  	   //reference coordinates w.r.t E1 (x)
	double E1_ed_gpts_y[nlgpts+1];     	   //reference coordinates w.r.t E1 (y)
	double E2_ed_gpts_x[nlgpts+1]; 	  	   //reference coordinates w.r.t E2 (x)
	double E2_ed_gpts_y[nlgpts+1]; 	  	   //reference coordinates w.r.t E2 (y)
	double ed_phys_coords_x[nlgpts+1]; 	   //physical coordinates (x)
	double ed_phys_coords_y[nlgpts+1]; 	   //physical coordinates (y)
	double ed_gpts_w[nlgpts+1];	           //quadrature weights
	double ed_gpts_basis_E1[Nlocmax*(nlgpts+1)];//values of basis functions w.r.t E1
	double ed_gpts_basis_E2[Nlocmax*(nlgpts+1)];//values of basis functions w.r.t E2
    //limiting gauss points
	double limit_E1_ed_gpts_x[nlgpts+1]; 	  	   //reference coordinates w.r.t E1 (x)
	double limit_E1_ed_gpts_y[nlgpts+1];     	   //reference coordinates w.r.t E1 (y)
	double limit_E2_ed_gpts_x[nlgpts+1]; 	  	   //reference coordinates w.r.t E2 (x)
	double limit_E2_ed_gpts_y[nlgpts+1]; 	  	   //reference coordinates w.r.t E2 (y)
	double limit_ed_phys_coords_x[nlgpts+1]; 	   //physical coordinates (x)
	double limit_ed_phys_coords_y[nlgpts+1]; 	   //physical coordinates (y)
	double limit_ed_gpts_w[nlgpts+1];	           //quadrature weights
	double limit_ed_gpts_basis_E1[Nlocmax*(nlgpts+1)];//values of basis functions w.r.t E1
	double limit_ed_gpts_basis_E2[Nlocmax*(nlgpts+1)];//values of basis functions w.r.t E2
    int marked;
}edge;
typedef struct{
	int closure;
	int test_case;
	double start_time;
	double Ftime;
	int flag_backwards_t_step;
	int flag_store_solution;
	int flag_store_solution_line;
	int flag_store_paraview_solution_cells;
	int flag_store_paraview_solution_points;
	int no_out;
	double nu;
	double mdx2;
	char* meshfilename;
	char* outputfilename;
	char* refsolfilename;
	char* startsolfilename;
	int ref_sol_flag;
	double alpha;
	int p_degree;
	char* densityfilename;
	double min_density;
	int flag_limiting;
	int flag_trans;
	int start_solution_file_flag;
	int custom_output_flag;
	double custom_out_a[2];
	double custom_out_b[2];
	double custom_out_length;
	int periodic;
	double periodic_domain_size_x;
	double periodic_domain_size_y;
	double voxel_domain_size[2];
	int reflective_bc;
	int num_cores;
	double CFL;
        ofstream* point_collector_of;
        ofstream* cell_collector_of;
	double R[nfmax][nphimax][3][3];
	double inv_R[nfmax][nphimax][3][3];
	int nf;
	int nphi;
	double fmax;
	double f[nfmax];
	double phi[nphimax];
	int flag_ensure_positivity;
	double delta_pos;
    int store_grid_solution;
    double mesh_min_diam;
    double mesh_max_edge_length;
    double mesh_min_edge_length;
    int symmetric_solution;
}parameters;
typedef struct{
	double n1[2];
	double n2[2];
	double n3[2];
	double n4[2];
	int num_elements;
	int element_list[VOXEL_MAX_NUM_ELTS];
}voxel;

//node data structure
typedef struct{
	double coord[2];
	int stokes_velocity_index;
	int node_type;
	int stokes_pressure_index;
	int darcy_pressure_index;
	int free_node_index;
}node;

typedef struct{
	double left_node_coord[3];
	double right_node_coord[3];
	int  neighbours[3];
	int stokes_local_index[3];
	int darcy_local_index[3];
	int stokes_edge_index;
	int darcy_edge_index;
	int global_node_index[3];
}interface;

typedef struct{
	int adjacent_nodes[50];
	int connecting_edge[50];
	int adjacent_count;
}new_edge_check;

//tec_plot data structures
typedef struct{
	int vertex[NUM_NODES+1];  
	int edge[NUM_EDGES+1];
	int refined;
}local_element;


typedef struct{
	double coord[2];
}local_node;

typedef struct{
	int vertex[3];
	int neighbour[3];
	int midpoint;
	int refined;
	int child[3];
}local_face;

#endif
