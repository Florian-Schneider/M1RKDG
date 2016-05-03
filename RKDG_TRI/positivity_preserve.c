#include"positivity_preserve.h"
#include"limiting_gauss_pts.h"
#include<math.h>
#include<complex.h>
void init_limiting_gpts(int Nelts,
                        int global_edge_count,
                        element* mesh_element,
                        edge* mesh_edge,
                        node* mesh_node)
{
    int E=0;
    int edge1,edge2,edge3=0;
    int k=0;
    int E2_edge1;
    int E2_edge2;
    int E2_edge3;
    int E1_edge1;
    int E1_edge2;
    int E1_edge3;



    //interior gpts
    double int_gpts_x[9];
    double int_gpts_y[9];
    double int_gpts_w[9];


    //edge gpts
    double edge1_gpts_x[3];
    double edge1_gpts_y[3];
    
    double edge2_gpts_x[3];
    double edge2_gpts_y[3];

    double edge3_gpts_x[3];
    double edge3_gpts_y[3];

    double edge_w[3];

    double phys_coords[2];
    double ref_coords[2];
    double ref_coords2[2];

    double phi_vec[NLOC];
    int deg=0;
    int idofs=0;
    int j=0;



    init_zero_d(int_gpts_x,9);
    init_zero_d(int_gpts_y,9);

    init_zero_d(edge1_gpts_x,3);
    init_zero_d(edge1_gpts_y,3);

    init_zero_d(edge2_gpts_x,3);
    init_zero_d(edge2_gpts_y,3);

    init_zero_d(edge3_gpts_x,3);
    init_zero_d(edge3_gpts_y,3);

    init_zero_d(phi_vec,NLOC);

    deg = mesh_element[1].degree;

    limiting_gpts(int_gpts_x,int_gpts_y,int_gpts_w,edge1_gpts_x,edge1_gpts_y,edge2_gpts_x,edge2_gpts_y,edge3_gpts_x,edge3_gpts_y,edge_w);

    for(E=1;E<Nelts+1;E++)
    {

        //interior points
        for(k=0;k<9;k++)
        {
            init_zero_d(phys_coords,2);
            map_to_physical_element(mesh_element,mesh_node,E,phys_coords,int_gpts_x[k],int_gpts_y[k]);

            init_zero_d(phi_vec,NLOC);
            init_monomial_basis(E,mesh_element,mesh_node,deg,int_gpts_x[k],int_gpts_y[k],phi_vec);
                       
            //physical coordinates
            mesh_element[E].interior_gpts_x[k] = phys_coords[0];
            mesh_element[E].interior_gpts_y[k] = phys_coords[1];
            //reference coordinates
            mesh_element[E].interior_gpts_ref_x[k] = int_gpts_x[k];
            mesh_element[E].interior_gpts_ref_y[k] = int_gpts_y[k];
            mesh_element[E].interior_gpts_w[k] = 0.5*int_gpts_w[k];
           
            for(j=0;j<NLOC;j++)
            {
                mesh_element[E].interior_phi_vals[k*NLOC+j] = phi_vec[j];
               
            }
        }//loop over gauss points


        //volume pts on edges
        //edge1
        for(k=0;k<3;k++)
        {
            init_zero_d(phys_coords,2);
            map_to_physical_element(mesh_element,mesh_node,E,phys_coords,edge1_gpts_x[k],edge1_gpts_y[k]);

            init_zero_d(phi_vec,NLOC);
            init_monomial_basis(E,mesh_element,mesh_node,deg,edge1_gpts_x[k],edge1_gpts_y[k],phi_vec);
                        
            mesh_element[E].vol_edge_gpts_ref_x[k] = edge1_gpts_x[k];
            mesh_element[E].vol_edge_gpts_ref_y[k] = edge1_gpts_y[k];

            mesh_element[E].vol_edge_gpts_x[k] = phys_coords[0];
            mesh_element[E].vol_edge_gpts_y[k] = phys_coords[1];
            mesh_element[E].vol_edge_w[k] = 0.5*edge_w[k];
            for(j=0;j<NLOC;j++)
            {
                mesh_element[E].vol_edge_phi_vals[k*NLOC+j] = phi_vec[j];
               
            }
        }
       //edge2
        for(k=0;k<3;k++)
        {
            init_zero_d(phys_coords,2);
            map_to_physical_element(mesh_element,mesh_node,E,phys_coords,edge2_gpts_x[k],edge2_gpts_y[k]);

            init_zero_d(phi_vec,NLOC);
            init_monomial_basis(E,mesh_element,mesh_node,deg,edge2_gpts_x[k],edge2_gpts_y[k],phi_vec);
                        
            mesh_element[E].vol_edge_gpts_ref_x[3+k] = edge2_gpts_x[k];
            mesh_element[E].vol_edge_gpts_ref_y[3+k] = edge2_gpts_y[k];

            mesh_element[E].vol_edge_gpts_x[3+k] = phys_coords[0];
            mesh_element[E].vol_edge_gpts_y[3+k] = phys_coords[1];
            mesh_element[E].vol_edge_w[3+k] = 0.5*edge_w[k];
           
            for(j=0;j<NLOC;j++)
            {
                mesh_element[E].vol_edge_phi_vals[(3+k)*NLOC+j] = phi_vec[j];
            }        
        }
        //edge3
        for(k=0;k<3;k++)
        {
            init_zero_d(phys_coords,2);
            map_to_physical_element(mesh_element,mesh_node,E,phys_coords,edge3_gpts_x[k],edge3_gpts_y[k]);

           
            init_zero_d(phi_vec,NLOC);
            init_monomial_basis(E,mesh_element,mesh_node,deg,edge3_gpts_x[k],edge3_gpts_y[k],phi_vec);
                        
            mesh_element[E].vol_edge_gpts_ref_x[6+k] = edge3_gpts_x[k];
            mesh_element[E].vol_edge_gpts_ref_y[6+k] = edge3_gpts_y[k];

            mesh_element[E].vol_edge_gpts_x[6+k] = phys_coords[0];
            mesh_element[E].vol_edge_gpts_y[6+k] = phys_coords[1];
            mesh_element[E].vol_edge_w[6+k] = 0.5*edge_w[k];
           
            for(j=0;j<NLOC;j++)
            {
                mesh_element[E].vol_edge_phi_vals[(6+k)*NLOC+j] = phi_vec[j];
            }        
        }





        edge1 = mesh_element[E].edge[1];  
        E1_edge1 = mesh_edge[edge1].neighbour[1];
        E2_edge1 = mesh_edge[edge1].neighbour[2];

        edge2 = mesh_element[E].edge[2];
        E1_edge2 = mesh_edge[edge2].neighbour[1];
        E2_edge2 = mesh_edge[edge2].neighbour[2];

        edge3 = mesh_element[E].edge[3];
        E1_edge3 = mesh_edge[edge3].neighbour[1];
        E2_edge3 = mesh_edge[edge3].neighbour[2];


        for(k=0;k <3;k++)
        {

            if(mesh_edge[edge1].marked ==0)
            {

                assert(E1_edge1==E);
                //evaluate basis functions at gauss points
                init_zero_d(phi_vec,NLOC);
                init_monomial_basis(E,mesh_element,mesh_node,deg,edge1_gpts_x[k],edge1_gpts_y[k],phi_vec);
                //store values of basis functions
                for(j=0;j<NLOC;j++)
                {
                    mesh_edge[edge1].limit_ed_gpts_basis_E1[k*NLOC+j] = phi_vec[j];
                }
                //edge 1 
                init_zero_d(phys_coords,2);
                map_to_physical_element(mesh_element,mesh_node,E,phys_coords,edge1_gpts_x[k],edge1_gpts_y[k]);

                //reference coords edge1 wrt E1
                mesh_edge[edge1].limit_E1_ed_gpts_x[k] = edge1_gpts_x[k];
                mesh_edge[edge1].limit_E1_ed_gpts_y[k] = edge1_gpts_y[k];
                //physical coordinates edge1
                mesh_edge[edge1].limit_ed_phys_coords_x[k] = phys_coords[0];
                mesh_edge[edge1].limit_ed_phys_coords_y[k] = phys_coords[1];
                //weights
                mesh_edge[edge1].limit_ed_gpts_w[k] = edge_w[k];
                
                if(mesh_edge[edge1].edge_type == INTERIOR)
                {
                    assert(E2_edge1);
                    //map to reference
                    init_zero_d(ref_coords2,2);
                    map_to_reference_element(mesh_element,mesh_node,E2_edge1,ref_coords2,phys_coords[0],phys_coords[1]); 
                    mesh_edge[edge1].limit_E2_ed_gpts_x[k] = ref_coords2[0];
                    mesh_edge[edge1].limit_E2_ed_gpts_y[k] = ref_coords2[1];
                    
                    init_zero_d(phi_vec,NLOC);
                    init_monomial_basis(E2_edge1,mesh_element,mesh_node,deg,ref_coords2[0],ref_coords2[1],phi_vec);
                    //store values of basis functions
                    for(j=0;j<NLOC;j++)
                    {
                        mesh_edge[edge1].limit_ed_gpts_basis_E2[k*NLOC+j] = phi_vec[j];
                    }
                }
            }//edge 1 unmarked

            if(mesh_edge[edge2].marked ==0)
            {
                assert(E1_edge2==E);
                //edge 2
                //evaluate basis functions at gauss points
                init_zero_d(phi_vec,NLOC);
                init_monomial_basis(E,mesh_element,mesh_node,deg,edge2_gpts_x[k],edge2_gpts_y[k],phi_vec);
                
                //store values of basis functions
                for(j=0;j<NLOC;j++)
                {
                    mesh_edge[edge2].limit_ed_gpts_basis_E1[k*NLOC+j] = phi_vec[j];
                }
                
                init_zero_d(phys_coords,2);
                map_to_physical_element(mesh_element,mesh_node,E,phys_coords,edge2_gpts_x[k],edge2_gpts_y[k]);
               
                //reference coords edge2 wrt E1
                mesh_edge[edge2].limit_E1_ed_gpts_x[k] = edge2_gpts_x[k];
                mesh_edge[edge2].limit_E1_ed_gpts_y[k] = edge2_gpts_y[k];
                //physical coordinates edge2
                mesh_edge[edge2].limit_ed_phys_coords_x[k] = phys_coords[0];
                mesh_edge[edge2].limit_ed_phys_coords_y[k] = phys_coords[1];

                
                //weights
                mesh_edge[edge2].limit_ed_gpts_w[k] = edge_w[k];
                
                if(mesh_edge[edge2].edge_type == INTERIOR)
                {
                    assert(E2_edge2);
                    //map to reference
                    init_zero_d(ref_coords2,2);
                    map_to_reference_element(mesh_element,mesh_node,E2_edge2,ref_coords2,phys_coords[0],phys_coords[1]); 
                    mesh_edge[edge2].limit_E2_ed_gpts_x[k] = ref_coords2[0];
                    mesh_edge[edge2].limit_E2_ed_gpts_y[k] = ref_coords2[1];

                    init_zero_d(phi_vec,NLOC);
                    init_monomial_basis(E2_edge2,mesh_element,mesh_node,deg,ref_coords2[0],ref_coords2[1],phi_vec);
                    //store values of basis functions
                    for(j=0;j<NLOC;j++)
                    {
                        mesh_edge[edge2].limit_ed_gpts_basis_E2[k*NLOC+j] = phi_vec[j];
                    }
                }
            }//edge 2 unmarked
            if(mesh_edge[edge3].marked==0)
            {

                assert(E1_edge3==E);
                //edge 3 
                //evaluate basis functions at gauss points
                init_zero_d(phi_vec,NLOC);
                init_monomial_basis(E,mesh_element,mesh_node,deg,edge3_gpts_x[k],edge3_gpts_y[k],phi_vec);
                
                //store values of basis functions
                for(j=0;j<NLOC;j++)
                {
                    mesh_edge[edge3].limit_ed_gpts_basis_E1[k*NLOC+j] = phi_vec[j];
                }
                init_zero_d(phys_coords,2);
                map_to_physical_element(mesh_element,mesh_node,E,phys_coords,edge3_gpts_x[k],edge3_gpts_y[k]);

                //reference coords edge1 wrt E1
                mesh_edge[edge3].limit_E1_ed_gpts_x[k] = edge3_gpts_x[k];
                mesh_edge[edge3].limit_E1_ed_gpts_y[k] = edge3_gpts_y[k];
                //physical coordinates edge1
                mesh_edge[edge3].limit_ed_phys_coords_x[k] = phys_coords[0];
                mesh_edge[edge3].limit_ed_phys_coords_y[k] = phys_coords[1];
                
                //weight
                mesh_edge[edge3].limit_ed_gpts_w[k] = edge_w[k];
                
                if(mesh_edge[edge3].edge_type == INTERIOR)
                {
                    assert(E1_edge3);
                    //map to reference
                    init_zero_d(ref_coords2,2);
                    map_to_reference_element(mesh_element,mesh_node,E2_edge3,ref_coords2,phys_coords[0],phys_coords[1]); 
                    mesh_edge[edge3].limit_E2_ed_gpts_x[k] = ref_coords2[0];
                    mesh_edge[edge3].limit_E2_ed_gpts_y[k] = ref_coords2[1];
                    
                    init_zero_d(phi_vec,NLOC);
                    init_monomial_basis(E2_edge3,mesh_element,mesh_node,deg,ref_coords2[0],ref_coords2[1],phi_vec);
                    //store values of basis functions
                    for(j=0;j<NLOC;j++)
                    {
                        mesh_edge[edge3].limit_ed_gpts_basis_E2[k*NLOC+j] = phi_vec[j];
                    }
                }
            }
        }
        
        mesh_edge[edge1].marked =1;
        mesh_edge[edge2].marked =1;
        mesh_edge[edge3].marked =1;

    }//loop over elements


}//init_limiting_gpts


/* void clocal_bdry_edge_loc 
 * computes the edge integral int_e(f(uh)\cdot n_{E,K}v(x)) for iedge on the boundary
 * in: index to triangle -- E
 * in: degrees of freedom	   -- Nloc
 * in: edge 			   -- iedge
 * in: mesh_element data structure -- mesh_element
 * in: mesh_node data structure    -- mesh_node
 * in: mesh_edge data structure    -- mesh_edge
 * in: parameters 		   -- user_param
 * in: solution at time tn         -- uh_tn_sol_vec
 * in: s_dim_count 	           -- dimension count
 * in: time tn 			   -- t_val
 * out: local_edge rhs vec         -- floc_vec2
 */
 void pos_clocal_bdry_edge_loc(int E,
			               int Nloc,
                           int iedge,
                           element* mesh_element,
                           node* mesh_node,
                           edge* mesh_edge,
			               parameters* user_param,
                           double* uh_tn_sol_vec,
                           int s_dim_count,
			               double t_val,
                           double* floc_vec2)
{
    double* xg = NULL;                  //quadrature points
    double* w  = NULL;                  //weights
    int E1=0;                           
    int E2=0;
    double* phi_vec1 = NULL;            //basis function values at quad points  
    double ref_coords1[2];              //reference coordinates 
    double normal_e[2];                 //unit normal at edge e
    double iedge_length=0.0;            //length of iedge
    double end_points[2];               //end points of physical edge
    double* approx_sol_E1= NULL;        //approximate solution on E1
    double* approx_sol_E2= NULL;        //approximate solution on E2
    double* reflec_vec = NULL;
    double* bdry_val= NULL;             //boundary value
    double numerical_flux=0.0;          //numerical flux value
    int deg=0;                          //polynomial degree of approximation of E
    int k=0;            
    int idofs=0;
    int normal_direction_flag=0;
    double temp_vec[Nloc];
    int p_edge=0;
    double p_edge_phys_coords[2];
    double p_edge_ref_coords[2];

    //allocate space for quad points and weights
    xg = (double*)malloc((nlgpts+1)*sizeof(double));
    w = (double*)malloc((nlgpts+1)*sizeof(double));
    approx_sol_E1 = (double*)malloc(SYSTEM_DIM*sizeof(double));
    approx_sol_E2 = (double*)malloc(SYSTEM_DIM*sizeof(double));
    reflec_vec = (double*)malloc(SYSTEM_DIM*sizeof(double));

    bdry_val = (double*)malloc(SYSTEM_DIM*sizeof(double));

    //get neigbors of edge
    E1 = mesh_edge[iedge].neighbour[1]; 

    //periodic boundary case
    if((*user_param).periodic==1)
    {
        E2 = mesh_edge[iedge].neighbour[2];
        assert(E2>0);
    }

    //attributes of iedge
    iedge_length =mesh_edge[iedge].length; //calc_length_face(mesh_edge,iedge,mesh_node);
    normal_direction_flag=0;
    calc_normal_face(iedge,mesh_node,mesh_edge,normal_direction_flag,normal_e);

    //get polynomial degree
    deg = mesh_element[E1].degree;

    //allocate memory for basis function and derivative
    phi_vec1 = (double*)malloc((Nloc)*sizeof(double));

    init_zero_d(ref_coords1,2);
    init_zero_d(phi_vec1, Nloc);
    init_zero_d(end_points,2);
    init_zero_d(p_edge_phys_coords,2);
    init_zero_d(p_edge_ref_coords,2);

    for(k=0;k<3; k++) 
    {
        get_approx_solution(uh_tn_sol_vec,E1,mesh_element,mesh_node,mesh_edge[iedge].limit_E1_ed_gpts_x[k],mesh_edge[iedge].limit_E1_ed_gpts_y[k],approx_sol_E1);

        if((*user_param).periodic)
        {
            p_edge = mesh_edge[iedge].periodic_nedge;
            
            if(mesh_edge[iedge].boundary_side == BOTTOM)
            {
                p_edge_phys_coords[0] = mesh_edge[iedge].limit_ed_phys_coords_x[k];
                p_edge_phys_coords[1] = mesh_edge[iedge].limit_ed_phys_coords_y[k] + (*user_param).periodic_domain_size_y;
            }
            else if(mesh_edge[iedge].boundary_side == TOP)
            {
                p_edge_phys_coords[0] = mesh_edge[iedge].limit_ed_phys_coords_x[k];
                p_edge_phys_coords[1] = mesh_edge[iedge].limit_ed_phys_coords_y[k] - (*user_param).periodic_domain_size_y;
            }
            else if(mesh_edge[iedge].boundary_side == LEFT)
            {
                p_edge_phys_coords[0] = mesh_edge[iedge].limit_ed_phys_coords_x[k] + (*user_param).periodic_domain_size_x;
                p_edge_phys_coords[1] = mesh_edge[iedge].limit_ed_phys_coords_y[k];
            }
            else if(mesh_edge[iedge].boundary_side == RIGHT)
            {
                p_edge_phys_coords[0] = mesh_edge[iedge].limit_ed_phys_coords_x[k] - (*user_param).periodic_domain_size_x;
                p_edge_phys_coords[1] = mesh_edge[iedge].limit_ed_phys_coords_y[k];
            }
            else
            {
                fprintf(stderr,"Error in phys coords translation \n");
                exit(1);
            }

            map_to_reference_element(mesh_element,mesh_node,E2,p_edge_ref_coords,p_edge_phys_coords[0],p_edge_phys_coords[1]);
            get_approx_solution(uh_tn_sol_vec,E2,mesh_element,mesh_node,p_edge_ref_coords[0],p_edge_ref_coords[1],approx_sol_E2);

            lax_friedrichs_flux(mesh_edge,iedge,mesh_element[E].el_density[k],approx_sol_E1,approx_sol_E2,normal_e,mesh_edge[iedge].limit_E1_ed_gpts_x[k],
                     mesh_edge[iedge].limit_E1_ed_gpts_y[k],user_param,s_dim_count,t_val,p_edge_phys_coords[0],p_edge_phys_coords[1],
                     mesh_edge[iedge].limit_ed_phys_coords_x[k],mesh_edge[iedge].limit_ed_phys_coords_y[k],&numerical_flux);
        }
        else if((*user_param).reflective_bc)
        {
            if(mesh_edge[iedge].reflective)
            {
                //use first element
                reflec_vec[0] = approx_sol_E1[0];
                reflec_vec[1] = approx_sol_E1[1] - 2.0*(normal_e[0]*approx_sol_E1[1] + normal_e[1]*approx_sol_E1[2])*normal_e[0];
                reflec_vec[2] = approx_sol_E1[2] - 2.0*(normal_e[0]*approx_sol_E1[1] + normal_e[1]*approx_sol_E1[2])*normal_e[1];
                lax_friedrichs_flux(mesh_edge,iedge,mesh_element[E].el_density[k],approx_sol_E1,reflec_vec,normal_e,mesh_edge[iedge].limit_E1_ed_gpts_x[k],
                    mesh_edge[iedge].limit_E1_ed_gpts_y[k],user_param,s_dim_count,t_val,0.0,0.0,mesh_edge[iedge].limit_ed_phys_coords_x[k],mesh_edge[iedge].limit_ed_phys_coords_y[k],
                    &numerical_flux);
            }
            else
            {
                //bdry_function_val(mesh_edge[iedge].limit_ed_phys_coords_x[k],mesh_edge[iedge].limit_ed_phys_coords_y[k],mesh_edge[iedge].limit_E1_ed_gpts_x[k], 
                //        mesh_edge[iedge].limit_E1_ed_gpts_y[k],t_val,E,mesh_element,mesh_node,uh_tn_sol_vec,user_param,bdry_val);

		        //use first element
                lax_friedrichs_flux(mesh_edge,iedge,mesh_element[E].el_density[k],approx_sol_E1,bdry_val,normal_e,mesh_edge[iedge].limit_E1_ed_gpts_x[k],mesh_edge[iedge].limit_E1_ed_gpts_y[k],
                        user_param,s_dim_count,t_val,0.0,0.0,mesh_edge[iedge].limit_ed_phys_coords_x[k],mesh_edge[iedge].limit_ed_phys_coords_y[k], &numerical_flux);
            }
        }
        else
        {
            //bdry_function_val(mesh_edge[iedge].limit_ed_phys_coords_x[k],mesh_edge[iedge].limit_ed_phys_coords_y[k],mesh_edge[iedge].limit_E1_ed_gpts_x[k],mesh_edge[iedge].limit_E1_ed_gpts_y[k],
            //        t_val,E,mesh_element,mesh_node,uh_tn_sol_vec,user_param,bdry_val);
            //use first element
            lax_friedrichs_flux(mesh_edge,iedge,mesh_element[E].el_density[k],approx_sol_E1,bdry_val,normal_e,mesh_edge[iedge].limit_E1_ed_gpts_x[k],mesh_edge[iedge].limit_E1_ed_gpts_y[k],
                    user_param,s_dim_count,t_val,0.0,0.0,mesh_edge[iedge].limit_ed_phys_coords_x[k],mesh_edge[iedge].limit_ed_phys_coords_y[k],&numerical_flux);

        }

        for(idofs=0;idofs<Nloc;idofs++)
        {
            floc_vec2[idofs] = floc_vec2[idofs] + mesh_edge[iedge].ed_gpts_w[k+1]*mesh_edge[iedge].limit_ed_gpts_basis_E1[k*Nloc+idofs]*numerical_flux*iedge_length;
        }//loop over idofs

    }//quadrature points
    
    free(phi_vec1);
    free(xg);
    free(w);
    free(approx_sol_E1);
    free(approx_sol_E2);
    free(bdry_val);
    free(reflec_vec);
}//clocal_bdry_edge_loc

/* void clocal_edge_loc 
 * computes the edge integral int_e(f(uh)\cdot n_{E,K}v(x)) for iedge on the interior
 * in: index to triangle -- E
 * in: degrees of freedom	   -- Nloc
 * in: edge -- iedge
 * in: mesh_element data structure -- mesh_element
 * in: mesh_node data structure    -- mesh_node
 * in: mesh_edge data structure    -- mesh_edge
 * in: parameters 		   -- user_param
 * in: solution at time tn         -- uh_tn_sol_vec
 * in: s_dim_count                 -- dimension of system counter
 * in: time tn 			   -- t_val
 * out: local_edge rhs vec         -- floc_vec2
 */
void pos_clocal_edge_loc(int E,
		     int Nloc,
                     int iedge,
                     element* mesh_element,
                     node* mesh_node,
                     edge* mesh_edge,
		     parameters* user_param,
                     double* uh_tn_sol_vec,
                     int s_dim_count,
		     double t_val,
                     double* floc_vec2)
{
    double* xg = NULL;                  //gauss points on [-1,1]
    double* w  = NULL;                  //weights 
    int E1=0;                           //first neighbour of E1
    int E2=0;                           //second neighbour of E2

    double* phi_vec1 = NULL;            //basis function values on edge from E1
    double* phi_vec2 = NULL;            //basis function values on edge from E2

    double ref_coords1[2];              //reference coordinates from E1
    double ref_coords2[2];              //reference coordinates from E2


    double normal_e[2];                 //normal vector points from E1 to E2
    double iedge_length=0.0;            //length of iedge
    double end_points[2];               //end points of iedge

    double* approx_sol_E1= NULL;           //local solution from E1
    double* approx_sol_E2= NULL;           //local solution from E2
    double numerical_flux=0.0;          //numerical flux

    int deg=0;                           //polynomial degree of approx
    int k=0;                            
    int idofs=0;
    int normal_direction_flag=0;

    //allocate space for quad points and weights
    xg = (double*)malloc((nlgpts+1)*sizeof(double));
    w = (double*)malloc((nlgpts+1)*sizeof(double));

    approx_sol_E1 = (double*)malloc((SYSTEM_DIM)*sizeof(double));
    approx_sol_E2 = (double*)malloc((SYSTEM_DIM)*sizeof(double));


    //get neigbors of edge
    E1 = mesh_edge[iedge].neighbour[1];
    E2 = mesh_edge[iedge].neighbour[2];

    assert(E1); 
    assert(E2);

    //mesh_edge[iedge].computed=1;

    if(E==E1) 
    {
	normal_direction_flag =0;
    }
    else if(E==E2) 
    {
	normal_direction_flag =1;
    }

    //get polynomial degree
    deg = mesh_element[E1].degree;

    //allocate memory for basis function and derivative for E1 and E2
    phi_vec1 = (double*)malloc((Nloc)*sizeof(double));

    phi_vec2 = (double*)malloc((Nloc)*sizeof(double));

    init_zero_d(ref_coords1,2);
    init_zero_d(ref_coords2,2);
    init_zero_d(phi_vec2, Nloc);
    init_zero_d(end_points,2);

    //attributes of iedge
    iedge_length = mesh_edge[iedge].length;//calc_length_face(mesh_edge,iedge,mesh_node);
    calc_normal_face(iedge,mesh_node,mesh_edge,normal_direction_flag,normal_e);


    //printf("iedge %d \n",iedge);
    for(k=0;k<3; k++)
    {
        get_approx_solution(uh_tn_sol_vec,E1,mesh_element,mesh_node,mesh_edge[iedge].limit_E1_ed_gpts_x[k],mesh_edge[iedge].limit_E1_ed_gpts_y[k],approx_sol_E1);
        get_approx_solution(uh_tn_sol_vec,E2,mesh_element,mesh_node,mesh_edge[iedge].limit_E2_ed_gpts_x[k],mesh_edge[iedge].limit_E2_ed_gpts_y[k],approx_sol_E2);

        /*printf("phys_coords %lf %lf \n",mesh_edge[iedge].limit_ed_phys_coords_x[k],mesh_edge[iedge].limit_ed_phys_coords_y[k]); 
        printf("k %d,mesh_edge[iedge].E1_ed_gpts_x[k] %lf mesh_edge[iedge].E1_ed_gpts_y[k] %lf \n",k,mesh_edge[iedge].limit_E1_ed_gpts_x[k],mesh_edge[iedge].limit_E1_ed_gpts_y[k]);

         printf("k %d,mesh_edge[iedge].E2_ed_gpts_x[k] %lf mesh_edge[iedge].E2_ed_gpts_y[k] %lf \n",k,mesh_edge[iedge].limit_E2_ed_gpts_x[k],mesh_edge[iedge].limit_E2_ed_gpts_y[k]);
         printf("approx_sol_E1 %lf %lf %lf \n",approx_sol_E1[0],approx_sol_E1[1],approx_sol_E1[2]);
         printf("approx_sol_E2 %lf %lf %lf \n",approx_sol_E2[0],approx_sol_E2[1],approx_sol_E2[2]);*/
         


        if(E==E1)
        {
            lax_friedrichs_flux(mesh_edge,iedge,mesh_element[E].el_density[k],approx_sol_E1,approx_sol_E2,normal_e,mesh_edge[iedge].limit_E1_ed_gpts_x[k],
                    mesh_edge[iedge].limit_E1_ed_gpts_y[k],user_param,s_dim_count,t_val,0.0,0.0,mesh_edge[iedge].limit_ed_phys_coords_x[k],mesh_edge[iedge].limit_ed_phys_coords_y[k],
                    &numerical_flux);
        }
        else if(E==E2)
        {
            lax_friedrichs_flux(mesh_edge,iedge,mesh_element[E].el_density[k],approx_sol_E2,approx_sol_E1,normal_e,mesh_edge[iedge].limit_E2_ed_gpts_x[k],
                    mesh_edge[iedge].limit_E2_ed_gpts_y[k],user_param,s_dim_count,t_val,0.0,0.0,mesh_edge[iedge].limit_ed_phys_coords_x[k],mesh_edge[iedge].limit_ed_phys_coords_y[k],
                    &numerical_flux);
        }

        for(idofs=0;idofs<Nloc;idofs++)
        {
            if(E==E1)
            {
                floc_vec2[idofs] = floc_vec2[idofs] + mesh_edge[iedge].ed_gpts_w[k+1]*mesh_edge[iedge].limit_ed_gpts_basis_E1[k*Nloc+idofs]*numerical_flux*iedge_length;
              // printf("k %d floc_vec2[idofs] %lf w[k] %lf basis value %lf numericalf %lf \n",k,floc_vec2[idofs],mesh_edge[iedge].ed_gpts_w[k+1],mesh_edge[iedge].limit_ed_gpts_basis_E1[k*Nloc+idofs],numerical_flux);
                
            }
            else if(E==E2)
            {
                floc_vec2[idofs] = floc_vec2[idofs] + mesh_edge[iedge].ed_gpts_w[k+1]*mesh_edge[iedge].limit_ed_gpts_basis_E2[k*Nloc+idofs]*numerical_flux*iedge_length;

              // printf("k %d floc_vec2[idofs] %lf w[k] %lf basis value %lf numericalf %lf \n",k,floc_vec2[idofs],mesh_edge[iedge].ed_gpts_w[k+1],mesh_edge[iedge].limit_ed_gpts_basis_E2[k*Nloc+idofs],numerical_flux);

            }

            //getchar();
        }//Nloc
        //getchar();

    }//quadrature points
    //getchar();

    free(xg);
    free(w);
    free(phi_vec1);
    free(phi_vec2);
    free(approx_sol_E1);
    free(approx_sol_E2);
}//clocal_edge_loc

/* void pos_flocal_vec1 
 * computes the volume integeral f(u_h)\cdot \nabla v 
 * in: index to triangle 	   -- E
 * in: time 			   -- t
 * in: degrees of freedom	   -- Nloc
 * in: mesh_element data structure -- mesh_element
 * in: mesh_node data structure    -- mesh_node
 * in: parameters 		   -- user_param
 * in: solution vector at time tn  -- uh_tn_sol_vec
 * in: s_dim_count 	           -- dimension count
 * out: value of the integral      -- floc_vec1
 */
void pos_flocal_vec1(int E,
		             double t,
		             int Nloc,
                     element* mesh_element,
                     node* mesh_node,
                     edge* mesh_edge,
                     parameters* user_param,
                     double* uh_tn_sol_vec,
                     int s_dim_counter,
                     double* floc_vec1)
{
    //local vars
    int k=0;                        //gauss point count
    int deg=0;                      //polynomail degree
    int idofs =0;                   //loop over degrees of freedom
    double* gpx = NULL;             //x coords of gauss points
    double* gpy = NULL;             //y coords of gauss points
    double* w = NULL;               // weights of quadrature rules
    double* phi_vec = NULL;         // values of basis function on reference element
    double* grad_phi_mat = NULL;    // gradient of basis function on reference element
    double det_BE_T=0.0;            // det of refcoord to physical element map also 2|E|
    //scalar case               
    double local_uh_tn_sol[SYSTEM_DIM];    // local solution at time tn
    
    double local_uh_tn_sol1[SYSTEM_DIM];
    double local_uh_tn_sol2[SYSTEM_DIM];

    double flux_local_uh_tn_sol[2];

    gpx = (double*)malloc((ngpts+1)*sizeof(double));
    gpy = (double*)malloc((ngpts+1)*sizeof(double));
    w = (double*)malloc((ngpts+1)*sizeof(double));

    init_zero_d(gpx,ngpts+1);
    init_zero_d(gpy,ngpts+1);
    init_zero_d(w,ngpts+1);

    //gauss points on reference element
    gauss_pt(gpx,gpy,w);

    //get polynomial degree
    deg = mesh_element[E].degree;
    double norm_psi_one =0.0;
    double psi_zero = 0.0;
    double diff =0.0;
    int count1=0;
    int count2=0;


    int edge1,edge2,edge3;
    int E1_edge1,E2_edge1;
    int E1_edge2,E2_edge2;
    int E1_edge3,E2_edge3;

    double numerical_flux =0.0;
    double normal_e[2];

    double test_int =0.0;
    double test_int1=0.0;
    double test_int2=0.0;
    double test_int3=0.0;





    //allocate memory for basis functions
    phi_vec = (double*)malloc(Nloc*sizeof(double));
    grad_phi_mat = (double*)malloc(2*Nloc*sizeof(double));

    init_zero_d(local_uh_tn_sol,SYSTEM_DIM);

    //get mappint from reference element to physical element
    det_BE_T = mesh_element[E].det;

    
    //interior gauss points
    for(k=0;k<9;k++) 
    {
        init_zero_d(grad_phi_mat,2*Nloc);
        init_zero_d(phi_vec,Nloc);
        init_zero_d(flux_local_uh_tn_sol,2);

        init_monomial_basis_deriv(E,mesh_element,mesh_node,deg,mesh_element[E].interior_gpts_ref_x[k],mesh_element[E].interior_gpts_ref_y[k],phi_vec,grad_phi_mat);

        //solution on E at gauss points 
        init_zero_d(local_uh_tn_sol,SYSTEM_DIM);
        get_approx_solution(uh_tn_sol_vec,E,mesh_element,mesh_node,mesh_element[E].interior_gpts_ref_x[k],mesh_element[E].interior_gpts_y[k],local_uh_tn_sol);
        
        //compute flux function

        flux_function(local_uh_tn_sol,mesh_element[E].interior_gpts_x[k],mesh_element[E].interior_gpts_y[k],t,mesh_element[E].el_density[k],user_param,
               s_dim_counter,flux_local_uh_tn_sol);

        for(idofs =1;idofs<Nloc+1;idofs++)
        {
            floc_vec1[idofs-1] = floc_vec1[idofs-1] +mesh_element[E].interior_gpts_w[k]*det_BE_T*((flux_local_uh_tn_sol[0]*grad_phi_mat(1,idofs,2))+
                         	 (flux_local_uh_tn_sol[1]*grad_phi_mat(2,idofs,2)));
               

        }//loop over idofs
    }//loop over quadrature nodes

     
    
    for(k=0;k<9;k++)
    {
        init_zero_d(grad_phi_mat,2*Nloc);
        init_zero_d(phi_vec,Nloc);
        init_zero_d(flux_local_uh_tn_sol,2);

        init_monomial_basis_deriv(E,mesh_element,mesh_node,deg,mesh_element[E].vol_edge_gpts_ref_x[k],mesh_element[E].vol_edge_gpts_ref_y[k],phi_vec,grad_phi_mat);

        //solution on E at gauss points 
        init_zero_d(local_uh_tn_sol,SYSTEM_DIM);
        get_approx_solution(uh_tn_sol_vec,E,mesh_element,mesh_node,mesh_element[E].vol_edge_gpts_ref_x[k],mesh_element[E].vol_edge_gpts_ref_y[k],local_uh_tn_sol);
        
        //compute flux function

        flux_function(local_uh_tn_sol,mesh_element[E].vol_edge_gpts_x[k],mesh_element[E].vol_edge_gpts_y[k],t,mesh_element[E].el_density[k],user_param,
               s_dim_counter,flux_local_uh_tn_sol);


        for(idofs =1;idofs<Nloc+1;idofs++)
        {
            floc_vec1[idofs-1] = floc_vec1[idofs-1] +mesh_element[E].vol_edge_w[k]*det_BE_T*((flux_local_uh_tn_sol[0]*grad_phi_mat(1,idofs,2))+
                         	 (flux_local_uh_tn_sol[1]*grad_phi_mat(2,idofs,2)));
               

        }//loop over idofs
    }

    free(gpx);
    free(gpy);
    free(w);
    free(grad_phi_mat);
    free(phi_vec);
}//flocal_vec1

/* void assemble_source
 * assemble source coefficients for the basis functions
 * on the Gauss points of element E
 * in: element	    		      -- E
 * in: mesh element data structure    -- mesh_element
 * in: mesh_node data structure       -- mesh_node
 * in: parameter structure	      -- user_param
 * in: solution vector                -- uh_tn_sol
 * in: index of moment		      -- s_dim_counter
 * in: time  			      -- t_val
 * out: local source coefficients     -- source_loc_vec
 */
void pos_assemble_source(int E,
		     element* mesh_element,
		     node* mesh_node,
		     parameters* user_param,
		     double* uh_tn_sol,
		     int s_dim_counter,
		     double t_val,
		     double* source_loc_vec)
{
    //local vars
    int k=0;                        //gauss point count
    int deg=0;                      //polynomail degree
    int Nloc=0;                     //degrees of freedoma per element
    int idofs =0;                   //loop over degrees of freedom
    double det_BE_T=0.0;            // det of refcoord to physical element map also 2|E|

    double source_fval[SYSTEM_DIM];
    double phys_coords[2];
    double factor;

    //get polynomial degree
    deg = mesh_element[E].degree;

    Nloc = (deg+1)*(deg+2)/2;
    //allocate memory for basis functions

    init_zero_d(source_loc_vec,Nloc);
    init_zero_d(source_fval,SYSTEM_DIM);

    //get mapping from reference element to physical element
    det_BE_T = mesh_element[E].det;

    for(k=0;k<9;k++) 
    {
        phys_coords[0] = mesh_element[E].interior_gpts_x[k];
        phys_coords[1] = mesh_element[E].interior_gpts_y[k];
        
        source_function(t_val,phys_coords,user_param,source_fval);
        
        factor = mesh_element[E].interior_gpts_w[k]*det_BE_T*source_fval[s_dim_counter];
        for(idofs=0;idofs<Nloc;idofs++) 
        {
            source_loc_vec[idofs] = source_loc_vec[idofs] + factor*mesh_element[E].interior_phi_vals[k*Nloc+idofs];
            /*source_loc_vec[idofs] = source_loc_vec[idofs]+ w[k]*det_BE_T*source_fval[s_dim_counter]*phi_vec[idofs];*/
	        //printf("idofs %d source_loc_val %lf fval %lf sdim %d \n",idofs,source_loc_vec[idofs],source_fval[s_dim_counter],s_dim_counter);
        }//loop over idofs
    }//loop over quadrature nodes

    
    //gauss points on edges
    for(k=0;k<9;k++)
    {
        phys_coords[0] = mesh_element[E].vol_edge_gpts_x[k];
        phys_coords[1] = mesh_element[E].vol_edge_gpts_y[k];
        
        source_function(t_val,phys_coords,user_param,source_fval);
        
        factor = mesh_element[E].vol_edge_w[k]*det_BE_T*source_fval[s_dim_counter];
        for(idofs=0;idofs<Nloc;idofs++) 
        {
            source_loc_vec[idofs] = source_loc_vec[idofs] + factor*mesh_element[E].vol_edge_phi_vals[k*Nloc+idofs];
            /*source_loc_vec[idofs] = source_loc_vec[idofs]+ w[k]*det_BE_T*source_fval[s_dim_counter]*phi_vec[idofs];*/
	        //printf("idofs %d source_loc_val %lf fval %lf sdim %d \n",idofs,source_loc_vec[idofs],source_fval[s_dim_counter],s_dim_counter);
        }//loop over idofs


    }
}

/* void scattering
 * assemble scattering coefficients for the basis functions
 * on the Gauss points of element E
 * in: element	    		      -- E
 * in: solution vector                -- uh_tn_sol_vec
 * in: mesh element data structure    -- mesh_element
 * in: mesh_node data structure       -- mesh_node
 * in: parameter structure	      -- user_param
 * in: index of moment		      -- s_dim_counter
 * in: time  			      -- t_val
 * out: local sigma_t coefficients    -- sigmat_vec
 * out: local sigma_s coefficients    -- sigmas_vec
 */
void pos_scattering(int E,
		double* uh_tn_sol_vec,
		element* mesh_element,
		node* mesh_node,
		parameters* user_param,
		int s_dim_counter,
		double t_val,
		double* sigmat_vec,
		double* sigmas_vec)
{
    //local vars
    int k = 0;                        //gauss point count
    int deg = 0;                      //polynomail degree
    int Nloc = 0;                     //degrees of freedoma per element
    int idofs = 0;                   //loop over degrees of freedom
    double det_BE_T=0.0;            // det of refcoord to physical element map also 2|E|
    //scalar case               
    double loc_sigma_s_vec[SYSTEM_DIM];
    double loc_sigma_t_vec[SYSTEM_DIM];
    double phys_coords[2];

    double sigma_s = 0.0;
    double sigma_t = 0.0;

    double factors, factort;

    //get polynomial degree
    deg = mesh_element[E].degree;

    Nloc = (deg+1)*(deg+2)/2;

    init_zero_d(loc_sigma_s_vec,SYSTEM_DIM);
    init_zero_d(loc_sigma_t_vec,SYSTEM_DIM);
    init_zero_d(sigmas_vec,Nloc);
    init_zero_d(sigmat_vec,Nloc);

    //get mapping from reference element to physical element
    det_BE_T = mesh_element[E].det;

    //interior points
    for(k=0;k<9;k++) 
    {
        //initialize vectors
        init_sigma_vectors(E,uh_tn_sol_vec,mesh_element,mesh_node,user_param,mesh_element[E].interior_gpts_ref_x[k],mesh_element[E].interior_gpts_ref_y[k],loc_sigma_s_vec,loc_sigma_t_vec);

        //evaluate scattering coeff at the gauss points
        sigma_s = eval_sigma_s(mesh_element[E].interior_gpts_x[k],mesh_element[E].interior_gpts_y[k],t_val,user_param); // Scattering cross section
        sigma_t = eval_sigma_t(mesh_element[E].interior_gpts_x[k],mesh_element[E].interior_gpts_y[k],t_val,user_param); // Total (scattering + absorption) cross section

        // Henyey-Greenstein model
        factors = 1.0/(4.0*PI)*sigma_s*det_BE_T*loc_sigma_s_vec[s_dim_counter]*mesh_element[E].el_density[k];
        factort = -1.0        *sigma_t*det_BE_T*loc_sigma_t_vec[s_dim_counter]*mesh_element[E].el_density[k];

        for(idofs=0;idofs<Nloc;idofs++) 
        {
            sigmas_vec[idofs] = sigmas_vec[idofs] + factors*mesh_element[E].interior_gpts_w[k]*mesh_element[E].interior_phi_vals[k*Nloc+idofs];
            sigmat_vec[idofs] = sigmat_vec[idofs] + factort*mesh_element[E].interior_gpts_w[k]*mesh_element[E].interior_phi_vals[k*Nloc+idofs];
            /*[>sigmas_vec[idofs] = sigmas_vec[idofs] + 1.0/(2.0)*sigma_s*w[k]*det_BE_T*loc_sigma_s_vec[s_dim_counter]*phi_vec[idofs];<]*/
        }//loop over idofs
    }//loop over quadrature nodes

    //edges
    for(k=0;k<9;k++)
    {
        //initialize vectors
        init_sigma_vectors(E,uh_tn_sol_vec,mesh_element,mesh_node,user_param,mesh_element[E].vol_edge_gpts_ref_x[k],mesh_element[E].vol_edge_gpts_ref_y[k],loc_sigma_s_vec,loc_sigma_t_vec);

        //evaluate scattering coeff at the gauss points
        sigma_s = eval_sigma_s(mesh_element[E].vol_edge_gpts_x[k],mesh_element[E].vol_edge_gpts_y[k],t_val,user_param); // Scattering cross section
        sigma_t = eval_sigma_t(mesh_element[E].vol_edge_gpts_x[k],mesh_element[E].vol_edge_gpts_y[k],t_val,user_param); // Total (scattering + absorption) cross section

        // Henyey-Greenstein model
        factors = 1.0/(4.0*PI)*sigma_s*det_BE_T*loc_sigma_s_vec[s_dim_counter]*mesh_element[E].el_density[k];
        factort = -1.0        *sigma_t*det_BE_T*loc_sigma_t_vec[s_dim_counter]*mesh_element[E].el_density[k];

        for(idofs=0;idofs<Nloc;idofs++) 
        {
            sigmas_vec[idofs] = sigmas_vec[idofs] + factors*mesh_element[E].vol_edge_w[k]*mesh_element[E].vol_edge_phi_vals[k*Nloc+idofs];
            sigmat_vec[idofs] = sigmat_vec[idofs] + factort*mesh_element[E].vol_edge_w[k]*mesh_element[E].vol_edge_phi_vals[k*Nloc+idofs];
            /*[>sigmas_vec[idofs] = sigmas_vec[idofs] + 1.0/(2.0)*sigma_s*w[k]*det_BE_T*loc_sigma_s_vec[s_dim_counter]*phi_vec[idofs];<]*/
        }//loop over idofs
    }


}

/* void project_true_solution_to_dg_space
 * sets up initial conditions
 * in: index to triangle -- E
 * in: mesh_node data structure    -- mesh_node
 * in: mesh_element data structure -- mesh_element
 * in: mesh_edge data structure    -- mesh_edge
 * in: parameters 		   -- user_param
 * in: s_dim_count 	           -- dimension count
 * out: local initional conditions -- u_zero_loc
 */
void  pos_project_true_sol_to_dg_space(int E,
                             node* mesh_node,
                             element* mesh_element,
                             edge* mesh_edge,
                             parameters* user_param,
                             int s_dim_count,
                             double* uh_tn_sol_vec,
                             double* vec_uval)
{
    int k=0;
    int idofs=0;

    int deg=0;
    int Nloc=0;

    double det_BE_T =0.0;

    double u_val[SYSTEM_DIM];

    det_BE_T = mesh_element[E].det;

    deg = mesh_element[E].degree;
    Nloc = ((deg+1)*(deg+2))/2;

    //interior points
    for(k=0;k<9;k++) 
    {
        //bdry_function_val(mesh_element[E].interior_gpts_x[k],mesh_element[E].interior_gpts_y[k],mesh_element[E].interior_gpts_ref_x[k],mesh_element[E].interior_gpts_ref_y[k],
        //        (*user_param).Ftime,E,mesh_element,mesh_node,uh_tn_sol_vec,user_param,u_val);
       // u_val[s_dim_count] =1.0;

        for(idofs=0;idofs<Nloc;idofs++) 
        {
            vec_uval[idofs] = vec_uval[idofs] + u_val[s_dim_count]*mesh_element[E].interior_phi_vals[k*Nloc+idofs]*det_BE_T*mesh_element[E].interior_gpts_w[k];

            //vec_uval[idofs] = vec_uval[idofs] + u_val[s_dim_count]*det_BE_T*mesh_element[E].interior_gpts_w[k];
            //printf("u_val %lf phi_vals %lf det %lf w %lf \n",u_val[s_dim_count],mesh_element[E].interior_phi_vals[k*Nloc+idofs],det_BE_T,mesh_element[E].interior_gpts_w[k]); 
        }//idofs
    }//ngpts
    getchar();

    //edges
    for(k=0;k<9;k++)
    {
        //bdry_function_val(mesh_element[E].vol_edge_gpts_x[k],mesh_element[E].vol_edge_gpts_y[k],mesh_element[E].vol_edge_gpts_ref_x[k],mesh_element[E].vol_edge_gpts_ref_y[k],
        //        (*user_param).Ftime,E,mesh_element,mesh_node,uh_tn_sol_vec,user_param,u_val);

        //u_val[s_dim_count] =1.0;
        for(idofs=0;idofs<Nloc;idofs++) 
        {
            vec_uval[idofs] = vec_uval[idofs] + u_val[s_dim_count]*mesh_element[E].vol_edge_phi_vals[k*Nloc+idofs]*det_BE_T*mesh_element[E].vol_edge_w[k];

            //vec_uval[idofs] = vec_uval[idofs] + u_val[s_dim_count]*det_BE_T*mesh_element[E].vol_edge_w[k];
            //printf("uval %lf edge phi_val %lf det %lf w %lf \n",u_val[s_dim_count],mesh_element[E].vol_edge_phi_vals[k*Nloc+idofs],det_BE_T,mesh_element[E].vol_edge_w[k]);
            //getchar();
        }//idofs


    }
}//initial conditions

void check_realizability_on_gpts(int E,
                            double* U_bar,
                            double* uh_tn_sol_vec,
                            node* mesh_node,
                            element* mesh_element,
                            edge* mesh_edge,
                            double epsilon,
                            double* U_hat)
{
    int idofs=0;
    int k=0;
    double approx_sol[SYSTEM_DIM];
    double U_s[SYSTEM_DIM];
    int rel_flag =0;
    double max_theta_s =0.0;
    double theta_s =0.0;

    int s_dim=0;
    int i=0;

    
    
    //interior points
   for(k=0;k<9;k++)
   { 
       init_zero_d(U_s,SYSTEM_DIM);
        get_approx_solution(uh_tn_sol_vec,E,mesh_element,mesh_node,mesh_element[E].interior_gpts_ref_x[k],mesh_element[E].interior_gpts_ref_y[k],U_s);
        rel_flag =is_realizable(U_s); 
        //printf("rel_flag %d \n",rel_flag); 
        //printf("U_bar %lf %lf %lf >>>>>>>>>>>>>>>>>>>>>> \n",U_bar[0],U_bar[1],U_bar[2]);
        if(rel_flag ==0)
        {
            //printf("rel_limit on interior \n");
            ensure_realizability(E,mesh_node,mesh_element,mesh_edge,U_bar,U_s,U_hat,epsilon,&theta_s);
            if(theta_s > max_theta_s)
            {
                max_theta_s = theta_s;
            }
        }
            
    }

    //edges 
    for(k=0;k<9;k++)
    {
        init_zero_d(U_s,SYSTEM_DIM);
        get_approx_solution(uh_tn_sol_vec,E,mesh_element,mesh_node,mesh_element[E].vol_edge_gpts_ref_x[k],mesh_element[E].vol_edge_gpts_ref_y[k],U_s);
 
        rel_flag =is_realizable(U_s);
        //printf("rel_flag in edge %d \n",rel_flag);
        //printf("U_s %10.16e %10.16e %10.16e >>>>>>>>>>>>>U_s flag  %d \n",U_s[0],U_s[1],U_s[2],rel_flag);
        //printf("U_bar %lf %lf %lf >>>>>>>>>>>>>>>>>>>>>> \n",U_bar[0],U_bar[1],U_bar[2]);
        if(rel_flag ==0)
        {
            //printf("rel_limit on edges \n");
            //printf("U_s %10.16e %10.16e %10.16e >>>>>>>>>>>>>U_s flag  %d \n",U_s[0],U_s[1],U_s[2],rel_flag);
            ensure_realizability(E,mesh_node,mesh_element,mesh_edge,U_bar,U_s,U_hat,epsilon,&theta_s);
            if(theta_s > max_theta_s)
            {
                max_theta_s = theta_s;
            }
        }       
    }

    //printf("max_theta_s  %10.16e \n",max_theta_s);

    if(max_theta_s > 0.0)
    {
       // printf("E limited %d \n",E);
        for(s_dim =0;s_dim < SYSTEM_DIM;s_dim++)
        {
            for(i=0;i<NLOC;i++)
            {
                if(i==0)
                {
                   uh_tn_sol_vec[SYSTEM_DIM*NLOC*(E-1) + (s_dim*NLOC)+i] = max_theta_s*U_bar[s_dim] + (1.0- max_theta_s)*uh_tn_sol_vec[SYSTEM_DIM*NLOC*(E-1) + (s_dim*NLOC)+i];
                }
                else
                {

                   uh_tn_sol_vec[SYSTEM_DIM*NLOC*(E-1) + (s_dim*NLOC)+i] = (1.0- max_theta_s)*uh_tn_sol_vec[SYSTEM_DIM*NLOC*(E-1) + (s_dim*NLOC)+i];
                }
            }
        }
    }
    
    
    for(k=0;k<9;k++)
   { 
       init_zero_d(U_s,SYSTEM_DIM);
        get_approx_solution(uh_tn_sol_vec,E,mesh_element,mesh_node,mesh_element[E].interior_gpts_ref_x[k],mesh_element[E].interior_gpts_ref_y[k],U_s);
        rel_flag =is_realizable(U_s);

        
        //printf("U_bar %lf %lf %lf >>>>>>>>>>>>>>>>>>>>>> \n",U_bar[0],U_bar[1],U_bar[2]);
        if(rel_flag ==0)
        {
             if(fabs(U_s[0]) > EPSILON)
             {

            /*printf("rel not flxed!!!!!!!!!!!!!!!!!!!!\n");
            printf("U_s %10.16e %10.16e %10.16e >>>>>>>>>>>>>U_s flag  %d \n",U_s[0],U_s[1],U_s[2],rel_flag);
            printf("max theta_s %10.16e \n",max_theta_s);
            exit(1);*/
             }
            
        }
            
    }
    //edges 
    for(k=0;k<9;k++)
    {
        init_zero_d(U_s,SYSTEM_DIM);
        get_approx_solution(uh_tn_sol_vec,E,mesh_element,mesh_node,mesh_element[E].vol_edge_gpts_ref_x[k],mesh_element[E].vol_edge_gpts_ref_y[k],U_s);
 
        rel_flag =is_realizable(U_s);
       //printf("rel_flag %d \n",rel_flag); 
        //printf("U_s %lf %lf %lf >>>>>>>>>>>>>U_s flag  %d \n",U_s[0],U_s[1],U_s[2],rel_flag);
        //printf("U_bar %lf %lf %lf >>>>>>>>>>>>>>>>>>>>>> \n",U_bar[0],U_bar[1],U_bar[2]);
        if(rel_flag ==0)
        {
            if(fabs(U_s[0]) > EPSILON)
             {

                //printf("rel flag not fixed in edges \n");
                // printf("U_s %10.16e %10.16e %10.16e >>>>>>>>>>>>>U_s flag  %d \n",U_s[0],U_s[1],U_s[2],rel_flag);
                // printf("max theta_s %10.16e \n",max_theta_s);
                 // exit(1);
             }
        }       
    }



}


void ensure_realizability(int E,
                     node* mesh_node,
                     element* mesh_element,
                     edge* mesh_edge,
                     double* U_bar,
                     double* U_s,
                     double* U_hat,
                     double epsilon,
                     double* theta_s
        )
{
    double E_bar=0.0;
    double F_bar[2];
    double E_s =0.0;
    double F_s[2];

    double a=0.0;
    double b=0.0;
    double c=0.0;
    double q=0.0;

    double sign_b = 0.0;

    double theta_1=0.0;
    double theta_2=0.0;
    

 


    E_bar = U_bar[0];
    F_bar[0] = U_bar[1];
    F_bar[1] = U_bar[2];

    E_s = U_s[0];
    F_s[0] = U_s[1];
    F_s[1] = U_s[2];



    a = (pow((F_bar[0]-F_s[0]),2.0) + pow((F_bar[1]-F_s[1]),2.0)) - pow((E_bar-E_s),2.0);
    //printf("E_bar %lf \n",E_bar);
   // printf("(pow((F_bar[0]-F_s[0]),2.0) + pow((F_bar[1]-F_s[1]),2.0)) = %10.16e \n",(pow((F_bar[0]-F_s[0]),2.0) + pow((F_bar[1]-F_s[1]),2.0)));
    //printf("a = %10.16e \n",a);
    b = 2.0*(pow(E_s,2.0) - E_bar*E_s - (pow(F_s[0],2.0) + pow(F_s[1],2.0))  + (F_bar[0]*F_s[0] + F_bar[1]*F_s[1]));
    //printf("b = %10.16e \n",b);

    if(b < 0.0)
    {
        sign_b = -1.0;
    }
    else if(b> 0.0)
    {
        sign_b = 1.0;
    }
    else
    {
        sign_b =0.0;
    }


    c = pow(F_s[0],2.0) + pow(F_s[1],2.0) - pow(E_s,2.0);
    //printf("c = %10.16e \n",c);
    //printf(" pow(b,2.0) - 4*a*c) = %10.16e \n",pow(b,2.0) - 4.0*a*c);

    if((pow(b,2.0) - 4.0*a*c) < 0.0)
    {
        q= -0.5*b;
    }
    else
    {
        q = -0.5*(b+ sign_b*sqrt(pow(b,2.0) - 4.0*a*c));
    }

    //printf("q = %10.16e \n",q);
    theta_1 = q/a;
    theta_2 = c/q;

    if((theta_1 >=0.0) && (theta_1 <=1.0))
    {
        *theta_s = theta_1 + epsilon;
    }
    else if((theta_2 >=0.0) && (theta_2 <=1.0))
    {
        *theta_s = theta_2 + epsilon;
    }
    else
    {
        *theta_s = 0.0;
    }
    //printf("theta_s %10.16e \n",*theta_s);getchar();


    U_hat[0] = *theta_s*U_bar[0] + (1.0-*theta_s)*U_s[0];
    U_hat[1] = *theta_s*U_bar[1] + (1.0-*theta_s)*U_s[1];
    U_hat[2] = *theta_s*U_bar[2] + (1.0-*theta_s)*U_s[2];
    //printf("U_s %lf %lf %lf \n",U_s[0],U_s[1],U_s[2]);
    //printf("U_bar %lf %lf %lf \n",U_bar[0],U_bar[1],U_bar[2]);
    //printf("U_hat %10.16e %10.16e %10.16e \n",U_hat[0],U_hat[1],U_hat[2]);getchar();

}

int is_realizable(double* U_sol)
{
    double psi_zero =0.0;
    double norm_psi_one =0.0;
    int rflag=0;

    double ratio =0.0;

    psi_zero = U_sol[0];
    norm_psi_one = sqrt(pow(U_sol[1],2.0) + pow(U_sol[2],2.0));
    ratio = norm_psi_one/psi_zero;


    if((psi_zero >= norm_psi_one) && (psi_zero >=EPSILON) ) //|| ((fabs(psi_zero-norm_psi_one) <EPSILON) && (psi_zero >=0.0)))
    {
        rflag =1;
        if(ratio >1.0)
            printf("ratio %10.16e psi_zero %10.16e norm_psi_one %10.16e \n",ratio,psi_zero,norm_psi_one);
        assert(ratio <=1.0);

    }
    else
    {
        rflag =0;
    }
    
    return rflag;
}
void compute_solution_average(int E,
                              node* mesh_node,
                              element* mesh_element,
                              edge* mesh_edge,
                              double* uh_tn_sol,
                              double* uh_quad_average
        )
{
    int idofs=0;
    int k=0;
    int s_dim=0;
    double approx_sol[SYSTEM_DIM];

    init_zero_d(uh_quad_average,SYSTEM_DIM);
    
    
    //interior points
    for(k=0;k<9;k++)
    {
        init_zero_d(approx_sol,SYSTEM_DIM);
        get_approx_solution(uh_tn_sol,E,mesh_element,mesh_node,mesh_element[E].interior_gpts_ref_x[k],mesh_element[E].interior_gpts_ref_y[k],approx_sol);
        
        for(s_dim=0;s_dim<SYSTEM_DIM;s_dim++)
        {
            uh_quad_average[s_dim] += (1.0/(0.5*mesh_element[E].det))*(mesh_element[E].det*approx_sol[s_dim]*mesh_element[E].interior_gpts_w[k]);
        }
    }

    //edges 
    for(k=0;k<9;k++)
    {
        init_zero_d(approx_sol,SYSTEM_DIM);
        get_approx_solution(uh_tn_sol,E,mesh_element,mesh_node,mesh_element[E].vol_edge_gpts_ref_x[k],mesh_element[E].vol_edge_gpts_ref_y[k],approx_sol);
        
        for(s_dim=0;s_dim<SYSTEM_DIM;s_dim++)
        {
            uh_quad_average[s_dim] += (1.0/(0.5*mesh_element[E].det))*(mesh_element[E].det*approx_sol[s_dim]*mesh_element[E].vol_edge_w[k]);
        }
    }
}
