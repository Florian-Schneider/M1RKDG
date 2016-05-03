/*
l1-error estimate
Prince Chidyagwai
Philipp Monreal
*/

/* double l1_dg_error
 * computes l1 error
 * in: number of elements             -- Nelts
 * in: mesh element data structure    -- mesh_element
 * in: mesh_node data structure       -- mesh_node
 * in: solution vector 		      -- uh_tn_sol_vec
 * in: test case 		      -- test_case
 * in: current time 		      -- current_time
 * in: number of basis function	      -- Nloc
 * out: l1 error
 */
#include"dg_error.h"
void l1_dg_error(int Nelts,
                 element* mesh_element,
                 node* mesh_node,
		         parameters* user_param,
                 double* uh_tn_sol_vec,
                 double current_time,
		         int s_dim_count,
                 int Nloc,
		         double* l1_error,
		         double* linf_error)
{
    int k=0;
    int E=0;
    double* gpx = NULL;
    double* gpy = NULL;
    double* w = NULL;
    double  phys_coords [2];
    double* phi_vec = NULL;
    double* grad_phi_mat = NULL;
    double det_BE_T;
    double local_error =0.0;
    double global_error=0.0;
    double exact_solution[SYSTEM_DIM];
    double approx_solution[SYSTEM_DIM];
    double sup_error =0.0;
    double new_sup_error =0.0;

    gpx = (double*)malloc((ngpts+1)*sizeof(double));
    gpy = (double*)malloc((ngpts+1)*sizeof(double));
    w = (double*)malloc((ngpts+1)*sizeof(double));

    //gauss points on reference element
    gauss_pt(gpx,gpy,w);

    //allocate memory for basis functions
    phi_vec = (double*)malloc(Nloc*sizeof(double));
    grad_phi_mat = (double*)malloc(2*Nloc*sizeof(double));

    for(E=1;E<Nelts+1;E++) 
    {
        det_BE_T = mesh_element[E].det;
        local_error = 0.0;
        //loop over gauss points
        for(k=1;k<ngpts+1;k++) 
        {
            init_zero_d(grad_phi_mat,2*Nloc);
            init_zero_d(phi_vec,Nloc);

            //map gauss points to physical element
            map_to_physical_element(mesh_element,mesh_node,E,phys_coords,gpx[k],gpy[k]);

            //exact solution
	        true_solution_val(phys_coords[0],phys_coords[1],gpx[k],gpy[k],current_time,E,mesh_element,mesh_node,uh_tn_sol_vec,user_param,exact_solution);

            get_approx_solution(uh_tn_sol_vec,E,mesh_element,mesh_node,gpx[k],gpy[k],approx_solution);

            sup_error = fabs(exact_solution[s_dim_count]-approx_solution[s_dim_count]);

            if(sup_error > new_sup_error)
            {
                new_sup_error = sup_error;
		// Debugging to find where largest errors occur
		/*printf("highest error found at: \tx = %lf\t y = %lf\n",phys_coords[0],phys_coords[1]);*/
            }

	    //approximate solution
	    //printf("approx %lf true %lf \n",approx_solution,true_solution);
	    // assumes that the area of the domain is 1 !!!
	    local_error = local_error+det_BE_T*fabs(exact_solution[s_dim_count]-approx_solution[s_dim_count])*w[k];
	    //printf("local_error %lf det_BET %lf w[k] %lf \n",local_error,det_BE_T,w[k]);
	}
	global_error = global_error + local_error;
    }
    *l1_error = global_error;
    *linf_error = new_sup_error;

    free(phi_vec);
    free(grad_phi_mat);
    free(w);
    free(gpy);
    free(gpx);
}


void check_dg_sol_realizability(element* mesh_element,
                                int Nelts,
                                node* mesh_node,
                                parameters* user_param,
                                double* uh_tn_sol_vec,
                                double current_time,
                                double* min_ratio_psi_one_psi_zero,
                                double* max_ratio_psi_one_psi_zero
                                )
{

    int E=0;
    int k=0;
    double approx_solution[SYSTEM_DIM];
    double psi_zero =0.0;
    double norm_psi_one =0.0;

    double min_ratio=-1.0;
    double max_ratio=0.0;
    double ratio =0.0;
    
    double* gpx = NULL;
    double* gpy = NULL;
    double* w = NULL;
    int rel_flag=0;
  
    gpx = (double*)malloc((ngpts+1)*sizeof(double));
    gpy = (double*)malloc((ngpts+1)*sizeof(double));
    w = (double*)malloc((ngpts+1)*sizeof(double));

    //gauss points on reference element
    gauss_pt(gpx,gpy,w);

    //allocate memory for basis functions
    int count=0;

     
    for(E=1;E<Nelts+1;E++) 
    {
        //loop over gauss points
        for(k=0;k<9;k++) 
        {
            get_approx_solution(uh_tn_sol_vec,E,mesh_element,mesh_node,mesh_element[E].interior_gpts_ref_x[k],mesh_element[E].interior_gpts_ref_y[k],approx_solution);
            rel_flag =is_realizable(approx_solution);
            psi_zero = approx_solution[0];
            norm_psi_one = sqrt(pow(approx_solution[1],2.0) + pow(approx_solution[2],2.0));

            if((psi_zero >EPSILON) && (norm_psi_one > EPSILON))
            {

                if(rel_flag ==0)
                {
                    //printf("approx_sol %10.16e %10.16e %10.16e \n",approx_solution[0],approx_solution[1],approx_solution[2]);
                    //exit(1);
                }
                ratio = norm_psi_one/psi_zero;
                //assert(psi_zero >0.0);
                //printf("ration %10.16e \n",ratio);
            }


            if(count ==0)
            {
                min_ratio = ratio;
                max_ratio = ratio;
                count =1;
            }
            else
            {
                if(ratio > max_ratio)
                {
                    max_ratio = ratio;
                }
                if(ratio < min_ratio)
                {
                    min_ratio = ratio;
                }
            }
        }
    }
    *min_ratio_psi_one_psi_zero = min_ratio;
    *max_ratio_psi_one_psi_zero = max_ratio;

    free(gpx);
    free(gpy);
    free(w);
}

void check_symmetry(element* mesh_element,
                   int Nelts,
                   node* mesh_node,
                   edge* mesh_edge,
                   double* uh_tn_sol_vec,
                   int s_dim_count,
                   int Nloc,
                   double current_time,
                   double* symmetry_error)
{
    int E=0;
    int Es =0; //symmetric element
    double* gpx = NULL;
    double* gpy = NULL;
    double* w = NULL;
    double sup_error =0.0;
    double new_sup_error=0.0;
    double local_error =0.0;
    double global_error =0.0;
    double Uh_E[SYSTEM_DIM];
    double Uh_Es[SYSTEM_DIM];
    int k=0;
  
    gpx = (double*)malloc((ngpts+1)*sizeof(double));
    gpy = (double*)malloc((ngpts+1)*sizeof(double));
    w = (double*)malloc((ngpts+1)*sizeof(double));

    
    double approx_solutionE[SYSTEM_DIM];
    double approx_solutionEs[SYSTEM_DIM];

    //gauss points on reference element
    gauss_pt(gpx,gpy,w);

    int i=0;
    int j=0;

    int s_dim=0;
    int idofs=0;
    int jdofs=0;
    double E_coef_sol=0.0;
    double Es_coef_sol=0.0;
    double diff_E_Es_coefs =0.0;
    double min_diff_E_Es_coefs =0.5;





    global_error =0.0;
    for(E=1;E<Nelts+1;E++)
    {
        Es = find_symmetric_element(mesh_element,Nelts,mesh_node,E); 

        min_diff_E_Es_coefs=0.5;
        /*for(idofs=0;idofs<NLOC;idofs++)
        {
            E_coef_sol= uh_tn_sol_vec[SYSTEM_DIM*(E-1)*NLOC + (NLOC*s_dim)+idofs];
           // printf("idofs = %d E_coef_sol %10.32lf \n",idofs,E_coef_sol);

            
            for(jdofs =0;jdofs < NLOC;jdofs++)
            {
                Es_coef_sol= uh_tn_sol_vec[SYSTEM_DIM*(Es-1)*NLOC + (NLOC*s_dim)+jdofs];

             //   printf("jdofs = %d Es_coef_sol %10.32lf \n",jdofs,Es_coef_sol);
                
                diff_E_Es_coefs = fabs(E_coef_sol - Es_coef_sol);
                if(diff_E_Es_coefs < min_diff_E_Es_coefs)
                {
                    min_diff_E_Es_coefs = diff_E_Es_coefs;
                    i = idofs;
                    j = jdofs;
                }
            }
            //printf("-----------------------------------------------------\n");
            sup_error = min_diff_E_Es_coefs;
        //    if(sup_error > 0)
        //        printf("mid_diff -------------------------------------------- %10.32lf \n",sup_error);
            if(sup_error > 1.0e-32)
            {
                if(Es > E)
                {
                    if(s_dim ==0)
                    {
                        uh_tn_sol_vec[SYSTEM_DIM*(Es-1)*NLOC + (NLOC*s_dim)+j]  = uh_tn_sol_vec[SYSTEM_DIM*(E-1)*NLOC + (NLOC*s_dim)+i];
                    }  
                }
            }
        }*/
   


       
       
        //printf("E %d Es %d \n",E,Es);
        
        local_error =0.0;
        init_zero_d(Uh_E,SYSTEM_DIM);
        init_zero_d(Uh_Es,SYSTEM_DIM);
        //Uh_E =0.0;
        //Uh_Es=0.0;
        if(Es==0)
        {
            exit(1);
        }
        else{

        compute_solution_average(E,mesh_node,mesh_element,mesh_edge,uh_tn_sol_vec,Uh_E);
        compute_solution_average(Es,mesh_node,mesh_element,mesh_edge,uh_tn_sol_vec,Uh_Es);
        sup_error = fabs(Uh_E[0] -Uh_Es[0]);

        if(sup_error > 1.0e-16)
        {
            //printf("s_dim_count %d symmetry Error :: E %d Es %d   %10.16e \n",s_dim_count,E,Es,global_error);
        }
        if(sup_error > new_sup_error)
        {
            new_sup_error = sup_error;
        }

        }
        //printf("E %d Es %d sup error %lf \n",E,Es,sup_error);
        global_error = new_sup_error;
    }
   printf("%lf %10.12e \n",current_time,global_error);

    *symmetry_error = global_error;
}

int find_symmetric_edge(edge* mesh_edge,
                        element* mesh_element,
                        node* mesh_node,
                        int E,
                        int Es,
                        int E_iedge
        )
{
    int n1,n2;
    int n1s,n2s;
    double x1,y1;
    double x2,y2;
    double x1s,y1s;
    double x2s,y2s;

    int found =0;
    int i=1;
    int Es_iedge=0;
    int verbose=0;

    n1 = mesh_edge[E_iedge].vertex[1];
    n2 = mesh_edge[E_iedge].vertex[2];


    x1 = mesh_node[n1].coord[0];
    y1 = mesh_node[n1].coord[1];
    x2 = mesh_node[n2].coord[0];
    y2 = mesh_node[n2].coord[1];


    while((!found) && (i <4))
    {
        Es_iedge = mesh_element[Es].edge[i];
        n1s = mesh_edge[Es_iedge].vertex[1];
        n2s = mesh_edge[Es_iedge].vertex[2];

      

        x1s = mesh_node[n1s].coord[0];
        y1s = mesh_node[n1s].coord[1];
        x2s = mesh_node[n2s].coord[0];
        y2s = mesh_node[n2s].coord[1];

        if(verbose)
        {
            printf("i =  %d \n",i);
            printf(" (n1 %d n2 %d ) and (n1s %d  n2s %d) \n",n1,n2,n1s,n2s);
            printf("(x1 %lf y1 %lf x2 %lf y2 %lf \n",x1,y1,x2,y2);
            printf("(x1s %lf y1s %lf x2s %lf y2s %lf \n",x1s,y1s,x2s,y2s);
        }

        

        if((((fabs(x1s-(1.0-y1)) < EPSILON) && (fabs(y1s-(1.0-x1)) <EPSILON)) ||
          ((fabs(x1s-(1.0-y2)) < EPSILON) && (fabs(y1s-(1.0-x2)) <EPSILON))) &&     
          (((fabs(x2s-(1.0-y1)) < EPSILON) && (fabs(y2s-(1.0-x1)) <EPSILON)) ||
          ((fabs(x2s-(1.0-y2)) < EPSILON) && (fabs(y2s-(1.0-x2)) <EPSILON))))    
        {
            found = i;
        }
       i++;
    }
    
    if(verbose)
        printf("-----found = %d ----------------------------\n", found);
   return found; 
}


int find_symmetric_element(element* mesh_element,
                           int Nelts,
                           node* mesh_node,
                           int E)
{
    int n1,n2,n3;
    int n1s,n2s,n3s;
    double x1,y1,x2,y2,x3,y3;
    double x1s,y1s,x2s,y2s,x3s,y3s;
    int Es=0;
    int found =0;

	n1 = mesh_element[E].vertex[1];
	n2 = mesh_element[E].vertex[2];
	n3 = mesh_element[E].vertex[3];

    x1 = mesh_node[n1].coord[0];
    y1 = mesh_node[n1].coord[1];
    x2 = mesh_node[n2].coord[0];
    y2 = mesh_node[n2].coord[1];
    x3 = mesh_node[n3].coord[0];
    y3 = mesh_node[n3].coord[1];

    //printf("(x1 %lf, y1 %lf), (x2 %lf, y2 %lf ) , (x3 %lf, y3 %lf ) \n",x1,y1,x2,y2,x3,y3);

    while((!found) && (Es <(Nelts+1)))
    {
        n1s = mesh_element[Es].vertex[1];
	    n2s = mesh_element[Es].vertex[2];
    	n3s = mesh_element[Es].vertex[3];


        x1s = mesh_node[n1s].coord[0];
        y1s = mesh_node[n1s].coord[1];
        x2s = mesh_node[n2s].coord[0];
        y2s = mesh_node[n2s].coord[1];
        x3s = mesh_node[n3s].coord[0];
        y3s = mesh_node[n3s].coord[1];


      
        /*printf("E %d  Es %d \n",E,Es);
        printf("(x1 %lf, y1 %lf), (x2 %lf, y2 %lf ) , (x3 %lf, y3 %lf ) \n",x1,y1,x2,y2,x3,y3);
        printf("(x1s %lf, y1s %lf), (x2s %lf, y2s %lf ) , (x3s %lf, y3s %lf ) \n",x1s,y1s,x2s,y2s,x3s,y3s);
        printf("(%lf %lf ) (%lf %lf ) (%lf %lf ) \n", fabs(x1s-(1.0-y1)), 
                                                      fabs(y1s-(1.0-x1)),
                                                      fabs(x2s-(1.0-y2)),
                                                      fabs(y2s-(1.0-x2)),
                                                      fabs(x3s-(1.0-y3)),
                                                      fabs(y3s-(1.0-x3)));*/
        //getchar();*/

        //symmetry about y=1-x

        if((((fabs(x1s-(1.0-y1)) < EPSILON) && (fabs(y1s-(1.0-x1)) <EPSILON)) ||
            ((fabs(x1s-(1.0-y2)) < EPSILON) && (fabs(y1s-(1.0-x2)) <EPSILON)) ||
            ((fabs(x1s-(1.0-y3)) < EPSILON) && (fabs(y1s-(1.0-x3)) <EPSILON)) ) &&     
           (((fabs(x2s-(1.0-y1)) < EPSILON) && (fabs(y2s-(1.0-x1)) <EPSILON)) ||
            ((fabs(x2s-(1.0-y2)) < EPSILON) && (fabs(y2s-(1.0-x2)) <EPSILON)) ||
            ((fabs(x2s-(1.0-y3)) < EPSILON) && (fabs(y2s-(1.0-x3)) <EPSILON))) &&
           (((fabs(x3s-(1.0-y1)) < EPSILON) && (fabs(y3s-(1.0-x1)) <EPSILON)) ||
            ((fabs(x3s-(1.0-y2)) < EPSILON) && (fabs(y3s-(1.0-x2)) <EPSILON))  ||
            ((fabs(x3s-(1.0-y3)) < EPSILON) && (fabs(y3s-(1.0-x3)) <EPSILON))))
           /*if((((fabs(x1s-(y1)) < EPSILON) && (fabs(y1s-(x1)) <EPSILON)) ||
            ((fabs(x1s-(y2)) < EPSILON) && (fabs(y1s-(x2)) <EPSILON)) ||
            ((fabs(x1s-(y3)) < EPSILON) && (fabs(y1s-(x3)) <EPSILON)) ) &&     
           (((fabs(x2s-(y1)) < EPSILON) && (fabs(y2s-(x1)) <EPSILON)) ||
            ((fabs(x2s-(y2)) < EPSILON) && (fabs(y2s-(x2)) <EPSILON)) ||
            ((fabs(x2s-(y3)) < EPSILON) && (fabs(y2s-(x3)) <EPSILON))) &&
           (((fabs(x3s-(y1)) < EPSILON) && (fabs(y3s-(x1)) <EPSILON)) ||
            ((fabs(x3s-(y2)) < EPSILON) && (fabs(y3s-(x2)) <EPSILON))  ||
            ((fabs(x3s-(y3)) < EPSILON) && (fabs(y3s-(x3)) <EPSILON))))*/


        {
            found =Es;
            mesh_element[E].symmetric_element = found;
           //printf("E %d Es %d \n",E, mesh_element[E].symmetric_element);
        }
        Es++;

    }
    if(found ==0)
    {
        printf("Syemmtric element not found: Mesh might be irregular \n");
        exit(1);
    }
    return found;
}


/* void true_solution
 * defines the boundary conditions
 * for NEUMANN_BC we need reference coordinates 
 * instead of physical ones
 * in: x coordinate physical -- x
 * in: y coordinate physical -- y
 * in: x coordinate reference -- x_ref
 * in: y coordinate reference -- y_ref
 * in: time	    -- t
 * in: index to triangle -- E
 * in: mesh_element data structure -- mesh_element
 * in: mesh_node data structure    -- mesh_node
 * in: parameters 		   -- user_param
 * in: solution at time tn         -- uh_tn_sol_vec
 * in: parameters   -- user_param
 * out: boundary value at x,y
 */
void  true_solution_val(double x,
                        double y,
                        double x_ref,
                        double y_ref,
                        double t,
                        int E,
                        element* mesh_element,
                        node* mesh_node,
                        double* uh_tn_sol_vec,
                        parameters* user_param,
                        double* bdry_f_val)
{
    int s_dim;
    double temp;


    switch ((*user_param).test_case) {
        case(1): {
	    for (s_dim=0;s_dim<SYSTEM_DIM;s_dim++) bdry_f_val[s_dim] = 1.0; 
	    break;
	}

    	case(4): {
	    if(SYSTEM_DIM ==3) {
		bdry_f_val[0] = exp(-t)*sin(2*PI*(x+y));
		bdry_f_val[1] = -1.0*exp(-t)*t*(2.0/3.0)*PI*cos(2*PI*(x+y));
		bdry_f_val[2] = -1.0*exp(-t)*t*(2.0/3.0)*PI*cos(2*PI*(x+y));
	    }
	    break;
	}

	case(7): {
	    if(SYSTEM_DIM ==3)
	    {
		bdry_f_val[0] = exp(-t)*sin(4*PI*((x-0.5)*(x-0.5)+(y-0.5)*(y-0.5)));
		bdry_f_val[1] = -4./3.*PI*t*exp(-t)*(2.0*x-1.0)*cos(4*PI*((x-0.5)*(x-0.5)+(y-0.5)*(y-0.5)));
		bdry_f_val[2] = -4./3.*PI*t*exp(-t)*(2.0*y-1.0)*cos(4*PI*((x-0.5)*(x-0.5)+(y-0.5)*(y-0.5)));
	    }
	    break;
	}
	
	case(8): {
	    if(SYSTEM_DIM ==3)
	    {
	        temp = -4./3.*PI*t*exp(-t)*cos(4*PI*((x-0.5)*(x-0.5)+(y-0.5)*(y-0.5)));
		bdry_f_val[0] = exp(-t)*sin(4*PI*((x-0.5)*(x-0.5)+(y-0.5)*(y-0.5)));
		bdry_f_val[1] = (2.0*x-1.0)*temp + 1.0/2.0*(sin(x-t) - cos(x-t));
		bdry_f_val[2] = (2.0*y-1.0)*temp + 1.0/2.0*(sin(y-t) + cos(y-t));
	    }
	    break;
	}

	case(11): {
	    if(SYSTEM_DIM ==3)
	    {
		if (x <= 0.00) {
		    bdry_f_val[0] = 1.0;
		    bdry_f_val[1] = 0.99;
		    bdry_f_val[2] = 0.0;
		}
		if (x >= 1.0) {
		    bdry_f_val[0] = 0.5;
		    bdry_f_val[1] = 0.0;
		    bdry_f_val[2] = 0.0;
		}
	    }
	    break;
	}

    	case(12): {
	    if(SYSTEM_DIM ==3) {
		bdry_f_val[0] = exp(-t)*sin(2.*PI*(x+2.*y)) + 3.;
		bdry_f_val[1] = exp(-t)*t*cos(2.*PI*(x+2.*y));
		bdry_f_val[2] = exp(-t)*t*cos(2.*PI*(x+2.*y));
	    }
	    break;
	}

    	case(14): {
	    if(SYSTEM_DIM ==3) {
		bdry_f_val[0] = exp(-t)*sin(2*PI*(x+2.*y));
		bdry_f_val[1] = -1.0*exp(-t)*t*(2.0/3.0)*PI*cos(2*PI*(x+2.*y));
		bdry_f_val[2] = -1.0*exp(-t)*t*(4.0/3.0)*PI*cos(2*PI*(x+2.*y));
	    }
	    break;
	}

    	case(16): {
	    if(SYSTEM_DIM ==3) {
		bdry_f_val[0] = exp(-t)*sin(2.*PI*(x+2.*y)) + 3.;
		bdry_f_val[1] = exp(-t)*t*cos(2.*PI*(x+2.*y));
		bdry_f_val[2] = exp(-t)*t*cos(2.*PI*(x+2.*y));
	    }
	    break;
	}

    	case(18): {
	    if(SYSTEM_DIM ==3) {
		bdry_f_val[0] = exp(-t)*sin(2.*PI*(x+2.*y)) + 3.;
		bdry_f_val[1] = exp(-t)*t*cos(2.*PI*(x+2.*y));
		bdry_f_val[2] = exp(-t)*t*cos(2.*PI*(x+2.*y));
	    }
	    break;
	}

    	case(19): {
	    if(SYSTEM_DIM ==3) {
		bdry_f_val[0] = exp(-t)*sin(2.*PI*(x+2.*y)) + 3.;
		bdry_f_val[1] = exp(-t)*t*cos(2.*PI*(x+2.*y));
		bdry_f_val[2] = exp(-t)*t*cos(2.*PI*(x+2.*y));
	    }
	    break;
	}

	case(23): {
	    if(SYSTEM_DIM ==3)
	    {
		bdry_f_val[0] = atan(200.*(y-0.51)) + 1;;
		bdry_f_val[1] = 0.0;
		bdry_f_val[2] = 0.0;
	    }
	    break;
	}
    case(24):{
                 bdry_f_val[0] = 0.0;
                 bdry_f_val[1] = 0.0;
                 bdry_f_val[2] = 0.0;
                 break;
             }
    case(25):{
                 bdry_f_val[0] = 0.0;
                 bdry_f_val[1] = 0.0;
                 bdry_f_val[2] = 0.0;
                 break;
             }
    case(26):
             {
                 bdry_f_val[0] = 0.0;
                 bdry_f_val[1] = 0.0;
                 bdry_f_val[2] = 0.0;
                 break;
             }
	// Zero boundary conditions by default
	default: {
	    for (s_dim=0;s_dim<SYSTEM_DIM;s_dim++) bdry_f_val[s_dim] = 0.0;
	    break;
	}
    }
}
void symmetry_correction(element* mesh_element,
                         node* mesh_node,
                         double* uh_tn_sol_vec,
                         int num_elts,
                         double current_time,
                         parameters* user_param
        )
{
    int E=0;
    int Es=0;
    int s_dim=0;
    int idofs=0;
    int j=0;
    int found =0;
    int k=0;
    double approx_solE[SYSTEM_DIM];
    double approx_solEs[SYSTEM_DIM];
    double phys_coordsE[2];
    double phys_coordsEs[2];
    double ref_coordsEs[2];
    double ref_coordsE[2];
    double error0=0.0;
    double error1=0.0;
    double error2=0.0;
    double max_error0 =0.0;
    double max_error1 =0.0;
    double max_error2 =0.0;
    double* Aloc_matrix = NULL;
    double* usolE_rhs = NULL;
    double* usolE_rhs1 = NULL;
    double* usolE_rhs2 = NULL;
    
    double** A_loc_QR;
    double* Q_coeff;
    double* R_coeff;
    int sing_decomp;
    int i=0;
    double diff_coefs=0.0;

    Aloc_matrix = (double*)malloc((NLOC*NLOC)*sizeof(double));
    Q_coeff = (double*)malloc((NLOC+1)*sizeof(double));
    R_coeff = (double*)malloc((NLOC+1)*sizeof(double));
    A_loc_QR = (double**)malloc((NLOC+1)*sizeof(double*));
    usolE_rhs = (double*)malloc(NLOC*sizeof(double));
    usolE_rhs1 = (double*)malloc(NLOC*sizeof(double));
    usolE_rhs2 = (double*)malloc(NLOC*sizeof(double));
     for(i=0;i<NLOC+1;i++)
         A_loc_QR[i] = (double*)malloc((NLOC+1)*sizeof(double));







    //printf("current_time %lf \n",current_time);


    for(E=1;E<num_elts+1;E++)
    {
        Es = find_symmetric_element(mesh_element,num_elts,mesh_node,E);
        
        max_error0 =0.0;
        max_error1 =0.0;
        max_error2 =0.0;

        for(k=1;k<ngpts+1;k++)
        {
            get_approx_solution(uh_tn_sol_vec,E,mesh_element,mesh_node,mesh_element[E].el_gpts_ref_x[k],mesh_element[E].el_gpts_ref_y[k],approx_solE);
            phys_coordsEs[0] = 1.0 - mesh_element[E].el_gpts_y[k];
            phys_coordsEs[1] = 1.0 - mesh_element[E].el_gpts_x[k];
            map_to_reference_element(mesh_element,mesh_node,Es,ref_coordsEs,phys_coordsEs[0],phys_coordsEs[1]);
            get_approx_solution(uh_tn_sol_vec,Es,mesh_element,mesh_node,ref_coordsEs[0],ref_coordsEs[1],approx_solEs);
            error0 = approx_solEs[0] - approx_solE[0];
            error1 = approx_solEs[1] + approx_solE[2];
            error2 = approx_solEs[2] + approx_solE[1];
            //printf("error0 =  %.12e \n",error0);
            if(error0 > max_error0)
            {
                max_error0 = error0;
            }
            if(error1 > max_error1)
            {
                max_error1 = error1;
            }
            if(error2 > max_error2)
            {
                max_error2 = error2;
            }

        }

        //if((max_error0 > EPSILON) ||(max_error1 > EPSILON) || (max_error2 > EPSILON))
         //   printf("%.6e %.6e %.6e \n",max_error0,max_error1,max_error2);

        if((max_error0 > EPSILON_S) ||(max_error1 > EPSILON_S) || (max_error2 > EPSILON_S))
        {
            init_zero_d(Aloc_matrix,NLOC*NLOC);
            init_zero_d(usolE_rhs,NLOC);
            init_zero_d(usolE_rhs1,NLOC);
            init_zero_d(usolE_rhs2,NLOC);
            Alocal_mat(Es,mesh_element,mesh_node,user_param,Aloc_matrix);
            assemble_qr(Aloc_matrix,NLOC,Q_coeff,R_coeff,A_loc_QR,&sing_decomp);

            for(k=1;k<ngpts+1;k++)
            {
                //get_approx_solution(uh_tn_sol_vec,E,mesh_element,mesh_node,mesh_element[E].el_gpts_ref_x[k],mesh_element[E].el_gpts_ref_y[k],approx_solE);
                phys_coordsE[0] = 1.0 - mesh_element[Es].el_gpts_y[k];
                phys_coordsE[1] = 1.0 - mesh_element[Es].el_gpts_x[k];
                map_to_reference_element(mesh_element,mesh_node,E,ref_coordsE,phys_coordsE[0],phys_coordsE[1]);
                get_approx_solution(uh_tn_sol_vec,E,mesh_element,mesh_node,ref_coordsE[0],ref_coordsE[1],approx_solE);

                for(idofs=0;idofs<NLOC;idofs++)
                {
                    usolE_rhs[idofs] += approx_solE[0]*mesh_element[Es].el_gpts_basis[(k-1)*NLOC+idofs]*mesh_element[Es].det*mesh_element[Es].el_gpts_w[k];
                    usolE_rhs1[idofs] += (-1.0*approx_solE[2]*mesh_element[Es].el_gpts_basis[(k-1)*NLOC+idofs]*mesh_element[Es].det*mesh_element[Es].el_gpts_w[k]);
                    usolE_rhs2[idofs] += (-1.0*approx_solE[1]*mesh_element[Es].el_gpts_basis[(k-1)*NLOC+idofs]*mesh_element[Es].det*mesh_element[Es].el_gpts_w[k]);
                }
            }
            if(max_error0 > EPSILON_S)
            {
                qrsolv(A_loc_QR,NLOC,Q_coeff,R_coeff,usolE_rhs);
                for(idofs = 0;idofs < NLOC;idofs++)
                {
                    uh_tn_sol_vec[SYSTEM_DIM*(Es-1)*NLOC+(0*NLOC+idofs)] = usolE_rhs[idofs];

                    //diff_coefs = fabs(uh_tn_sol_vec[SYSTEM_DIM*(Es-1)*NLOC+(0*NLOC+idofs)]-usolE_rhs[idofs]);
                    //printf("Es = %10.18lf correction %10.18lf diff = %.8e\n", uh_tn_sol_vec[SYSTEM_DIM*(Es-1)*NLOC+(0*NLOC+idofs)],usolE_rhs[idofs],diff_coefs);
                }
            }
            //printf("-------------------------\n");
             init_zero_d(Aloc_matrix,NLOC*NLOC);
             init_zero_m(A_loc_QR,NLOC+1,NLOC+1);
             Alocal_mat(Es,mesh_element,mesh_node,user_param,Aloc_matrix);
             assemble_qr(Aloc_matrix,NLOC,Q_coeff,R_coeff,A_loc_QR,&sing_decomp);
             if(max_error1 > EPSILON_S)
             {
                 qrsolv(A_loc_QR,NLOC,Q_coeff,R_coeff,usolE_rhs1);
                for(idofs = 0;idofs < NLOC;idofs++)
                {

                    uh_tn_sol_vec[SYSTEM_DIM*(Es-1)*NLOC+(1*NLOC+idofs)] =usolE_rhs1[idofs];
                    //diff_coefs = fabs(uh_tn_sol_vec[SYSTEM_DIM*(Es-1)*NLOC+(1*NLOC+idofs)]-usolE_rhs1[idofs]);
                    //printf("Es = %10.18lf correction %10.18lf diff = %.8e\n", uh_tn_sol_vec[SYSTEM_DIM*(Es-1)*NLOC+(1*NLOC+idofs)],usolE_rhs1[idofs],diff_coefs);
                }
             }
             //printf("-------------------------\n");
             init_zero_d(Aloc_matrix,NLOC*NLOC);
             init_zero_m(A_loc_QR,NLOC+1,NLOC+1);
             Alocal_mat(Es,mesh_element,mesh_node,user_param,Aloc_matrix);
             assemble_qr(Aloc_matrix,NLOC,Q_coeff,R_coeff,A_loc_QR,&sing_decomp);
             if(max_error2 > EPSILON_S)
             {
                 qrsolv(A_loc_QR,NLOC,Q_coeff,R_coeff,usolE_rhs2);
                for(idofs = 0;idofs < NLOC;idofs++)
                {

                    uh_tn_sol_vec[SYSTEM_DIM*(Es-1)*NLOC+(2*NLOC+idofs)]=usolE_rhs2[idofs];
                    //diff_coefs = fabs(uh_tn_sol_vec[SYSTEM_DIM*(Es-1)*NLOC+(2*NLOC+idofs)]-usolE_rhs2[idofs]);
                    //printf("Es = %10.18lf correction %10.18lf diff = %.8e\n", uh_tn_sol_vec[SYSTEM_DIM*(Es-1)*NLOC+(2*NLOC+idofs)],usolE_rhs2[idofs],diff_coefs);
                }
             }

            //getchar();
        }
    }
}
