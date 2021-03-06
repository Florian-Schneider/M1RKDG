/*
Hyperbolic solver
Prince Chidyagwai
Philipp Monreal
*/
#include"hyperbolic_solver.h"

/* void rk_dg_time_step
 * Runge-Kutta time stepping with Discrete-Galerkin space discretisation
 * main solver-routine, assembling the ODE system and calling all auxiliary programs
 * in: mesh element data structure    -- mesh_element
 * in: mesh_node data structure       -- mesh_node
 * in: mesh edge data structure       -- mesh_edge
 * in: number of time steps	      -- N_time_steps
 * in: number of basis functions      -- Nloc
 * in: number of nodes  	      -- num_nodes
 * in: number of elements             -- Nelts
 * in: size of time step	      -- delta_t
 * in: parameter structure	      -- user_param
 * out: solution vector at final time -- uh_tn1_sol_vec
 */
void rk_dg_time_step(element* mesh_element,
                     node* mesh_node,
                     edge* mesh_edge,
		             FILE* init_coefs_file,
		             double start_time,
                     int N_time_steps,
                     int Nloc,
                     int num_nodes,
                     int Nelts,
                     double delta_t,
                     parameters* user_param,
                     double* uh_tn1_sol_vec)
{
    // Loop variable
    int E=0;
    int idofs=0;
    int i=0;
    int s_dim=0;
    int tsteps =0;
    int size =0;
    int rel_flag =0;
    double* uzero_loc = NULL;
    double* U_bar = NULL;
    double* U_bar_post = NULL;
    double* uval = NULL;
    double* Aloc_matrix = NULL;
    double* temp_Aloc = NULL;
    double* loc_pde_rhs_vec = NULL;
    double* loc_source_rhs_vec = NULL;
    double* sigmas_vec = NULL;
    double* sigmat_vec = NULL;

    // QR-decomposition and solving
    double** A_loc_QR; 
    double* Q_coeff; 
    double* R_coeff; 
    int sing_decomp; 

    double current_time  =0.0;
    double current_time1 =0.0;
    double current_time2 =0.0;
    double current_time3 =0.0;
    double current_time4 =0.0;
    double current_time5 =0.0;
    double current_time6 =0.0;
    double current_time7 =0.0;
    double current_time8 =0.0;
    double current_time9 =0.0;
    double current_time10 =0.0;

    // solution vector
    double* uh_tn_sol_vec = NULL;
    double* uh_pre_limit = NULL;
    double* uh_limit = NULL;

    double* U_hat = NULL;
    double s_epsilon = 1.0e-12;

    // intermediate solution vectors
    double* uh_aux1 = NULL;
    double* uh_aux2 = NULL;
    double* uh_aux3 = NULL;
    double* uh_aux4 = NULL;
    double* uh_aux5 = NULL;
    double* uh_aux6 = NULL;
    double* uh_aux7 = NULL;
    double* uh_aux8 = NULL;
    double* uh_aux9 = NULL;

    double* projected_vec = NULL;

    // zero vector 
    double* zero_vec = NULL;

    // variables for intermediate output
    char file_out[80];

    char meshname_str[80];
    int result =0;

    file_out[0] = '\0';
    int ind_out = 0;

    // Information about when a thread finishes a timestep
    int verbose_parallel = 0;
    double symm_error=0.0;

    // Allocate variables
    size =SYSTEM_DIM*Nloc*Nelts;
    uh_tn_sol_vec = (double*)malloc(SYSTEM_DIM*(Nelts*Nloc)*sizeof(double));
    projected_vec = (double*)malloc(SYSTEM_DIM*(Nelts*Nloc)*sizeof(double));
    uh_aux1 = (double*)malloc(SYSTEM_DIM*(Nelts*Nloc)*sizeof(double));
    uh_aux2 = (double*)malloc(SYSTEM_DIM*(Nelts*Nloc)*sizeof(double));
    uh_aux3 = (double*)malloc(SYSTEM_DIM*(Nelts*Nloc)*sizeof(double));
    uh_aux4 = (double*)malloc(SYSTEM_DIM*(Nelts*Nloc)*sizeof(double));
    uh_aux5 = (double*)malloc(SYSTEM_DIM*(Nelts*Nloc)*sizeof(double));
    uh_aux6 = (double*)malloc(SYSTEM_DIM*(Nelts*Nloc)*sizeof(double));
    uh_aux7 = (double*)malloc(SYSTEM_DIM*(Nelts*Nloc)*sizeof(double));
    uh_aux8 = (double*)malloc(SYSTEM_DIM*(Nelts*Nloc)*sizeof(double));
    uh_aux9 = (double*)malloc(SYSTEM_DIM*(Nelts*Nloc)*sizeof(double));
        
    zero_vec = (double*)malloc(SYSTEM_DIM*(Nelts*Nloc)*sizeof(double));
    uh_pre_limit = (double*)malloc(SYSTEM_DIM*(Nelts*3)*sizeof(double));
    uh_limit = (double*)malloc(SYSTEM_DIM*(Nelts*Nloc)*sizeof(double));


    // Initializations
    init_zero_d(uh_tn_sol_vec,size);
    init_zero_d(uh_tn1_sol_vec,size);
    init_zero_d(uh_aux1,size);
    init_zero_d(uh_aux2,size);
    init_zero_d(uh_aux3,size);
    init_zero_d(uh_aux4,size);
    init_zero_d(uh_aux5,size);
    init_zero_d(uh_aux6,size);
    init_zero_d(uh_aux7,size);
    init_zero_d(uh_aux8,size);
    init_zero_d(uh_aux9,size);
    init_zero_d(zero_vec,size);
    init_zero_d(uh_pre_limit,SYSTEM_DIM*Nelts*3);

    // Start parallel region with declaration of variables used
    #pragma omp parallel \
    default(none) \
    shared(uh_tn_sol_vec,uh_tn1_sol_vec,uh_aux1,uh_aux2,uh_aux3,uh_aux4,uh_aux5,uh_aux6,uh_aux7,uh_aux8,uh_aux9,zero_vec,uh_pre_limit,uh_limit,user_param,mesh_edge,\
            mesh_element,mesh_node,current_time,current_time1,current_time2,current_time3,current_time4,current_time5,current_time6,current_time7,current_time8,current_time9,\
            current_time10,projected_vec,symm_error) \
    private(stderr,result,meshname_str,i,tsteps,E,s_dim,idofs,temp_Aloc,Aloc_matrix,loc_pde_rhs_vec,sigmas_vec,sigmat_vec,loc_source_rhs_vec,Q_coeff,R_coeff,sing_decomp,\
            uzero_loc,A_loc_QR,uval,U_bar,U_bar_post,rel_flag,U_hat,s_epsilon) \
    firstprivate(start_time,init_coefs_file,Nloc,Nelts,size,file_out,ind_out,delta_t,N_time_steps,num_nodes,verbose_parallel)
    {
        // Allocate variables
        uzero_loc = (double*)malloc(Nloc*sizeof(double));
        uval = (double*)malloc(Nloc*sizeof(double));
        Aloc_matrix = (double*)malloc((Nloc*Nloc)*sizeof(double));
        temp_Aloc = (double*)malloc((Nloc*Nloc)*sizeof(double));
        loc_pde_rhs_vec = (double*)malloc(Nloc*sizeof(double));
        sigmas_vec = (double*)malloc(Nloc*sizeof(double));
        sigmat_vec = (double*)malloc(Nloc*sizeof(double));
        loc_source_rhs_vec = (double*)malloc(Nloc*sizeof(double));

        Q_coeff = (double*)malloc((Nloc+1)*sizeof(double));
        R_coeff = (double*)malloc((Nloc+1)*sizeof(double));
        A_loc_QR = (double**)malloc((Nloc+1)*sizeof(double*));
        U_bar = (double*)malloc(SYSTEM_DIM*sizeof(double));
        U_bar_post = (double*)malloc(SYSTEM_DIM*sizeof(double));
        U_hat = (double*)malloc(SYSTEM_DIM*sizeof(double));
        s_epsilon = 1.0e-12;

        for(i=0;i<Nloc+1;i++)
	    A_loc_QR[i] = (double*)malloc((Nloc+1)*sizeof(double));

        if(!((*user_param).start_solution_file_flag))
        {
            // set up the initial solution vector coefficients
            
#pragma omp for schedule(dynamic)
            for(E=1;E<Nelts+1;E++)
            {
                init_zero_d(Aloc_matrix,Nloc*Nloc);
                //init_zero_m(A_loc_QR,Nloc+1,Nloc+1);
                Alocal_mat(E,mesh_element,mesh_node,user_param,Aloc_matrix);

                // Perform QR-decomposition
                assemble_qr(Aloc_matrix,Nloc,Q_coeff,R_coeff,A_loc_QR,&sing_decomp);

                for(s_dim=0;s_dim < SYSTEM_DIM;s_dim++)
                {
                    init_zero_d(uzero_loc,Nloc);
                    init_initial_conditions(E,mesh_node,mesh_element,mesh_edge,user_param,s_dim,uzero_loc);

                    // Solve the linear system
                    qrsolv(A_loc_QR,Nloc,Q_coeff,R_coeff,uzero_loc);

                    init_zero_d(Aloc_matrix,Nloc*Nloc);
                    //init_zero_m(A_loc_QR,Nloc+1,Nloc+1);
                    Alocal_mat(E,mesh_element,mesh_node,user_param,Aloc_matrix);

                    // Perform QR-decomposition
                    assemble_qr(Aloc_matrix,Nloc,Q_coeff,R_coeff,A_loc_QR,&sing_decomp);
                    init_zero_d(uval,Nloc);
                    project_true_sol_to_dg_space(E,mesh_node,mesh_element,mesh_edge,user_param,s_dim,zero_vec,uval);

                     // Solve the linear system
                    qrsolv(A_loc_QR,Nloc,Q_coeff,R_coeff,uval);

                    for(idofs=0;idofs<Nloc;idofs++)
                    {
                        uh_tn_sol_vec[SYSTEM_DIM*(E-1)*Nloc+(s_dim*Nloc+idofs)] = uzero_loc[idofs];
                        uh_tn1_sol_vec[SYSTEM_DIM*(E-1)*Nloc+(s_dim*Nloc+idofs)] = uzero_loc[idofs];
                        projected_vec[SYSTEM_DIM*(E-1)*Nloc+(s_dim*Nloc+idofs)] = uval[idofs];
                    }                    
                }//sys_dim
            }//loop over E
#pragma omp barrier  //sync threads (initial vector)
        }
        else //read initial solution vector from file
        {
#pragma omp master
            {
                load_sol_coeffs(uh_tn_sol_vec,mesh_element,mesh_node,init_coefs_file,Nloc,Nelts);
                for(E=1;E<Nelts+1;E++)
                {
                    for(s_dim=0;s_dim<SYSTEM_DIM;s_dim++)
                    {
                        for(idofs=0;idofs<Nloc;idofs++)
                        {
                                uh_tn1_sol_vec[SYSTEM_DIM*(E-1)*Nloc+(s_dim*Nloc+idofs)] = uh_tn_sol_vec[SYSTEM_DIM*(E-1)*Nloc+(s_dim*Nloc+idofs)];
                        }
                    }
                } 
                printf("Read initial state from file starting at time %lf \n",start_time);
                fclose(init_coefs_file);
            }
#pragma omp barrier
        }

#pragma omp for schedule(dynamic)
        for(E=1;E<Nelts+1;E++)
        {
            compute_solution_average(E,mesh_node,mesh_element,mesh_edge,uh_tn_sol_vec,U_bar);

        }
#pragma omp barrier
        

        /*store initial solution up to the first moment
          Make sure only one processor writes output*/
        #pragma omp master
        {
            if((*user_param).flag_store_solution)
            {
                save_solution_coefficients(uh_tn_sol_vec,mesh_element,mesh_node,user_param,Nloc,Nelts,0.0,"_0");
                plot_realizability(uh_tn_sol_vec,mesh_element,mesh_node,user_param,Nloc,Nelts,0.0,"_0");
            }

            //init collector files for paraview output
            if((*user_param).flag_store_paraview_solution_points)
            {
                init_point_collector(user_param);
            } 
            if((*user_param).flag_store_paraview_solution_cells)
            {
                init_cell_collector(user_param);
            }
        }//master thread write solution
        #pragma omp barrier


	#pragma omp for schedule(dynamic)
	for(E=1;E<Nelts+1;E++) {
	    
        if((*user_param).p_degree ==2)
        {
            project_quadr_sol_vec_to_lin(E, uh_tn_sol_vec, uh_pre_limit);
        }
        else if((*user_param).p_degree ==1)
        {
            project_solution_to_linear(E,mesh_element,mesh_node,user_param,uh_tn_sol_vec,uh_pre_limit);

        }
	}
	#pragma omp barrier

        //do limiting in parallel
        #pragma omp for schedule(dynamic)
        for (E=1;E<Nelts+1;E++)
        {
	    shu_osher_limit(E,mesh_element,mesh_node,mesh_edge,user_param,Nelts,Nloc,uh_pre_limit,uh_tn_sol_vec);
        }
        #pragma omp barrier

        if((*user_param).flag_ensure_positivity ==1)
        {
#pragma omp for schedule(dynamic)
            for(E=1;E<Nelts+1;E++)
            {
                compute_solution_average(E,mesh_node,mesh_element,mesh_edge,uh_tn_sol_vec,U_bar);
                rel_flag = is_realizable(U_bar);
                
                if(rel_flag ==0)
                {
                    //printf("ERROR @ E = %d uh_tn_sol_vec UBar %10.16e %10.16e %10.16e >>>>>>flag %d \n",E,U_bar[0],U_bar[1],U_bar[2],rel_flag);

                    //realizability_hack(U_bar);
                    //exit(1);
                }            

                check_realizability_on_gpts(E,U_bar,uh_tn_sol_vec,mesh_node,mesh_element,mesh_edge,s_epsilon,U_hat);

                compute_solution_average(E,mesh_node,mesh_element,mesh_edge,uh_tn_sol_vec,U_bar_post);
                if((fabs(U_bar[0] - U_bar_post[0]) > EPSILON) ||(fabs(U_bar[1] - U_bar_post[1]) > EPSILON)||(fabs(U_bar[2] - U_bar_post[2]) > EPSILON))
                {
                     printf("Error limiting changed average value @ E= %d diff = %10.12lf \n",E,fabs(U_bar[0] - U_bar_post[0]));
                     //exit(1);
                }
            }
#pragma omp barrier
        }

        //thread 0 updates the time
        #pragma omp master
        {
            current_time = start_time;
            printf("Computation starting at time t = %lf \n",current_time);
        }
        #pragma omp barrier


        //time stepping loop
        for(tsteps=1;tsteps< N_time_steps+1;tsteps++)
        {
#pragma omp master
            {
                if(((*user_param).p_degree ==2) || ((*user_param).p_degree==1)||((*user_param).p_degree==0))
                {
                    // backwards time-stepping
                    if ( (*user_param).flag_backwards_t_step) 
                    {
                        current_time1 = current_time - delta_t;
                        current_time2 = current_time - 0.5*delta_t;
                    }
                    else 
                    {
                        current_time1 = current_time + delta_t;
                        current_time2 = current_time + 0.5*delta_t;
                    }
                }
            }
#pragma omp barrier
            if(((*user_param).p_degree == 2)||((*user_param).p_degree == 1) ||((*user_param).p_degree == 0))

            {
#pragma omp for schedule(dynamic)
            for(E=1;E<Nelts+1;E++)
            {
                fractional_step(E,mesh_node,mesh_element,mesh_edge,user_param,current_time,delta_t,
                        uh_aux1,zero_vec,0.0,uh_tn_sol_vec,1.0,1.0,Nloc);     
            }//uh_aux1 (E)
#pragma omp barrier  //sync threads stage 1
#pragma omp master
            {
                //check_symmetry(mesh_element,Nelts,mesh_node,mesh_edge,uh_tn_sol_vec,0,Nloc,current_time,&symm_error);
 
                if((*user_param).symmetric_solution)
                {
                    //symmetry_correction(mesh_element,mesh_node,uh_aux1,Nelts,current_time,user_param);
                }
            }
#pragma omp barrier
//shu osher limiter
#pragma omp for schedule(dynamic)
	    for(E=1;E<Nelts+1;E++) 
        {
            if((*user_param).p_degree ==2)
            {
                project_quadr_sol_vec_to_lin(E, uh_aux1, uh_pre_limit);
            }
            else if((*user_param).p_degree ==1)
            {
                project_solution_to_linear(E,mesh_element,mesh_node,user_param,uh_tn_sol_vec,uh_pre_limit);
            }
	    }
#pragma omp barrier
#pragma omp for schedule(dynamic)
            for (E=1;E<Nelts+1;E++)
            {
                shu_osher_limit(E,mesh_element,mesh_node,mesh_edge,user_param,Nelts,Nloc,uh_pre_limit,uh_aux1);
            }
#pragma omp barrier
//realizability limiter
            if((*user_param).flag_ensure_positivity ==1)
            {
#pragma omp for schedule(dynamic)
                for(E=1;E<Nelts+1;E++)
                {
                    compute_solution_average(E,mesh_node,mesh_element,mesh_edge,uh_aux1,U_bar); 
                    rel_flag = is_realizable(U_bar);
                    if(rel_flag ==0)
                    {
                        //printf("ERROR @ E = %d uh_aux1 UBar %10.16e %10.16e %10.16e >>>>>>flag %d \n",E,U_bar[0],U_bar[1],U_bar[2],rel_flag);
                        //exit(1);
                    }

                    check_realizability_on_gpts(E,U_bar,uh_aux1,mesh_node,mesh_element,mesh_edge,s_epsilon,U_hat);
                    //printf("uh_aux1 check \n");
                    
                    compute_solution_average(E,mesh_node,mesh_element,mesh_edge,uh_aux1,U_bar_post);
                    if((fabs(U_bar[0] - U_bar_post[0]) > EPSILON) ||(fabs(U_bar[1] - U_bar_post[1]) > EPSILON)||(fabs(U_bar[2] - U_bar_post[2]) > EPSILON))
                    {
                        printf("Error limiting changed average value @ E= %d diff = %10.12lf \n",E,fabs(U_bar[0] - U_bar_post[0]));
                        //exit(1);
                    }
                }
#pragma omp barrier
            }

            /********************************stage 1**********************************************************/

#pragma omp for schedule(dynamic)
            for(E=1;E<Nelts+1;E++)
            {
               fractional_step(E,mesh_node,mesh_element,mesh_edge,user_param,current_time1,delta_t,
                        uh_aux2,uh_tn_sol_vec,0.75,uh_aux1,0.25,0.25,Nloc);    
            }//loop over elements
#pragma omp barrier //sync threads stage 2
#pragma omp master
            {

                //check_symmetry(mesh_element,Nelts,mesh_node,mesh_edge,uh_tn_sol_vec,0,Nloc,current_time,&symm_error);
                if((*user_param).symmetric_solution)
                {
                    //symmetry_correction(mesh_element,mesh_node,uh_aux2,Nelts,current_time,user_param);
                }
            }
#pragma omp barrier
//shu_osher limiter
#pragma omp for schedule(dynamic)
            for(E=1;E<Nelts+1;E++) 
            {
                if((*user_param).p_degree ==2)
                {
                    project_quadr_sol_vec_to_lin(E, uh_aux2, uh_pre_limit);
                }
                else if((*user_param).p_degree ==1)
                {
                    project_solution_to_linear(E,mesh_element,mesh_node,user_param,uh_tn_sol_vec,uh_pre_limit);
                }
            }
#pragma omp barrier
#pragma omp for schedule(dynamic)
            for (E=1;E<Nelts+1;E++)
            {
                shu_osher_limit(E,mesh_element,mesh_node,mesh_edge,user_param,Nelts,Nloc,uh_pre_limit,uh_aux2);
            }
#pragma omp barrier
            if((*user_param).flag_ensure_positivity ==1)
            {
#pragma omp for schedule(dynamic)
                for(E=1;E<Nelts+1;E++)
                {
                    compute_solution_average(E,mesh_node,mesh_element,mesh_edge,uh_aux2,U_bar);
                    rel_flag = is_realizable(U_bar);
                    if(rel_flag ==0)
                    {
                        //printf("ERROR @ E = %d  uh_aux2 UBar %10.16e %10.16e %10.16e >>>>>>flag %d \n",E,U_bar[0],U_bar[1],U_bar[2],rel_flag);
                        //exit(1);
                    }             
                    //printf("UBar %10.16e %10.16e %10.16e >>>>>>flag %d \n",U_bar[0],U_bar[1],U_bar[2],rel_flag);
                    //printf("i::::::::::::::::::::::::::::::::UBar %10.16e %10.16e %10.16e >>>>>>flag %d \n",U_bar[0],U_bar[1],U_bar[2],rel_flag);

                    check_realizability_on_gpts(E,U_bar,uh_aux2,mesh_node,mesh_element,mesh_edge,s_epsilon,U_hat);

                    //printf("uh_aux2 check \n");
                    
                    compute_solution_average(E,mesh_node,mesh_element,mesh_edge,uh_aux2,U_bar_post);
                    if((fabs(U_bar[0] - U_bar_post[0]) > EPSILON) ||(fabs(U_bar[1] - U_bar_post[1]) > EPSILON)||(fabs(U_bar[2] - U_bar_post[2]) > EPSILON))
                    {
                        printf("Error limiting changed average value @ E= %d diff = %10.12lf \n",E,fabs(U_bar[0] - U_bar_post[0]));
                        //exit(1);
                    }       
                } 
#pragma omp barrier
            }
 

            /*******************************stage 2************************************************************/

#pragma omp for schedule(dynamic)
            for(E=1;E< Nelts+1;E++)
            {
                fractional_step(E,mesh_node,mesh_element,mesh_edge,user_param,current_time2,delta_t,
                        uh_tn1_sol_vec,uh_tn_sol_vec,(1.0/3.0),uh_aux2,(2.0/3.0),(2.0/3.0),Nloc);
            }//loop over elements
#pragma omp barrier //sync threads stage 3
#pragma omp master
            {
                //check_symmetry(mesh_element,Nelts,mesh_node,mesh_edge,uh_tn_sol_vec,0,Nloc,current_time,&symm_error);
                if((*user_param).symmetric_solution)
                {
                    //symmetry_correction(mesh_element,mesh_node,uh_tn1_sol_vec,Nelts,current_time,user_param);
                }
            }
#pragma omp barrier
#pragma omp for schedule(dynamic)
	    for(E=1;E<Nelts+1;E++) 
        {
            if((*user_param).p_degree ==2)
            {
                project_quadr_sol_vec_to_lin(E, uh_tn1_sol_vec, uh_pre_limit);
            }
            else if((*user_param).p_degree ==1)
            {
                project_solution_to_linear(E,mesh_element,mesh_node,user_param,uh_tn_sol_vec,uh_pre_limit);
            }	    
        }
#pragma omp barrier
#pragma omp for schedule(dynamic)
            for (E=1;E<Nelts+1;E++)
            {
                shu_osher_limit(E,mesh_element,mesh_node,mesh_edge,user_param,Nelts,Nloc,uh_pre_limit,uh_tn1_sol_vec);
            }
#pragma omp barrier

#pragma omp barrier
            if((*user_param).flag_ensure_positivity ==1)
            {
#pragma omp for schedule(dynamic)
                for(E=1;E<Nelts+1;E++)
                {
                    compute_solution_average(E,mesh_node,mesh_element,mesh_edge,uh_tn1_sol_vec,U_bar);
                    rel_flag = is_realizable(U_bar);
                    //printf("UBar %10.16e %10.16e %10.16e >>>>>>flag %d \n",U_bar[0],U_bar[1],U_bar[2],rel_flag);
                    cell_average(uh_tn1_sol_vec,E,mesh_element,mesh_node,U_bar);

                    //printf("i::::::::::::::::::::::::::::::::UBar %10.16e %10.16e %10.16e >>>>>>flag %d \n",U_bar[0],U_bar[1],U_bar[2],rel_flag);
                    if(rel_flag ==0)
                    {
                        //printf("ERROR @ E = %d  uh_tn1_sol_vec UBar %10.16e %10.16e %10.16e >>>>>>flag %d \n",E,U_bar[0],U_bar[1],U_bar[2],rel_flag);
                        //exit(1);
                    }             

                    check_realizability_on_gpts(E,U_bar,uh_tn1_sol_vec,mesh_node,mesh_element,mesh_edge,s_epsilon,U_hat);

                //    printf("uh_tn1_sol_check  \n");
         
                    compute_solution_average(E,mesh_node,mesh_element,mesh_edge,uh_tn1_sol_vec,U_bar_post);
                    if((fabs(U_bar[0] - U_bar_post[0]) > EPSILON) ||(fabs(U_bar[1] - U_bar_post[1]) > EPSILON)||(fabs(U_bar[2] - U_bar_post[2]) > EPSILON))
                    {
                        printf("Error limiting changed average value @ E= %d diff = %10.12lf \n",E,fabs(U_bar[0] - U_bar_post[0]));
                        //exit(1);
                    }   
                }     
#pragma omp barrier
            }
 
            /********************************stage 3***********************************************************/
        }
        // Make sure only one processor writes output
        // store intermediate solution
            #pragma omp master
            {
                if (tsteps >= (ind_out*N_time_steps)/((*user_param).no_out-1))
                {
                    sprintf(file_out,"_%d",ind_out+1);
                   // printf("Stage %d \n",ind_out+1);
                    if((*user_param).symmetric_solution)
                    {
                        //printf("current_time %lf \n",current_time);
                        check_symmetry(mesh_element,Nelts,mesh_node,mesh_edge,uh_tn_sol_vec,0,Nloc,current_time,&symm_error);
                    }
                    // printf("symm_error = %10.12e \n",symm_error);
                    if((*user_param).flag_store_solution)
                    {
                         save_solution_coefficients(uh_tn_sol_vec,mesh_element,mesh_node,user_param,Nloc,Nelts,current_time,file_out);
                         plot_realizability(uh_tn_sol_vec,mesh_element,mesh_node,user_param,Nloc,Nelts,current_time,file_out);
                    }
                    if((*user_param).flag_store_paraview_solution_points)
                    {
                         save_point_solution_vtp(uh_tn_sol_vec,mesh_element,mesh_node,user_param,Nelts,num_nodes,file_out);
                         // error compared to analytical solution (not always available)
                         save_point_solution_error_vtp(uh_tn_sol_vec,mesh_element,mesh_node,user_param,Nelts,num_nodes,current_time,file_out);
                         write_point_collector(user_param,current_time,file_out);
                    }
	            if((*user_param).flag_store_paraview_solution_cells)
                    {
                         save_cell_solution_vtp(uh_tn_sol_vec,mesh_element,mesh_node,user_param,Nelts,num_nodes,file_out);
                         write_cell_collector(user_param,current_time,file_out);
                    }
                    if(ind_out==(*user_param).no_out-1)
                    {
                         if((*user_param).flag_store_paraview_solution_points)
                         {
                              end_point_collector(user_param);
                         } 
                         if((*user_param).flag_store_paraview_solution_cells)
                         {
                              end_cell_collector(user_param);
                         } 
                    }

                    ind_out += 1;
                }
                for(i=0;i<size;i++)
                {
                    uh_tn_sol_vec[i] = uh_tn1_sol_vec[i];
                    //printf("uh_tn_sol %lf \n", uh_tn_sol_vec[i]);
                }
                if(fabs((*user_param).Ftime - current_time) < delta_t)
                {
                    delta_t = fabs((*user_param).Ftime-current_time);
                }

                // backwards time-stepping
                if ( (*user_param).flag_backwards_t_step) 
                {
                    current_time = current_time - delta_t;
                }
                else 
                {
                    current_time = current_time + delta_t;
                }
            }
            #pragma omp barrier

            if (verbose_parallel == 1)
                printf("thread %d reached end of timestep %d\n",omp_get_thread_num(),tsteps);

        }//end loop timesteps
        #pragma omp barrier

        #pragma omp master
        {
            save_error_coeffs(uh_tn_sol_vec,projected_vec,user_param,Nelts,Nloc);
        }

	// Free allocated memory
	for(i=0;i<Nloc+1;i++)
	{
	    free(A_loc_QR[i]);
	}
	free(A_loc_QR);
	free(R_coeff);
	free(Q_coeff);
	free(loc_source_rhs_vec);
	free(sigmat_vec);
	free(sigmas_vec);
	free(loc_pde_rhs_vec);
	free(temp_Aloc);
	free(Aloc_matrix);
	free(uzero_loc);
        free(uval);
    } // end parallelized region

    free(uh_pre_limit);
    free(uh_aux1);
    free(uh_aux2);
    free(uh_aux3);
    free(uh_aux4);
    free(uh_aux5);
    free(uh_aux6);
    free(uh_aux7);
    free(uh_aux8);
    free(uh_aux9);
    free(projected_vec);
    free(zero_vec);
    //free(uh_limit);
    free(uh_tn_sol_vec);
}//rk_dg

/* void discrete_pde_rhs
 * assembles right-hand side Lh locally
 * in: index to triangle           -- E
 * in: degrees of freedom	   -- Nloc
 * in: mesh_node data structure    -- mesh_node
 * in: mesh element data structure -- mesh_element
 * in: mesh edge data structure    -- mesh_edge
 * in: parameters 		   -- user_param
 * in: solution vector at time tn  -- uh_tn_sol_vec
 * in: index to moment of system   -- s_dim_count
 * in: time tn 			   -- t_val
 * out: rhs of discretized system  -- pde_rhs_vec
 */
void discrete_pde_rhs(int E,
		              int Nloc,
                      node* mesh_node,
                      element* mesh_element,
                      edge* mesh_edge,
                      parameters* user_param,
                      double* uh_tn_sol_vec,
                      int s_dim_count,
		              double t_val,
                      double* pde_rhs_vec_E)
{
    double* vol_fvec1 = NULL;               //vector for int_E((f(uh)\cdot\nabla v(x))dE
    double* edge_fvec2 = NULL;              //vector for sum over all edges int_e( f(uh)\cdot \nabla v(x))dx 
    double* iedge_fvec2 = NULL;             //vector for interior and boundary edge int_e( f(uh)\cdot \nabla v(x))dx
    int i=0;                                //counter over edges
    int iedge=0;                            //edge counter
    int idofs=0;                            //degrees of freedom counter

    //allocate memory for local rhs functions on each element E
    vol_fvec1 = (double*)malloc(Nloc*sizeof(double));
    edge_fvec2 = (double*)malloc(Nloc*sizeof(double));
    iedge_fvec2 = (double*)malloc(Nloc*sizeof(double));

    init_zero_d(vol_fvec1,Nloc);
    init_zero_d(edge_fvec2,Nloc);
    init_zero_d(iedge_fvec2,Nloc);
    flocal_vec1(E,t_val,Nloc,mesh_element,mesh_node,user_param,uh_tn_sol_vec,s_dim_count,vol_fvec1);

    //loop over the edges
    for(i=1;i<4;i++) 
    {
        iedge = mesh_element[E].edge[i];
        if(mesh_edge[iedge].edge_type == INTERIOR) 
        {
            init_zero_d(iedge_fvec2,Nloc);
            clocal_edge_loc(E,Nloc,iedge,mesh_element,mesh_node,mesh_edge,user_param,uh_tn_sol_vec,s_dim_count,t_val,iedge_fvec2);
        }

        else if(mesh_edge[iedge].edge_type == EXTERIOR) 
        {
            init_zero_d(iedge_fvec2,Nloc);
            clocal_bdry_edge_loc(E,Nloc,iedge,mesh_element,mesh_node,mesh_edge,user_param,uh_tn_sol_vec,s_dim_count,t_val,iedge_fvec2);
        }

        for(idofs=0;idofs<Nloc;idofs++) 
        {
            edge_fvec2[idofs] = edge_fvec2[idofs] + iedge_fvec2[idofs];
            
        }

    }//loop over edges

    //printf("E %d \n",E);
    for(idofs=0;idofs<Nloc;idofs++) 
    {
        pde_rhs_vec_E[idofs] = vol_fvec1[idofs]-edge_fvec2[idofs];


        //if((E==2) || (E==24))
        //{
       
        
        
        //printf("vol %lf edge %lf \n",vol_fvec1[idofs],edge_fvec2[idofs]);
    }
    //getchar();

    free(vol_fvec1);
    free(edge_fvec2);
    free(iedge_fvec2);
}

/* void Alocal_mat
 * sets up local stiffness matrix
 * in: index to triangle -- E
 * in: mesh_element data structure -- mesh_element
 * in: mesh_node data structure    -- mesh_node
 * in: parameters 		   -- user_param
 * out: stiffness matrix	   -- Aloc_matrix
 */
void Alocal_mat(int E, 
		element* mesh_element,
		node* mesh_node,
		parameters* user_param,
		double* Aloc_matrix)
{
    //local vars
    int k=0;
    int idofs=0;
    int jdofs=0;
    int deg=0;
    int Nloc=0;

    //get polynomial degree
    deg = mesh_element[E].degree;
    Nloc = (deg+1)*(deg+2)/2;

    //loop over the gauss points
    for(k=1;k<ngpts+1;k++) 
    {
        for(idofs =1; idofs < Nloc+1;idofs++) 
        {
            for(jdofs =1;jdofs < Nloc+1;jdofs++) 
            {
                Aloc_matrix(idofs,jdofs,Nloc) = Aloc_matrix(idofs,jdofs,Nloc)+
                mesh_element[E].el_gpts_basis[(k-1)*Nloc+(idofs-1)]
                *mesh_element[E].el_gpts_basis[(k-1)*Nloc+(jdofs-1)]
                *mesh_element[E].el_gpts_w[k]*mesh_element[E].det;
            }//j
        }//i
    }//loop over quad points
}//Alocal_mat

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
void  project_true_sol_to_dg_space(int E,
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

    /*init_zero_d(u0_val,SYSTEM_DIM);*/

    for(k=1;k<ngpts+1;k++) 
    {
        //bdry_function_val(mesh_element[E].el_gpts_x[k],mesh_element[E].el_gpts_y[k],mesh_element[E].el_gpts_ref_x[k],mesh_element[E].el_gpts_ref_y[k],
         //       (*user_param).Ftime,E,mesh_element,mesh_node,uh_tn_sol_vec,user_param,u_val);

        for(idofs=0;idofs<Nloc;idofs++) 
        {
	    vec_uval[idofs] = vec_uval[idofs] + u_val[s_dim_count]*mesh_element[E].el_gpts_basis[(k-1)*Nloc+idofs]*det_BE_T*mesh_element[E].el_gpts_w[k];
        }//idofs
    }//ngpts
}//initial conditions


/* void init_initial_conditions
 * sets up initial conditions
 * in: index to triangle -- E
 * in: mesh_node data structure    -- mesh_node
 * in: mesh_element data structure -- mesh_element
 * in: mesh_edge data structure    -- mesh_edge
 * in: parameters 		   -- user_param
 * in: s_dim_count 	           -- dimension count
 * out: local initional conditions -- u_zero_loc
 */
void  init_initial_conditions(int E,
                             node* mesh_node,
                             element* mesh_element,
                             edge* mesh_edge,
                             parameters* user_param,
                             int s_dim_count,
                             double* u_zero_loc)
{
    int k=0;
    int idofs=0;

    int deg=0;
    int Nloc=0;

    double det_BE_T =0.0;

    double u0_val[SYSTEM_DIM];

    det_BE_T = mesh_element[E].det;

    deg = mesh_element[E].degree;
    Nloc = ((deg+1)*(deg+2))/2;

    /*init_zero_d(u0_val,SYSTEM_DIM);*/

    for(k=1;k<ngpts+1;k++) 
    {
        u_zero_val(mesh_element[E].el_gpts_x[k],mesh_element[E].el_gpts_y[k],user_param,u0_val);

        for(idofs=0;idofs<Nloc;idofs++) 
        {
            u_zero_loc[idofs] = u_zero_loc[idofs] + u0_val[s_dim_count]*mesh_element[E].el_gpts_basis[(k-1)*Nloc+idofs]*det_BE_T*mesh_element[E].el_gpts_w[k];
            //printf("u_zero_loc %lf u0_val %lf basis %lf det %lf w %lf \n",u_zero_loc[idofs],u0_val[s_dim_count],mesh_element[E].el_gpts_basis[(k-1)*Nloc+idofs],det_BE_T,mesh_element[E].el_gpts_w[k]);
            //getchar();
            
        }//idofs
    }//ngpts

}//initial conditions

/* void flocal_vec1 
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
void flocal_vec1(int E,
		 double t,
		 int Nloc,
                 element* mesh_element,
                 node* mesh_node,
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
    double flux_local_uh_tn_sol[2];

    gpx = (double*)malloc((ngpts+1)*sizeof(double));
    gpy = (double*)malloc((ngpts+1)*sizeof(double));
    w = (double*)malloc((ngpts+1)*sizeof(double));

    init_zero_d(gpx,ngpts+1);
    init_zero_d(gpy,ngpts+1);
    init_zero_d(w,ngpts+1);
    int rel_flag=0;

    //gauss points on reference element
    gauss_pt(gpx,gpy,w);

    //get polynomial degree
    deg = mesh_element[E].degree;

    //allocate memory for basis functions
    phi_vec = (double*)malloc(Nloc*sizeof(double));
    grad_phi_mat = (double*)malloc(2*Nloc*sizeof(double));

    init_zero_d(local_uh_tn_sol,SYSTEM_DIM);

    //get mappint from reference element to physical element
    det_BE_T = mesh_element[E].det;

    for(k=1;k<ngpts+1;k++) 
    {
	init_zero_d(grad_phi_mat,2*Nloc);
	init_zero_d(phi_vec,Nloc);
	init_zero_d(flux_local_uh_tn_sol,2);

	//for each gauss point get values of phi and grad_phi
	init_monomial_basis_deriv(E,mesh_element,mesh_node,deg,gpx[k],gpy[k],phi_vec,grad_phi_mat);

	//map reference coordinates to physical coordinates
	//solution on E at gauss points 
	get_approx_solution(uh_tn_sol_vec,E,mesh_element,mesh_node,gpx[k],gpy[k],local_uh_tn_sol);

	// Calculate flux

    rel_flag = is_realizable(local_uh_tn_sol);
    double ratio =0;
    
    /*if(rel_flag ==0)
    {
        printf("local_uh_tn_sol %10.16e %10.16e %10.16e \n",local_uh_tn_sol[0],local_uh_tn_sol[1],local_uh_tn_sol[2]);
        ratio = sqrt(pow(local_uh_tn_sol[1],2.0) + pow(local_uh_tn_sol[2],2.0))/(local_uh_tn_sol[0]);
        printf("norm psi one %10.16e s_dim %d \n",sqrt(pow(local_uh_tn_sol[1],2.0) + pow(local_uh_tn_sol[2],2.0)),s_dim_counter);
        
        printf("ratio %lf \n",ratio);
        printf("Element %d E, k %d  @ %lf %lf \n",E,k,gpx[k],gpx[k]);
    }*/
    //assert(rel_flag);
        
    flux_function(local_uh_tn_sol,mesh_element[E].el_gpts_x[k],mesh_element[E].el_gpts_y[k],t,mesh_element[E].el_density[k],user_param,s_dim_counter,flux_local_uh_tn_sol);

        
  
 
        for(idofs =1;idofs<Nloc+1;idofs++)
        {
            floc_vec1[idofs-1] = floc_vec1[idofs-1]
                + w[k]*det_BE_T*((flux_local_uh_tn_sol[0]*grad_phi_mat(1,idofs,2))+
                         	 (flux_local_uh_tn_sol[1]*grad_phi_mat(2,idofs,2)));

        }//loop over idofs
    }//loop over quadrature nodes
    
    free(gpx);
    free(gpy);
    free(w);
    free(grad_phi_mat);
    free(phi_vec);
}//flocal_vec1

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
 void clocal_bdry_edge_loc(int E,
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

    for(k=1;k<nlgpts+1; k++) 
    {
        get_approx_solution(uh_tn_sol_vec,E1,mesh_element,mesh_node,mesh_edge[iedge].E1_ed_gpts_x[k],mesh_edge[iedge].E1_ed_gpts_y[k],approx_sol_E1);

        if((*user_param).periodic)
        {
            p_edge = mesh_edge[iedge].periodic_nedge;
            
            if(mesh_edge[iedge].boundary_side == BOTTOM)
            {
                p_edge_phys_coords[0] = mesh_edge[iedge].ed_phys_coords_x[k];
                p_edge_phys_coords[1] = mesh_edge[iedge].ed_phys_coords_y[k] + (*user_param).periodic_domain_size_y;
            }
            else if(mesh_edge[iedge].boundary_side == TOP)
            {
                p_edge_phys_coords[0] = mesh_edge[iedge].ed_phys_coords_x[k];
                p_edge_phys_coords[1] = mesh_edge[iedge].ed_phys_coords_y[k] - (*user_param).periodic_domain_size_y;
            }
            else if(mesh_edge[iedge].boundary_side == LEFT)
            {
                p_edge_phys_coords[0] = mesh_edge[iedge].ed_phys_coords_x[k] + (*user_param).periodic_domain_size_x;
                p_edge_phys_coords[1] = mesh_edge[iedge].ed_phys_coords_y[k];
            }
            else if(mesh_edge[iedge].boundary_side == RIGHT)
            {
                p_edge_phys_coords[0] = mesh_edge[iedge].ed_phys_coords_x[k] - (*user_param).periodic_domain_size_x;
                p_edge_phys_coords[1] = mesh_edge[iedge].ed_phys_coords_y[k];
            }
            else
            {
                fprintf(stderr,"Error in phys coords translation \n");
                exit(1);
            }

            map_to_reference_element(mesh_element,mesh_node,E2,p_edge_ref_coords,p_edge_phys_coords[0],p_edge_phys_coords[1]);
            get_approx_solution(uh_tn_sol_vec,E2,mesh_element,mesh_node,p_edge_ref_coords[0],p_edge_ref_coords[1],approx_sol_E2);

            lax_friedrichs_flux(mesh_edge,iedge,mesh_element[E].el_density[k],approx_sol_E1,approx_sol_E2,normal_e,mesh_edge[iedge].E1_ed_gpts_x[k],
                     mesh_edge[iedge].E1_ed_gpts_y[k],user_param,s_dim_count,t_val,p_edge_phys_coords[0],p_edge_phys_coords[1],
                     mesh_edge[iedge].ed_phys_coords_x[k],mesh_edge[iedge].ed_phys_coords_y[k],&numerical_flux);
        }
        else if((*user_param).reflective_bc)
        {
            if(mesh_edge[iedge].reflective)
            {
                //use first element
		/*printf("normal = [ %lf, %lf ]\n",normal_e[0],normal_e[1]);*/
                reflec_vec[0] = approx_sol_E1[0];
                reflec_vec[1] = approx_sol_E1[1] - 2.0*(normal_e[0]*approx_sol_E1[1] + normal_e[1]*approx_sol_E1[2])*normal_e[0];
                reflec_vec[2] = approx_sol_E1[2] - 2.0*(normal_e[0]*approx_sol_E1[1] + normal_e[1]*approx_sol_E1[2])*normal_e[1];
                lax_friedrichs_flux(mesh_edge,iedge,mesh_element[E].el_density[k],approx_sol_E1,reflec_vec,normal_e,mesh_edge[iedge].E1_ed_gpts_x[k],
                    mesh_edge[iedge].E1_ed_gpts_y[k],user_param,s_dim_count,t_val,0.0,0.0,mesh_edge[iedge].ed_phys_coords_x[k],mesh_edge[iedge].ed_phys_coords_y[k],
                    &numerical_flux);
            }
            else
            {
                bdry_function_val(mesh_edge[iedge].ed_phys_coords_x[k],mesh_edge[iedge].ed_phys_coords_y[k],mesh_edge[iedge].E1_ed_gpts_x[k], 
                        mesh_edge[iedge].E1_ed_gpts_y[k],t_val,E,iedge,mesh_element,mesh_node,mesh_edge,uh_tn_sol_vec,user_param,bdry_val);
                
                lax_friedrichs_flux(mesh_edge,iedge,mesh_element[E].el_density[k],approx_sol_E1,bdry_val,normal_e,mesh_edge[iedge].E1_ed_gpts_x[k], 
                        mesh_edge[iedge].E1_ed_gpts_y[k],user_param,s_dim_count,t_val,0.0,0.0,mesh_edge[iedge].ed_phys_coords_x[k],
                        mesh_edge[iedge].ed_phys_coords_y[k], &numerical_flux);
            }
        }
        else
        {
            bdry_function_val(mesh_edge[iedge].ed_phys_coords_x[k],mesh_edge[iedge].ed_phys_coords_y[k],mesh_edge[iedge].E1_ed_gpts_x[k],mesh_edge[iedge].E1_ed_gpts_y[k],
                    t_val,E,iedge,mesh_element,mesh_node,mesh_edge,uh_tn_sol_vec,user_param,bdry_val);

            if((E==2) || (E==24))
            {
                if(s_dim_count == 0)
                {
                    //printf("k %d bdry_val =  %10.32lf @ E = %d \n",k,bdry_val[0],E);
                }
            }

            //use first element
            lax_friedrichs_flux(mesh_edge,iedge,mesh_element[E].el_density[k],approx_sol_E1,bdry_val,normal_e,mesh_edge[iedge].E1_ed_gpts_x[k],
                    mesh_edge[iedge].E1_ed_gpts_y[k],user_param,s_dim_count,t_val,0.0,0.0,mesh_edge[iedge].ed_phys_coords_x[k],
                    mesh_edge[iedge].ed_phys_coords_y[k],&numerical_flux);
            //printf("numerical_flux %lf \n",numerical_flux);getchar();


        }

        for(idofs=0;idofs<Nloc;idofs++)
        {
            floc_vec2[idofs] = floc_vec2[idofs] + mesh_edge[iedge].ed_gpts_w[k]*mesh_edge[iedge].ed_gpts_basis_E1[(k-1)*Nloc+idofs]*numerical_flux*iedge_length;
            if((E == 2) || (E==24))
            {
                if(s_dim_count ==0)
                {
                    //printf("k= %d E %d numerical flux  %10.32lf floc_vec2[idofs] =  %10.32lf \n",k,E,numerical_flux,floc_vec2[idofs]);
                  
                  //  printf("E %d numerical flux  %10.32lf basis =  %10.32lf \n",E,numerical_flux,mesh_edge[iedge].ed_gpts_basis_E1[(k-1)*Nloc+idofs]);
 
                }
            }
        }//loop over idofs
       
        /*if((E==2)||(E==24))
        {
           if(s_dim_count==0)
                //getchar();
        }*/

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
void clocal_edge_loc(int E,
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


    for(k=1;k<nlgpts+1; k++)
    {
        get_approx_solution(uh_tn_sol_vec,E1,mesh_element,mesh_node,mesh_edge[iedge].E1_ed_gpts_x[k],mesh_edge[iedge].E1_ed_gpts_y[k],approx_sol_E1);
        get_approx_solution(uh_tn_sol_vec,E2,mesh_element,mesh_node,mesh_edge[iedge].E2_ed_gpts_x[k],mesh_edge[iedge].E2_ed_gpts_y[k],approx_sol_E2);

        if(E==E1)
        {
            lax_friedrichs_flux(mesh_edge,iedge,mesh_element[E].el_density[k],approx_sol_E1,approx_sol_E2,normal_e,mesh_edge[iedge].E1_ed_gpts_x[k],
                    mesh_edge[iedge].E1_ed_gpts_y[k],user_param,s_dim_count,t_val,0.0,0.0,mesh_edge[iedge].ed_phys_coords_x[k],mesh_edge[iedge].ed_phys_coords_y[k],
                    &numerical_flux);
        }
        else if(E==E2)
        {
            lax_friedrichs_flux(mesh_edge,iedge,mesh_element[E].el_density[k],approx_sol_E2,approx_sol_E1,normal_e,mesh_edge[iedge].E2_ed_gpts_x[k],
                    mesh_edge[iedge].E2_ed_gpts_y[k],user_param,s_dim_count,t_val,0.0,0.0,mesh_edge[iedge].ed_phys_coords_x[k],mesh_edge[iedge].ed_phys_coords_y[k],
                    &numerical_flux);
        }

        for(idofs=0;idofs<Nloc;idofs++)
        {
            if(E==E1)
            {
                floc_vec2[idofs] = floc_vec2[idofs] + mesh_edge[iedge].ed_gpts_w[k]*mesh_edge[iedge].ed_gpts_basis_E1[(k-1)*Nloc+idofs]*numerical_flux*iedge_length;
                
            }
            else if(E==E2)
            {
                floc_vec2[idofs] = floc_vec2[idofs] + mesh_edge[iedge].ed_gpts_w[k]*mesh_edge[iedge].ed_gpts_basis_E2[(k-1)*Nloc+idofs]*numerical_flux*iedge_length;

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

/* void lax_friedrichs_flux
 * evaluates Lax-Friedrichs flux
 * in: edge over which the flux is evaluated -- edge
 * in: density value  				 -- density
 * in: values of flux-function in first  element -- a
 * in: values of flux-function in second element -- b
 * in: normal vector between the two cells -- normal_e
 * in: reference x coordinate -- xhat
 * in: reference y coordinate -- yhat
 * in: parameters 	      -- user_param
 * in: s_dim_count            -- dimension of system counter
 * in: time tn 	 	      -- t_val
 * in: p_edge_phys_coords x --- for use with periodic bc
 * in: p_edge_phys_coords y --- for use with periodic bc
 * in: physical x coordinate  -- phys_coords_x
 * in: physical y coordinate  -- phys_coords_y
 * out: numerical flux -- numerical flux
 */
void lax_friedrichs_flux(edge* mesh_edge,
                         int edge,
			             double density,
                         double* approx_sol_element_a,
                         double* approx_sol_element_b,
                         double* normal_e,
                         double xhat,
                         double yhat,
                         parameters* user_param,
                         int s_dim_counter,
			             double t_val,
                         double p_edge_phys_coords_x,
                         double p_edge_phys_coords_y,
			             double phys_coords_x,
			             double phys_coords_y,
                         double* numerical_flux)
{
    double flux_val_a[2];
    double flux_val_b[2];
    double adv_vector[2];

    double num_flux_val=0.0;
    init_zero_d(flux_val_a,2);
    init_zero_d(flux_val_b,2);

    double ALPHA_EK=0.0;

    ALPHA_EK=1.0;

    //for closure 7
    /*ALPHA_EK = max(fabs(approx_sol_element_a[0]),fabs(approx_sol_element_b[0]));*/

    flux_function(approx_sol_element_a,phys_coords_x,phys_coords_y,t_val,density,user_param,s_dim_counter,flux_val_a);
    if(((*user_param).periodic) && (mesh_edge[edge].edge_type == EXTERIOR))
    {
        flux_function(approx_sol_element_b,p_edge_phys_coords_x,p_edge_phys_coords_y,t_val,density,user_param,s_dim_counter,flux_val_b);
    }
    else if(((*user_param).reflective_bc) && (mesh_edge[edge].edge_type == EXTERIOR) &&(mesh_edge[edge].reflective))
    {
        flux_function(approx_sol_element_b,phys_coords_x,phys_coords_y,t_val,density,user_param,s_dim_counter,flux_val_b);
    }
    else
    {
        flux_function(approx_sol_element_b,phys_coords_x,phys_coords_y,t_val,density,user_param,s_dim_counter,flux_val_b);
    }

       
    num_flux_val = 0.5*((flux_val_a[0]*normal_e[0] + flux_val_a[1]*normal_e[1]) +
                (flux_val_b[0]*normal_e[0] + flux_val_b[1]*normal_e[1]) - ALPHA_EK*(approx_sol_element_b[s_dim_counter]-approx_sol_element_a[s_dim_counter]));
 


    *numerical_flux = num_flux_val;
}

/* void get_approx_solution
 * returns approximate solution value
 * in: solution at time tn -- uh_tn_sol
 * in: index of triangle -- E
 * in: mesh element data structure -- mesh_element
 * in: mesh node data structure -- mesh_node
 * in: reference x coordinate -- xhat 
 * in: reference y coordinate -- yhat
 * out: the solution at tn,xhat,yhat -- local_approx_sol
 */
void get_approx_solution(double* uh_tn_sol,
                         int E,
			 element* mesh_element, 
			 node* mesh_node,
                         double xhat, 
                         double yhat,
                         double* local_approx_sol)
{
    int Nloc=0;                     //dimension of local system Nloc = ((k+1)(k+2))/2
    int i =0;                       //loop over Nloc
    double* phi_vec = NULL;         //value of basis functions at quad points
    int deg=0;                      //polynomial degree of approx
    double temp =0.0;               //temp var to store sol during computation
    int s_dim=0;

    deg = mesh_element[E].degree;
    Nloc =((deg+1)*(deg+2))/2;

    phi_vec = (double*)malloc(Nloc*sizeof(double));

    // initialise vectors to zero
    init_zero_d(phi_vec,Nloc);

    init_monomial_basis(E,mesh_element,mesh_node,deg,xhat,yhat,phi_vec);	

    for(s_dim=0;s_dim< SYSTEM_DIM;s_dim++)
    {
	temp=0.0;
	//loop over monomials
	for(i=0;i< Nloc;i++) 
	{
	    temp = temp+ uh_tn_sol[SYSTEM_DIM*Nloc*(E-1)+(s_dim*Nloc)+i]*phi_vec[i];
	}
	local_approx_sol[s_dim] = temp;
    }//loop over dimensions

    free(phi_vec);
}//get_approx_sol

/* void get_approx_solution_lin
 * returns approximate solution value of a linear function
 * in: solution at time tn -- uh_tn_sol
 * in: index of triangle -- E
 * in: mesh element data structure -- mesh_element
 * in: mesh node data structure -- mesh_node
 * in: reference x coordinate -- xhat 
 * in: reference y coordinate -- yhat
 * out: the solution at tn,xhat,yhat -- local_approx_sol
 */
void get_approx_solution_lin(double* uh_tn_sol,
                             int E,
			     element* mesh_element, 
			     node* mesh_node,
			     double xhat, 
			     double yhat,
			     double* local_approx_sol)
{
    int i =0;                       //loop over Nloc
    double* phi_vec = NULL;         //value of basis functions at quad points
    double temp =0.0;               //temp var to store sol during computation
    int s_dim=0;

    phi_vec = (double*)malloc(3*sizeof(double));

    // initialise vectors to zero
    init_zero_d(phi_vec,3);

    init_monomial_basis(E,mesh_element,mesh_node,1,xhat,yhat,phi_vec);	

    for(s_dim=0;s_dim< SYSTEM_DIM;s_dim++)
    {
	temp=0.0;
	//loop over monomials
	for(i=0;i< 3;i++) 
	{
	    temp = temp+ uh_tn_sol[SYSTEM_DIM*3*(E-1)+(s_dim*3)+i]*phi_vec[i];
	}
	local_approx_sol[s_dim] = temp;
    }//loop over dimensions

    free(phi_vec);
}//get_approx_sol

/* void get_approx_sol_one_mom
 * returns approximate solution value
 * in: solution at time tn -- uh_tn_sol
 * in: index of triangle -- E
 * in: mesh element data structure -- mesh_element
 * in: mesh node data structure -- mesh_node
 * in: reference x coordinate -- xhat 
 * in: reference y coordinate -- yhat
 * out: the solution at tn,xhat,yhat -- local_approx_sol
 */
void get_approx_solution_one_mom(double* uh_tn_sol,
                                 int E,
				 element* mesh_element, 
				 node* mesh_node,
                                 double xhat, 
                                 double yhat,
                                 double* local_approx_sol)
{
    int Nloc=0;                     //dimension of local system Nloc = ((k+1)(k+2))/2
    int i =0;                       //loop over Nloc
    double* phi_vec = NULL;         //value of basis functions at quad points
    int deg=0;                      //polynomial degree of approx
    int sol_index=0;                //position of solution in global solution vector
    double temp =0.0;               //temp var to store sol during computation

    deg = mesh_element[E].degree;
    Nloc =((deg+1)*(deg+2))/2;

    phi_vec = (double*)malloc(Nloc*sizeof(double));

    // initialise vectors to zero
    init_zero_d(phi_vec,Nloc);

    init_monomial_basis(E,mesh_element,mesh_node,deg,xhat,yhat,phi_vec);	

    for(i=1;i< Nloc+1;i++) 
    {
	sol_index = Nloc*(E-1)+(i-1);
	temp = temp+ uh_tn_sol[sol_index]*phi_vec[i-1];
    }

    *local_approx_sol = temp;
    free(phi_vec);
}//get_approx_sol

/* void calc_normal_face
 * computes the unit normal vector pointing from E1 to E2
 * in: current edge -- iedge        
 * in: node data structure -- mesh_node    
 * in: edge data structure -- mesh_edge    
 * in: normal pointing in right direction? -- normal_direction_flag
 * out: normal vector -- normal_vec
 */
void calc_normal_face(int iedge,
		      node* mesh_node, 
		      edge* mesh_edge,
		      int normal_direction_flag,
		      double* normal_vec)
{
    double xc1,yc1,xc2,yc2;
    int node_a,node_b;
    double area =0.0;
    area = calc_length_face(mesh_edge,iedge,mesh_node);
    node_a = mesh_edge[iedge].vertex[1];
    node_b = mesh_edge[iedge].vertex[2];
    xc1 = mesh_node[node_a].coord[0];
    yc1 = mesh_node[node_a].coord[1];
    xc2 = mesh_node[node_b].coord[0];
    yc2 = mesh_node[node_b].coord[1];

    //unit normal vector
    if(normal_direction_flag == 0) 
    {
        normal_vec[0] = (yc2 -yc1)/area;
        normal_vec[1] = -(xc2-xc1)/area;
    }
    else if(normal_direction_flag==1) 
    {
        normal_vec[0] = -(yc2 -yc1)/area;
        normal_vec[1] = (xc2-xc1)/area;
    }

    /*double norm = sqrt(normal_vec[0]*normal_vec[0] + normal_vec[1]*normal_vec[1]);*/
    /*if (fabs(norm - 1.0) > 1e-12) {*/
	/*printf("Norm of normal vector not one!");*/
	/*getchar();*/
    /*}*/
}

/* void get_end_points
 * get 1D-interval on which integration is being performed
 * return x coordinates unless vertical line
 * in: index to triangle -- E
 * in: current edge -- iedge        
 * in: node data structure -- mesh_node    
 * in: element data structure -- mesh_element 
 * in: edge data structure -- mesh_edge    
 * out: end points -- end_points
 */
void get_end_points(int E,
		    int iedge, 
		    node* mesh_node, 
		    element* mesh_element,
		    edge* mesh_edge,
		    double* end_points)
{
    int node_a,node_b;
    double a1,a2,b1,b2;

    // Get nodes of iedge
    node_a = mesh_edge[iedge].vertex[1];
    node_b = mesh_edge[iedge].vertex[2];

    // Get ordering relative to current element
    if((mesh_element[E].vertex[1]==node_a ||mesh_element[E].vertex[1]==node_b) &&
		(mesh_element[E].vertex[2]==node_a ||mesh_element[E].vertex[2]==node_b)) 
    {
        node_a = mesh_element[E].vertex[1];
        node_b = mesh_element[E].vertex[2];
    }
    else if((mesh_element[E].vertex[2]==node_a || mesh_element[E].vertex[2] ==node_b)&&
		(mesh_element[E].vertex[3] ==node_a || mesh_element[E].vertex[3] == node_b)) 
    {
        node_a = mesh_element[E].vertex[2];
        node_b = mesh_element[E].vertex[3];
    }
    else if((mesh_element[E].vertex[3] ==node_a || mesh_element[E].vertex[3] ==node_b)&&
		(mesh_element[E].vertex[1] ==node_a || mesh_element[E].vertex[1] == node_b)) 
    {
        node_a = mesh_element[E].vertex[3];
        node_b = mesh_element[E].vertex[1];
    }

    // Get coordinates of each node
    a1 = mesh_node[node_a].coord[0];
    b1 = mesh_node[node_a].coord[1];

    a2 = mesh_node[node_b].coord[0];
    b2 = mesh_node[node_b].coord[1];

    // Case1:: vertical line
    if (fabs(a1-a2) < EPSILON) 
    {
        end_points[0] = b1;
        end_points[1] = b2;
    }
    else 
    {
        end_points[0] = a1;
        end_points[1] = a2;
    }
}//get_end_points

/* void get_phys_coords_edge
 * transform reference point on edge to physical coordinates
 * in: current edge -- iedge        
 * in: node data structure -- mesh_node    
 * in: element data structure -- mesh_element 
 * in: edge data structure -- mesh_edge    
 * in: reference point on the edge -- lgpt    
 * out: physical coordintes -- phys_coords
 */
void get_phys_coords_edge(int iedge,
			  node* mesh_node,
                          element* mesh_element,
                          edge* mesh_edge, 
                          double lgpt,
                          double* phys_coords)
{
    int node_a, node_b;
    double a1,a2,b1,b2;

    node_a = mesh_edge[iedge].vertex[1];
    node_b = mesh_edge[iedge].vertex[2];

    a1 = mesh_node[node_a].coord[0];
    b1 = mesh_node[node_a].coord[1];

    a2 = mesh_node[node_b].coord[0];
    b2 = mesh_node[node_b].coord[1];

    //case 1 horizontal line
    if(fabs(b1- b2) <EPSILON) 
    {
        phys_coords[0] = lgpt;
        phys_coords[1] = b2;
    }
    //case 2 vertical line
    else if(fabs(a1- a2) < EPSILON) 
    {
        phys_coords[0] = a1;
        phys_coords[1] = lgpt;
    }
    else 
    {
        phys_coords[0] = lgpt;
        phys_coords[1] = b1 + ((b2-b1)/(a2-a1))*(lgpt - a1);
    }
}

/* double calc_length_face
 * calculates length of an edge
 * in: edge data structure -- mesh_edge
 * in: current edge -- iedge
 * in: node data structure -- mesh_node
 * out: length of edge
 */
double calc_length_face(edge* mesh_edge, 
	                int iedge, 
		        node* mesh_node)
{
    double length;
    int node_a, node_b;
    double xc1,yc1,xc2,yc2;
    node_a = mesh_edge[iedge].vertex[1];
    node_b = mesh_edge[iedge].vertex[2];

    xc1 = mesh_node[node_a].coord[0];
    yc1 = mesh_node[node_a].coord[1];
    xc2 = mesh_node[node_b].coord[0];
    yc2 = mesh_node[node_b].coord[1];
    length = (xc1 - xc2)*(xc1-xc2)+(yc1-yc2)*(yc1-yc2);
    length = sqrt(length);
    return length;
}

/* void copy_matrix
 * copies matrix
 * in: matrix of size Nloc*Nloc	 -- temp_Aloc
 * in: size of square matrix     -- Nloc
 * out: copy of matrix -- Aloc_matrix
 */
void copy_matrix(double* temp_Aloc,
		 int Nloc,
                 double* Aloc_matrix) 
{
    int idofs=0;
    int jdofs=0;

    for(idofs=1;idofs< Nloc+1;idofs++) 
    {
        for(jdofs=1;jdofs < Nloc+1;jdofs++) 
        {
            Aloc_matrix(idofs,jdofs,Nloc) = temp_Aloc(idofs,jdofs,Nloc);
        }
    }
}

/* void shu_osher_limit_scalar
 * performs flux limiting on the solution vector 
 * in: mesh_element data structure   -- mesh_element
 * in: mesh_node data structure      -- mesh_node
 * in: mesh_edge data structure      -- mesh_edge
 * in: number of triangles           -- Nelts
 * in: number of degrees of freedom  -- Nloc
 * in: limiting parameter nu	     -- nu
 * in: estimate of limiting paramter -- mdx2
 * in: solution at time tn           -- uh_tn_sol_vec
 * out: flux limited solution        -- uh_limit
 */
void shu_osher_limit_scalar(element* mesh_element,
                            node* mesh_node,
                            edge* mesh_edge,
                            int Nelts,
                            int Nloc,
                            double nu, 
                            double mdx2,
                            double* uh_tn_sol_vec,
                            double* uh_limit)
{
    int E;				// Index of triangle E
    int edge1,edge2,edge3;		// Index of edges of E
    int E1,E2,E3;	                // Neighbours of E
    int n1,n2,n3;			// Nodes of E
    int idofs;				// Index of parameter for basis fcts
    double baryctr0[2]; 		// Barycenters 
    double baryctr1[2];
    double baryctr2[2];
    double baryctr3[2];
    double uh_E,uh_E1,uh_E2,uh_E3; 	// Values at barycenters
    double m1[2]; 			// Midpoints of edges
    double m2[2];
    double m3[2];
    double uh_m1,uh_m2,uh_m3; 		// Values at edge-midpoints
    double alpha1[2]; 			// Needed to compare slopes
    double alpha2[2]; 			
    double alpha3[2]; 		
    double D_uh_m1,D_uh_m2,D_uh_m3; 	// 'Slope' between neighbouring triangles
    double Delta1,Delta2,Delta3; 	// Output of the minmod for all edges
    double Deltah1,Deltah2,Deltah3; 	// Modified output of the minmod for all edges
    double pos,neg,sigma_p,sigma_m; 	// Temporary variables to process Deltas
    double uh_limit_local[3];		// Limited solution local
    double uh_lin[3];		// Transformed linear approximation 
    double temp;
    int verbose = 0;

    // Loop over all triangles
    for(E=1;E<Nelts+1;E++) {
        // Initializations
	Deltah1 = 6.66; Deltah2 = 6.66; Deltah3 = 6.66;
	uh_E1 = 6.66; uh_E2 = 6.66; uh_E3 = 6.66;
	E1 = E; E2 = E; E3 = E;

        // Edges of E
	edge1 = mesh_element[E].edge[1];
	edge2 = mesh_element[E].edge[2];
	edge3 = mesh_element[E].edge[3];

	// Coordinates of midpoints of edges
	n1 = mesh_element[E].vertex[1];
	n2 = mesh_element[E].vertex[2];
	n3 = mesh_element[E].vertex[3];

	m1[0] = (mesh_node[n1].coord[0]+mesh_node[n2].coord[0])/2.0;
	m1[1] = (mesh_node[n1].coord[1]+mesh_node[n2].coord[1])/2.0;
	m2[0] = (mesh_node[n2].coord[0]+mesh_node[n3].coord[0])/2.0;
	m2[1] = (mesh_node[n2].coord[1]+mesh_node[n3].coord[1])/2.0;
	m3[0] = (mesh_node[n3].coord[0]+mesh_node[n1].coord[0])/2.0;
	m3[1] = (mesh_node[n3].coord[1]+mesh_node[n1].coord[1])/2.0;

	// Value of uh in E at midpoints of edges
	get_approx_solution_one_mom(uh_tn_sol_vec,E,mesh_element,mesh_node,0.5,0.0,&uh_m1);
	get_approx_solution_one_mom(uh_tn_sol_vec,E,mesh_element,mesh_node,0.5,0.5,&uh_m2);
	get_approx_solution_one_mom(uh_tn_sol_vec,E,mesh_element,mesh_node,0.0,0.5,&uh_m3);

	// Get cell average
	uh_E = (uh_m1 + uh_m2 + uh_m3)/3.0;

	// Transform to linear part
	uh_lin[0] = uh_m1 - uh_m2 + uh_m3;
	uh_lin[1] = 2.0*(uh_m2 - uh_m3);
	uh_lin[2] = 2.0*(uh_m2 - uh_m1);

	// Compute barycenter
	find_barycenter(E,mesh_element,mesh_node,baryctr0);

	// Repeat for all neighbours of E
        if(mesh_edge[edge1].edge_type == INTERIOR){
	    if (mesh_edge[edge1].neighbour[1] == E) {
		E1 = mesh_edge[edge1].neighbour[2];
	    }
	    else {
		E1 = mesh_edge[edge1].neighbour[1];
	    }
	    find_barycenter(E1,mesh_element,mesh_node,baryctr1);
	    uh_E1 = cell_average_one_mom(uh_tn_sol_vec,E1,mesh_element,mesh_node);
	}
	// If E is on the boundary, mirror the triangle along its boundary edge to find barycenter and use 
	// the value of the midpoint as the cell-average on the ghost triangle
	else if (mesh_edge[edge1].edge_type == EXTERIOR) {
	    E1 = E;
	    baryctr1[0] = baryctr0[0] + 2*(m1[0] - baryctr0[0]); 
	    baryctr1[1] = baryctr0[1] + 2*(m1[1] - baryctr0[1]); 
	    uh_E1 = uh_m1;
	}

        if(mesh_edge[edge2].edge_type == INTERIOR){
	    if (mesh_edge[edge2].neighbour[1] == E) {
		E2 = mesh_edge[edge2].neighbour[2];
	    }
	    else {
		E2 = mesh_edge[edge2].neighbour[1];
	    }
	    find_barycenter(E2,mesh_element,mesh_node,baryctr2);
	    uh_E2 = cell_average_one_mom(uh_tn_sol_vec,E2,mesh_element,mesh_node);
	}
	else if (mesh_edge[edge2].edge_type == EXTERIOR){
	    E2 = E;
	    baryctr2[0] = baryctr0[0] + 2*(m2[0] - baryctr0[0]); 
	    baryctr2[1] = baryctr0[1] + 2*(m2[1] - baryctr0[1]); 
	    uh_E2 = uh_m2;
	}

        if(mesh_edge[edge3].edge_type == INTERIOR){
	    if (mesh_edge[edge3].neighbour[1] == E) {
		E3 = mesh_edge[edge3].neighbour[2];
	    }
	    else {
		E3 = mesh_edge[edge3].neighbour[1];
	    }
	    find_barycenter(E3,mesh_element,mesh_node,baryctr3);
	    uh_E3 = cell_average_one_mom(uh_tn_sol_vec,E3,mesh_element,mesh_node);
	}
	else if (mesh_edge[edge3].edge_type == EXTERIOR){
	    E3 = E;
	    baryctr3[0] = baryctr0[0] + 2*(m3[0] - baryctr0[0]); 
	    baryctr3[1] = baryctr0[1] + 2*(m3[1] - baryctr0[1]); 
	    uh_E3 = uh_m3;
	}

	// Compute alphas and Delta_uh_bar(mi,E)
	compute_alphas(m1,baryctr0,baryctr1,baryctr2,alpha1);
	if ((alpha1[0] < -1e-12) || (alpha1[1] < -1e-12)) {
	    compute_alphas(m1,baryctr0,baryctr2,baryctr3,alpha1);
	    if ((alpha1[0] < -1e-12) || (alpha1[1] < -1e-12)) {
		compute_alphas(m1,baryctr0,baryctr1,baryctr3,alpha1);
		D_uh_m1 = alpha1[0]*(uh_E1-uh_E) + alpha1[1]*(uh_E3-uh_E);
		if ((alpha1[0] < -1e-12) || (alpha1[1] < -1e-12)) {
		    printf("couldn't find non-negative alphas!\n");
		}
	    }
	    else {
		D_uh_m1 = alpha1[0]*(uh_E2-uh_E) + alpha1[1]*(uh_E3-uh_E);
	    }
	}
	else {
	    D_uh_m1 = alpha1[0]*(uh_E1-uh_E) + alpha1[1]*(uh_E2-uh_E);
	}

	compute_alphas(m2,baryctr0,baryctr1,baryctr2,alpha2);
	if ((alpha2[0] < -1e-12) || (alpha2[1] < -1e-12)) {
	    compute_alphas(m2,baryctr0,baryctr2,baryctr3,alpha2);
	    if ((alpha2[0] < -1e-12) || (alpha2[1] < -1e-12)) {
		compute_alphas(m2,baryctr0,baryctr1,baryctr3,alpha2);
		D_uh_m2 = alpha2[0]*(uh_E1-uh_E) + alpha2[1]*(uh_E3-uh_E);
		if ((alpha2[0] < -1e-12) || (alpha2[1] < -1e-12)) {
		    printf("couldn't find non-negative alphas!\n");
		}
	    }
	    else {
		D_uh_m2 = alpha2[0]*(uh_E2-uh_E) + alpha2[1]*(uh_E3-uh_E);
	    }
	}
	else {
	    D_uh_m2 = alpha2[0]*(uh_E1-uh_E) + alpha2[1]*(uh_E2-uh_E);
	}

	compute_alphas(m3,baryctr0,baryctr1,baryctr2,alpha3);
	if ((alpha3[0] < -1e-12) || (alpha3[1] < -1e-12)) {
	    compute_alphas(m3,baryctr0,baryctr2,baryctr3,alpha3);
	    if ((alpha3[0] < -1e-12) || (alpha3[1] < -1e-12)) {
		compute_alphas(m3,baryctr0,baryctr1,baryctr3,alpha3);
		D_uh_m3 = alpha3[0]*(uh_E1-uh_E) + alpha3[1]*(uh_E3-uh_E);
		if ((alpha3[0] < -1e-12) || (alpha3[1] < -1e-12)) {
		    printf("couldn't find non-negative alphas!\n");
		}
	    }
	    else {
		D_uh_m3 = alpha3[0]*(uh_E2-uh_E) + alpha3[1]*(uh_E3-uh_E);
	    }
	}
	else {
	    D_uh_m3 = alpha3[0]*(uh_E1-uh_E) + alpha3[1]*(uh_E2-uh_E);
	}

	// Compute Deltas
	Delta1 = minmod(uh_m1-uh_E,nu*D_uh_m1,mdx2);
	Delta2 = minmod(uh_m2-uh_E,nu*D_uh_m2,mdx2);
	Delta3 = minmod(uh_m3-uh_E,nu*D_uh_m3,mdx2);

	// Limit only linear part, i.e. up to third coefficient...
	/*uh_limit_local[0] = uh_tn_sol_vec[(E-1)*Nloc]+Delta1-Delta2+Delta3;*/
	if (fabs(Delta1+Delta2+Delta3) > 1e-12) {
	    pos =  (0.0<Delta1)?Delta1:0.0;
	    pos += (0.0<Delta2)?Delta2:0.0;
	    pos += (0.0<Delta3)?Delta3:0.0;
	    neg =  (0.0<-Delta1)?-Delta1:0.0;
	    neg += (0.0<-Delta2)?-Delta2:0.0;
	    neg += (0.0<-Delta3)?-Delta3:0.0;
	    sigma_p = (1.0<(neg/pos))?1.0:(neg/pos);
	    sigma_m = (1.0<(pos/neg))?1.0:(pos/neg);
	    Deltah1 = sigma_p*((0.0<Delta1)?Delta1:0.0);
	    Deltah1 -= sigma_m*((0.0<-Delta1)?-Delta1:0.0);
	    Deltah2 = sigma_p*((0.0<Delta2)?Delta2:0.0);
	    Deltah2 -= sigma_m*((0.0<-Delta2)?-Delta2:0.0);
	    Deltah3 = sigma_p*((0.0<Delta3)?Delta3:0.0);
	    Deltah3 -= sigma_m*((0.0<-Delta3)?-Delta3:0.0);
	    uh_limit_local[0] = uh_E+Deltah1-Deltah2+Deltah3;
	    uh_limit_local[1] = 2.0*(Deltah2-Deltah3);
	    uh_limit_local[2] = 2.0*(Deltah2-Deltah1);
	}
	else {
	    uh_limit_local[0] = uh_E+Delta1-Delta2+Delta3;
	    uh_limit_local[1] = 2.0*(Delta2-Delta3);
	    uh_limit_local[2] = 2.0*(Delta2-Delta1);
	}

	// Check whether limiting changes the solution
	if ( (fabs(uh_limit_local[0] - uh_lin[0]) > 1e-12)
	  || (fabs(uh_limit_local[1] - uh_lin[1]) > 1e-12)
	  || (fabs(uh_limit_local[2] - uh_lin[2]) > 1e-12)) {
 	    // ... if it does, copy the linear part...
	    uh_limit[(E-1)*Nloc  ] = uh_limit_local[0];
	    uh_limit[(E-1)*Nloc+1] = uh_limit_local[1];
	    uh_limit[(E-1)*Nloc+2] = uh_limit_local[2];
	    // and cut the non-linear part
	    for(idofs=3;idofs<Nloc;idofs++) {
		uh_limit[(E-1)*Nloc + idofs] = 0.0;
	    }
	    if (verbose == 1) {
		printf("m1 %lf,%lf :: uh_m1 %lf\n",m1[0],m1[1],uh_m1);
		printf("m2 %lf,%lf :: uh_m2 %lf\n",m2[0],m2[1],uh_m2);
		printf("m3 %lf,%lf :: uh_m3 %lf\n",m3[0],m3[1],uh_m3);
		printf("u_transf: alpha1 %lf :: alpha2 %lf :: alpha3 %lf\n",uh_lin[0],uh_lin[1],uh_lin[2]);
		printf("E %d :: uh_E %lf :: bary0 %lf,%lf\nE1 %d :: uh_E1 %lf :: bary1 %lf,%lf\nE2 %d :: uh_E2 %lf :: bary2 %lf,%lf\nE3 %d :: uh_E3 %lf :: bary3 %lf,%lf\n",E,uh_E,baryctr0[0],baryctr0[1],E1,uh_E1,baryctr1[0],baryctr1[1],E2,uh_E2,baryctr2[0],baryctr2[1],E3,uh_E3,baryctr3[0],baryctr3[1]);
		printf("alpha1 %lf,%lf :: alpha2 %lf,%lf :: alpha3 %lf,%lf\n",alpha1[0],alpha1[1],alpha2[0],alpha2[1],alpha3[0],alpha3[1]);
		printf("D_uh_m1 %lf :: D_uh_m2 %lf :: D_uh_m3 %lf\n",D_uh_m1,D_uh_m2,D_uh_m3);
		printf("uh_m1-uh_e %lf :: nu*D_uh_m1 %lf\n",uh_m1-uh_E,nu*D_uh_m1);
		printf("uh_m2-uh_e %lf :: nu*D_uh_m2 %lf\n",uh_m2-uh_E,nu*D_uh_m2);
		printf("uh_m3-uh_e %lf :: nu*D_uh_m3 %lf\n",uh_m3-uh_E,nu*D_uh_m3);
		printf("pos %lf :: neg %lf :: sigma_p %lf :: sigma_m %lf\n",pos,neg,sigma_p,sigma_m);
		printf("        Delta1 %lf :: Delta2 %lf :: Delta3 :: %lf\n",Delta1,Delta2,Delta3);
		printf("Hatted: Delta1 %lf :: Delta2 %lf :: Delta3 :: %lf\n",Deltah1,Deltah2,Deltah3);
		temp = cell_average_one_mom(uh_limit,E,mesh_element,mesh_node);
		if (mesh_element[E].degree == 2) {
		    printf("Before limiting: alpha_1 %lf :: alpha_2 %lf :: alpha_3 %lf :: alpha_4 %lf :: alpha_5 %lf :: alpha_6 %lf\n",uh_tn_sol_vec[(E-1)*Nloc],uh_tn_sol_vec[(E-1)*Nloc+1],uh_tn_sol_vec[(E-1)*Nloc+2],uh_tn_sol_vec[(E-1)*Nloc+3],uh_tn_sol_vec[(E-1)*Nloc+4],uh_tn_sol_vec[(E-1)*Nloc+5]);
		    printf("After  limiting: alpha_1 %lf :: alpha_2 %lf :: alpha_3 %lf :: alpha_4 %lf :: alpha_5 %lf :: alpha_6 %lf\n",uh_limit[(E-1)*Nloc],uh_limit[(E-1)*Nloc+1],uh_limit[(E-1)*Nloc+2],uh_limit[(E-1)*Nloc+3],uh_limit[(E-1)*Nloc+4],uh_limit[(E-1)*Nloc+5]);
		}
		else if (mesh_element[E].degree == 1) {
		    printf("Before limiting: alpha_1 %lf :: alpha_2 %lf :: alpha_3 %lf\n",uh_tn_sol_vec[(E-1)*Nloc],uh_tn_sol_vec[(E-1)*Nloc+1],uh_tn_sol_vec[(E-1)*Nloc+2]);
		    printf("After  limiting: alpha_1 %lf :: alpha_2 %lf :: alpha_3 %lf\n",uh_limit_local[0],uh_limit_local[1],uh_limit_local[2]);
		}
		printf("New value at barycenter :: %lf\n",temp);
		printf("==============================\n");
		getchar();
	    }
	}
	// otherwise perform no limiting
	else {
	    for(idofs=0;idofs<Nloc;idofs++) {
		uh_limit[(E-1)*Nloc+idofs] = uh_tn_sol_vec[(E-1)*Nloc+idofs];
	    }
	}
    }
}

/*void find_barycenter
 * computes the coordinates of the barycenter of a triangle
 * in: index to triangle -- E
 * in: mesh_element data structure -- mesh_element
 * in: mesh_node data structure    -- mesh_node
 * out: coordinates of barycenter  -- bary_coords
 */
void find_barycenter(int E,
		element* mesh_element,
		node* mesh_node,
		double* bary_coords)
{
    int n1,n2,n3;
    double x1,y1,x2,y2,x3,y3;

    n1 = mesh_element[E].vertex[1];
    n2 = mesh_element[E].vertex[2];
    n3 = mesh_element[E].vertex[3];

    x1 = mesh_node[n1].coord[0];
    y1 = mesh_node[n1].coord[1];
    x2 = mesh_node[n2].coord[0];
    y2 = mesh_node[n2].coord[1];
    x3 = mesh_node[n3].coord[0];
    y3 = mesh_node[n3].coord[1];

    bary_coords[0] = (x1+x2+x3)/3.0;
    bary_coords[1] = (y1+y2+y3)/3.0;

    /*printf("c1 %lf,%lf\n",x1,y1);*/
    /*printf("c2 %lf,%lf\n",x2,y2);*/
    /*printf("c3 %lf,%lf\n",x3,y3);*/
    /*printf("bary %lf,%lf\n",bary_coords[0],bary_coords[1]);*/
    /*getchar();*/
}

/*void compute_alphas
 * computes alpha1, alpha2 s.t.
 * mi-bo = alpha1*(b1-b0) + alpha2*(b2-b0)
 * in: midpoint of edge        -- mi
 * in: barycentric coordinates -- b0,b1,b2
 * out: parameters             -- alphai
 */
void compute_alphas(double* mi,
		    double* b0,
		    double* b1,
		    double* b2,
		    double* alphai)
{
    alphai[0] = 0.0; alphai[1] = 0.0;

    double a1n = -(b0[0]*b2[1] - b0[1]*b2[0] - b0[0]*mi[1] + b0[1]*mi[0] + b2[0]*mi[1] - b2[1]*mi[0]);
    double a2n =   b0[0]*b1[1] - b0[1]*b1[0] - b0[0]*mi[1] + b0[1]*mi[0] + b1[0]*mi[1] - b1[1]*mi[0];

    double ad = b0[0]*b1[1] - b0[1]*b1[0] - b0[0]*b2[1] + b0[1]*b2[0] + b1[0]*b2[1] - b1[1]*b2[0];

    alphai[0] = a1n/ad;
    alphai[1] = a2n/ad;

    /*double diffx = alphai[0]*(b1[0] - b0[0]) + alphai[1]*(b2[0] - b0[0]) - mi[0] + b0[0];*/
    /*double diffy = alphai[0]*(b1[1] - b0[1]) + alphai[1]*(b2[1] - b0[1]) - mi[1] + b0[1];*/

    /*if ( ( fabs(diffx) > 1e-12 ) || ( fabs(diffy) > 1e-12 ) ){*/
	/*double diffx_ch = alphai[1]*(b1[0] - b0[0]) + alphai[0]*(b2[0] - b0[0]) - mi[0] + b0[0];*/
	/*double diffy_ch = alphai[1]*(b1[1] - b0[1]) + alphai[0]*(b2[1] - b0[1]) - mi[1] + b0[1];*/

	/*printf("diff = %.12lf ,%.12lf\n",diffx,diffy);*/
	/*printf("diff_ch = %.12lf ,%.12lf\n",diffx_ch,diffy_ch);*/
	/*getchar();*/
    /*}*/

    /*double diff_x = mi[0] - b0[0] - alphai[0]*( b1[0] - b0[0] ) - alphai[1]*(b2[0] - b0[0]);*/
    /*double diff_y = mi[1] - b0[1] - alphai[0]*( b1[1] - b0[1] ) - alphai[1]*(b2[1] - b0[1]);*/
    
    /*if ( ( fabs(diff_x) > 1e-12 ) || ( fabs(diff_y) > 1e-12 ) ){*/
	/*printf("diff_x = %lf \n diff_y = %lf\n !!!\n",diff_x,diff_y);*/
	/*getchar();*/
    /*}*/

    /*if (isnan(alphai[0]) || isnan(alphai[1])){*/
	/*printf("Alpha %lf :: %lf!\t",alphai[0],alphai[1]);*/
	/*printf("mi %lf,%lf :: bary0 %lf,%lf :: bary1 %lf,%lf :: bary2 %lf,%lf \n",mi[0],mi[1],b0[0],b0[1],b1[0],b1[1],b2[0],b2[1]);*/
	/*printf("t1 %lf :: t2 %lf",t1,t2);*/
	/*getchar();*/
    /*}*/

    /*double t1, t2;*/
    /*if ((b2[1] != b0[1])) {*/
	/*t1 = ((b1[0]-b0[0]) - (b2[0]-b0[0])/(b2[1]-b0[1])*(b1[1]-b0[1]));*/
	/*t2 = ((mi[0]-b0[0]) - (b2[0]-b0[0])/(b2[1]-b0[1])*(mi[1]-b0[1]));*/
	/*alphai[0] = t2/t1;*/
	/*alphai[1] = ((mi[1]-b0[1]) - alphai[0]*(b1[1]-b0[1]))/(b2[1]-b0[1]);*/
    /*}*/
    /*else if ((b2[0] != b0[0])) {*/
	/*t1 = ((b1[1]-b0[1]) - (b2[1]-b0[1])/(b2[0]-b0[0])*(b1[0]-b0[0]));*/
	/*t2 = ((mi[1]-b0[1]) - (b2[1]-b0[1])/(b2[0]-b0[0])*(mi[0]-b0[0]));*/
	/*alphai[0] = t2/t1;*/
	/*alphai[1] = ((mi[0]-b0[0]) - alphai[0]*(b1[0]-b0[0]))/(b2[0]-b0[0]);*/
    /*}*/
    /*else {*/
	/*printf("Cannot determine coefficients alpha for limiting, since two barycenters coincide!");*/
	/*alphai[0] = -1.0;*/
    /*}*/
}

/* double minmod
 * modified minmod function
 * mdx2=M*(dx^2); where M is an upper bound of the absolute of value 
 * of the second-order derivative
 * in: two real values   -- x,y
 * in: parameters	 -- mdx2
 * out: minmod value     -- minmod
 */
double minmod(double x, double y, double mdx2)
{
 
 
    
    /*return(y);*/
    if (fabs(x) < mdx2){
	return(x);
    }
    else if ((x>0.0) && (y>0.0)){
        if ( x < y ) {
	    return(x);
	    }
	else {
	    return(y);
	}
	/*return(x<y)?x:y;*/
    }
    else if ((x<0.0) && (y<0.0)){
        if ( x < y ) {
	    return(y);
	    }
	else {
	    return(x);
	}
	/*return(x<y)?y:x;*/
    }
    else{
        return(0.0);
    }
}

/* double cell_average_one_mom
 * returns cell average
 * in: index to triangle -- E
 * in: mesh_element data structure -- mesh_element
 * in: mesh_node data structure    -- mesh_node
 * in: solution vector -- uh_sol
 * out: cell average -- 
 */
double cell_average_one_mom(double* uh_sol,
		    int E,
		    element* mesh_element,
		    node* mesh_node)
{
     double uh_m1,uh_m2,uh_m3;
    // Value of uh in E at midpoints of edges
    get_approx_solution_one_mom(uh_sol,E,mesh_element,mesh_node,0.5,0.0,&uh_m1);
    get_approx_solution_one_mom(uh_sol,E,mesh_element,mesh_node,0.5,0.5,&uh_m2);
    get_approx_solution_one_mom(uh_sol,E,mesh_element,mesh_node,0.0,0.5,&uh_m3);

    return((uh_m1 + uh_m2 + uh_m3)/3.0);
}

/* void cell_average
 * returns cell average for each system component
 * in: solution vector -- uh_sol
 * in: index to triangle -- E
 * in: mesh_element data structure -- mesh_element
 * in: mesh_node data structure    -- mesh_node
 * out: cell average -- uh_E
 */
void cell_average(double* uh_sol,
		  int E,
		  element* mesh_element,
		  node* mesh_node,
		  double* uh_E)
{
    double uh_m1[SYSTEM_DIM]; 		// Values at edge-midpoints
    double uh_m2[SYSTEM_DIM];
    double uh_m3[SYSTEM_DIM];
    int s_dim;

    // Value of uh in E at midpoints of edges
    get_approx_solution(uh_sol,E,mesh_element,mesh_node,0.5,0.0,uh_m1);
    get_approx_solution(uh_sol,E,mesh_element,mesh_node,0.5,0.5,uh_m2);
    get_approx_solution(uh_sol,E,mesh_element,mesh_node,0.0,0.5,uh_m3);

    for (s_dim=0;s_dim<SYSTEM_DIM;s_dim++) {
	    uh_E[s_dim] = (uh_m1[s_dim] + uh_m2[s_dim] + uh_m3[s_dim])/3.0;
	/*uh_E[s_dim] = (uh_m1[s_dim] + uh_m2[s_dim] + uh_m3[s_dim])/3.0*mesh_element[E].det;*/
    }
    /*printf("det = %lf \n",mesh_element[E].det);*/
    /*printf("uh_E = [%lf,\t%lf,\t%lf] \n", uh_E[0], uh_E[1], uh_E[2]);*/
    /*getchar();*/
    return;
}

/* void cell_average_lin
 * returns cell average for each system component of a linear function
 * in: solution vector -- uh_sol
 * in: index to triangle -- E
 * in: mesh_element data structure -- mesh_element
 * in: mesh_node data structure    -- mesh_node
 * out: cell average -- uh_E
 */
void cell_average_lin(double* uh_sol,
		      int E,
		      element* mesh_element,
		      node* mesh_node,
		      double* uh_E)
{
    double uh_m1[SYSTEM_DIM]; 		// Values at edge-midpoints
    double uh_m2[SYSTEM_DIM];
    double uh_m3[SYSTEM_DIM];
    int s_dim;

    /*// Value of uh in E at midpoints of edges*/
    /*get_approx_solution_lin(uh_sol,E,mesh_element,mesh_node,0.5,0.0,uh_m1);*/
    /*get_approx_solution_lin(uh_sol,E,mesh_element,mesh_node,0.5,0.5,uh_m2);*/
    /*get_approx_solution_lin(uh_sol,E,mesh_element,mesh_node,0.0,0.5,uh_m3);*/

    /*for (s_dim=0;s_dim<SYSTEM_DIM;s_dim++) {*/
	/*uh_E[s_dim] = (uh_m1[s_dim] + uh_m2[s_dim] + uh_m3[s_dim])/3.0;*/
    /*}*/

    /*double bla[3];*/
    get_approx_solution_lin(uh_sol,E,mesh_element,mesh_node,1.0/3.0,1.0/3.0,uh_E);

    /*printf("uh_E - bla = [%lf,\t%lf,\t%lf] \n", uh_E[0] - bla[0], uh_E[1] - bla[1], uh_E[2] - bla[2]);*/
    /*getchar();*/
    return;
}

/* void qrdcmp
 * constructs the QR decomposition of a
 * Q is represented as the product of Q1*...*Q_n-1 where 
 * Qj = 1 - u_j * u_j / c_j
 * in: matrix -- a
 * in: dimension of a -- n
 * out: coefficients in Q  -- c
 * out: diagonal elements -- d
 * out: flag whether singularity was encountered -- sing
 * out: upper triangular matrix R -- a
 */
void qrdcmp(double **a, int n, double *c, double *d, int *sing)
{
    int i,j,k; 
    double scale,sigma,sum,tau;
    *sing=0; 
    for (k=1;k<n;k++) {
	scale=0.0; 
	for (i=k;i<=n;i++) {
	    scale=(scale<fabs(a[i][k]))?fabs(a[i][k]):scale; 
	}
	if (scale == 0.0) {	// Singular case.
	    *sing=1;
	    c[k]=d[k]=0.0; 
	} 
	else {			// Form Qk and Qk · A.
	    for (i=k;i<=n;i++) {
		a[i][k] /= scale; 
	    }
	    sum = 0.0;
	    for (i=k;i<=n;i++) { 
		sum += a[i][k]*a[i][k]; 
	    }
	    sigma=((a[k][k])>=0.0?fabs(sqrt(sum)):-fabs(sqrt(sum)));
	    a[k][k] += sigma; 
	    c[k]=sigma*a[k][k]; 
	    d[k] = -scale*sigma; 
	    for (j=k+1;j<=n;j++) {
		for (sum=0.0,i=k;i<=n;i++) sum += a[i][k]*a[i][j]; 
		tau=sum/c[k]; 
		for (i=k;i<=n;i++) a[i][j] -= tau*a[i][k];
	    }
	}
    } 
    d[n]=a[n][n]; 
    if (fabs(d[n]) < 1e-12) *sing=1;
}

/* void qrsolv
 * solves the set of n linear equations Ax = b
 * uses input from qrdcmp
 * in: matrix -- a
 * in: dimension of a -- n
 * in: coefficients in Q  -- c
 * in: diagonal elements -- d
 * in: rhs vector -- b
 * out: solution vector -- b
 */
void qrsolv(double **a, int n, double c[], double d[], double b[])
{
    /*void rsolv(double **a, int n, double d[], double b[]);*/
    int i,j; 
    double sum,tau;
    for (j=1;j<n;j++) {		// Form QT · b. 
	for (sum=0.0,i=j;i<=n;i++) sum += a[i][j]*b[i-1];
	tau=sum/c[j]; 
	for (i=j;i<=n;i++) b[i-1] -= tau*a[i][j];
    } 
    rsolv(a,n,d,b);		// Solve R · x = QT · b.
}

/* void rsolv
 * solves the set of n linear equations Rx = b
 * in: matrix -- a
 * in: dimension of a -- n
 * in: diagonal elements -- d
 * in: rhs vector -- b
 * out: solution vector -- b
 */
void rsolv(double **a, int n, double d[], double b[])
{
    int i,j; 
    double sum;
    b[n-1] /= d[n]; 
    for (i=n-1;i>=1;i--) {
	for (sum=0.0,j=i+1;j<=n;j++) sum += a[i][j]*b[j-1]; 
	b[i-1]=(b[i-1]-sum)/d[i];
    }
}

/* void assemble_solve_qr
 * driver routine solving the system Ax=b using 
 * QR-decomposition
 * in: matrix A -- Aloc_matrix
 * in: rhs b -- rhs
 * in: dimension -- Nloc
 * out: coefficients in Q -- Q_coeff
 * out: coefficients in R -- R_coeff
 * out: rhs from QR -- rhs_qr
 * out: R -- A_loc_QR
 * out: flag whether singularity was encountered -- sing
 */
void assemble_solve_qr(double *Aloc_matrix, double *rhs, int Nloc, double *Q_coeff, double *R_coeff, double *rhs_qr, double **A_loc_QR, int *sing)
{
    int i,j;
    init_zero_d(Q_coeff,(Nloc+1));
    init_zero_d(R_coeff,(Nloc+1));
    init_zero_d(rhs_qr,(Nloc+1));

    for(i=1;i<Nloc+1;i++) {
	for(j=1;j<Nloc+1;j++) {
	    A_loc_QR[i][j] = Aloc_matrix(i,j,Nloc);
	}
	rhs_qr[i] = rhs[i-1];
    }

    qrdcmp(A_loc_QR,Nloc,Q_coeff,R_coeff,sing);
    qrsolv(A_loc_QR,Nloc,Q_coeff,R_coeff,rhs_qr);

    for(i=1;i<Nloc+1;i++) {
	rhs[i-1] = rhs_qr[i];
    }
}

/* void assemble_qr
 * Computes QR-decomposition
 * in: matrix A -- Aloc_matrix
 * in: dimension -- Nloc
 * out: coefficients in Q -- Q_coeff
 * out: coefficients in R -- R_coeff
 * out: R -- A_loc_QR
 * out: flag whether singularity was encountered -- sing
 */
void assemble_qr(double *Aloc_matrix, int Nloc, double *Q_coeff, double *R_coeff, double **A_loc_QR, int *sing)
{
    int i,j;
    init_zero_d(Q_coeff,(Nloc+1));
    init_zero_d(R_coeff,(Nloc+1));

    for(i=1;i<Nloc+1;i++) {
	for(j=1;j<Nloc+1;j++) {
	    A_loc_QR[i][j] = Aloc_matrix(i,j,Nloc);		
	}
    }
    qrdcmp(A_loc_QR,Nloc,Q_coeff,R_coeff,sing);
}

/* void upwind_flux
 * evaluates Upwind flux
 * in: index to element and it's neighbour -- E, E_neighbour
 * in: solution vector			   -- uh_tn_sol_vec
 * in: mesh_element data structure   -- mesh_element
 * in: mesh_node data structure      -- mesh_node
 * in: parameters 		     -- user_param
 * in: normal vector between the two cells -- normal_e
 * in: s_dim_counter          -- dimension of system counter
 * in: time tn 	 	      -- t_val
 * in: physical coordinates   -- phys_coords
 * in: reference coordinates on E -- ref_coords_E 
 * in: reference coordinates on En-- ref_coords_En
 * out: numerical flux -- numerical flux
 */
void upwind_flux(int E,
                 int E_neighbour,
                 double* uh_tn_sol_vec,
                 element* mesh_element,
                 node* mesh_node,
		 parameters* user_param,
                 double* normal_e,
		 int s_dim_counter,
                 double t_val,
                 double* phys_coords,
                 double* ref_coords_E,
                 double* ref_coords_En,
                 double* numerical_flux)
{
	double udotn=0.0;
	double adv_vector[2];
	double approx_sol_E[SYSTEM_DIM];
	init_zero_d(approx_sol_E,SYSTEM_DIM);
	advection_v(t_val,phys_coords,user_param,s_dim_counter,adv_vector);
	udotn = adv_vector[0]*normal_e[0] + adv_vector[1]*normal_e[1];
	if(udotn>0)
	{
		get_approx_solution(uh_tn_sol_vec,E,mesh_element,mesh_node,ref_coords_E[0],ref_coords_E[1],approx_sol_E);
	}
	else
	{
		get_approx_solution(uh_tn_sol_vec,E_neighbour,mesh_element,mesh_node,ref_coords_En[0],ref_coords_En[1],approx_sol_E);
	}

	*numerical_flux = approx_sol_E[s_dim_counter]*(adv_vector[0]*normal_e[0] + adv_vector[1]*normal_e[1]);
}

/* void shu_osher_limit
 * performs flux limiting on the solution vector 
 * ASSUMES THAT UH_PRE_LIMIT IS LINEAR
 * controls transformation to characteristic fields and limiting
 * in: index of current triangle     -- E
 * in: mesh_element data structure   -- mesh_element
 * in: mesh_node data structure      -- mesh_node
 * in: mesh_edge data structure      -- mesh_edge
 * in: parameters 		     -- user_param
 * in: number of triangles           -- Nelts
 * in: number of basis functions     -- Nloc
 * in: solution before limiting      -- uh_pre_limit
 * out: flux limited solution        -- uh_limit
 */
/*void shu_osher_limit(int E,
                     element* mesh_element,
                     node* mesh_node,
                     edge* mesh_edge,
                     parameters* user_param,
                     int Nelts,
                     int Nloc,
                     double* uh_pre_limit,
                     double* uh_limit)
{
    int edge1,edge2,edge3;		// Index of edges of E
    int E1,E2,E3;	            // Neighbours of E
    int n1,n2,n3;			    // Nodes of E
    int s_dim;				    // Index for system component
    int idofs;				    // Index of parameter for basis fcts

    int i;				
    int size = SYSTEM_DIM*Nelts*Nloc;

    double uh_E[SYSTEM_DIM]; 		// Values at barycenters
    double uh_E1[SYSTEM_DIM];
    double uh_E2[SYSTEM_DIM];
    double uh_E3[SYSTEM_DIM];

    double m1[2]; 			// Midpoints of edges
    double m2[2];
    double m3[2];
    
    double uh_m1[SYSTEM_DIM]; 		// Values at edge-midpoints
    double uh_m2[SYSTEM_DIM];
    double uh_m3[SYSTEM_DIM];

    double baryctr0[2]; 		//b0
    double baryctr1[2];         //b1
    double baryctr2[2];         //b2
    double baryctr3[2];         //b3

    // Normal vectors pointing form the barycenter to the i-th edge midpoint
    double vec_b0_m1[2]; 			
    double vec_b0_m2[2]; 			
    double vec_b0_m3[2]; 			

    double alpha1[2]; 			// Needed to compare slopes
    double alpha2[2]; 			
    double alpha3[2]; 		

    int alpha1_elements[2];
    int alpha2_elements[2];
    int alpha3_elements[2];

    double diff_E1_E[SYSTEM_DIM]; 	// Difference in call averages
    double diff_E2_E[SYSTEM_DIM];	// between neighbouring triangles
    double diff_E3_E[SYSTEM_DIM];

    // Quantities needed for limiting, before and after transformation
    double D_uh_m1[SYSTEM_DIM]; 	// 'Slope' between neighbouring triangles
    double D_uh_m2[SYSTEM_DIM];
    double D_uh_m3[SYSTEM_DIM];
    
    double diff_m1_E[SYSTEM_DIM];	// 'Slope' within one triangle
    double diff_m2_E[SYSTEM_DIM];
    double diff_m3_E[SYSTEM_DIM];

    double trans_D_uh_m1[SYSTEM_DIM]; 
    double trans_D_uh_m2[SYSTEM_DIM]; 
    double trans_D_uh_m3[SYSTEM_DIM]; 
    
    double trans_diff_m1_E[SYSTEM_DIM];
    double trans_diff_m2_E[SYSTEM_DIM];
    double trans_diff_m3_E[SYSTEM_DIM];

    // Deltas, before and after transformation
    double Delta1_trans[SYSTEM_DIM];
    double Delta2_trans[SYSTEM_DIM];
    double Delta3_trans[SYSTEM_DIM];
    
    double Delta1[SYSTEM_DIM];
    double Delta2[SYSTEM_DIM];
    double Delta3[SYSTEM_DIM];

    // Temporary variables to process Deltas
    double pos = 0; 
    double neg = 0;
    double sigma_p = 0;
    double sigma_m = 0; 	

    double Deltah1[SYSTEM_DIM];
    double Deltah2[SYSTEM_DIM];
    double Deltah3[SYSTEM_DIM];

    double uh_lin_limit[3];

    double uh_local[SYSTEM_DIM];	// Local values of the moments

    int verbose = 0;

    int flag_close_to_real_bound = 0;	
    double n1_norm;			// Needed to check whether the previous flag needs to be set
    //double trans_tol_bound = 1e-4;
    double trans_tol_bound = 0.0;

    double alpha_tol = -1e-10;

    if ((*user_param).flag_limiting == 1) {
	// Initializations
	E1 = E; E2 = E; E3 = E;

	// Edges of E
	edge1 = mesh_element[E].edge[1];
	edge2 = mesh_element[E].edge[2];
	edge3 = mesh_element[E].edge[3];

	// Coordinates of midpoints of edges
	n1 = mesh_element[E].vertex[1];
	n2 = mesh_element[E].vertex[2];
	n3 = mesh_element[E].vertex[3];

	m1[0] = (mesh_node[n1].coord[0]+mesh_node[n2].coord[0])/2.0;
	m1[1] = (mesh_node[n1].coord[1]+mesh_node[n2].coord[1])/2.0;
	m2[0] = (mesh_node[n2].coord[0]+mesh_node[n3].coord[0])/2.0;
	m2[1] = (mesh_node[n2].coord[1]+mesh_node[n3].coord[1])/2.0;
	m3[0] = (mesh_node[n3].coord[0]+mesh_node[n1].coord[0])/2.0;
	m3[1] = (mesh_node[n3].coord[1]+mesh_node[n1].coord[1])/2.0;

	// Compute barycenter
	find_barycenter(E,mesh_element,mesh_node,baryctr0);


	// Value of uh in E at midpoints of edges
	get_approx_solution_lin(uh_pre_limit,E,mesh_element,mesh_node,0.5,0.0,uh_m1);
	get_approx_solution_lin(uh_pre_limit,E,mesh_element,mesh_node,0.5,0.5,uh_m2);
	get_approx_solution_lin(uh_pre_limit,E,mesh_element,mesh_node,0.0,0.5,uh_m3);

	//s_dim = 0;i


	cell_average_lin(uh_pre_limit,E,mesh_element,mesh_node,uh_E);


	// Repeat for all neighbours of E
	if(mesh_edge[edge1].edge_type == INTERIOR){
	    if (mesh_edge[edge1].neighbour[1] == E) E1 = mesh_edge[edge1].neighbour[2]; 
	    else E1 = mesh_edge[edge1].neighbour[1]; 

	    find_barycenter(E1,mesh_element,mesh_node,baryctr1);
	    cell_average_lin(uh_pre_limit,E1,mesh_element,mesh_node,uh_E1);
	}
	// If E is on the boundary, mirror the triangle along its boundary edge to find barycenter and use
	// the value of the midpoint as the cell-average on the ghost triangle
	else if (mesh_edge[edge1].edge_type == EXTERIOR) {
	    E1 = E;
	    baryctr1[0] = baryctr0[0] + 2*(m1[0] - baryctr0[0]); 
	    baryctr1[1] = baryctr0[1] + 2*(m1[1] - baryctr0[1]); 
	    for (s_dim=0;s_dim < SYSTEM_DIM; s_dim++) {
		uh_E1[s_dim] = uh_m1[s_dim];
	    }
	}

	if(mesh_edge[edge2].edge_type == INTERIOR){
	    if (mesh_edge[edge2].neighbour[1] == E)  E2 = mesh_edge[edge2].neighbour[2]; 
	    else E2 = mesh_edge[edge2].neighbour[1];

	    find_barycenter(E2,mesh_element,mesh_node,baryctr2);
	    cell_average_lin(uh_pre_limit,E2,mesh_element,mesh_node,uh_E2);
	}
	else if (mesh_edge[edge2].edge_type == EXTERIOR){
	    E2 = E;
	    baryctr2[0] = baryctr0[0] + 2*(m2[0] - baryctr0[0]); 
	    baryctr2[1] = baryctr0[1] + 2*(m2[1] - baryctr0[1]); 
	    for (s_dim=0;s_dim < SYSTEM_DIM; s_dim++) {
	        uh_E2[s_dim] = uh_m2[s_dim];
	    }
	}

	if(mesh_edge[edge3].edge_type == INTERIOR){
	    if (mesh_edge[edge3].neighbour[1] == E) E3 = mesh_edge[edge3].neighbour[2]; 
	    else E3 = mesh_edge[edge3].neighbour[1]; 

	    find_barycenter(E3,mesh_element,mesh_node,baryctr3);
	    cell_average_lin(uh_pre_limit,E3,mesh_element,mesh_node,uh_E3);
	}
	else if (mesh_edge[edge3].edge_type == EXTERIOR){
	    E3 = E;
	    baryctr3[0] = baryctr0[0] + 2*(m3[0] - baryctr0[0]); 
	    baryctr3[1] = baryctr0[1] + 2*(m3[1] - baryctr0[1]); 
	    for (s_dim=0;s_dim < SYSTEM_DIM; s_dim++) {
	        uh_E3[s_dim] = uh_m3[s_dim];
	    }
	}

	for (s_dim=0;s_dim < SYSTEM_DIM; s_dim++){ 

	    // Get cell average
	    //uh_E[s_dim] = (uh_m1[s_dim] + uh_m2[s_dim] + uh_m3[s_dim])/3.0;

	    diff_E1_E[s_dim] = uh_E1[s_dim] - uh_E[s_dim];
	    diff_E2_E[s_dim] = uh_E2[s_dim] - uh_E[s_dim];
	    diff_E3_E[s_dim] = uh_E3[s_dim] - uh_E[s_dim];
	}

	// Compute alphas and Delta_uh_bar(mi,E)
	compute_alphas(m1,baryctr0,baryctr1,baryctr2,alpha1);
	if ((alpha1[0] < alpha_tol) || (alpha1[1] < alpha_tol)) {
	    compute_alphas(m1,baryctr0,baryctr2,baryctr3,alpha1);
	    if ((alpha1[0] < alpha_tol) || (alpha1[1] < alpha_tol)) {
		compute_alphas(m1,baryctr0,baryctr1,baryctr3,alpha1);
		alpha1_elements[0] = E1;
		alpha1_elements[1] = E3;
		for (s_dim=0;s_dim < SYSTEM_DIM; s_dim++) {
		    D_uh_m1[s_dim] = alpha1[0]*(diff_E1_E[s_dim]) + alpha1[1]*(diff_E3_E[s_dim]);
		}
		if ((alpha1[0] < alpha_tol) || (alpha1[1] < alpha_tol)) {
		    printf("couldn't find non-negative alphas!\n");
		}
	    }
	    else {
		alpha1_elements[0] = E2;
		alpha1_elements[1] = E3;
		for (s_dim=0;s_dim < SYSTEM_DIM; s_dim++) {
		    D_uh_m1[s_dim] = alpha1[0]*(diff_E2_E[s_dim]) + alpha1[1]*(diff_E3_E[s_dim]);
		}
	    }
	}
	else {
	    alpha1_elements[0] = E1;
	    alpha1_elements[1] = E2;
	    for (s_dim=0;s_dim < SYSTEM_DIM; s_dim++) {
		D_uh_m1[s_dim] = alpha1[0]*(diff_E1_E[s_dim]) + alpha1[1]*(diff_E2_E[s_dim]);
	    }
	}

	compute_alphas(m2,baryctr0,baryctr1,baryctr2,alpha2);
	if ((alpha2[0] < alpha_tol) || (alpha2[1] < alpha_tol)) {
	    compute_alphas(m2,baryctr0,baryctr2,baryctr3,alpha2);
	    if ((alpha2[0] < alpha_tol) || (alpha2[1] < alpha_tol)) {
		compute_alphas(m2,baryctr0,baryctr1,baryctr3,alpha2);
		alpha2_elements[0] = E1;
		alpha2_elements[1] = E3;
		for (s_dim=0;s_dim < SYSTEM_DIM; s_dim++) {
		    D_uh_m2[s_dim] = alpha2[0]*(diff_E1_E[s_dim]) + alpha2[1]*(diff_E3_E[s_dim]);
		}
		if ((alpha2[0] < alpha_tol) || (alpha2[1] < alpha_tol)) {
		    printf("couldn't find non-negative alphas!\n");
		}
	    }
	    else {
		alpha2_elements[0] = E2;
		alpha2_elements[1] = E3;
		for (s_dim=0;s_dim < SYSTEM_DIM; s_dim++) {
		    D_uh_m2[s_dim] = alpha2[0]*(diff_E2_E[s_dim]) + alpha2[1]*(diff_E3_E[s_dim]);
		}
	    }
	}
	else {
	    alpha2_elements[0] = E1;
	    alpha2_elements[1] = E2;
	    for (s_dim=0;s_dim < SYSTEM_DIM; s_dim++) {
		D_uh_m2[s_dim] = alpha2[0]*(diff_E1_E[s_dim]) + alpha2[1]*(diff_E2_E[s_dim]);
	    }
	}

	compute_alphas(m3,baryctr0,baryctr1,baryctr2,alpha3);
	if ((alpha3[0] < alpha_tol) || (alpha3[1] < alpha_tol)) {
	    compute_alphas(m3,baryctr0,baryctr2,baryctr3,alpha3);
	    if ((alpha3[0] < alpha_tol) || (alpha3[1] < alpha_tol)) {
		compute_alphas(m3,baryctr0,baryctr1,baryctr3,alpha3);
		alpha3_elements[0] = E1;
		alpha3_elements[1] = E3;
		for (s_dim=0;s_dim < SYSTEM_DIM; s_dim++) {
		    D_uh_m3[s_dim] = alpha3[0]*(diff_E1_E[s_dim]) + alpha3[1]*(diff_E3_E[s_dim]);
		}
		if ((alpha3[0] < alpha_tol) || (alpha3[1] < alpha_tol)) {
		    printf("couldn't find non-negative alphas!\n");
		}
	    }
	    else {
		alpha3_elements[0] = E2;
		alpha3_elements[1] = E3;
		for (s_dim=0;s_dim < SYSTEM_DIM; s_dim++) {
		    D_uh_m3[s_dim] = alpha3[0]*(diff_E2_E[s_dim]) + alpha3[1]*(diff_E3_E[s_dim]);
		}
	    }
	}
	else {
	    alpha3_elements[0] = E1;
	    alpha3_elements[1] = E2;
	    for (s_dim=0;s_dim < SYSTEM_DIM; s_dim++) {
		D_uh_m3[s_dim] = alpha3[0]*(diff_E1_E[s_dim]) + alpha3[1]*(diff_E2_E[s_dim]);
	    }
	}

	// Calculate quantities needed for limiting
	for (s_dim=0;s_dim < SYSTEM_DIM; s_dim++) { 
	    diff_m1_E[s_dim] = uh_m1[s_dim] - uh_E[s_dim];
	    diff_m2_E[s_dim] = uh_m2[s_dim] - uh_E[s_dim];
	    diff_m3_E[s_dim] = uh_m3[s_dim] - uh_E[s_dim];
	}

	// Local solution vector
	get_approx_solution_lin(uh_pre_limit,E,mesh_element,mesh_node,0.33,0.33,uh_local);

	// Check whether we are close to the boundary of the realizability domain
	if (SYSTEM_DIM > 1) {
	    n1_norm = sqrt(uh_local[1]*uh_local[1] + uh_local[2]*uh_local[2]);

	    if (uh_local[0] < trans_tol_bound) {
		flag_close_to_real_bound = 1;
	    }
	    else if (n1_norm > uh_local[0] - trans_tol_bound) {
		flag_close_to_real_bound = 1;
	    }
	}
	flag_close_to_real_bound = 0;

	vec_b0_m1[0] = m1[0] - baryctr0[0];
	vec_b0_m1[1] = m1[1] - baryctr0[1];
	vec_b0_m2[0] = m2[0] - baryctr0[0];
	vec_b0_m2[1] = m2[1] - baryctr0[1];
	vec_b0_m3[0] = m3[0] - baryctr0[0];
	vec_b0_m3[1] = m3[1] - baryctr0[1];

 	// N needs to be normal for the transformation
	normalize_vector(2,vec_b0_m1);
	normalize_vector(2,vec_b0_m2);
	normalize_vector(2,vec_b0_m3);

	if (( (*user_param).flag_trans == 1) && (flag_close_to_real_bound == 0) )  {
	    // Transformation to characteristic fields
	   

        transform_to_char_var(uh_local,diff_m1_E,vec_b0_m1,user_param,trans_diff_m1_E,TRANS_TO_CHAR); 
        transform_to_char_var(uh_local,diff_m2_E,vec_b0_m2,user_param,trans_diff_m2_E,TRANS_TO_CHAR); 
	    transform_to_char_var(uh_local,diff_m3_E,vec_b0_m3,user_param,trans_diff_m3_E,TRANS_TO_CHAR); 

        //transf_to_char_var(uh_local,diff_m1_E,vec_b0_m1,user_param,trans_diff_m1_E); 
	    //transf_to_char_var(uh_local,diff_m2_E,vec_b0_m2,user_param,trans_diff_m2_E); 
	    //transf_to_char_var(uh_local,diff_m3_E,vec_b0_m3,user_param,trans_diff_m3_E);

	    transform_to_char_var(uh_local,D_uh_m1,vec_b0_m1,user_param,trans_D_uh_m1,TRANS_TO_CHAR); 
	    transform_to_char_var(uh_local,D_uh_m2,vec_b0_m2,user_param,trans_D_uh_m2,TRANS_TO_CHAR); 
	    transform_to_char_var(uh_local,D_uh_m3,vec_b0_m3,user_param,trans_D_uh_m3,TRANS_TO_CHAR); 
	}
	else {
	    for (s_dim=0;s_dim < SYSTEM_DIM; s_dim++) { 
		trans_diff_m1_E[s_dim] = diff_m1_E[s_dim];
		trans_diff_m2_E[s_dim] = diff_m2_E[s_dim];
		trans_diff_m3_E[s_dim] = diff_m3_E[s_dim];

		trans_D_uh_m1[s_dim] = D_uh_m1[s_dim];
		trans_D_uh_m2[s_dim] = D_uh_m2[s_dim];
		trans_D_uh_m3[s_dim] = D_uh_m3[s_dim];
	    }
	}

	// Perform actual limiting on each component of the (transformed) system
	for (s_dim=0;s_dim < SYSTEM_DIM; s_dim++){ 
	    // Compute Deltas
	    Delta1_trans[s_dim] = minmod(trans_diff_m1_E[s_dim],(*user_param).nu*trans_D_uh_m1[s_dim],(*user_param).mdx2);
	    Delta2_trans[s_dim] = minmod(trans_diff_m2_E[s_dim],(*user_param).nu*trans_D_uh_m2[s_dim],(*user_param).mdx2);
	    Delta3_trans[s_dim] = minmod(trans_diff_m3_E[s_dim],(*user_param).nu*trans_D_uh_m3[s_dim],(*user_param).mdx2);
	}

	if (( (*user_param).flag_trans == 1) && (flag_close_to_real_bound == 0) )  {
	    // transform back
	    //transf_from_char_var(uh_local,Delta1_trans,vec_b0_m1,user_param,Delta1);
	    //transf_from_char_var(uh_local,Delta2_trans,vec_b0_m2,user_param,Delta2);
	    //transf_from_char_var(uh_local,Delta3_trans,vec_b0_m3,user_param,Delta3);

        transform_to_char_var(uh_local,Delta1_trans,vec_b0_m1,user_param,Delta1,TRANS_FROM_CHAR);
	    transform_to_char_var(uh_local,Delta2_trans,vec_b0_m2,user_param,Delta2,TRANS_FROM_CHAR);
	    transform_to_char_var(uh_local,Delta3_trans,vec_b0_m3,user_param,Delta3,TRANS_FROM_CHAR);	
    }
	else {
	    for (s_dim=0;s_dim < SYSTEM_DIM; s_dim++){ 
		Delta1[s_dim] = Delta1_trans[s_dim];
		Delta2[s_dim] = Delta2_trans[s_dim];
		Delta3[s_dim] = Delta3_trans[s_dim];
	    }
	}

	// Perform actual limiting on each component of the transformed system
	for (s_dim=0;s_dim < SYSTEM_DIM; s_dim++){
	    // Limit only linear part, i.e. up to third coefficient...
	    if (fabs(Delta1[s_dim]+Delta2[s_dim]+Delta3[s_dim]) > 1e-12) {
		pos =  max(Delta1[s_dim],0.0);
		pos += max(Delta2[s_dim],0.0);
		pos += max(Delta3[s_dim],0.0);

		neg =  max(-Delta1[s_dim],0.0);
		neg += max(-Delta2[s_dim],0.0);
		neg += max(-Delta3[s_dim],0.0);

		sigma_p = min(1.0,neg/pos);
		sigma_m = min(1.0,pos/neg);

		Deltah1[s_dim] = sigma_p*max(0.0,Delta1[s_dim]) - sigma_m*max(0.0,-Delta1[s_dim]);
		Deltah2[s_dim] = sigma_p*max(0.0,Delta2[s_dim]) - sigma_m*max(0.0,-Delta2[s_dim]);
		Deltah3[s_dim] = sigma_p*max(0.0,Delta3[s_dim]) - sigma_m*max(0.0,-Delta3[s_dim]);

		uh_lin_limit[0] = uh_E[s_dim]+Deltah1[s_dim]-Deltah2[s_dim]+Deltah3[s_dim];
		uh_lin_limit[1] = 2.0*(Deltah2[s_dim]-Deltah3[s_dim]);
		uh_lin_limit[2] = 2.0*(Deltah2[s_dim]-Deltah1[s_dim]);
	    }
	    else {
		uh_lin_limit[0] = uh_E[s_dim]+Delta1[s_dim]-Delta2[s_dim]+Delta3[s_dim];
		uh_lin_limit[1] = 2.0*(Delta2[s_dim]-Delta3[s_dim]);
		uh_lin_limit[2] = 2.0*(Delta2[s_dim]-Delta1[s_dim]);
	    }

	    // Check whether limiting changes the solution
	    if ( (fabs(uh_lin_limit[0] - uh_pre_limit[SYSTEM_DIM*(E-1)*3 + s_dim*3  ]) > 1e-16)
	      || (fabs(uh_lin_limit[1] - uh_pre_limit[SYSTEM_DIM*(E-1)*3 + s_dim*3+1]) > 1e-16)
	      || (fabs(uh_lin_limit[2] - uh_pre_limit[SYSTEM_DIM*(E-1)*3 + s_dim*3+2]) > 1e-16)) {
        //    printf("limiting changes!...\n");
            

		// ... if it does, copy the linear part...
		uh_limit[SYSTEM_DIM*(E-1)*Nloc+s_dim*Nloc]     = uh_lin_limit[0];
		uh_limit[SYSTEM_DIM*(E-1)*Nloc+s_dim*Nloc + 1] = uh_lin_limit[1];
		uh_limit[SYSTEM_DIM*(E-1)*Nloc+s_dim*Nloc + 2] = uh_lin_limit[2];

		// ... and cut off the higher orders
		if(Nloc == 6)
        {
            for(idofs=3;idofs<Nloc;idofs++) 
            {
                uh_limit[SYSTEM_DIM*(E-1)*Nloc+s_dim*Nloc + idofs] = 0.0;
            }
        }

		if (verbose == 1) {
		    printf("m1 %lf,%lf :: uh_m1 %lf\n",m1[0],m1[1],uh_m1[s_dim]);
		    printf("m2 %lf,%lf :: uh_m2 %lf\n",m2[0],m2[1],uh_m2[s_dim]);
		    printf("m3 %lf,%lf :: uh_m3 %lf\n",m3[0],m3[1],uh_m3[s_dim]);
		    printf("E %d :: uh_E %lf :: bary0 %lf,%lf\nE1 %d :: uh_E1 %lf :: bary1 %lf,%lf\nE2 %d :: uh_E2 %lf :: bary2 %lf,%lf\nE3 %d :: uh_E3 %lf :: bary3 %lf,%lf\n",E,uh_E[s_dim],baryctr0[0],baryctr0[1],E1,uh_E1[s_dim],baryctr1[0],baryctr1[1],E2,uh_E2[s_dim],baryctr2[0],baryctr2[1],E3,uh_E3[s_dim],baryctr3[0],baryctr3[1]);
		    printf("alpha1 :: elements = [%d, %d] \t vals = %.16lf,%.16lf\n", alpha1_elements[0], alpha1_elements[1], alpha1[0],alpha1[1]);
		    printf("alpha2 :: elements = [%d, %d] \t vals = %.16lf,%.16lf\n", alpha2_elements[0], alpha2_elements[1], alpha2[0],alpha2[1]);
		    printf("alpha3 :: elements = [%d, %d] \t vals = %.16lf,%.16lf\n", alpha3_elements[0], alpha3_elements[1], alpha3[0],alpha3[1]);
		    printf("D_uh_m1 %.16lf :: D_uh_m2 %.16lf :: D_uh_m3 %.16lf\n",D_uh_m1[s_dim],D_uh_m2[s_dim],D_uh_m3[s_dim]);
		    printf("diff_m1_E %.16lf :: diff_m2_E %.16lf :: diff_m3_E %.16lf\n",diff_m1_E[s_dim],diff_m2_E[s_dim],diff_m3_E[s_dim]);
		    printf("uh_m1-uh_e %lf :: nu*D_uh_m1 %lf\n",uh_m1[s_dim]-uh_E[s_dim],(*user_param).nu*D_uh_m1[s_dim]);
		    printf("uh_m2-uh_e %lf :: nu*D_uh_m2 %lf\n",uh_m2[s_dim]-uh_E[s_dim],(*user_param).nu*D_uh_m2[s_dim]);
		    printf("uh_m3-uh_e %lf :: nu*D_uh_m3 %lf\n",uh_m3[s_dim]-uh_E[s_dim],(*user_param).nu*D_uh_m3[s_dim]);
		    printf("pos %lf :: neg %lf :: sigma_p %lf :: sigma_m %lf\n",pos,neg,sigma_p,sigma_m);
		    printf("        Delta1 %.16lf :: Delta2 %.16lf :: Delta3 :: %.16lf\n",Delta1[s_dim],Delta2[s_dim],Delta3[s_dim]);
		    printf("Hatted: Delta1 %.16lf :: Delta2 %.16lf :: Delta3 :: %.16lf\n",Deltah1[s_dim],Deltah2[s_dim],Deltah3[s_dim]);
		    printf("Sum Delta hat: %.16lf\n",Deltah1[s_dim]+Deltah2[s_dim]+Deltah3[s_dim]);

		    printf("Linear part before limiting: alpha_1 %lf :: alpha_2 %lf :: alpha_3 %lf\n",uh_pre_limit[SYSTEM_DIM*(E-1)*3 + s_dim*3],uh_pre_limit[SYSTEM_DIM*(E-1)*3 + s_dim*3+1],uh_pre_limit[SYSTEM_DIM*(E-1)*3 + s_dim*3+2]);
		    printf("Linear part after limiting: alpha_1 %lf :: alpha_2 %lf :: alpha_3 %lf\n",uh_lin_limit[0],uh_lin_limit[1],uh_lin_limit[2]);
		    printf("Difference of both: alpha_1 %1.12lf :: alpha_2 %1.12lf :: alpha_3 %1.12lf\n",uh_pre_limit[SYSTEM_DIM*(E-1)*3 + s_dim*3] - uh_lin_limit[0],uh_pre_limit[SYSTEM_DIM*(E-1)*3 + s_dim*3+1] - uh_lin_limit[1],uh_pre_limit[SYSTEM_DIM*(E-1)*3 + s_dim*3+2] - uh_lin_limit[2]);
		    printf("After limiting: alpha_1 %lf :: alpha_2 %lf :: alpha_3 %lf\n",uh_limit[SYSTEM_DIM*(E-1)*Nloc + s_dim*Nloc],uh_limit[SYSTEM_DIM*(E-1)*Nloc + s_dim*Nloc+1],uh_limit[SYSTEM_DIM*(E-1)*Nloc + s_dim*Nloc+2]);
		    printf("==============================\n");
		    getchar();
		}
	    }
	}
    }
}*/

/* void normalize_vector
 * normalizes a vector of length n
 * in: length of the vector -- n
 * in: vector -- vec
 * out: normalized vector -- vec
 */
void normalize_vector(int n, 
		      double* vec)
{
    double norm;
    int i;

    norm = 0;
    for (i=0;i<n;i++) norm += vec[i]*vec[i];
    norm = sqrt(norm);

    for (i=0;i<n;i++) vec[i] = vec[i]/norm;
}

void sum_vectors(int E,
                 double* output_vec,
                 double* input_vec1,
                 double input_vec1_coef,
                 double* input_vec2,
                 double input_vec2_coef,
                 int Nloc)
{
    int s_dim =0;
    int idofs =0;

    for(s_dim=0;s_dim < SYSTEM_DIM;s_dim++)
    {
        //printf("s_dim %d \n",s_dim);
        //getchar();
        for(idofs=0;idofs<Nloc;idofs++)
        {
            output_vec[SYSTEM_DIM*(E-1)*Nloc + (Nloc*s_dim)+idofs] = input_vec1_coef*input_vec1[SYSTEM_DIM*(E-1)*Nloc + (Nloc*s_dim)+idofs] 
                                                                   + input_vec2_coef*input_vec2[SYSTEM_DIM*(E-1)*Nloc + (Nloc*s_dim)+idofs];
        }
        //getchar();
    }
}

void fractional_step(int E,
                     node* mesh_node,
                     element* mesh_element,
                     edge* mesh_edge,
                     parameters* user_param,
                     double current_time,
                     double delta_t,
                     double* frac_step_output_vec,
                     double* uh_tn_sol_vec,
                     double uh_tn_sol_vec_coef,
                     double* frac_input_vec,
                     double frac_input_vec_coef,
                     double pde_rhs_coef,
                     int Nloc)
{
    int s_dim =0;
    int i=0;
    int idofs=0;
    double* uzero_loc = NULL;
    double* Aloc_matrix = NULL;
    double* temp_Aloc = NULL;
    double* loc_pde_rhs_vec = NULL;
    double* loc_source_rhs_vec = NULL;
    double* sigmas_vec = NULL;
    double* sigmat_vec = NULL;

    // QR-decomposition and solving
    double** A_loc_QR;
    double* Q_coeff;
    double* R_coeff;
    int sing_decomp;

    // solution vector
    double* uh_pre_limit = NULL;
    double* uh_limit = NULL;

    // intermediate solution vectors
    double* uh_aux1 = NULL;
    double* uh_aux2 = NULL;

    //allocate memory
    uzero_loc = (double*)malloc(Nloc*sizeof(double));
    Aloc_matrix = (double*)malloc((Nloc*Nloc)*sizeof(double));
    temp_Aloc = (double*)malloc((Nloc*Nloc)*sizeof(double));
    loc_pde_rhs_vec = (double*)malloc(Nloc*sizeof(double));
    sigmas_vec = (double*)malloc(Nloc*sizeof(double));
    sigmat_vec = (double*)malloc(Nloc*sizeof(double));
    loc_source_rhs_vec = (double*)malloc(Nloc*sizeof(double));

    Q_coeff = (double*)malloc((Nloc+1)*sizeof(double));
    R_coeff = (double*)malloc((Nloc+1)*sizeof(double));
    A_loc_QR = (double**)malloc((Nloc+1)*sizeof(double*));

    for(i=0;i<Nloc+1;i++)
    A_loc_QR[i] = (double*)malloc((Nloc+1)*sizeof(double));

    // local rhs vector for element E
    init_zero_d(Aloc_matrix,Nloc*Nloc);
    init_zero_m(A_loc_QR,Nloc,Nloc);
    Alocal_mat(E,mesh_element,mesh_node,user_param,Aloc_matrix);
    // Perform QR-decomposition
    assemble_qr(Aloc_matrix,Nloc,Q_coeff,R_coeff,A_loc_QR,&sing_decomp);

    for(s_dim=0;s_dim < SYSTEM_DIM;s_dim++)
    {
        init_zero_d(loc_pde_rhs_vec,Nloc);
        init_zero_d(loc_source_rhs_vec,Nloc);
        init_zero_d(sigmas_vec,Nloc);
        init_zero_d(sigmat_vec,Nloc);

        discrete_pde_rhs(E,Nloc,mesh_node,mesh_element,mesh_edge,user_param,frac_input_vec,s_dim,current_time,loc_pde_rhs_vec);
        //add source term
        assemble_source(E,mesh_element,mesh_node,user_param,frac_input_vec,s_dim,current_time,loc_source_rhs_vec);
        //add scattering
        scattering(E,frac_input_vec,mesh_element,mesh_node,user_param,s_dim,current_time,sigmas_vec,sigmat_vec);

         // Assemble DG RHS vector
        for(idofs=0;idofs<Nloc;idofs++)
        {
            loc_pde_rhs_vec[idofs] += loc_source_rhs_vec[idofs]+ sigmas_vec[idofs] + sigmat_vec[idofs];
        }

        // we have an ODE system of the form Aloc_matrix d/dt uh = L_h(uh)
        // Solve the linear system
        qrsolv(A_loc_QR,Nloc,Q_coeff,R_coeff,loc_pde_rhs_vec);
        for(idofs=0;idofs<Nloc;idofs++)
        {
            frac_step_output_vec[SYSTEM_DIM*(E-1)*Nloc + (Nloc*s_dim)+idofs] = uh_tn_sol_vec_coef*uh_tn_sol_vec[SYSTEM_DIM*(E-1)*Nloc+(Nloc*s_dim)+idofs] 
                                                                            +frac_input_vec_coef*frac_input_vec[SYSTEM_DIM*(E-1)*Nloc+(Nloc*s_dim)+idofs]
                                                                            +pde_rhs_coef*delta_t*loc_pde_rhs_vec[idofs];

        }//Nloc
    }//sys_dim

    for(i=0;i<Nloc+1;i++)
    {
        free(A_loc_QR[i]);
    }
    free(uzero_loc);
    free(Aloc_matrix);
    free(temp_Aloc);
    free(loc_pde_rhs_vec);
    free(sigmas_vec);
    free(sigmat_vec);
    free(loc_source_rhs_vec);
    free(Q_coeff);
    free(R_coeff);
    free(A_loc_QR);
}

/* void project_quadr_sol_vec_to_lin
 * orthogonal L2 projection of pieveise quadratic solution 
 * to the space of piecewise linear functions
 * in: current element 			   -- E
 * in: coefficients of p.w. quadr. function -- uh_sol
 * out: coefficients of p.w. lin. function -- vec
 */
void project_quadr_sol_vec_to_lin(int E,
				  double* uh_sol,
			    	  double* uh_lin)
{
    int idofs;
    int s_dim;
    double alpha[6]; // Assumes p.w. quadr. solution!
    double gamma[3];

    for(s_dim=0;s_dim < SYSTEM_DIM;s_dim++) 
    {
        for(idofs=0;idofs<6;idofs++) 
        {
            alpha[idofs] = uh_sol[SYSTEM_DIM*(E-1)*6+s_dim*6+idofs];
        }
        gamma[0] = alpha[0] - 1.0/10.0*(alpha[3] + alpha[5]) - 1.0/20.0*alpha[4];
        gamma[1] = alpha[1] + 4.0/5.0*alpha[3] + 1.0/5.0*alpha[4];
        gamma[2] = alpha[2] + 1.0/5.0*alpha[4] + 4.0/5.0*alpha[5];

        for(idofs=0;idofs<3;idofs++) 
        {
            uh_lin[SYSTEM_DIM*(E-1)*3+(s_dim*3+idofs)] = gamma[idofs];
        }
    }

}
void project_solution_to_linear(int E,
                            element* mesh_element,
                            node* mesh_node,
                            parameters* user_param,
                            double* uh_tn_sol,
                            double* uh_linear)
{
    double* linear_Aloc_matrix = NULL;
    double quad_sol_linear_phi_vec[linear_NLOC];
    int s_dim=0;
    // QR-decomposition and solving
    double** A_loc_QR; 
    double* Q_coeff; 
    double* R_coeff; 
    int sing_decomp;
    int i,j;
    int idofs=0;

    Q_coeff = (double*)malloc((linear_NLOC+1)*sizeof(double));
    R_coeff = (double*)malloc((linear_NLOC+1)*sizeof(double));
    A_loc_QR = (double**)malloc((linear_NLOC+1)*sizeof(double*));
    linear_Aloc_matrix = (double*)malloc((linear_NLOC*linear_NLOC)*sizeof(double));

    for(i=0;i<linear_NLOC+1;i++)
    A_loc_QR[i] = (double*)malloc((linear_NLOC+1)*sizeof(double));


    for(s_dim =0;s_dim <SYSTEM_DIM;s_dim++)
    {

        init_zero_d(linear_Aloc_matrix,linear_NLOC*linear_NLOC);
        init_zero_d(quad_sol_linear_phi_vec,linear_NLOC);

        linear_Aloc_mat(E,mesh_element,mesh_node,user_param,linear_Aloc_matrix);
        assemble_qr(linear_Aloc_matrix,linear_NLOC,Q_coeff,R_coeff,A_loc_QR,&sing_decomp);
        init_zero_d(quad_sol_linear_phi_vec,linear_NLOC);
        dg_sol_linear_phi(E,mesh_element,mesh_node,user_param,uh_tn_sol,s_dim,quad_sol_linear_phi_vec);
        //solve linear system
        qrsolv(A_loc_QR,linear_NLOC,Q_coeff,R_coeff,quad_sol_linear_phi_vec);

        //printf("E = %d s_dim = %d \n",E,s_dim);
        for(idofs=0;idofs<linear_NLOC;idofs++)
        {
            uh_linear[SYSTEM_DIM*(E-1)*linear_NLOC + (s_dim*linear_NLOC) +idofs] = uh_tn_sol[SYSTEM_DIM*(E-1)*linear_NLOC + (s_dim*linear_NLOC) +idofs]; //quad_sol_linear_phi_vec[idofs];
            //if(quad_sol_linear_phi_vec[idofs] > EPSILON)
               // printf("%10.7e \n",uh_linear[SYSTEM_DIM*(E-1)*linear_NLOC + (s_dim*linear_NLOC) +idofs]);
            
            //assert(fabs(quad_sol_linear_phi_vec[idofs]-uh_tn_sol[SYSTEM_DIM*(E-1)*linear_NLOC + (s_dim*linear_NLOC) +idofs]) < EPSILON);
        }
        //printf("----------------------------\n");


    }

    free(Q_coeff);
    free(R_coeff);
    free(linear_Aloc_matrix);
    for(i=0;i<linear_NLOC+1;i++)
        free(A_loc_QR[i]);
}
void dg_sol_linear_phi(int E,
                         element* mesh_element,
                         node* mesh_node,
                         parameters* user_param,
                         double* uh_tn_sol,
                         int s_dim,
                         double* quad_sol_linear_phi_vec)
{
    int idofs=0;
    int k=0;
    double linear_phi_vec[3];
    double quad_sol[SYSTEM_DIM];

    init_zero_d(quad_sol_linear_phi_vec,linear_NLOC);
    for(k=1;k<ngpts+1;k++) 
    {
        //get quadratic solution at gauss points
        init_zero_d(quad_sol,SYSTEM_DIM);
        get_approx_solution(uh_tn_sol,E,mesh_element,mesh_node,mesh_element[E].el_gpts_ref_x[k],mesh_element[E].el_gpts_ref_y[k],quad_sol);
        limiting_linear_basis_functions(mesh_element[E].el_gpts_ref_x[k],mesh_element[E].el_gpts_ref_y[k],linear_phi_vec);
        for(idofs =0; idofs < linear_NLOC;idofs++) 
        {
            quad_sol_linear_phi_vec[idofs]+=mesh_element[E].det*(quad_sol[s_dim]*linear_phi_vec[idofs])*mesh_element[E].el_gpts_w[k];

        }//idofs
    }//loop over quad points
}
void limiting_linear_basis_functions(double x,
                                    double y,
                                    double* phi_linear_vec)
{
   
    phi_linear_vec[0] = 1.0 - 2.0*y;
    phi_linear_vec[1] = 2.0*x+2.0*y-1.0;
    phi_linear_vec[2] = 1.0-2.0*x;
    
}

void linear_Aloc_mat(int E,
                     element* mesh_element,
                     node* mesh_node,
                     parameters* user_param,
                     double* linear_Aloc_matrix)
{
    int idofs=0;
    int jdofs=0;
    int k=0;
    double linear_phi_vec[linear_NLOC];

    init_zero_d(linear_Aloc_matrix,linear_NLOC*linear_NLOC);
    for(k=1;k<ngpts+1;k++) 
    {
        //evaluate linear basis functions at gauss points
        init_zero_d(linear_phi_vec,linear_NLOC);
        limiting_linear_basis_functions(mesh_element[E].el_gpts_ref_x[k],mesh_element[E].el_gpts_ref_y[k],linear_phi_vec);
        for(idofs =1; idofs < linear_NLOC+1;idofs++) 
        {
            for(jdofs =1;jdofs < linear_NLOC+1;jdofs++) 
            {
                linear_Aloc_matrix(idofs,jdofs,linear_NLOC) +=mesh_element[E].det*(linear_phi_vec[idofs-1]*linear_phi_vec[jdofs-1])*mesh_element[E].el_gpts_w[k];
            }//j
        }//
    }//loop over quad points
}

/* void shu_osher_limit
 * performs flux limiting on the solution vector 
 * ASSUMES THAT UH_PRE_LIMIT IS LINEAR
 * controls transformation to characteristic fields and limiting
 * in: index of current triangle     -- E
 * in: mesh_element data structure   -- mesh_element
 * in: mesh_node data structure      -- mesh_node
 * in: mesh_edge data structure      -- mesh_edge
 * in: parameters 		     -- user_param
 * in: number of triangles           -- Nelts
 * in: number of basis functions     -- Nloc
 * in: solution before limiting      -- uh_pre_limit
 * out: flux limited solution        -- uh_limit
 */
void shu_osher_limit(int E,
                     element* mesh_element,
                     node* mesh_node,
                     edge* mesh_edge,
                     parameters* user_param,
                     int Nelts,
                     int Nloc,
                     double* uh_pre_limit,
                     double* uh_limit)
{
    int edge1,edge2,edge3;		// Index of edges of E
    int E1,E2,E3;	            // Neighbours of E
    int n1,n2,n3;			    // Nodes of E
    int s_dim;				    // Index for system component
    int idofs;				    // Index of parameter for basis fcts

    int i;				
    int size = SYSTEM_DIM*Nelts*Nloc;

    double uh_E[SYSTEM_DIM]; 		// Values at barycenters
    double uh_E1[SYSTEM_DIM];
    double uh_E2[SYSTEM_DIM];
    double uh_E3[SYSTEM_DIM];

    double m1[2]; 			// Midpoints of edges
    double m2[2];
    double m3[2];
    
    double uh_m1[SYSTEM_DIM]; 		// Values at edge-midpoints
    double uh_m2[SYSTEM_DIM];
    double uh_m3[SYSTEM_DIM];

    double baryctr0[2]; 		//b0
    double baryctr1[2];         //b1
    double baryctr2[2];         //b2
    double baryctr3[2];         //b3

    // Normal vectors pointing form the barycenter to the i-th edge midpoint
    double vec_b0_m1[2]; 			
    double vec_b0_m2[2]; 			
    double vec_b0_m3[2]; 			

    double alpha1[2]; 			// Needed to compare slopes
    double alpha2[2]; 			
    double alpha3[2]; 		

    int alpha1_elements[2];
    int alpha2_elements[2];
    int alpha3_elements[2];

    double diff_E1_E[SYSTEM_DIM]; 	// Difference in call averages
    double diff_E2_E[SYSTEM_DIM];	// between neighbouring triangles
    double diff_E3_E[SYSTEM_DIM];

    // Quantities needed for limiting, before and after transformation
    double D_uh_m1[SYSTEM_DIM]; 	// 'Slope' between neighbouring triangles
    double D_uh_m2[SYSTEM_DIM];
    double D_uh_m3[SYSTEM_DIM];
    
    double diff_m1_E[SYSTEM_DIM];	// 'Slope' within one triangle
    double diff_m2_E[SYSTEM_DIM];
    double diff_m3_E[SYSTEM_DIM];

    double trans_D_uh_m1[SYSTEM_DIM]; 
    double trans_D_uh_m2[SYSTEM_DIM]; 
    double trans_D_uh_m3[SYSTEM_DIM]; 
    
    double trans_diff_m1_E[SYSTEM_DIM];
    double trans_diff_m2_E[SYSTEM_DIM];
    double trans_diff_m3_E[SYSTEM_DIM];

    // Deltas, before and after transformation
    double Delta1_trans[SYSTEM_DIM];
    double Delta2_trans[SYSTEM_DIM];
    double Delta3_trans[SYSTEM_DIM];
    
    double Delta1[SYSTEM_DIM];
    double Delta2[SYSTEM_DIM];
    double Delta3[SYSTEM_DIM];

    // Temporary variables to process Deltas
    double pos = 0; 
    double neg = 0;
    double sigma_p = 0;
    double sigma_m = 0; 	

    double Deltah1[SYSTEM_DIM];
    double Deltah2[SYSTEM_DIM];
    double Deltah3[SYSTEM_DIM];

    double uh_lin_limit[3];

    double uh_local[SYSTEM_DIM];	// Local values of the moments

    int verbose = 0;

    int flag_close_to_real_bound = 0;	
    double n1_norm;			// Needed to check whether the previous flag needs to be set
    double trans_tol_bound = 0.0;

    double alpha_tol = -1e-10;

    if ((*user_param).flag_limiting == 1) 
    {
        E1 = E; 
        E2 = E; 
        E3 = E;
        // Edges of E
        edge1 = mesh_element[E].edge[1];
        edge2 = mesh_element[E].edge[2];
        edge3 = mesh_element[E].edge[3];

        // Coordinates of midpoints of edges
        n1 = mesh_element[E].vertex[1];
        n2 = mesh_element[E].vertex[2];
        n3 = mesh_element[E].vertex[3];

        m1[0] = (mesh_node[n1].coord[0]+mesh_node[n2].coord[0])/2.0;
        m1[1] = (mesh_node[n1].coord[1]+mesh_node[n2].coord[1])/2.0;
        m2[0] = (mesh_node[n2].coord[0]+mesh_node[n3].coord[0])/2.0;
        m2[1] = (mesh_node[n2].coord[1]+mesh_node[n3].coord[1])/2.0;
        m3[0] = (mesh_node[n3].coord[0]+mesh_node[n1].coord[0])/2.0;
        m3[1] = (mesh_node[n3].coord[1]+mesh_node[n1].coord[1])/2.0;
        
        find_barycenter(E,mesh_element,mesh_node,baryctr0);


        // Value of uh in E at midpoints of edges
        get_approx_solution_lin(uh_pre_limit,E,mesh_element,mesh_node,0.5,0.0,uh_m1);
        get_approx_solution_lin(uh_pre_limit,E,mesh_element,mesh_node,0.5,0.5,uh_m2);
        get_approx_solution_lin(uh_pre_limit,E,mesh_element,mesh_node,0.0,0.5,uh_m3);
        
        cell_average_lin(uh_pre_limit,E,mesh_element,mesh_node,uh_E);
        
        
        if(mesh_edge[edge1].edge_type == INTERIOR)
        {
            if (mesh_edge[edge1].neighbour[1] == E) 
                E1 = mesh_edge[edge1].neighbour[2]; 
            else 
                E1 = mesh_edge[edge1].neighbour[1]; 
            find_barycenter(E1,mesh_element,mesh_node,baryctr1);
            cell_average_lin(uh_pre_limit,E1,mesh_element,mesh_node,uh_E1);
        }
        else if (mesh_edge[edge1].edge_type == EXTERIOR) 
        {
            // If E is on the boundary, mirror the triangle along its boundary edge to find barycenter and use
	        // the value of the midpoint as the cell-average on the ghost triangle
            E1 = E;
            baryctr1[0] = baryctr0[0] + 2*(m1[0] - baryctr0[0]); 
            baryctr1[1] = baryctr0[1] + 2*(m1[1] - baryctr0[1]); 
            for (s_dim=0;s_dim < SYSTEM_DIM; s_dim++) 
            {
                uh_E1[s_dim] = uh_m1[s_dim];
            }
        }
        if(mesh_edge[edge2].edge_type == INTERIOR)
        {
            if (mesh_edge[edge2].neighbour[1] == E)  
                E2 = mesh_edge[edge2].neighbour[2]; 
            else 
                E2 = mesh_edge[edge2].neighbour[1];
            find_barycenter(E2,mesh_element,mesh_node,baryctr2);
            cell_average_lin(uh_pre_limit,E2,mesh_element,mesh_node,uh_E2);
        }
        else if (mesh_edge[edge2].edge_type == EXTERIOR)
        {
            E2 = E;
            baryctr2[0] = baryctr0[0] + 2*(m2[0] - baryctr0[0]); 
            baryctr2[1] = baryctr0[1] + 2*(m2[1] - baryctr0[1]); 
            for (s_dim=0;s_dim < SYSTEM_DIM; s_dim++) 
            {
                uh_E2[s_dim] = uh_m2[s_dim];
            }
        }
        if(mesh_edge[edge3].edge_type == INTERIOR)
        {
            if (mesh_edge[edge3].neighbour[1] == E) 
                E3 = mesh_edge[edge3].neighbour[2]; 
            else 
                E3 = mesh_edge[edge3].neighbour[1]; 
            find_barycenter(E3,mesh_element,mesh_node,baryctr3);
            cell_average_lin(uh_pre_limit,E3,mesh_element,mesh_node,uh_E3);
        }
        else if (mesh_edge[edge3].edge_type == EXTERIOR)
        {
            E3 = E;
            baryctr3[0] = baryctr0[0] + 2*(m3[0] - baryctr0[0]); 
            baryctr3[1] = baryctr0[1] + 2*(m3[1] - baryctr0[1]); 
            for (s_dim=0;s_dim < SYSTEM_DIM; s_dim++) 
            {
                uh_E3[s_dim] = uh_m3[s_dim];
            }
        }
        for (s_dim=0;s_dim < SYSTEM_DIM; s_dim++)
        { 
            diff_E1_E[s_dim] = uh_E1[s_dim] - uh_E[s_dim];
            diff_E2_E[s_dim] = uh_E2[s_dim] - uh_E[s_dim];
            diff_E3_E[s_dim] = uh_E3[s_dim] - uh_E[s_dim];
        }
	// Compute alphas and Delta_uh_bar(mi,E)
       compute_alphas(m1,baryctr0,baryctr1,baryctr2,alpha1);
       if ((alpha1[0] < alpha_tol) || (alpha1[1] < alpha_tol)) 
       {
           compute_alphas(m1,baryctr0,baryctr2,baryctr3,alpha1);
           if ((alpha1[0] < alpha_tol) || (alpha1[1] < alpha_tol)) 
           {
                compute_alphas(m1,baryctr0,baryctr1,baryctr3,alpha1);
                alpha1_elements[0] = E1;
                alpha1_elements[1] = E3;
                for (s_dim=0;s_dim < SYSTEM_DIM; s_dim++) 
                {
                    D_uh_m1[s_dim] = alpha1[0]*(diff_E1_E[s_dim]) + alpha1[1]*(diff_E3_E[s_dim]);
                }
                if ((alpha1[0] < alpha_tol) || (alpha1[1] < alpha_tol)) 
                {
                    printf("couldn't find non-negative alphas!\n");
                }
            }
            else 
            {
                alpha1_elements[0] = E2;
                alpha1_elements[1] = E3;
                for (s_dim=0;s_dim < SYSTEM_DIM; s_dim++) 
                {
                    D_uh_m1[s_dim] = alpha1[0]*(diff_E2_E[s_dim]) + alpha1[1]*(diff_E3_E[s_dim]);
                }
            }
        }
        else 
        {
            alpha1_elements[0] = E1;
            alpha1_elements[1] = E2;
            for (s_dim=0;s_dim < SYSTEM_DIM; s_dim++) 
            {
                D_uh_m1[s_dim] = alpha1[0]*(diff_E1_E[s_dim]) + alpha1[1]*(diff_E2_E[s_dim]);
            }
        }
        compute_alphas(m2,baryctr0,baryctr1,baryctr2,alpha2);
        if ((alpha2[0] < alpha_tol) || (alpha2[1] < alpha_tol)) 
        {
            compute_alphas(m2,baryctr0,baryctr2,baryctr3,alpha2);
            if ((alpha2[0] < alpha_tol) || (alpha2[1] < alpha_tol)) 
            {
                compute_alphas(m2,baryctr0,baryctr1,baryctr3,alpha2);
                alpha2_elements[0] = E1;
                alpha2_elements[1] = E3;
                for (s_dim=0;s_dim < SYSTEM_DIM; s_dim++) 
                {
                    D_uh_m2[s_dim] = alpha2[0]*(diff_E1_E[s_dim]) + alpha2[1]*(diff_E3_E[s_dim]);
                }
                if ((alpha2[0] < alpha_tol) || (alpha2[1] < alpha_tol)) 
                {
                    printf("couldn't find non-negative alphas!\n");
                }
            }
            else 
            {
                alpha2_elements[0] = E2;
                alpha2_elements[1] = E3;
                for (s_dim=0;s_dim < SYSTEM_DIM; s_dim++) 
                {
                    D_uh_m2[s_dim] = alpha2[0]*(diff_E2_E[s_dim]) + alpha2[1]*(diff_E3_E[s_dim]);
                }
            }
        }
        else 
        {
            alpha2_elements[0] = E1;
            alpha2_elements[1] = E2;
            for (s_dim=0;s_dim < SYSTEM_DIM; s_dim++) 
            {
                D_uh_m2[s_dim] = alpha2[0]*(diff_E1_E[s_dim]) + alpha2[1]*(diff_E2_E[s_dim]);
            }
        }
        compute_alphas(m3,baryctr0,baryctr1,baryctr2,alpha3);
        if ((alpha3[0] < alpha_tol) || (alpha3[1] < alpha_tol)) 
        {
            compute_alphas(m3,baryctr0,baryctr2,baryctr3,alpha3);
            if ((alpha3[0] < alpha_tol) || (alpha3[1] < alpha_tol)) 
            {
                compute_alphas(m3,baryctr0,baryctr1,baryctr3,alpha3);
                alpha3_elements[0] = E1;
                alpha3_elements[1] = E3;
                for (s_dim=0;s_dim < SYSTEM_DIM; s_dim++) 
                {
                    D_uh_m3[s_dim] = alpha3[0]*(diff_E1_E[s_dim]) + alpha3[1]*(diff_E3_E[s_dim]);
                }
                if ((alpha3[0] < alpha_tol) || (alpha3[1] < alpha_tol)) 
                {
                    printf("couldn't find non-negative alphas!\n");
                }
            }
            else 
            {
                alpha3_elements[0] = E2;
                alpha3_elements[1] = E3;
                for (s_dim=0;s_dim < SYSTEM_DIM; s_dim++) 
                {
                    D_uh_m3[s_dim] = alpha3[0]*(diff_E2_E[s_dim]) + alpha3[1]*(diff_E3_E[s_dim]);
                }
            }
        }
        else 
        {
            alpha3_elements[0] = E1;
            alpha3_elements[1] = E2;
            for (s_dim=0;s_dim < SYSTEM_DIM; s_dim++) 
            {
                D_uh_m3[s_dim] = alpha3[0]*(diff_E1_E[s_dim]) + alpha3[1]*(diff_E2_E[s_dim]);
            }
        }

        // Calculate quantities needed for limiting
        for (s_dim=0;s_dim < SYSTEM_DIM; s_dim++) 
        { 
            diff_m1_E[s_dim] = uh_m1[s_dim] - uh_E[s_dim];
            diff_m2_E[s_dim] = uh_m2[s_dim] - uh_E[s_dim];
            diff_m3_E[s_dim] = uh_m3[s_dim] - uh_E[s_dim];

        }

        // Local solution vector
        get_approx_solution_lin(uh_pre_limit,E,mesh_element,mesh_node,(1.0/3.0),(1.0/3.0),uh_local);
        
        vec_b0_m1[0] = m1[0] - baryctr0[0];
        vec_b0_m1[1] = m1[1] - baryctr0[1];
        vec_b0_m2[0] = m2[0] - baryctr0[0];
        vec_b0_m2[1] = m2[1] - baryctr0[1];
        vec_b0_m3[0] = m3[0] - baryctr0[0];
        vec_b0_m3[1] = m3[1] - baryctr0[1];

        // N needs to be normal for the transformation
        normalize_vector(2,vec_b0_m1);
        normalize_vector(2,vec_b0_m2);
        normalize_vector(2,vec_b0_m3);
        if (( (*user_param).flag_trans == 1) && (flag_close_to_real_bound == 0)  && (fabs(uh_local[0] > EPSILON)))  
        {
            // Transformation to characteristic field
            transform_to_char_var(uh_local,diff_m1_E,vec_b0_m1,user_param,trans_diff_m1_E,TRANS_TO_CHAR); 
            transform_to_char_var(uh_local,diff_m2_E,vec_b0_m2,user_param,trans_diff_m2_E,TRANS_TO_CHAR); 
            transform_to_char_var(uh_local,diff_m3_E,vec_b0_m3,user_param,trans_diff_m3_E,TRANS_TO_CHAR); 

            //transf_to_char_var(uh_local,diff_m1_E,vec_b0_m1,user_param,trans_diff_m1_E); 
            //transf_to_char_var(uh_local,diff_m2_E,vec_b0_m2,user_param,trans_diff_m2_E); 
            //transf_to_char_var(uh_local,diff_m3_E,vec_b0_m3,user_param,trans_diff_m3_E); 

            transform_to_char_var(uh_local,D_uh_m1,vec_b0_m1,user_param,trans_D_uh_m1,TRANS_TO_CHAR); 
            transform_to_char_var(uh_local,D_uh_m2,vec_b0_m2,user_param,trans_D_uh_m2,TRANS_TO_CHAR); 
            transform_to_char_var(uh_local,D_uh_m3,vec_b0_m3,user_param,trans_D_uh_m3,TRANS_TO_CHAR); 
        }
        else 
        {
            for (s_dim=0;s_dim < SYSTEM_DIM; s_dim++) 
            { 
                trans_diff_m1_E[s_dim] = diff_m1_E[s_dim];
                trans_diff_m2_E[s_dim] = diff_m2_E[s_dim];
                trans_diff_m3_E[s_dim] = diff_m3_E[s_dim];
                
                trans_D_uh_m1[s_dim] = D_uh_m1[s_dim];
                trans_D_uh_m2[s_dim] = D_uh_m2[s_dim];
                trans_D_uh_m3[s_dim] = D_uh_m3[s_dim];
            }
        }

        // Perform actual limiting on each component of the (transformed) system
        for (s_dim=0;s_dim < SYSTEM_DIM; s_dim++)
        { 
            // Compute Deltas
            Delta1_trans[s_dim] = minmod(trans_diff_m1_E[s_dim],(*user_param).nu*trans_D_uh_m1[s_dim],(*user_param).mdx2);
            Delta2_trans[s_dim] = minmod(trans_diff_m2_E[s_dim],(*user_param).nu*trans_D_uh_m2[s_dim],(*user_param).mdx2);
            Delta3_trans[s_dim] = minmod(trans_diff_m3_E[s_dim],(*user_param).nu*trans_D_uh_m3[s_dim],(*user_param).mdx2);
        }



        if (( (*user_param).flag_trans == 1) && (flag_close_to_real_bound == 0) && (fabs(uh_local[0] > EPSILON)) )  
        {
            // transform back
            //transf_from_char_var(uh_local,Delta1_trans,vec_b0_m1,user_param,Delta1);
            //transf_from_char_var(uh_local,Delta2_trans,vec_b0_m2,user_param,Delta2);
            //transf_from_char_var(uh_local,Delta3_trans,vec_b0_m3,user_param,Delta3);

            transform_to_char_var(uh_local,Delta1_trans,vec_b0_m1,user_param,Delta1,TRANS_FROM_CHAR);
            transform_to_char_var(uh_local,Delta2_trans,vec_b0_m2,user_param,Delta2,TRANS_FROM_CHAR);
            transform_to_char_var(uh_local,Delta3_trans,vec_b0_m3,user_param,Delta3,TRANS_FROM_CHAR);	
        }
        else 
        {
            for (s_dim=0;s_dim < SYSTEM_DIM; s_dim++)
            { 
                Delta1[s_dim] = Delta1_trans[s_dim];
                Delta2[s_dim] = Delta2_trans[s_dim];
                Delta3[s_dim] = Delta3_trans[s_dim];
            }
        }

        // Perform actual limiting on each component of the transformed system
        for (s_dim=0;s_dim < SYSTEM_DIM; s_dim++)
        {
            // Limit only linear part, i.e. up to third coefficient...
            if (fabs(Delta1[s_dim]+Delta2[s_dim]+Delta3[s_dim]) > 1e-12) 
            {
                pos =  max(Delta1[s_dim],0.0);
                pos += max(Delta2[s_dim],0.0);
                pos += max(Delta3[s_dim],0.0);

                neg =  max(-Delta1[s_dim],0.0);
                neg += max(-Delta2[s_dim],0.0);
                neg += max(-Delta3[s_dim],0.0);

                sigma_p = min(1.0,neg/pos);
                sigma_m = min(1.0,pos/neg);

                Deltah1[s_dim] = sigma_p*max(0.0,Delta1[s_dim]) - sigma_m*max(0.0,-Delta1[s_dim]);
                Deltah2[s_dim] = sigma_p*max(0.0,Delta2[s_dim]) - sigma_m*max(0.0,-Delta2[s_dim]);
                Deltah3[s_dim] = sigma_p*max(0.0,Delta3[s_dim]) - sigma_m*max(0.0,-Delta3[s_dim]);

                uh_lin_limit[0] = uh_E[s_dim]+Deltah1[s_dim]-Deltah2[s_dim]+Deltah3[s_dim];
                uh_lin_limit[1] = 2.0*(Deltah2[s_dim]-Deltah3[s_dim]);
                uh_lin_limit[2] = 2.0*(Deltah2[s_dim]-Deltah1[s_dim]);
            }
            else 
            {
                uh_lin_limit[0] = uh_E[s_dim]+Delta1[s_dim]-Delta2[s_dim]+Delta3[s_dim];
                uh_lin_limit[1] = 2.0*(Delta2[s_dim]-Delta3[s_dim]);
                uh_lin_limit[2] = 2.0*(Delta2[s_dim]-Delta1[s_dim]);
            }

           // Check whether limiting changes the solution
            if ( (fabs(uh_lin_limit[0] - uh_pre_limit[SYSTEM_DIM*(E-1)*3 + s_dim*3  ]) > 1e-12)
              || (fabs(uh_lin_limit[1] - uh_pre_limit[SYSTEM_DIM*(E-1)*3 + s_dim*3+1]) > 1e-12)
              || (fabs(uh_lin_limit[2] - uh_pre_limit[SYSTEM_DIM*(E-1)*3 + s_dim*3+2]) > 1e-12)) 
            {
                
                // ... if it does, copy the linear part...
                uh_limit[SYSTEM_DIM*(E-1)*Nloc+s_dim*Nloc]     = uh_lin_limit[0];
                uh_limit[SYSTEM_DIM*(E-1)*Nloc+s_dim*Nloc + 1] = uh_lin_limit[1];
                uh_limit[SYSTEM_DIM*(E-1)*Nloc+s_dim*Nloc + 2] = uh_lin_limit[2];
                //if((s_dim==0) && ((E==5) || (E==17)))
                //{
                //    printf("E %d uh_limit %1.12lf %1.12lf %1.12lf \n",E,uh_lin_limit[0],uh_lin_limit[1],uh_lin_limit[2]);
               // }


                // ... and cut off the higher orders
                for(idofs=3;idofs<Nloc;idofs++) 
                {
                    uh_limit[SYSTEM_DIM*(E-1)*Nloc+s_dim*Nloc + idofs] = 0.0;
                }
            }
        }//limiting
    }//flag
}//shu_osher


