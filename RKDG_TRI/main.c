/* main.c
 * Driver routine for 2d hyperbolic solver
 * This code sets up the mesh data structure = edges,elements,nodes in 2d
 * Prince Chidyagwai
 * Philipp Monreal
*/
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>
#include<string.h>
#include"data_structure.h"
#include"refine.h"
#include"mesh_utils.h"
#include"hyperbolic_solver.h"
#include"output.h"
#include"compute_reference_mesh_error.h"
#include"dg_error.h"
#include"parameters.h"
#include"positivity_preserve.h"
int main(int argc, char **args)
{
    //mesh data structures	
    element *mesh_element;
    node *mesh_node;
    edge* mesh_edge;
    element *ref_mesh_element;
    node* ref_mesh_node;
    voxel* mesh_voxel;

    int num_elts=0;
    int ref_num_elts=0;
    int num_nodes=0;
    int ref_num_nodes=0;
    int i =0;
    int j =0;
    int result=0;
    int E=0;
    int s_dim_count=0;
    int size=0;

    // Read in mesh from given file	 
    FILE *meshfile = NULL;
    FILE *ref_meshfile_coefs = NULL;
    FILE *ref_meshfile = NULL;
    FILE *init_coefs_meshfile = NULL;
    char ref_filename[80];
    char mesh[80];
    char ref_data_str_mesh[80];
    char ref_sol_mesh[80];
    char ref_sol_data_str[80];
    char mesh_filename[80];
    char scenario_temp[80];
    char init_coefs_data_str[80];
    int global_edge_count=0;
    int active_elt_count=0;

    double* uh_tn_sol = NULL;
    double* uh_tn1_sol = NULL;
    double* ref_uh_tn_sol = NULL;
    double det =0.0;

    int Nloc = 0;
    double l1_error[SYSTEM_DIM];
    double linf_error[SYSTEM_DIM];

    double ref_l1_error[SYSTEM_DIM];
    double ref_linf_error[SYSTEM_DIM];

    double l1_e =0.0;
    double linf_e =0.0;
    double min_diam =10;
    double diam=0.0;

    double delta_t =0.0;
    int Ntsteps =0;
    int iedge=0;
    double iedge_length=0.0;
    double E_edge_length[3];

    double alpha = 0.0;
    int num_voxels=0;
    int added_elements=0;
    double min_edge_length =10.0;
    double max_edge_length =0.0;
    double min_ratio =0.0;
    double max_ratio =0.0;

    double saved_solution_time =0.0;

    // Time measurement
    double elapsed_t;
    double x0,y0,x1,y1;
    int N_grid =0;
    double mid_xc, mid_yc=0.0;
    double mid_point_vec[2];

    int found =0;
    int V =0;
    int V_found =0;
    int found_e =0;
    int e_in_v =0;
    int num_e_in_V=0;
    int E_v =0;
    int e_found =0;
    double local_uh_tn_sol[3];
    double ref_coords[2];

    double delta_x=0.0;
    double delta_y=0.0;


    FILE* coef_outfile;
    char out_filename[80];
    int grid_count=0;
    double norm_f_over_e =0.0;




    printf("=======================================================================\n");
    printf("== 2D RKDG code for hyperbolic equations on unstructured grids v 0.2 ==\n");
    printf("=======================================================================\n");

    // Parameter struct
    parameters *user_param=NULL;

    user_param = (parameters*)malloc(sizeof(parameters));

    // Check for available cores for parallelization
   
    #pragma omp parallel default(none) shared(user_param) 
    {
	(*user_param).num_cores = omp_get_num_threads();
    }
	
    //Reading parameters from inputfile
	cout<<"Init parameters"<<endl;
	ifstream* p_infile;
	string in;

	if(argc==1)
	{
	    ifstream infile("input/user_parameters.dat");
	    if(infile.good())
	    {	
		p_infile = &infile;
		init_parameters(user_param,p_infile);
	    }
	    else
	    {
		cout<<"Required inputfile "<<"\"user_parameters.dat\""<<"is not okay"<<endl;
		exit(1);			
	    }
	}
	else
	{	
	    char* filename=new char[50];
	    filename=args[1];
	    ifstream infile(filename);

	    if(infile.good())
	    {
		p_infile = &infile;
		init_parameters(user_param,p_infile);
	    }
	    else
	    {
		cout<<"Inputfile "<<"\""<<filename<<"\""<<"is not okay"<<endl;
		exit(1);
	    }
	}
	cout<<"--------------------------------------------------------------"<<endl;
	
     // Use set number of threads
    omp_set_num_threads((*user_param).num_cores);
  

    printf("Using %d CPU core(s) for parallel computing.\n",(*user_param).num_cores);


    //open file to read the mesh into element and node data structure
    if((*user_param).start_solution_file_flag)
    {
        init_coefs_meshfile = fopen((*user_param).startsolfilename,"r");
        if(init_coefs_meshfile == NULL)
        {
            fprintf(stderr,"In main.c :: Cannot open initial state file \n");
            exit(1);
        }
        //first line is reference data structure file
        result = fscanf(init_coefs_meshfile,"%s",init_coefs_data_str);
        if(!(result))
        {
            fprintf(stderr,"Wrong file format first line should be meshfilename \n");
            exit(1);
        }
        mesh_filename[0] = '\0';
        strcat(mesh_filename,"meshes/");
        strcat(mesh_filename,init_coefs_data_str);
        meshfile = fopen(mesh_filename,"r");
        if(meshfile == NULL)
        {
            fprintf(stderr,"Cannot open mesh data struture file \n");
            exit(1);
        }
        result = fscanf(init_coefs_meshfile,"%lf",&saved_solution_time);

        if(result ==0)
        {
            fprintf(stderr,"Error in file input, expecting stored time \n");
        }
        (*user_param).start_time = saved_solution_time;
    }
    else
    {
        printf("--------------------------------------------------------------\n");	
        meshfile = fopen((*user_param).meshfilename, "r");

        if(meshfile == NULL)
        {
            fprintf(stderr, "Can't Open Mesh File\n");
            exit(1);
        }
        printf("Using mesh: %s\n",(*user_param).meshfilename);

        if (!((*user_param).ref_sol_flag))
        {
	    printf("--------------------------------------------------------------\n");	
            printf("Storing solution coefficients\n");
        }
    }

    if((*user_param).ref_sol_flag)
    {
        ref_meshfile_coefs = fopen((*user_param).refsolfilename,"r");
        if(ref_meshfile_coefs == NULL)
        {
            fprintf(stderr,"Cannot open ref coefs mesh file \n");
            exit(1);
        }
        //first line is reference data structure file
        result = fscanf(ref_meshfile_coefs,"%s",ref_sol_data_str);
        if(!(result))
        {
            fprintf(stderr,"Wrong file format first line should be meshfilename \n");
            exit(1);
        }
        ref_data_str_mesh[0] = '\0';
        strcat(ref_data_str_mesh,"meshes/");
        strcat(ref_data_str_mesh,ref_sol_data_str);
        ref_meshfile = fopen(ref_data_str_mesh,"r");
        if(ref_meshfile == NULL)
        {
            fprintf(stderr,"Cannot open reference mesh data struture file \n");
            exit(1);
        }
    }

    Nloc = (((*user_param).p_degree+1)*((*user_param).p_degree+2))/2;

    //mesh data structures
    mesh_element = (element*)malloc(NUM_ELMTS_MAX*sizeof(element));
    mesh_node = (node*)malloc(NUM_NODES_MAX*sizeof(node));
    mesh_edge = (edge*)malloc(NUM_EDGES_MAX*sizeof(edge));

    mesh_voxel = (voxel*)malloc(NUM_VOXEL_MAX*sizeof(voxel)); 
    ref_mesh_element = (element*)malloc(NUM_ELMTS_MAX*sizeof(element));
    ref_mesh_node = (node*)malloc(NUM_NODES_MAX*sizeof(node));

    load_mesh_data_structure(meshfile,mesh_node,mesh_element,(*user_param).p_degree,&num_nodes,&num_elts,user_param);
    fclose(meshfile);

    if((*user_param).ref_sol_flag)
    {
        load_mesh_data_structure(ref_meshfile,ref_mesh_node,ref_mesh_element,(*user_param).p_degree,&ref_num_nodes,&ref_num_elts,user_param);
        fclose(ref_meshfile);
    }
    size = SYSTEM_DIM*num_elts*Nloc;

    uh_tn_sol = (double*)malloc(SYSTEM_DIM*(num_elts*Nloc)*sizeof(double));
    uh_tn1_sol = (double*)malloc(SYSTEM_DIM*(num_elts*Nloc)*sizeof(double));

    if((*user_param).ref_sol_flag)
    {
        ref_uh_tn_sol = (double*)malloc(SYSTEM_DIM*(ref_num_elts*Nloc)*sizeof(double));
        fscanf(ref_meshfile_coefs,"%lf",&saved_solution_time);
        if(!(saved_solution_time ==(*user_param).Ftime))
        {
            printf("WARNGING:::reference solution saved time does not match (*user_param).Ftime \n");
        }

        load_sol_coeffs(ref_uh_tn_sol,ref_mesh_element,ref_mesh_node,ref_meshfile_coefs,Nloc,ref_num_elts);
    }

    if((*user_param).ref_sol_flag)
    {
	printf("Read reference mesh file with %d elements and %d nodes \n",ref_num_elts,ref_num_nodes);
    }

    init_mesh_data_structure(mesh_element,mesh_node,mesh_edge,num_elts,num_nodes,&global_edge_count,&active_elt_count);
    printf("--------------------------------------------------------------\n");	
    printf("Initialised %d mesh edges\n",global_edge_count);
    printf("Current active element count is %d \n",active_elt_count);
    printf("Current node count is %d \n",num_nodes);
    printf("--------------------------------------------------------------\n");
    active_elt_count=0;


    if((*user_param).symmetric_solution)
    {
        printf("Checking symmetry in edge alignments\n");
        //check_symmetry_edges(mesh_element,mesh_node,mesh_edge,num_elts);
    }

    //label reflective boundary
    if((*user_param).reflective_bc)
    {
        label_reflective_boundary(global_edge_count,mesh_edge,mesh_node,user_param);
    }

    for(E=1;E<(num_elts+1);E++)
    {
        init_zero_d(E_edge_length,3);
        det = mesh_element[E].det;
        for(i=1;i<4;i++)
        {
            iedge = mesh_element[E].edge[i];
            iedge_length = calc_length_face(mesh_edge,iedge,mesh_node);
            mesh_edge[iedge].length = iedge_length;

            if(iedge_length < min_edge_length)
            {
                min_edge_length = iedge_length;
            }
            if(iedge_length > max_edge_length)
            {
                max_edge_length = iedge_length;
            }
            E_edge_length[i-1] = iedge_length;
        }
        diam = (E_edge_length[0]*E_edge_length[1]*E_edge_length[2])/(2.0*det);
        if(diam < min_diam)
        {
            min_diam=diam;
        }
    }

    if((*user_param).periodic)
    {
        make_boundary_periodic(mesh_edge,mesh_node,mesh_element,global_edge_count,user_param);
    }

    if(((*user_param).ref_sol_flag) || ((*user_param).store_grid_solution))
    {
        if((*user_param).test_case ==5)
        {
            x0 = 0.0;
            y0 = 0.0;

            x1 = 7.0;
            y1 = 7.0;
        }
        else if((*user_param).test_case == 6)
        {
            x0 = -0.5;
            y0 = -0.5;

            x1 = 0.5;
            y1 = 0.5;
        }
        else if((*user_param).test_case == 24)
        {
            x0 = -5.0;
            y0 = -5.0;

            x1 = 5.0;
            y1 = 5.0;
        }
        else if((*user_param).test_case == 25)
        {
            x0 = -10.0;
            y0 = -10.0;

            x1 = 10.0;
            y1 = 10.0;
        }
        else if((*user_param).test_case == 26)
        {
            x0 = 0.0;
            y0 = 0.0;

            x1 = 12.0;
            y1 = 6.0;
        }
        else if ((*user_param).test_case ==27)
        {
            x0 = 0.0;
            y0 = 0.0;
            
            x1 = 7.0;
            y1 = 7.0;
        }
        else
        {
            x0 = 0.0;
            y0 = 0.0;

            x1 = 1.0;
            y1 = 1.0;
        }
    }

    (*user_param).mesh_min_diam = min_diam;
    (*user_param).mesh_max_edge_length = max_edge_length;
    (*user_param).mesh_min_edge_length = min_edge_length;
    printf("Minimum triangle diameter = %lf \n",min_diam);
    printf("Minimum edge length = %lf \n",min_edge_length);
    printf("Maximum edge length = %lf \n",max_edge_length);
    printf("--------------------------------------------------------------\n");

    alpha = 5.0*max_edge_length;
    if((*user_param).ref_sol_flag)
    {
        create_voxel(alpha,ref_mesh_element,ref_mesh_node,x0,y0,x1,y1,mesh_voxel,&num_voxels);
        printf("Voxel structure with %d voxels created \n",num_voxels);

        add_elements_to_voxel(ref_num_elts,num_voxels,mesh_voxel,ref_mesh_element,ref_mesh_node,&added_elements);
        printf("%d Elements have been assigned to voxel structure \n",added_elements);
    }
    printf("Initialize Gauss points in all triangles... \n");
    init_gpts(mesh_element,mesh_node,mesh_edge,num_elts,global_edge_count);
    printf("done.\n");
    if((*user_param).symmetric_solution)
    {
        //check_symmetry_of_gpts(mesh_element,mesh_node,mesh_edge,num_elts,global_edge_count);
    }
    
    init_limiting_gpts(num_elts,global_edge_count,mesh_element,mesh_edge,mesh_node);
    printf("Initialize positivity preserving points \n");

    printf("Initialize density values on all Gauss points... ");
    init_density(mesh_element,mesh_node,mesh_edge,user_param,num_elts);
    printf("done.\n");

    printf("--------------------------------------------------------------\n");
    printf("Minimal density in the medium = %lf\n",(*user_param).min_density);

    // Calculate timestep
    delta_t = (*user_param).CFL*(min_edge_length)*(*user_param).min_density;
    /*delta_t = 0.01/4.0;*/
    /*delta_t = 0.02/4.0;*/
    Ntsteps = ceil(fabs((*user_param).Ftime- (*user_param).start_time)/delta_t);

    printf("--------------------------------------------------------------\n");
    printf("Running test-case ");
    switch ((*user_param).test_case) 
    {
	case(1): {
	    printf("1: Constant initial and boundary conditions\n");
	    break;
	}
	case(2): {
	    printf("2: Initial Gaussian\n");
	    break;
	}
	case(3): {
	    printf("3: Initial plane source\n");
	    break;
	}
	case(4): {
	    printf("Manufactured solutions: P1 closure; constant along along x+y = c\n");
	    break;
	}
	case(5): {
	    printf("5: Checkerboard\n");
	    break;
	}
	case(6): {
	    printf("6: Linesource\n");
	    break;
	}
	case(7): {
	    printf("Manufactured solutions: P1 closure; radially symmetric\n");
	    break;
	}
	case(8): {
	    printf("8: Manufactured solutions: radially symmetric with time- and space-dependent flux\n");
	    break;
	}
	case(9): {
	    printf("9: 3 beams space-dependent flux using CT-data\n");
	    break;
	}
	case(10): {
	    printf("10: Single beam\n");
	    break;
	}
	case(11): {
	    printf("11 : Quasi 1D Riemann-problem\n");
	    break;
	}
	case(12): {
	    printf("12 : Manufactured solutions: M1 closure; constant along x+2y=c\n");
	    break;
	}
	case(13): {
	    printf("13 : Two initial Gaussians\n");
	    break;
	}
	case(14): {
	    printf("14 : Manufactured solutions: P1 closure; constant along x+2y=c\n");
	    break;
	}
	case(15): {
	    printf("15 : Quasi 1D test-case with discontinuous coefficients\n");
	    break;
	}
	case(16): {
	    printf("16 : Manufactured solutions: K1 closure; constant along x+2y=c\n");
	    break;
	}
	case(17): {
	    printf("17 : Reflective arc\n");
	    break;
	}
	case(18): {
	    printf("18 : Manufactured solutions: P1 closure on checkerboard geometry\n");
	    break;
	}
	case(19): {
	    printf("19 : Manufactured solutions: P1 closure; discont. parameters in a box\n");
	    break;
	}
	case(20): {
	    printf("20 : Smooth checkerboard\n");
	    break;
	}
	case(21): {
	    printf("21 : Reed's problem\n");
	    break;
	}
	case(22): {
	    printf("22 : Beams enter vacuum ring\n");
	    break;
	}
	case(23): {
	    printf("23 : Initial atan\n");
	    break;
	}
    }

    printf("--------------------------------------------------------------\n");	
    printf("Order of the polynomial degree approximation is: %d \n",(*user_param).p_degree);

    printf("--------------------------------------------------------------\n");
    printf("Using closure ");
    switch ((*user_param).closure) 
    {
	case(1): {
	    printf("1: P1\n");
	    break;
	}
	case(2): {
	    printf("2: K1\n");
	    break;
	}
	case(3): {
	    printf("3: M1\n");
	    break;
	}
	case(4): {
	    printf("4: P2\n");
	    break;
	}
	case(5): {
	    printf("5: K2\n");
	    break;
	}
	case(6): {
	    printf("6: P3\n");
	    break;
	}
	case(7): {
	    printf("7: Simple nonlinear flux function\n");
	    break;
	}
    }
    printf("--------------------------------------------------------------\n");
    printf("Time interval is [%lf,%lf]\n",(*user_param).start_time,(*user_param).Ftime);
    printf("Timestep is dt = %lf \n",delta_t);
    printf("Number of timesteps is %d \n",Ntsteps);
    printf("--------------------------------------------------------------\n");
    printf("Number of intermediate outputs (=Stages) is %d \n",(*user_param).no_out);
    printf("--------------------------------------------------------------\n");
    printf("Computing...\n");


    // =============================================
    // pass initial vector into time stepping scheme
    // Actual computation
    elapsed_t = omp_get_wtime();
    rk_dg_time_step(mesh_element,mesh_node,mesh_edge,init_coefs_meshfile,(*user_param).start_time,Ntsteps,Nloc,num_nodes,num_elts,delta_t,user_param,uh_tn1_sol);

    elapsed_t = omp_get_wtime()-elapsed_t;

    if(!((*user_param).ref_sol_flag))
    {
        // calculate error
        for(s_dim_count=0;s_dim_count<SYSTEM_DIM;s_dim_count++)
        {
            l1_dg_error(num_elts,mesh_element,mesh_node,user_param,uh_tn1_sol,(*user_param).Ftime,s_dim_count,Nloc,&l1_e,&linf_e);
            l1_error[s_dim_count]=l1_e;
            linf_error[s_dim_count]=linf_e;
        }
        for(s_dim_count=0;s_dim_count<SYSTEM_DIM;s_dim_count++)
        {
            printf("l1_error moment %d == %10.7e \n",s_dim_count,l1_error[s_dim_count]);
            printf("linf_error moment %d == %10.7e \n",s_dim_count,linf_error[s_dim_count]);
        }
    }

    // Write solution
    save_solution_coefficients(uh_tn1_sol,mesh_element,mesh_node,user_param,Nloc,num_elts,(*user_param).Ftime,"_end");
    plot_realizability(uh_tn1_sol,mesh_element,mesh_node,user_param,Nloc,num_elts,0.0,"_end");

    // output solution along line specified in parameter-struct
    if((*user_param).flag_store_solution_line) 
    {
	    ctrl_output_line(uh_tn1_sol,mesh_node,mesh_element,user_param,num_elts);
    }

    // calculate relative error
    if((*user_param).ref_sol_flag)
    {
        printf("Computing Relative error \n");
        for(s_dim_count=0;s_dim_count<SYSTEM_DIM;s_dim_count++)
        {
            compute_relative_error(uh_tn1_sol,ref_uh_tn_sol,mesh_node,ref_mesh_node,mesh_element,ref_mesh_element,mesh_voxel,num_elts,num_voxels,s_dim_count,&l1_e,&linf_e);
            ref_l1_error[s_dim_count]=l1_e;
            ref_linf_error[s_dim_count]=linf_e;
        }	

        for(s_dim_count=0;s_dim_count<SYSTEM_DIM;s_dim_count++)
        {
            printf("ref_l1_error moment %d == %10.7e \n",s_dim_count,ref_l1_error[s_dim_count]);
            printf("ref_linf_error moment %d == %10.7e \n",s_dim_count,ref_linf_error[s_dim_count]);
        }
        //write ref_error outputfile 
        save_ref_error(uh_tn1_sol,ref_uh_tn_sol,mesh_node,ref_mesh_node,mesh_element,ref_mesh_element,mesh_voxel,num_nodes,num_elts,ref_num_elts,num_voxels,user_param);
    }

    store_solution(mesh_element,mesh_node,uh_tn1_sol,num_elts,num_nodes,1,"_end",user_param);
    check_dg_sol_realizability(mesh_element,num_elts,mesh_node,user_param,uh_tn1_sol,(*user_param).Ftime,&min_ratio,&max_ratio);
    printf("norm_psi_one/psi_zero max %10.16e min %10.16e \n",max_ratio,min_ratio);
    
    if((*user_param).store_grid_solution)
    {
        create_voxel(alpha,mesh_element,mesh_node,x0,y0,x1,y1,mesh_voxel,&num_voxels);
        printf("Voxel structure with %d voxels created \n",num_voxels);

        add_elements_to_voxel(num_elts,num_voxels,mesh_voxel,mesh_element,mesh_node,&added_elements);
        printf("%d Elements have been assigned to voxel structure \n",added_elements);

        out_filename[0] = '\0';
        strcat(out_filename,"matlab_unif_grid/");
        sprintf(scenario_temp, "s%d", (*user_param).test_case);
        strcat(out_filename,scenario_temp);
        strcat(out_filename,(*user_param).outputfilename);

        coef_outfile = fopen(out_filename,"w");

        N_grid = (*user_param).store_grid_solution;

        delta_x = fabs(x1-x0)/N_grid;
        delta_y = fabs(y1-y0)/N_grid;
        printf("N_grid %d \n",N_grid);


        //store final solution on a uniform grid.
        for(j=0;j<N_grid;j++)
        {
            for(i=0;i<N_grid;i++)
            {
                mid_xc = 0.5*((x0 + i*delta_x) + (x0+(i+1)*delta_x));
                mid_yc = 0.5*((y0 + j*delta_y) + (y0+(j+1)*delta_y));
                mid_point_vec[0] = mid_xc;
                mid_point_vec[1] = mid_yc;


                V=1;
                E_v=0;

                while((!found) && (V < num_voxels+1))
                {
                    found = check_point_in_voxel(mid_point_vec,mesh_voxel,V);
                    if(found)
                    {
                        V_found = V;
                    }
                    V++;
                }
                if((V == num_voxels) &&(found ==0))
                {
                    fprintf(stderr,"Error point outside of voxel domain \n");
                    exit(1);
                }

                num_e_in_V = mesh_voxel[V_found].num_elements;
                while(!(found_e) && (E_v < num_e_in_V))
                {
                    e_in_v = mesh_voxel[V_found].element_list[E_v];
                    found_e = check_point_in_element(e_in_v,mesh_element,mesh_node,mid_point_vec);
                    if(found_e)
                    {
                        e_found = e_in_v;
                    }
                    E_v++;
                }
                if(found_e == 0) 
                {
                    fprintf(stderr,"Error in voxel element list :: Cannot find element \n");
                    exit(1);
                }

                init_zero_d(local_uh_tn_sol,SYSTEM_DIM);
                map_to_reference_element(mesh_element,mesh_node,e_found,ref_coords,mid_xc,mid_yc);
                get_approx_solution(uh_tn1_sol,e_found,mesh_element,mesh_node,ref_coords[0],ref_coords[1],local_uh_tn_sol);
                norm_f_over_e = sqrt(pow(local_uh_tn_sol[1],2.0) + pow(local_uh_tn_sol[2],2.0))/local_uh_tn_sol[0];

                fprintf(coef_outfile,"%lf %lf %.16lf %.16lf %.16lf\n",mid_xc,mid_yc,local_uh_tn_sol[0],local_uh_tn_sol[1],local_uh_tn_sol[2]);
                grid_count++;

                V_found =0;
                found =0;
                found_e =0;
                
            
            }
        }
        fclose(coef_outfile);
    }//store on grid
    printf("grid_count = %d \n",grid_count);

    //elapsed_t = omp_get_wtime()-elapsed_t;
    printf("Elapsed Time = %lf s \n",elapsed_t);


    // free variables
    delete[] (*user_param).meshfilename;
    delete[] (*user_param).outputfilename;
    delete[] (*user_param).refsolfilename;
    delete[] (*user_param).startsolfilename;
    delete[] (*user_param).densityfilename;
    free(user_param);
    free(mesh_voxel);
    free(ref_mesh_element);
    free(ref_mesh_node);
    free(uh_tn_sol);
    free(uh_tn1_sol);
    free(mesh_element);
    free(mesh_node);
    free(mesh_edge);
    free(ref_meshfile_coefs);
    /*free(ref_meshfile);*/
    free(init_coefs_meshfile);
    free(ref_uh_tn_sol);

    return(0);
}
