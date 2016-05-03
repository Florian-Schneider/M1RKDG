/*
Parameters
Philipp Monreal
*/
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include"parameters.h"

/* void init_parameters
 * Initialize parameter structure
 * out: parameter structure -- user_param
 */
void init_parameters(parameters* user_param)
{
    /*Closures:*/
	/*1: P1*/
	/*2: K1*/
	/*3: M1*/
	/*4: P2*/
	/*5: K2*/
	/*6: P3*/
	/*7: Simple nonlinear flux function*/
    (*user_param).closure = 3;

    /*Test cases:*/
	/*1 : Constant initial and boundary conditions*/
	/*2 : Initial Gaussian */
	/*3 : Initial plane source*/
	/*4 : Manufactured solutions: P1 closure; constant along along x+y = c*/
	/*5 : Checkerboard*/
	/*6 : Linesource*/
	/*7 : Manufactured solutions: P1 closure; radially symmetric*/
	/*8 : Manufactured solutions: radially symmetric with time- and space-dependent flux*/
	/*9 : 3 beams using CT-data*/
	/*10 : 1 beam*/
	/*11 : Quasi 1D Riemann-problem*/
	/*12 : Manufactured solutions: M1 closure; constant along along x+2y = c*/
	/*13 : Two initial Gaussians*/
	/*14 : Manufactured solutions: P1 closure; constant along along x+2y = c*/
	/*15 : Quasi 1D test-case with discontinuous coefficients */
	/*16 : Manufactured solutions: K1 closure; constant along along x+2y = c*/
	/*17 : Reflective arc*/
	/*18 : Manufactured solutions: P1 closure on checkerboard geometry*/
	/*19 : Manufactured solutions: P1 closure; discont. parameters in a box */
	/*20 : Smooth checkerboard*/
	/*21 : Reed's problem*/
    //24 : Homogeneous circle test
    //25 : Flash test
    //26 : Shadow test
    //27 : 3 Beam 
    (*user_param).test_case = 24;

    //rectangular mesh
    /*(*user_param).meshfilename = "meshes/rect1";*/
    /*(*user_param).outputfilename = "rect1";*/
    /*(*user_param).meshfilename = "meshes/rect1_refine";*/
    /*(*user_param).outputfilename = "rect1_refine";*/
    /*(*user_param).meshfilename = "meshes/rect1_refine_refine";*/
    /*(*user_param).outputfilename = "rect1_refine_refine";*/

    /*(*user_param).meshfilename = "meshes/mesh32_2";*/
    /*(*user_param).outputfilename = "mesh32_2";*/
    /*(*user_param).meshfilename = "meshes/mesh128_2";*/
    /*(*user_param).outputfilename = "mesh128_2";*/
    /*(*user_param).meshfilename = "meshes/mesh512_2";*/
    /*(*user_param).outputfilename = "mesh512_2";*/
    /*(*user_param).meshfilename = "meshes/mesh2048_2";*/
    /*(*user_param).outputfilename = "mesh2048_2";*/
    //(*user_param).meshfilename = "meshes/mesh8192_2";
    //(*user_param).outputfilename = "mesh8192_2";
    (*user_param).meshfilename = "meshes/mesh32768_2";
    (*user_param).outputfilename = "mesh32768_2";

    /*(*user_param).meshfilename = "meshes/unif00025";*/
    /*(*user_param).outputfilename = "unif00025";*/
    /*(*user_param).meshfilename = "meshes/unif0005";*/
    /*(*user_param).outputfilename = "unif0005";*/
    /*(*user_param).meshfilename = "meshes/unif0.01";*/
    /*(*user_param).outputfilename = "unif0.01";*/
    //(*user_param).meshfilename = "meshes/unif0.02";
    //(*user_param).outputfilename = "unif0.02";
    /*(*user_param).meshfilename = "meshes/unif0.05";*/
    /*(*user_param).outputfilename = "unif0.05";*/
    /*(*user_param).meshfilename = "meshes/unif0.10";*/
    /*(*user_param).outputfilename = "unif0.10";*/
    /*(*user_param).meshfilename = "meshes/unif0.25";*/
    /*(*user_param).outputfilename = "unif0.25";*/
    /*(*user_param).meshfilename = "meshes/unif0.50";*/
    /*(*user_param).outputfilename = "unif0.50";*/

    /*(*user_param).meshfilename = "meshes/cunif0.01";*/
    /*(*user_param).outputfilename = "cunif0.01";*/
    /*(*user_param).meshfilename = "meshes/cunif0.015625";*/
    /*(*user_param).outputfilename = "cunif0.015625";*/
    /*(*user_param).meshfilename = "meshes/cunif0.0125";*/
    /*(*user_param).outputfilename = "cunif0.0125";*/
    /*(*user_param).meshfilename = "meshes/cunif0.02";*/
    /*(*user_param).outputfilename = "cunif0.02";*/
    /*(*user_param).meshfilename = "meshes/cunif0.025";*/
    /*(*user_param).outputfilename = "cunif0.025";*/
    /*(*user_param).meshfilename = "meshes/cunif0.05";*/
    /*(*user_param).outputfilename = "cunif0.05";*/
    /*(*user_param).meshfilename = "meshes/cunif0.10";*/
    /*(*user_param).outputfilename = "cunif0.10";*/
    /*(*user_param).meshfilename = "meshes/cunif0.15";*/
    /*(*user_param).outputfilename = "cunif0.15";*/
    /*(*user_param).meshfilename = "meshes/cunif0.20";*/
    /*(*user_param).outputfilename = "cunif0.20";*/
    /*(*user_param).meshfilename = "meshes/cunif0.25";*/
    /*(*user_param).outputfilename = "cunif0.25";*/
    /*(*user_param).meshfilename = "meshes/cunif0.40";*/
    /*(*user_param).outputfilename = "cunif0.40";*/
    /*(*user_param).meshfilename = "meshes/cunif0.50";*/
    /*(*user_param).outputfilename = "cunif0.50";*/
    /*(*user_param).meshfilename = "meshes/cunif0.80";*/
    /*(*user_param).outputfilename = "cunif0.80";*/

    /*(*user_param).meshfilename = "meshes/chkr0.00625";*/
    /*(*user_param).outputfilename = "chkr0.00625";*/
    /*(*user_param).meshfilename = "meshes/chkr0.01";*/
    /*(*user_param).outputfilename = "chkr0.01";*/
    /*(*user_param).meshfilename = "meshes/chkr0.015625";*/
    /*(*user_param).outputfilename = "chkr0.015625";*/
    /*(*user_param).meshfilename = "meshes/chkr0.0125";*/
    /*(*user_param).outputfilename = "chkr0.0125";*/
    /*(*user_param).meshfilename = "meshes/chkr0.02";*/
    /*(*user_param).outputfilename = "chkr0.02";*/
    /*(*user_param).meshfilename = "meshes/chkr0.025";*/
    /*(*user_param).outputfilename = "chkr0.025";*/
    /*(*user_param).meshfilename = "meshes/chkr0.05";*/
    /*(*user_param).outputfilename = "chkr0.05";*/
    /*(*user_param).meshfilename = "meshes/chkr0.10";*/
    /*(*user_param).outputfilename = "chkr0.10";*/
    /*(*user_param).meshfilename = "meshes/chkr0.20";*/
    /*(*user_param).outputfilename = "chkr0.20";*/
    /*(*user_param).meshfilename = "meshes/chkr0.25";*/
    /*(*user_param).outputfilename = "chkr0.25";*/
    /*(*user_param).meshfilename = "meshes/chkr0.40";*/
    /*(*user_param).outputfilename = "chkr0.40";*/
    /*(*user_param).meshfilename = "meshes/chkr0.50";*/
    /*(*user_param).outputfilename = "chkr0.50";*/

    /*(*user_param).meshfilename = "meshes/chkr_unif1.msh";*/
    /*(*user_param).outputfilename = "chkr_unif1.msh";*/
    /*(*user_param).meshfilename = "meshes/chkr_unif2.msh";*/
    /*(*user_param).outputfilename = "chkr_unif2.msh";*/
    /*(*user_param).meshfilename = "meshes/chkr_unif3.msh";*/
    /*(*user_param).outputfilename = "chkr_unif3.msh";*/
    /*(*user_param).meshfilename = "meshes/chkr_unif4.msh";*/
    /*(*user_param).outputfilename = "chkr_unif4.msh";*/
    /*(*user_param).meshfilename = "meshes/chkr_unif5.msh";*/
    /*(*user_param).outputfilename = "chkr_unif5.msh";*/
    /*(*user_param).meshfilename = "meshes/chkr_unif6.msh";*/
    /*(*user_param).outputfilename = "chkr_unif6.msh";*/

    /*(*user_param).meshfilename = "meshes/mesh_earfemale_0p08000";*/
    /*(*user_param).outputfilename = "mesh_earfemale_0p08000";*/

    /*(*user_param).meshfilename = "meshes/mesh_earfemale_0p18000";*/
    /*(*user_param).outputfilename = "mesh_earfemale_0p18000";*/

    /*(*user_param).meshfilename = "meshes/reflect37";*/
    /*(*user_param).outputfilename = "reflect37";*/
    /*(*user_param).meshfilename = "meshes/reflect128";*/
    /*(*user_param).outputfilename = "reflect128";*/
    /*(*user_param).meshfilename = "meshes/reflect472";*/
    /*(*user_param).outputfilename = "reflect472";*/
    /*(*user_param).meshfilename = "meshes/reflect1808";*/
    /*(*user_param).outputfilename = "reflect1808";*/
    /*(*user_param).meshfilename = "meshes/reflect7072";*/
    /*(*user_param).outputfilename = "reflect7072";*/
    /*(*user_param).meshfilename = "meshes/reflect27968";*/
    /*(*user_param).outputfilename = "reflect27968";*/
    /*(*user_param).meshfilename = "meshes/reflect294912";*/
    /*(*user_param).outputfilename = "reflect294912";*/

    (*user_param).start_time = 0.0;
    /*(*user_param).start_time = 21.0;*/
    /*(*user_param).start_time = 16.0;*/
    /*(*user_param).start_time = 4.8;*/

    /*(*user_param).Ftime = 30.0;*/
    /*(*user_param).Ftime = 20.0;*/
    /*(*user_param).Ftime = 16.;*/
    /*(*user_param).Ftime = 10.0;*/
    /*(*user_param).Ftime = 3.2;*/
    /*(*user_param).Ftime = 2.8;*/
    /*(*user_param).Ftime = 2.8;*/
    /*(*user_param).Ftime = 1.0;*/
    (*user_param).Ftime = 0.55;
    /*(*user_param).Ftime = 0.3;*/
    /*(*user_param).Ftime = 0.1;*/
    /*(*user_param).Ftime = 0.0;*/

    // backwards time-stepping 
    if ( (*user_param).start_time > (*user_param).Ftime) 
    {
	(*user_param).flag_backwards_t_step = 1;
    }
    else 
    {
	(*user_param).flag_backwards_t_step = 0;
    }

    // store intermediate solutions
    (*user_param).flag_store_solution = 1;

    // number of intermediate solutions that are stored
    (*user_param).no_out = 20;

    // Limiting
    (*user_param).flag_limiting = 0;

    // Transformation into the charactersitic fields during limiting
    (*user_param).flag_trans = 0;

    // nu and mdx2 determine how strict the limiting is; 
    // 1 (strict) <= nu <= 2 (less strict)
    // mdx2=M*(dx^2); where M is an upper bound of the absolute of value 
    // of the second-order derivative of the initial conditions
    // high values of mdx2 effectively switch of limiting
    (*user_param).nu = 1.5;

    switch ((*user_param).test_case) {
	case(3): {
	    (*user_param).mdx2 = 1e-4;
	    break;
	}
	case(5): {
	    (*user_param).mdx2 = 1e-10;
	    /*(*user_param).nu = 1.0;*/
	    break;
	}
	case(6): {
	    (*user_param).mdx2 = 50.0;
	    /*(*user_param).nu = 1.0;*/
	    break;
	}
	case(9): {
	    (*user_param).mdx2 = 1e-16;
	    break;
	}
	case(10): {
	    (*user_param).mdx2 = 1e-4;
	    /*(*user_param).nu = 1.0;*/
	    break;
	}
	case(11): {
	    (*user_param).mdx2 = 1e-18;
	    (*user_param).nu = 1.0;
	    break;
	}
	case(22): {
	    (*user_param).mdx2 = 1e-14;
	    (*user_param).nu = 1.0;
	    break;
	}
	default: {
	    (*user_param).mdx2 = 1e-3;
	    break;
	}
	break;
    }

    (*user_param).ref_sol_flag = 0;

    //specify the reference solution file
    if((*user_param).ref_sol_flag)
    {
	(*user_param).refsolfilename = "solution/s20chkr_unif6.msh_end";
	/*(*user_param).refsolfilename = "solution/s20cunif0.025_end";*/
    }

    // Polynomial degree of approximation
    (*user_param).p_degree = 2;

    //specify density file
    switch ((*user_param).test_case) 
    {
        case(9): 
            {
                (*user_param).densityfilename = "raddata_earfemale_density.dat";
                break;
            }
        default: 
            {
	        // !!! Does not allow for nodes with negative coordinates !!!
                (*user_param).densityfilename = "unit_density.dat";
		/*(*user_param).densityfilename = "density.dat";*/
                break;
            }
            break;
    }

    /*specify if the computation starts with a precomputed solution*/
    (*user_param).start_solution_file_flag = 0;

    if((*user_param).start_solution_file_flag)
    {
	(*user_param).startsolfilename = "solution/s7mesh128_2_0.5";
    }

    // produce 1D cut of the solution?
    (*user_param).custom_output_flag = 0;

    if((*user_param).custom_output_flag) 
    {
        // define start (a) and end (b) point of custom output vector
        (*user_param).custom_out_a[0] = 0.0;
        (*user_param).custom_out_a[1] = 0.50;

        (*user_param).custom_out_b[0] = 1.0;
	(*user_param).custom_out_b[1] = 0.50;

        // Number of points to along the line a to b
        (*user_param).custom_out_length = 10000;
    }

    if( ((*user_param).test_case ==5) || ((*user_param).test_case == 20))
    {
        (*user_param).voxel_domain_size[0] = 7.0;
        (*user_param).voxel_domain_size[1] = 7.0;
	(*user_param).periodic_domain_size_x = 7.0;
	(*user_param).periodic_domain_size_y = 7.0;
    }
    else if ( strcmp((*user_param).meshfilename,"meshes/reflective_arc") == 0 )
    {
        (*user_param).voxel_domain_size[0] = 11.0;
        (*user_param).voxel_domain_size[1] = 12.0;
    }
    else
    {
        (*user_param).voxel_domain_size[0] = 1.0;
        (*user_param).voxel_domain_size[1] = 1.0;
	(*user_param).periodic_domain_size_x = 1.0;
	(*user_param).periodic_domain_size_y = 1.0;
    }

    //Periodic boundary conditions
    (*user_param).periodic = 0;

    //reflective boundary conditions
    (*user_param).reflective_bc = 0;

    // Enforce reflective b.c. for reflective arc test-case
    if((*user_param).test_case == 17)
    {
	(*user_param).reflective_bc = 1;
    }
(*user_param).CFL=0.1;
(*user_param).num_cores=24;
}

/* void init_parameters
*  Initialize parameter structure from a inputfile, checks if all required parameters were set
*  and stores parameter configuration in text file
*  in: infilestream -- infile
*  out: parameter structure -- user_param
*/
void init_parameters(parameters* user_param, std::ifstream* infile)
{
using namespace std;
    //variables to save parameters from inputfile
	string variable;
	double double_value;
	int int_value;
	string string_value;
	char* meshfilename=new char[50];
	char* outputfilename=new char[50];
	char* refsolfilename=new char[50];
	char* densityfilename=new char[50];
	char* startsolfilename=new char[50];
	//array to check if required parameters were set
	bool  required_parameters[29];
	// ofstream to save parameters
	char filename[80];
	filename[0] = '\0';
	strcat(filename,"solution/");
	strcat(filename,"used_parameters.dat");

	ofstream parameters_out(filename);
	
    	for(int i=1;i<18;i++)
	{
	    required_parameters[i]=false;
	}
	//set default values for some parameters
	(*user_param).densityfilename="unit_density.dat";
	while(!(*infile).eof())
	{	
	    (*infile)>>variable;
	    if(variable=="closure")
	    {
		(*infile)>>int_value;
		(*user_param).closure=int_value;
		required_parameters[1]=true;
		parameters_out<<variable<<": "<<int_value<<endl;
	    }
	    else if(variable=="test_case")
	    {
		(*infile)>>int_value;
		(*user_param).test_case=int_value;
		required_parameters[2]=true;
		parameters_out<<variable<<": "<<int_value<<endl;
	    }
	    else if(variable=="meshfilename")
	    {
		(*infile)>>string_value;
		strcpy(meshfilename,string_value.c_str());
		(*user_param).meshfilename=meshfilename;
		required_parameters[3]=true;
		parameters_out<<variable<<": "<<string_value<<endl;
	    }
	    else if(variable=="outpufilename")
	    {
		(*infile)>>string_value;
		strcpy(outputfilename,string_value.c_str());
		(*user_param).outputfilename=outputfilename;
		required_parameters[4]=true;
		parameters_out<<variable<<": "<<string_value<<endl;
	    }
	    else if(variable=="start_time")
	    {
		(*infile)>>double_value;
		(*user_param).start_time=double_value;
		required_parameters[5]=true;
		parameters_out<<variable<<": "<<double_value<<endl;
	    }
	    else if(variable=="Ftime")
	    {
		(*infile)>>double_value;
		(*user_param).Ftime=double_value;
		required_parameters[6]=true;
		parameters_out<<variable<<": "<<double_value<<endl;
	    }
	    else if(variable=="flag_backwards_t_step")
	    {
		(*infile)>>int_value;
		(*user_param).flag_backwards_t_step=int_value;
		parameters_out<<variable<<": "<<int_value<<endl;
	    }
	    else if(variable=="flag_store_solution")
	    {
		(*infile)>>int_value;
		(*user_param).flag_store_solution=int_value;
		required_parameters[7]=true;
		parameters_out<<variable<<": "<<int_value<<endl;
	    }
	    else if(variable=="flag_store_solution_line")
	    {
		(*infile)>>int_value;
		(*user_param).flag_store_solution_line=int_value;
		required_parameters[8]=true;
		parameters_out<<variable<<": "<<int_value<<endl;
	    }
	    else if(variable=="flag_store_paraview_solution_cells")
	    {
		(*infile)>>int_value;
		(*user_param).flag_store_paraview_solution_cells=int_value;
		required_parameters[9]=true;
		parameters_out<<variable<<": "<<int_value<<endl;
	    }
	    else if(variable=="flag_store_paraview_solution_points")
	    {
		(*infile)>>int_value;
		(*user_param).flag_store_paraview_solution_points=int_value;
		required_parameters[10]=true;
		parameters_out<<variable<<": "<<int_value<<endl;
	    }
	    else if(variable=="no_out")
	    {
		(*infile)>>int_value;
		(*user_param).no_out=int_value;
		required_parameters[18]=true;
		parameters_out<<variable<<": "<<int_value<<endl;				
	    }
	    else if(variable=="flag_limiting")
	    {
		(*infile)>>int_value;
		(*user_param).flag_limiting=int_value;
		required_parameters[11]=true;
		parameters_out<<variable<<": "<<int_value<<endl;
	    }
	    else if(variable=="flag_trans")
	    {
		(*infile)>>int_value;
		(*user_param).flag_trans=int_value;
		required_parameters[26]=true;
		parameters_out<<variable<<": "<<int_value<<endl;
	    }
	    else if(variable=="nu")
	    {
		(*infile)>>double_value;
		(*user_param).nu=double_value;
		required_parameters[27]=true;
		parameters_out<<variable<<": "<<double_value<<endl;
	    }
	    else if(variable=="mdx2")
	    {
		(*infile)>>double_value;
		(*user_param).mdx2=double_value;
		required_parameters[28]=true;
		parameters_out<<variable<<": "<<double_value<<endl;
	    }
	    else if(variable=="ref_sol_flag")
	    {
		(*infile)>>int_value;
		(*user_param).ref_sol_flag=int_value;
		required_parameters[12]=true;
		parameters_out<<variable<<": "<<int_value<<endl;				
	    }
	    else if(variable=="refsolfilename")
	    {
		(*infile)>>string_value;
		strcpy(refsolfilename,string_value.c_str());
		(*user_param).refsolfilename=refsolfilename;
		required_parameters[19]=true;
		parameters_out<<variable<<": "<<string_value<<endl;
	    }
	    else if(variable=="p_degree")
	    {
		(*infile)>>int_value;
		(*user_param).p_degree=int_value;
		required_parameters[13]=true;
		parameters_out<<variable<<": "<<int_value<<endl;				
	    }
	    else if(variable=="densityfilename")
	    {
		(*infile)>>string_value;
		strcpy(densityfilename,string_value.c_str());
		(*user_param).densityfilename=densityfilename;
		parameters_out<<variable<<": "<<string_value<<endl;
	    }
	    else if(variable=="start_solution_file_flag")
	    {
		(*infile)>>int_value;
		(*user_param).start_solution_file_flag=int_value;
		required_parameters[14]=true;
		parameters_out<<variable<<": "<<int_value<<endl;
	    }
	    else if(variable=="startsolfilename")
	    {
		(*infile)>>string_value;
		strcpy(startsolfilename,string_value.c_str());
		(*user_param).startsolfilename=startsolfilename;
		required_parameters[20]=true;
		parameters_out<<variable<<": "<<string_value<<endl;				
	    }
	    else if(variable=="custom_out_a_0")
	    {
		(*infile)>>double_value;
		(*user_param).custom_out_a[0]=double_value;
		required_parameters[21]=true;
		parameters_out<<variable<<": "<<double_value<<endl;
	    }
	    else if(variable=="custom_out_a_1")
	    {
		(*infile)>>double_value;
		(*user_param).custom_out_a[1]=double_value;
		required_parameters[22]=true;
		parameters_out<<variable<<": "<<double_value<<endl;
	    }
	    else if(variable=="custom_out_b_0")
	    {
		(*infile)>>double_value;
		(*user_param).custom_out_b[0]=double_value;
		required_parameters[23]=true;
		parameters_out<<variable<<": "<<double_value<<endl;
	    }
	    else if(variable=="custom_out_b_1")
	    {
		(*infile)>>double_value;
		(*user_param).custom_out_b[1]=double_value;
		required_parameters[24]=true;
		parameters_out<<variable<<": "<<double_value<<endl;
	    }
	    else if(variable=="custom_out_length")
	    {
		(*infile)>>double_value;
		(*user_param).custom_out_length=double_value;
		required_parameters[25]=true;
		parameters_out<<variable<<": "<<double_value<<endl;
	    }
	    else if(variable=="voxel_domain_size_0")
	    {
		(*infile)>>double_value;
		(*user_param).voxel_domain_size[0]=double_value;
		parameters_out<<variable<<": "<<double_value<<endl;
	    }
	    else if(variable=="voxel_domain_size_1")
	    {
		(*infile)>>double_value;
		(*user_param).voxel_domain_size[1]=double_value;
		parameters_out<<variable<<": "<<double_value<<endl;
	    }
	    else if(variable=="periodic")
	    {
		(*infile)>>int_value;
		(*user_param).periodic=int_value;
		required_parameters[15]=true;
		parameters_out<<variable<<": "<<int_value<<endl;
	    }
	    else if(variable=="reflective_bc")
	    {
		(*infile)>>int_value;
		(*user_param).reflective_bc=int_value;
		required_parameters[16]=true;
		parameters_out<<variable<<": "<<int_value<<endl;
	    }
	    else if(variable=="CFL")
	    {
		(*infile)>>double_value;
		(*user_param).CFL=double_value;
		required_parameters[17]=true;
		parameters_out<<variable<<": "<<double_value<<endl;
	    }
	    else if(variable=="num_cores")
	    {
		(*infile)>>int_value;
		(*user_param).num_cores=int_value;
		parameters_out<<variable<<": "<<int_value<<endl;
	    }
	    else if(variable=="flag_ensure_positivity")
	    {
		(*infile)>>int_value;
		(*user_param).flag_ensure_positivity=int_value;
		parameters_out<<variable<<": "<<int_value<<endl;
	    }
	    else if(variable=="delta_pos")
	    {
		(*infile)>>double_value;
		(*user_param).delta_pos=double_value;
		parameters_out<<variable<<": "<<double_value<<endl;
	    }
        else if(variable =="store_grid_solution")
        {
            (*infile)>>int_value;
            (*user_param).store_grid_solution=int_value;
            parameters_out<<variable<<": "<<double_value<<endl;
        }
        else if(variable =="symmetric_solution")
        {
            (*infile)>>int_value;
            (*user_param).symmetric_solution=int_value;
            parameters_out<<variable<<": "<<int_value<<endl;
        }
        else if(variable!="-------------------------")
	    {
		cout<<"Invalid inputfile! Variable "<<variable<<" not known"<<endl;
		exit(1);
	    }
	}
	//Check if all required parameters have been set
	int i=1;
	while(required_parameters[i]==true&&i<18)
	{
		i++;
	}
	if(i!=18)
	{
		cout<<"Required parameter ";
		switch(i)
		{
			case(1):
			{
				cout<<"'closure'";
				break;
			}
			case(2):
			{
				cout<<"'test_case'";
				break;
			}
			case(3):
			{
				cout<<"'meshfilename'";
				break;
			}
			case(4):
			{
				cout<<"'outputfilename'";
				break;
			}
			case(5):
			{
				cout<<"'start_time'";
				break;
			}
			case(6):
			{
				cout<<"'Ftime'";
				break;
			}
			case(7):
			{
				cout<<"'flag_store_solution'";
				break;
			}
			case(8):
			{
				cout<<"'flag_store_solution_line'";
				break;
			}
			case(9):
			{
				cout<<"'flag_store_paraview_solution_cells'";
				break;
			}
			case(10):
			{
				cout<<"'flag_store_paraview_solution_points'";
				break;
			}
			case(11):
			{
				cout<<"'flag_limiting'";
				break;
			}
			case(12):
			{
				cout<<"'ref_sol_flag'";
				break;
			}
			case(13):
			{
				cout<<"'p_degree'";
				break;
			}
			case(14):
			{
				cout<<"'start_solution_file_flag'";
				break;
			}
			case(15):
			{
				cout<<"'periodic'";
				break;
			}
			case(16):
			{
				cout<<"'reflective_bc'";
				break;
			}
			case(17):
			{
				cout<<"'CFL'";
				break;
			}

		}
		cout<<" was not set in the inputfile"<<endl;
		(*infile).close();
		exit(1);
	}
	(*infile).close();

    // backwards time-stepping 
    if ( (*user_param).start_time > (*user_param).Ftime) 
    {
	(*user_param).flag_backwards_t_step = 1;
	parameters_out<<"flag_backwards_t_step: 1"<<endl;
    }
    else 
    {
	(*user_param).flag_backwards_t_step = 0;
    }
    //check parameters which are dependent from required parameters
    if((*user_param).flag_store_solution==1&&!(required_parameters[18]))
    {
        cout<<"Required parameter 'no_out' (due to flag_store_solution = 1) was not set in the inputfile"<<endl;
    }
    if((*user_param).ref_sol_flag==1&&!(required_parameters[19]))
    {
        cout<<"Required parameter 'refsolfilename' (due to ref_sol_flag = 1) was not set in the inputfile"<<endl;
    }
    if((*user_param).start_solution_file_flag==1&&!(required_parameters[20]))
    {
        cout<<"Required parameter 'startsolfilename' (due to start_solution_file_flag = 1) was not set in the inputfile"<<endl;
    }
	if((*user_param).flag_store_solution_line==1&&!(required_parameters[21]))
    {
        cout<<"Required parameter 'custom_out_a_0' (due to flag_store_solution_line = 1) was not set in the inputfile"<<endl;
    }
	if((*user_param).flag_store_solution_line==1&&!(required_parameters[22]))
    {
        cout<<"Required parameter 'custom_out_a_1' (due to flag_store_solution_line = 1) was not set in the inputfile"<<endl;
    }
	if((*user_param).flag_store_solution_line==1&&!(required_parameters[23]))
    {
        cout<<"Required parameter 'custom_out_b_0' (due to flag_store_solution_line = 1) was not set in the inputfile"<<endl;
    }
	if((*user_param).flag_store_solution_line==1&&!(required_parameters[24]))
    {
        cout<<"Required parameter 'custom_out_b_1' (due to flag_store_solution_line = 1) was not set in the inputfile"<<endl;
    }
	if((*user_param).flag_store_solution_line==1&&!(required_parameters[25]))
    {
        cout<<"Required parameter 'custom_out_length' (due to flag_store_solution_line = 1) was not set in the inputfile"<<endl;
    }
	if((*user_param).flag_limiting==1&&!(required_parameters[26]))
    {
        cout<<"Required parameter 'flag_trans' (due to flag_limiting = 1) was not set in the inputfile"<<endl;
    }
	if((*user_param).flag_limiting==1&&!(required_parameters[27]))
    {
        cout<<"Required parameter 'nu' (due to flag_limiting = 1) was not set in the inputfile"<<endl;
    }
	if((*user_param).flag_limiting==1&&!(required_parameters[28]))
    {
        cout<<"Required parameter 'mdx2' (due to flag_limiting = 1) was not set in the inputfile"<<endl;
    }
	
	if( ((*user_param).test_case ==5) || ((*user_param).test_case == 20) )
    {
        (*user_param).voxel_domain_size[0] = 7.0;
		parameters_out<<"voxel_domain_size_0: 7.0"<<endl;
        (*user_param).voxel_domain_size[1] = 7.0;
		parameters_out<<"voxel_domain_size_1: 7.0"<<endl;
		(*user_param).periodic_domain_size_x = 7.0;
		parameters_out<<"periodic_domain_size_x: 7.0"<<endl;
		(*user_param).periodic_domain_size_y = 7.0;
		parameters_out<<"periodic_domain_size_y: 7.0"<<endl;
    }
    else if ( strcmp((*user_param).meshfilename,"meshes/reflective_arc") == 0 )
    {
        (*user_param).voxel_domain_size[0] = 11.0;
		parameters_out<<"voxel_domain_size_0: 11.0"<<endl;
        (*user_param).voxel_domain_size[1] = 12.0;
		parameters_out<<"voxel_domain_size_1: 12.0"<<endl;
    }
    else
    {
        (*user_param).voxel_domain_size[0] = 1.0;
		parameters_out<<"voxel_domain_size_0: 1.0"<<endl;
        (*user_param).voxel_domain_size[1] = 1.0;
		parameters_out<<"voxel_domain_size_1: 1.0"<<endl;
		(*user_param).periodic_domain_size_x = 1.0;
		parameters_out<<"periodic_domain_size_x: 1.0"<<endl;
		(*user_param).periodic_domain_size_y = 1.0;
		parameters_out<<"periodic_domain_size_y: 7.0"<<endl;
    }
	parameters_out<<"densityfilename: unit_density.dat"<<endl;
}

/* void init_density
 * Initialize density values on all Gauss points
 * in: mesh element data structure    -- mesh_element
 * in: mesh_node data structure       -- mesh_node
 * in: mesh edge data structure       -- mesh_edge
 * in: parameter structure	      -- user_param
 * in: number of elements	      -- Nelts
 * out: parameter structure 	      -- mesh_element, mesh_edge
 */
void init_density(element* mesh_element,
		  node* mesh_node,
		  edge* mesh_edge,
		  parameters* user_param,
		  int Nelts)
{ 
    int i;
    int j;
    int k;
    int E;
    int iedge;			// Index to edge
    int result;
    double hx, hy;
    double left_b_x, left_b_y, right_b_x, right_b_y; // Define the domain
    int nrx, nry;
    int cx,cy;

    FILE *file;// *fopen();

    file = fopen((*user_param).densityfilename, "r");

    // Read information about the domain in the first three rows
    result = fscanf(file,"%lf %lf",&left_b_x,&right_b_x);
    if(!(result)) {
	fprintf(stderr,"Wrong density file format in first line!\n");
	exit(1);
    }

    result = fscanf(file,"%lf %lf",&left_b_y,&right_b_y);
    if(!(result)) {
	fprintf(stderr,"Wrong density file format in second line!\n");
	exit(1);
    }

    result = fscanf(file,"%d %d",&nrx,&nry);
    if(!(result)) {
	fprintf(stderr,"Wrong density file format in third line!\n");
	exit(1);
    }

    hx = (right_b_x - left_b_x)/nrx;
    hy = (right_b_y - left_b_y)/nry;

    /*printf("left_b_x = %lf\t right_b_x = %lf\n",left_b_x,right_b_x);*/
    /*printf("left_b_y = %lf\t right_b_y = %lf\n",left_b_y,right_b_y);*/
    /*printf("nrx = %d\t nry = %d\n",nrx,nry);*/

    double* density_data = NULL;
    density_data = (double*)malloc(nrx*nry*sizeof(double));

    for(j=nry-1;j>=0;j--) {
	for(i=0;i<nrx;i++) {
	    fscanf(file, "%lf ", &density_data[j*nrx+i]);
	    //printf("density_data[%d] = %lf",j*nrx+i,density_data[j*nrx+i]);
	}
    }

    // Initialize density values on elements
    for(E=1;E<Nelts+1;E++) {
	for(k=1;k<ngpts+1; k++) {
	    cx = floor(mesh_element[E].el_gpts_x[k] / hx);
	    cy = floor(mesh_element[E].el_gpts_y[k] / hy);
	    mesh_element[E].el_density[k] = 1.0;//density_data[cy*nrx+cx];
	    //printf("gpx = %lf\t gpy = %lf\t density = %lf\n", mesh_element[E].el_gpts_x[k], mesh_element[E].el_gpts_y[k], mesh_element[E].el_density[k]);
	}
	/*getchar();*/
    }

    // Find minimal density
    (*user_param).min_density = density_data[0];
    for (i=0; i < nrx*nry; i++){
	if ( density_data[i] < (*user_param).min_density)
	    (*user_param).min_density = 1.0;//density_data[i];
    }


    // Initialize density values on edges
    for(E=1;E<Nelts+1;E++) 
    {
	for(i=1;i<NUM_EDGES+1;i++) 
	{
	    iedge = mesh_element[E].edge[i];
	    for(k=1;k<nlgpts+1; k++) 
	    {
		cx = floor(mesh_edge[iedge].ed_phys_coords_x[k] / hx);
		if (cx >= nrx)
		    cx = nrx - 1;
		cy = floor(mesh_edge[iedge].ed_phys_coords_y[k] / hy);
		if (cy >= nry)
		    cy = nry - 1;
		mesh_edge[iedge].ed_density[k] = 1.0;//density_data[cy*nrx+cx];
	    }
	}
    }
    fclose (file);

    free(density_data);
}

/* void init_gpts
 * Initialize information of Gauss points on the mesh
 * in: mesh element data structure    -- mesh_element
 * in: mesh_node data structure       -- mesh_node
 * in: mesh edge data structure       -- mesh_edge
 * in: number of elements	      -- Nelts
 * out: parameter structure 	      -- mesh_element, mesh_edge
 */
void init_gpts(element* mesh_element,
	       node* mesh_node,
	       edge* mesh_edge,
	       int Nelts,
	       int global_edge_count)
{ 
    int i;
    int j;
    int k;
    int E;
    int deg;
    int Nloc;
    int iedge;			// Index to edge
    double end_points[2];       // End points of physical edge
    double* xg = NULL;          // Quadrature points on the edge
    double* w_edge  = NULL;     // Weights on edge
    double phys_coords[2];      //physical coordinates of gauss points

    double* gpx = NULL;
    double* gpy = NULL;
    double* w_el = NULL;

    double* phi_vec = NULL;
    int E1=0;
    int E2=0;
    double ref_coords1[2];
    double ref_coords2[2];

    xg = (double*)malloc((nlgpts+1)*sizeof(double));
    w_edge = (double*)malloc((nlgpts+1)*sizeof(double));

    gpx = (double*)malloc((ngpts+1)*sizeof(double));
    gpy = (double*)malloc((ngpts+1)*sizeof(double));
    w_el = (double*)malloc((ngpts+1)*sizeof(double));

    // Get polynomial degree
    deg = mesh_element[1].degree; // Assumes degree identical on all elements
    Nloc = (deg+1)*(deg+2)/2;

    // Allocate memory for basis functions
    phi_vec = (double*)malloc(Nloc*sizeof(double));

    // Gauss points on reference element
    gauss_pt(gpx,gpy,w_el);

    for(E=1;E<Nelts+1;E++) {
	for(k=1;k<ngpts+1;k++) {
	    // Map gauss points to physical element
	    map_to_physical_element(mesh_element,mesh_node,E,phys_coords,gpx[k],gpy[k]);

	    mesh_element[E].el_gpts_x[k] = phys_coords[0];
	    mesh_element[E].el_gpts_y[k] = phys_coords[1];
	    mesh_element[E].el_gpts_ref_x[k] = gpx[k];
	    mesh_element[E].el_gpts_ref_y[k] = gpy[k];
	    mesh_element[E].el_gpts_w[k] = w_el[k];

	    // Compute the values of basis functions
	    init_zero_d(phi_vec,Nloc);
	    init_monomial_basis(E,mesh_element,mesh_node,deg,gpx[k],gpy[k],phi_vec);

	    for (j=0;j < Nloc; j++)
		mesh_element[E].el_gpts_basis[(k-1)*Nloc+j] = phi_vec[j];
	}
    }

    // Initialize Gauss points on edges
    for(iedge=1;iedge<global_edge_count+1;iedge++) 
    {
	if(mesh_edge[iedge].edge_type == INTERIOR)
	{
	    E1 = mesh_edge[iedge].neighbour[1];
	    E2 = mesh_edge[iedge].neighbour[2];

	    // Get end points of iedge for 
	    get_end_points(E1,iedge,mesh_node,mesh_element,mesh_edge,end_points);
	    // Initialize quadrature points and weights for line [a,b]
	    // these are corresponding physical coordinates
	    gl_weight(nlgpts,end_points[0],end_points[1],xg,w_edge);
	    for(k=1;k<nlgpts+1;k++) 
	    {
		get_phys_coords_edge(iedge,mesh_node,mesh_element,mesh_edge,xg[k],phys_coords);

		mesh_edge[iedge].ed_phys_coords_x[k] = phys_coords[0];
		mesh_edge[iedge].ed_phys_coords_y[k] = phys_coords[1];
		mesh_edge[iedge].ed_gpts_w[k] = w_edge[k];

		//each physical gauss point is mapped to the reference element edge of E1
		init_zero_d(ref_coords1,2);
		map_to_reference_element(mesh_element,mesh_node,E1,ref_coords1,phys_coords[0],phys_coords[1]);
		mesh_edge[iedge].E1_ed_gpts_x[k] = ref_coords1[0];
		mesh_edge[iedge].E1_ed_gpts_y[k] = ref_coords1[1];

		// Compute the values of basis functions
		init_zero_d(phi_vec,Nloc);
		init_monomial_basis(E1,mesh_element,mesh_node,deg,ref_coords1[0],ref_coords1[1],phi_vec);

		//store the value of basis functions at gauss points
		for(j=0;j<Nloc;j++)
		{
		    mesh_edge[iedge].ed_gpts_basis_E1[(k-1)*Nloc+j] = phi_vec[j];
		}

		//each physical gauss point is mapped to reference element edge of E2
		init_zero_d(ref_coords2,2);
		map_to_reference_element(mesh_element,mesh_node,E2,ref_coords2,phys_coords[0],phys_coords[1]);
		mesh_edge[iedge].E2_ed_gpts_x[k] = ref_coords2[0];
		mesh_edge[iedge].E2_ed_gpts_y[k] = ref_coords2[1];

		// Compute the values of basis functions
		init_zero_d(phi_vec,Nloc);
		init_monomial_basis(E2,mesh_element,mesh_node,deg,ref_coords2[0],ref_coords2[1],phi_vec);

		//store the value of basis functions at gauss points
		for(j=0;j<Nloc;j++)
		{
		    mesh_edge[iedge].ed_gpts_basis_E2[(k-1)*Nloc+j] = phi_vec[j];
		}

	    }//gauss pts 
	}//loop over interior
	else
	{
	    E1 = mesh_edge[iedge].neighbour[1];
	    get_end_points(E1,iedge,mesh_node,mesh_element,mesh_edge,end_points);

	    // Initialize quadrature points and weights for line [a,b]
	    // these are corresponding physical coordinates
	    gl_weight(nlgpts,end_points[0],end_points[1],xg,w_edge);
	    for(k=1;k<nlgpts+1;k++)
	    {
            get_phys_coords_edge(iedge,mesh_node,mesh_element,mesh_edge,xg[k],phys_coords);
		    mesh_edge[iedge].ed_phys_coords_x[k] = phys_coords[0];
		    mesh_edge[iedge].ed_phys_coords_y[k] = phys_coords[1];
            mesh_edge[iedge].ed_gpts_w[k] = w_edge[k];

            init_zero_d(ref_coords1,2);
            map_to_reference_element(mesh_element,mesh_node,E1,ref_coords1,phys_coords[0],phys_coords[1]);
	
            mesh_edge[iedge].E1_ed_gpts_x[k] = ref_coords1[0];
            mesh_edge[iedge].E1_ed_gpts_y[k] = ref_coords1[1];

        /*    if((E1 == 2) || (E1 ==24))
            {
                printf("E %d idge %d k %d coords %lf %lf xg %lf \n",E1,iedge,k,phys_coords[0],phys_coords[1],xg[k]);
                printf("E %d idge %d k %d coords %lf %lf xg %lf \n",E1,iedge,k,ref_coords1[0],ref_coords1[1],xg[k]);

            getchar();
            }*/

		
            // Compute the values of basis functions
            init_zero_d(phi_vec,Nloc);
            init_monomial_basis(E1,mesh_element,mesh_node,deg,ref_coords1[0],ref_coords1[1],phi_vec);
            //store the value of basis functions at gauss points
            for(j=0;j<Nloc;j++)
            {
                        
                mesh_edge[iedge].ed_gpts_basis_E1[(k-1)*Nloc+j] = phi_vec[j];

               //printf("iedge %d phi_vec %lf @j %d \n",iedge,phi_vec[j],j);
               // getchar();
            }
	    }
	}//boundary edges
    }//loop over edges

    free(phi_vec);
    free(w_el);
    free(gpy);
    free(gpx);
    free(w_edge);
    free(xg);
}

/* void init_M1_trafo_matrices
 * Initialize M1 transformation matrices
 * these are stored in corresponding MATLAB rountine
 * inout: parameter structure 	      -- user_param
 */
void init_M1_trafo_matrices(parameters* user_param)
{
    int i, j, r, c;
    int nf, nphi; 		// Number of points in f and phi
    int result;

    FILE *file;// *fopen();

    file = fopen("trafo_M1.dat", "r");
    // Read information about the domain in the first three rows
    result = fscanf(file,"%d %d",&nf,&nphi);
    if(!(result)) {
	fprintf(stderr,"Wrong trafo file format in first line!\n");
	exit(1);
    }

    fscanf(file, "%lf ", &(*user_param).fmax);

    for (i=0;i<nf;i++) {
	fscanf(file, "%lf ", &(*user_param).f[i]);
    }

    for (j=0;j<nphi;j++) {
	fscanf(file, "%lf ", &(*user_param).phi[j]);
    }

    (*user_param).nf = nf;
    (*user_param).nphi = nphi;

    for (i=0;i<nf;i++) {
	for (j=0;j<nphi;j++) {
	    for (r=0;r<3;r++) {
		for (c=0;c<3;c++) {
		    fscanf(file, "%lf ", &(*user_param).R[i][j][r][c]);
		}
	    }
	}
    }

    for (i=0;i<nf;i++) {
	for (j=0;j<nphi;j++) {
	    for (r=0;r<3;r++) {
		for (c=0;c<3;c++) {
		    fscanf(file, "%lf ", &(*user_param).inv_R[i][j][r][c]);
		}
	    }
	}
    }

    fclose (file);
	    /*printf("gpx = %lf\t gpy = %lf\t density = %lf\n", mesh_element[E].el_gpts_x[k], mesh_element[E].el_gpts_y[k], mesh_element[E].el_density[k]);*/
	/*getchar();*/

}
