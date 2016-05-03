/* output.c
 * Philipp Monreal
*/
#include"output.h"

/* void store_solution
 * writes file "solution_'suffix'.dat" containing the solution value at all nodes
 * in: mesh_element
 * in: mesh_node
 * in: solution_vector
 * in: num_elt
 * in: num_nodes
 * in: suffix
*/
void store_solution(
	    element* mesh_element, 
	    node* mesh_node, 
	    double* solution_vector,
	    int num_elt,
	    int num_nodes,
	    int s_dim_count,
	    char* suffix,       
	    parameters* user_param)
{
    //triangle nodes
    int n1=0;  //node 1
    int n2=0;  //node 2
    int n3=0;  //node 3

    //coordinates of nodes 1-3
    double x1=0.0;
    double y1=0.0;
    double x2=0.0;
    double y2=0.0;
    double x3=0.0;
    double y3=0.0;
    int E=0;
    double exact_sol1[SYSTEM_DIM];
    double exact_sol2[SYSTEM_DIM];
    double exact_sol3[SYSTEM_DIM];

    double error1,error2,error3;

    //value of solution at each node
    double sol_val1[SYSTEM_DIM];
    double sol_val2[SYSTEM_DIM];
    double sol_val3[SYSTEM_DIM];
    char sdim[10];

    //output variables
    FILE* output_file;
    char filename[30];

    // Handle filename
    filename[0] = '\0';

    strcat(filename,"solution/");
    sprintf(sdim,"%d",s_dim_count);
    strcat(filename,sdim);
    strcat(filename,"solution_");
    strcat(filename,suffix);
    strcat(filename,".dat");

    // Store solution
    output_file = fopen(filename,"w");

    for(E=1;E<num_elt+1;E++) 
    {
        // get nodes of triangle E   
        n1 = mesh_element[E].vertex[1];
        n2 = mesh_element[E].vertex[2];
        n3 = mesh_element[E].vertex[3];

        //physical coordinates of node 1
        x1 = mesh_node[n1].coord[0];
        y1 = mesh_node[n1].coord[1];
        
        //physical coordinates of node 2
        x2 = mesh_node[n2].coord[0];
        y2 = mesh_node[n2].coord[1];
        
        //physical coordinates of node 3
        x3 = mesh_node[n3].coord[0];
        y3 = mesh_node[n3].coord[1];

        /*// get approximate solution*/
        init_zero_d(sol_val1,SYSTEM_DIM);
        init_zero_d(sol_val2,SYSTEM_DIM);
        init_zero_d(sol_val3,SYSTEM_DIM);

        get_approx_solution(solution_vector,E,mesh_element,mesh_node,0.0,0.0,sol_val1);
	//bdry_function_val(x1,y1,0.0,0.0,1.0,E,mesh_element,mesh_node,solution_vector,user_param,exact_sol1);
        get_approx_solution(solution_vector,E,mesh_element,mesh_node,1.0,0.0,sol_val2);
	//bdry_function_val(x2,y2,0.0,0.0,1.0,E,mesh_element,mesh_node,solution_vector,user_param,exact_sol2);
        get_approx_solution(solution_vector,E,mesh_element,mesh_node,0.0,1.0,sol_val3);
	//bdry_function_val(x3,y3,0.0,0.0,1.0,E,mesh_element,mesh_node,solution_vector,user_param,exact_sol3);
      
        error1 =fabs((sol_val1[s_dim_count]-exact_sol1[s_dim_count]));
        error2 =fabs((sol_val2[s_dim_count]-exact_sol2[s_dim_count]));
        error3 =fabs((sol_val3[s_dim_count]-exact_sol3[s_dim_count]));
            
        fprintf(output_file,"%.16lf %.16lf %.16lf \n",x1,y1,error1);
        fprintf(output_file,"%.16lf %.16lf %.16lf \n",x2,y2,error2);
        fprintf(output_file,"%.16lf %.16lf %.16lf \n",x3,y3,error3);
    }
    fclose(output_file);
}

void save_error_coeffs(double* uh_tn_sol_vec,
                       double* uzero_sol_vec,
                       parameters* user_param,
                       int Nelts,
                       int Nloc)
{
    int idofs=0;
    int s_dim=0;
    int E =0;
    FILE* coef_outfile;
    char out_filename[80];
    char scenario_temp[80];
    double error_coef =0.0;

    out_filename[0] = '\0';
    strcat(out_filename,"errors/");
    sprintf(scenario_temp, "s%d", (*user_param).test_case);
    strcat(out_filename,scenario_temp);
    strcat(out_filename,(*user_param).outputfilename);

    coef_outfile = fopen(out_filename,"w");

    fprintf(coef_outfile,"%s \n",(*user_param).outputfilename);
    fprintf(coef_outfile,"%.4f \n",0.0);
    for(E=1;E<Nelts+1;E++)
    {
        for(s_dim=0;s_dim < SYSTEM_DIM;s_dim++)
        {
            for(idofs=0;idofs<Nloc;idofs++)
            {
                error_coef = uh_tn_sol_vec[SYSTEM_DIM*(E-1)*Nloc + (Nloc*s_dim)+idofs]- uzero_sol_vec[SYSTEM_DIM*(E-1)*Nloc + (Nloc*s_dim)+idofs]; 

                fprintf(coef_outfile,"%d %d %.16lf \n",E,s_dim,error_coef);
            }
        }
    }
    fclose(coef_outfile);
}

void save_solution_coefficients(double* uh_tn_sol_vec,
                                element* mesh_element,
                                node* mesh_node,
                                parameters* user_param,
                                int Nloc,
                                int Nelts,
                                double time,
                                char* suffix)
{
    int idofs=0;
    int s_dim=0;
    int E =0;
    FILE* coef_outfile;
    char out_filename[80];
    char scenario_temp[80];

    out_filename[0] = '\0';
    strcat(out_filename,"solution/");
    sprintf(scenario_temp, "s%d", (*user_param).test_case);
    strcat(out_filename,scenario_temp);
    strcat(out_filename,(*user_param).outputfilename);
    strcat(out_filename,suffix);

    coef_outfile = fopen(out_filename,"w");

    fprintf(coef_outfile,"%s \n",(*user_param).outputfilename);
    fprintf(coef_outfile,"%.5f \n",time);
    for(E=1;E<Nelts+1;E++)
    {
        for(s_dim=0;s_dim < SYSTEM_DIM;s_dim++)
        {
            for(idofs=0;idofs<Nloc;idofs++)
            {
                fprintf(coef_outfile,"%d %d %.16lf \n",E,s_dim,uh_tn_sol_vec[SYSTEM_DIM*(E-1)*Nloc + (Nloc*s_dim)+idofs]);
            }
        }
    }
    fclose(coef_outfile);
}//save_solution_coeffients

void plot_realizability(double* uh_tn_sol_vec,
                       element* mesh_element,
                       node* mesh_node,
                       parameters* user_param,
                       int Nloc,
                       int Nelts,
                       double time,
                       char* suffix)
{
    int idofs=0;
    int s_dim=0;
    int E =0;
    FILE* coef_outfile;
    char out_filename[80];
    char scenario_temp[80];

    double psi_zero =0.0;
    double norm_psi_one =0.0;
    double approx_solution[SYSTEM_DIM];
    double ratio = 0.0;
    double max_ratio =0.0;
    int rel_flag=0;
    double plot_real=0.0;
    int k=0.0;

    int real_count=0;



    out_filename[0] = '\0';
    strcat(out_filename,"real/");
    sprintf(scenario_temp, "s%d", (*user_param).test_case);
    strcat(out_filename,scenario_temp);
    strcat(out_filename,(*user_param).outputfilename);
    strcat(out_filename,suffix);

    coef_outfile = fopen(out_filename,"w");

    fprintf(coef_outfile,"%s \n",(*user_param).outputfilename);
    fprintf(coef_outfile,"%.5f \n",time);
    for(E=1;E<Nelts+1;E++)
    {
       
        for(k=0;k<9;k++) 
        {
            init_zero_d(approx_solution,SYSTEM_DIM);
            get_approx_solution(uh_tn_sol_vec,E,mesh_element,mesh_node,mesh_element[E].interior_gpts_ref_x[k],mesh_element[E].interior_gpts_ref_y[k],approx_solution);
            rel_flag =is_realizable(approx_solution);

            psi_zero = approx_solution[0];
            norm_psi_one = sqrt(pow(approx_solution[1],2.0) + pow(approx_solution[2],2.0));

            ratio = norm_psi_one/psi_zero;
            if(ratio > max_ratio)
            {
                max_ratio = ratio;
            }
        }

        if(!(rel_flag))
        {
            real_count++;
        }

        for(s_dim=0;s_dim < SYSTEM_DIM;s_dim++)
        {
            for(idofs=0;idofs<Nloc;idofs++)
            {
                if(idofs==0)
                {
                    if(rel_flag ==1)
                    {
                        plot_real = 0.0;
                    }
                    else
                    {
                        plot_real = 1.0;
                    }
                    fprintf(coef_outfile,"%d %d %.16lf \n",E,s_dim,plot_real);
                }
                else
                {

                    fprintf(coef_outfile,"%d %d %.16lf \n",E,s_dim,0.0);
                }
            }
        }
        max_ratio=0.0;
    }
    //printf("----------------------------------------\n");
    //printf(" %d Elements out of %d are not realizable \n",real_count,Nelts);
    fclose(coef_outfile);
}//save_solution_coeffients




/*void output_sol_at_vec
 *writes solution evaluated at coordinates stored in x,y
 *in: uh_tn_sol_vec   	-- the solution of the current computation
 *in: mesh_node 	-- current mesh node data structure
 *in: mesh_element 	-- current mesh element data structure
 *in: nelts 		-- number of elements in the current mesh
 *in: x 		-- x-coordinates
 *in: y 		-- y-coordinates
 *in: length		-- length of x,y*/
void output_sol_at_vec(double* uh_tn_sol_vec,
		       node* mesh_node,
                       element* mesh_element,
                       int nelts,
                       double* x,
                       double* y,
                       int length)
{
    int E;	// index to element
    int i = 0;  // index for vector

    int found = 0; // determine whether point lies in element

    double coords[2];
    double ref_coords[2];
    double local_uh_tn_sol[SYSTEM_DIM];

    int s_dim_count = 0; // Moment that is being stored

    FILE* output_file;
    char filename[50];

    // Handle filename
    filename[0] = '\0';

    strcat(filename,"solution/sol_at_vec.dat");

    // Store solution
    output_file = fopen(filename,"w");

    for (i=0; i<length; i++)
    {
        found = 0;
        E = 1;
        coords[0] = x[i];
        coords[1] = y[i];

        // Find element containing physical coordinates
        while(!(found) && (E < nelts)) 
        {
            found = check_point_in_element(E,mesh_element,mesh_node,coords);

            if(found) 
            {
                 // find reference coordinates on the element
                map_to_reference_element(mesh_element,mesh_node,E,ref_coords,x[i],y[i]);

                // evaluate solution vector
                get_approx_solution(uh_tn_sol_vec,E,mesh_element,mesh_node,ref_coords[0],ref_coords[1],local_uh_tn_sol);

                // write output
                fprintf(output_file,"%.15lf %.15lf %.15lf \n",x[i],y[i],local_uh_tn_sol[s_dim_count]);
            }

            E++;
        }

        if ((E = nelts) && !(found))
        {
            printf("xcoords %lf %lf \n",coords[0],coords[1]);
            fprintf(stderr,"Error output coordinates not within domain! \n");
            exit(1);
        }
    }//i

    // tidy up
    fclose(output_file);
}

/*void ctrl_output_line
 *in :uh_tn_sol_vec 	-- the solution of the current computation
 *in :mesh_node 	-- current mesh node data structure
 *in: mesh_element 	-- current mesh element data structure
 *in: parameters   	-- user_param
 *in: num_elts 		-- number of elements in the current mesh*/
void ctrl_output_line(double* uh_tn_sol_vec,
                      node* mesh_node,
                      element* mesh_element,
                      parameters* user_param,
                      int num_elts)
{
    // Custom output at specified vector
    int length_out = (*user_param).custom_out_length;
    double out_x[length_out];
    double out_y[length_out];

    double x_temp;

    int j;

    for (j=0; j<length_out; j++)
    {
        out_x[j] = (*user_param).custom_out_a[0] + j*((*user_param).custom_out_b[0] - (*user_param).custom_out_a[0])/(length_out-1);
        out_y[j] = (*user_param).custom_out_a[1] + j*((*user_param).custom_out_b[1] - (*user_param).custom_out_a[1])/(length_out-1);
    }
    output_sol_at_vec(uh_tn_sol_vec,mesh_node,mesh_element,num_elts,out_x,out_y,length_out);
}

/*void save_cell_solution_vtp
 *in: uh_tn_sol_vec 	-- the solution of the current computation
 *in: mesh_element 	-- current mesh element data structure					   
 *in: mesh_node 	-- current mesh node data structure					   
 *in: user_param	-- user parameter set				   
 *in: num_elts	-- number of elements				   
 *in: num_nodes	-- number of nodes				   					  
 *in: suffix		-- suffix for filename	
 *Generates solutionfile in .vtp-format for paraview calculating the solution for each cell of a triangle	   
 */
void save_cell_solution_vtp(double* uh_tn_sol_vec,
					   element* mesh_element,
					   node* mesh_node,
					   parameters* user_param,
					   int num_elt,
					   int num_nodes,
					   char* suffix)
{
    using namespace std;

    int nPoints,nVerts,nLines,nStrips,nPolys;
    nPoints=num_nodes;
    nVerts=0;
    nLines=0;
    nStrips=0;
    nPolys=num_elt;

    ofstream outfile;
    
    char out_filename[80];
    char scenario_temp[80];

    out_filename[0] = '\0';
    strcat(out_filename,"solution/");
    sprintf(scenario_temp, "paraview_s%d", (*user_param).test_case);
    strcat(out_filename,scenario_temp);
    strcat(out_filename,(*user_param).outputfilename);
    strcat(out_filename,"_cells_");
    strcat(out_filename,suffix);
    strcat(out_filename,".vtp");

    outfile.open(out_filename);
  
    //Header
    outfile<<"<?xml version=\"1.0\"?>"<<endl;
    outfile<<"<VTKFile type=\"PolyData\" version=\"0.1\" byte_order=\"LittleEndian\">"<<endl;
    outfile<<"<PolyData>"<<endl;
    outfile<<"<Piece NumberOfPoints=\""<<nPoints<<"\"";
    outfile<<" NumberOfVerts=\""<<nVerts<<"\"";
    outfile<<" NumberOfLines=\""<<nLines<<"\"";
    outfile<<" NumberOfStrips=\""<<nStrips<<"\"";
    outfile<<" NumberOfPolys=\""<<nPolys<<"\">"<<endl;

    //'Points' section, defining the points by physical coordinates of the nodes
    outfile<<"<Points>"<<endl;
    outfile<<"<DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">"<<endl;
    double x,y;
    x=0;
    y=0;
    for(int i=1;i<=nPoints;i++)
    {
        x = mesh_node[i].coord[0];
        y = mesh_node[i].coord[1];
        outfile<<x<<" "<<y<<" "<<"0"<<endl;
    }
    outfile<<"</DataArray>"<<endl;
    outfile<<"</Points>"<<endl;

    //'CellData' section, calculating cell values by average value of each triangle
    outfile<<"<CellData Scalars=\"average_cell_data\">"<<endl;
    double average_cell_data[SYSTEM_DIM];
    for(int s_dim=0;s_dim<SYSTEM_DIM;s_dim++)
    {
        outfile<<"<DataArray type=\"Float64\" Name=\"moment "<<s_dim<<"\" format=\"ascii\">"<<endl;
        for(int E=1;E<=nPolys;E++)
        {
            cell_average(uh_tn_sol_vec,E,mesh_element,mesh_node,average_cell_data);
            outfile<<average_cell_data[s_dim]<<endl;
        }
        outfile<<"</DataArray>"<<endl;
    }
    outfile<<"</CellData>"<<endl;
    
    //'Polys' section, defining the triangles by connectivity of nodes
    outfile<<"<Polys>"<<endl;
    outfile<<"<DataArray type=\"Int64\" Name=\"connectivity\" format=\"ascii\">"<<endl;
    int n1,n2,n3;
    n1=0;
    n2=0;
    n3=0;
    for(int E=1;E<=nPolys;E++)
    {
        // get nodes of triangle E   
        n1 = mesh_element[E].vertex[1];
        n2 = mesh_element[E].vertex[2];
        n3 = mesh_element[E].vertex[3];
        outfile<<n1-1<<" "<<n2-1<<" "<<n3-1<<endl;
    }
    outfile<<"</DataArray>"<<endl;
    outfile<<"<DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\">"<<endl;
    for(int E=1;E<=nPolys;E++)
    {
        outfile<<3*E<<" ";
    }
    outfile<<endl;
    outfile<<"</DataArray>"<<endl;
    outfile<<"</Polys>"<<endl;

    //Footer
    outfile<<"</Piece>"<<endl;
    outfile<<"</PolyData>"<<endl;
    outfile<<"</VTKFile>";
    outfile.close();

}
/*void save_point_solution_vtp
 *in: uh_tn_sol_vec 	-- the solution of the current computation
 *in: mesh_element 	-- current mesh element data structure					   
 *in: mesh_node 	-- current mesh node data structure					   
 *in: user_param	-- user parameter set				   
 *in: num_elts	-- number of elements				   
 *in: num_nodes	-- number of nodes				   					  
 *in: suffix		-- suffix for filename
 *Generates solutionfile in .vtp-format for paraview calculating the solution at each node of a triangle
 */
void save_point_solution_vtp(double* uh_tn_sol_vec,
					   element* mesh_element,
					   node* mesh_node,
					   parameters* user_param,
					   int num_elt,
					   int num_nodes,
					   char* suffix)
{
    using namespace std;

    int nPoints,nVerts,nLines,nStrips,nPolys;
    nPoints=3*num_elt;
    nVerts=0;
    nLines=0;
    nStrips=0;
    nPolys=num_elt;

    ofstream outfile;
   	
    char out_filename[80];
    char scenario_temp[80];

    out_filename[0] = '\0';
    strcat(out_filename,"solution/");
    sprintf(scenario_temp, "paraview_s%d", (*user_param).test_case);
    strcat(out_filename,scenario_temp);
    strcat(out_filename,(*user_param).outputfilename);
    strcat(out_filename,"_points_");
    strcat(out_filename,suffix);
    strcat(out_filename,".vtp");
    outfile.open(out_filename);
    
    //Header
    outfile<<"<?xml version=\"1.0\"?>"<<endl;
    outfile<<"<VTKFile type=\"PolyData\" version=\"0.1\" byte_order=\"LittleEndian\">"<<endl;
    outfile<<"<PolyData>"<<endl;
    outfile<<"<Piece NumberOfPoints=\""<<nPoints<<"\"";
    outfile<<" NumberOfVerts=\""<<nVerts<<"\"";
    outfile<<" NumberOfLines=\""<<nLines<<"\"";
    outfile<<" NumberOfStrips=\""<<nStrips<<"\"";
    outfile<<" NumberOfPolys=\""<<nPolys<<"\">"<<endl;

    //'Points' section, defining the points by physical coordinates of the nodes
    outfile<<"<Points>"<<endl;
    outfile<<"<DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">"<<endl;
    double x,y;
    x=0;
    y=0;

	for(int E=1;E<=nPolys;E++)
	{
		for(int j=1;j<4;j++)
		{
			x = mesh_node[(mesh_element[E].vertex[j])].coord[0];
			y = mesh_node[(mesh_element[E].vertex[j])].coord[1];
			outfile<<x<<" "<<y<<" "<<"0"<<endl;
        
		}
	}
	
    outfile<<"</DataArray>"<<endl;
    outfile<<"</Points>"<<endl;

    //'PointData' section, map physical to reference coordinates and calculate approx solution for all nodes of each triangle
    outfile<<"<PointData Scalars=\"my_scalars\">"<<endl;
    double pointdata[SYSTEM_DIM];
    double ref_coords[2];
    for(int s_dim=0;s_dim<SYSTEM_DIM;s_dim++)
    {
        outfile<<"<DataArray type=\"Float64\" Name=\"moment "<<s_dim<<"\" format=\"ascii\">"<<endl;
	 for(int E=1;E<=nPolys;E++)
	 {
	     for(int j=1;j<4;j++)
	     {
	         x = mesh_node[(mesh_element[E].vertex[j])].coord[0];
		  y = mesh_node[(mesh_element[E].vertex[j])].coord[1];
		  map_to_reference_element(mesh_element,mesh_node,E,ref_coords,x,y);
		  get_approx_solution(uh_tn_sol_vec,E,mesh_element,mesh_node,ref_coords[0],ref_coords[1],pointdata);
		  outfile<<(pointdata[s_dim])<<endl;
	     }
	 }
        outfile<<"</DataArray>"<<endl;
    }
    outfile<<"</PointData>"<<endl;
    
    //'Polys' section, defining the triangles by connectivity of nodes
    outfile<<"<Polys>"<<endl;
    outfile<<"<DataArray type=\"Int64\" Name=\"connectivity\" format=\"ascii\">"<<endl;
	int n=0;
	for(int E=0;E<nPolys;E++)
	{
		n=3*E;
		outfile<<n<<" "<<n+1<<" "<<n+2<<endl;
	}
    outfile<<"</DataArray>"<<endl;
    outfile<<"<DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\">"<<endl;
    for(int i=1;i<=nPolys;i++)
    {
		outfile<<3*i<<" ";
    }
    outfile<<endl;
    outfile<<"</DataArray>"<<endl;
    outfile<<"</Polys>"<<endl;
    
    //Footer
    outfile<<"</Piece>"<<endl;
    outfile<<"</PolyData>"<<endl;
    outfile<<"</VTKFile>";
    outfile.close();
}

/*void save_point_solution_error_vtp
 *in: uh_tn_sol_vec 	-- the solution of the current computation
 *in: mesh_element 	-- current mesh element data structure					   
 *in: mesh_node 	-- current mesh node data structure					   
 *in: user_param	-- user parameter set				   
 *in: num_elts	-- number of elements				   
 *in: num_nodes	-- number of nodes				   					  
 *in: suffix		-- suffix for filename
 *Generates solutionfile in .vtp-format for paraview calculating the solution at each node of a triangle and its error compared to a given solution
 */
void save_point_solution_error_vtp(double* uh_tn_sol_vec,
					   element* mesh_element,
					   node* mesh_node,
					   parameters* user_param,
					   int num_elt,
					   int num_nodes,
                                      double time,
					   char* suffix)
{
    using namespace std;

    int nPoints,nVerts,nLines,nStrips,nPolys;
    nPoints=3*num_elt;
    nVerts=0;
    nLines=0;
    nStrips=0;
    nPolys=num_elt;

    ofstream outfile;
   	
    char out_filename[80];
    char scenario_temp[80];

    out_filename[0] = '\0';
    strcat(out_filename,"solution/");
    sprintf(scenario_temp, "paraview_s%d", (*user_param).test_case);
    strcat(out_filename,scenario_temp);
    strcat(out_filename,(*user_param).outputfilename);
    strcat(out_filename,"_points_error_");
    strcat(out_filename,suffix);
    strcat(out_filename,".vtp");
    outfile.open(out_filename);
    
    //Header
    outfile<<"<?xml version=\"1.0\"?>"<<endl;
    outfile<<"<VTKFile type=\"PolyData\" version=\"0.1\" byte_order=\"LittleEndian\">"<<endl;
    outfile<<"<PolyData>"<<endl;
    outfile<<"<Piece NumberOfPoints=\""<<nPoints<<"\"";
    outfile<<" NumberOfVerts=\""<<nVerts<<"\"";
    outfile<<" NumberOfLines=\""<<nLines<<"\"";
    outfile<<" NumberOfStrips=\""<<nStrips<<"\"";
    outfile<<" NumberOfPolys=\""<<nPolys<<"\">"<<endl;

    //'Points' section, defining the points by physical coordinates of the nodes
    outfile<<"<Points>"<<endl;
    outfile<<"<DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">"<<endl;
    double x,y;
    x=0;
    y=0;

	for(int E=1;E<=nPolys;E++)
	{
		for(int j=1;j<4;j++)
		{
			x = mesh_node[(mesh_element[E].vertex[j])].coord[0];
			y = mesh_node[(mesh_element[E].vertex[j])].coord[1];
			outfile<<x<<" "<<y<<" "<<"0"<<endl;
        
		}
	}
	
    outfile<<"</DataArray>"<<endl;
    outfile<<"</Points>"<<endl;

    //'PointData' section, map physical to reference coordinates and calculate approx solution for all nodes of each triangle and its error
    outfile<<"<PointData Scalars=\"my_scalars\">"<<endl;
    double pointdata[SYSTEM_DIM];
    double ref_coords[2];
    for(int s_dim=0;s_dim<SYSTEM_DIM;s_dim++)
    {
        outfile<<"<DataArray type=\"Float64\" Name=\"moment "<<s_dim<<"\" format=\"ascii\">"<<endl;
	 for(int E=1;E<=nPolys;E++)
	 {
	     for(int j=1;j<4;j++)
	     {
		  x = mesh_node[(mesh_element[E].vertex[j])].coord[0];
		  y = mesh_node[(mesh_element[E].vertex[j])].coord[1];
		  map_to_reference_element(mesh_element,mesh_node,E,ref_coords,x,y);
		  get_approx_solution(uh_tn_sol_vec,E,mesh_element,mesh_node,ref_coords[0],ref_coords[1],pointdata);
                outfile<<(pointdata[s_dim])<<endl;
	     }
        } 
        outfile<<"</DataArray>"<<endl;
    }
    double error[SYSTEM_DIM];
    double ref_data[SYSTEM_DIM];
    for(int s_dim=0;s_dim<SYSTEM_DIM;s_dim++)
    {
        outfile<<"<DataArray type=\"Float64\" Name=\"error moment "<<s_dim+1<<"\" format=\"ascii\">"<<endl;
        for(int E=1;E<=nPolys;E++)
        {
	     for(int j=1;j<4;j++)
	     {
                x = mesh_node[(mesh_element[E].vertex[j])].coord[0];
	         y = mesh_node[(mesh_element[E].vertex[j])].coord[1];
	         map_to_reference_element(mesh_element,mesh_node,E,ref_coords,x,y);
                get_approx_solution(uh_tn_sol_vec,E,mesh_element,mesh_node,ref_coords[0],ref_coords[1],pointdata);
                //bdry_function_val(x,y,ref_coords[0],ref_coords[1],time,E,mesh_element,mesh_node,uh_tn_sol_vec,user_param,ref_data);
                error[s_dim]=fabs(pointdata[s_dim]-ref_data[s_dim]);
                outfile<<(error[s_dim])<<endl;
	     }
	 }
        outfile<<"</DataArray>"<<endl;
    }
    outfile<<"</PointData>"<<endl;
    
    //'Polys' section, defining the triangles by connectivity of nodes
    outfile<<"<Polys>"<<endl;
    outfile<<"<DataArray type=\"Int64\" Name=\"connectivity\" format=\"ascii\">"<<endl;
	int n=0;
	for(int E=0;E<nPolys;E++)
	{
		n=3*E;
		outfile<<n<<" "<<n+1<<" "<<n+2<<endl;
	}
    outfile<<"</DataArray>"<<endl;
    outfile<<"<DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\">"<<endl;
    for(int i=1;i<=nPolys;i++)
    {
		outfile<<3*i<<" ";
    }
    outfile<<endl;
    outfile<<"</DataArray>"<<endl;
    outfile<<"</Polys>"<<endl;
    
    //Footer
    outfile<<"</Piece>"<<endl;
    outfile<<"</PolyData>"<<endl;
    outfile<<"</VTKFile>";
    outfile.close();
}

/* void save_ref_error
*  writes error of last timestep referring to reference solution
*  in: uh_tn_sol_vec 	-- the solution of the current computation
*  in: mesh_element 	-- current mesh element data structure					   
*  in: mesh_node 	-- current mesh node data structure					   
*  in: user_param	-- user parameter set				   
*  in: num_elts	-- number of elements				   
*  in: num_nodes	-- number of nodes				   					  
*  in: suffix		-- suffix for filename 
*/
void save_ref_error(double* uh_tn_sol_vec,
					   double* ref_uh_tn_sol_vec, 
                                     node* mesh_node,
                                     node* ref_mesh_node,                                   
					   element* mesh_element,
					   element* ref_mesh_element,
                                     voxel* mesh_voxel,
					   int num_nodes,
					   int num_elt,
                                     int ref_nelts,
                                     int num_voxels,
					   parameters* user_param)
{
    using namespace std;

    int nPoints,nVerts,nLines,nStrips,nPolys;
    nPoints=num_nodes;
    nVerts=0;
    nLines=0;
    nStrips=0;
    nPolys=num_elt;

    ofstream outfile;
    
    char out_filename[80];
    char scenario_temp[80];

    out_filename[0] = '\0';
    strcat(out_filename,"solution/");
    sprintf(scenario_temp, "paraview_s%d", (*user_param).test_case);
    strcat(out_filename,scenario_temp);
    strcat(out_filename,(*user_param).outputfilename);
    strcat(out_filename,"_ref_error");
    strcat(out_filename,".vtp");

    outfile.open(out_filename);
  
    //Header
    outfile<<"<?xml version=\"1.0\"?>"<<endl;
    outfile<<"<VTKFile type=\"PolyData\" version=\"0.1\" byte_order=\"LittleEndian\">"<<endl;
    outfile<<"<PolyData>"<<endl;
    outfile<<"<Piece NumberOfPoints=\""<<nPoints<<"\"";
    outfile<<" NumberOfVerts=\""<<nVerts<<"\"";
    outfile<<" NumberOfLines=\""<<nLines<<"\"";
    outfile<<" NumberOfStrips=\""<<nStrips<<"\"";
    outfile<<" NumberOfPolys=\""<<nPolys<<"\">"<<endl;

    //'Points' section, defining the points by physical coordinates of the nodes
    outfile<<"<Points>"<<endl;
    outfile<<"<DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">"<<endl;
    double x,y;
    x=0;
    y=0;
    for(int i=1;i<=nPoints;i++)
    {
        x = mesh_node[i].coord[0];
        y = mesh_node[i].coord[1];
        outfile<<x<<" "<<y<<" "<<"0"<<endl;;
    }
    outfile<<"</DataArray>"<<endl;
    outfile<<"</Points>"<<endl;

    //'CellData' section, calculating the error on each cell (average over error on all gauss points)
    outfile<<"<CellData Scalars=\"average_cell_data\">"<<endl;
    double average_cell_error[SYSTEM_DIM];
    for(int s_dim=0;s_dim<SYSTEM_DIM;s_dim++)
    {
        outfile<<"<DataArray type=\"Float64\" Name=\"ref_error moment "<<s_dim<<"\" format=\"ascii\">"<<endl;
        for(int E=1;E<=nPolys;E++)
        {
            compute_cell_error(uh_tn_sol_vec,ref_uh_tn_sol_vec,E,mesh_node,ref_mesh_node,mesh_element,ref_mesh_element,mesh_voxel,ref_nelts,num_voxels,average_cell_error);
            outfile<<average_cell_error[s_dim]<<endl;
        }
        outfile<<"</DataArray>"<<endl;
    }
    outfile<<"</CellData>"<<endl;
    
    //'Polys' section, defining the triangles by connectivity of nodes
    outfile<<"<Polys>"<<endl;
    outfile<<"<DataArray type=\"Int64\" Name=\"connectivity\" format=\"ascii\">"<<endl;
    int n1,n2,n3;
    n1=0;
    n2=0;
    n3=0;
    for(int E=1;E<=nPolys;E++)
    {
        // get nodes of triangle E   
        n1 = mesh_element[E].vertex[1];
        n2 = mesh_element[E].vertex[2];
        n3 = mesh_element[E].vertex[3];
        outfile<<n1-1<<" "<<n2-1<<" "<<n3-1<<endl;
    }
    outfile<<"</DataArray>"<<endl;
    outfile<<"<DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\">"<<endl;
    for(int E=1;E<=nPolys;E++)
    {
        outfile<<3*E<<" ";
    }
    outfile<<endl;
    outfile<<"</DataArray>"<<endl;
    outfile<<"</Polys>"<<endl;

    //Footer
    outfile<<"</Piece>"<<endl;
    outfile<<"</PolyData>"<<endl;
    outfile<<"</VTKFile>";
    outfile.close();
}

/* void init_point_collector
*  creates collector file for point solution files and writes header
*  in:  parameters* -- user_param
*/
void init_point_collector(parameters* user_param)
{
    char filename[80];
    char scenario_temp[80];
    filename[0] = '\0';
    strcat(filename,"solution/");
    sprintf(scenario_temp, "paraview_s%d", (*user_param).test_case);
    strcat(filename,scenario_temp);
    strcat(filename,(*user_param).outputfilename);
    strcat(filename,"_point_collector");
    strcat(filename,".pvd");
    (*user_param).point_collector_of = new ofstream(filename);
    (*((*user_param).point_collector_of))<<"<?xml version=\"1.0\"?>"<<endl;
    (*((*user_param).point_collector_of))<<"<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"LittleEndian\">"<<endl;
    (*((*user_param).point_collector_of))<<"<Collection>"<<endl;
}

/* void init_cell_collector
*  creates collector file for cell solution files and writes header
*  in:  parameters* -- user_param
*/
void init_cell_collector(parameters* user_param)
{
    char filename[80];
    char scenario_temp[80];
    filename[0] = '\0';
    strcat(filename,"solution/");
    sprintf(scenario_temp, "paraview_s%d", (*user_param).test_case);
    strcat(filename,scenario_temp);
    strcat(filename,(*user_param).outputfilename);
    strcat(filename,"_cell_collector");
    strcat(filename,".pvd");
    (*user_param).cell_collector_of = new ofstream(filename);
    (*((*user_param).cell_collector_of))<<"<?xml version=\"1.0\"?>"<<endl;
    (*((*user_param).cell_collector_of))<<"<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"LittleEndian\">"<<endl;
    (*((*user_param).cell_collector_of))<<"<Collection>"<<endl;    
}

/* void write_point_collector
*  collects all point solution files with time attribute
*  in:  parameters* -- user_param
*  in:  double      -- time
*  in:  char*       -- suffix
*/
void write_point_collector(parameters* user_param,double time, char* suffix)
{
    char filename[80];
    char scenario_temp[80];
    filename[0] = '\0';
    sprintf(scenario_temp, "paraview_s%d", (*user_param).test_case);
    strcat(filename,scenario_temp);
    strcat(filename,(*user_param).outputfilename);
    strcat(filename,"_points_");
    strcat(filename,suffix);
    strcat(filename,".vtp");
    (*((*user_param).point_collector_of))<<"<DataSet timestep=\""<<time<<"\" ";
    (*((*user_param).point_collector_of))<<"group=\"\" part=\"0\" file=\""<<filename<<"\"/>"<<endl;
    
}

/* void write_cell_collector
*  collects all cell solution files with time attribute
*  in:  parameters* -- user_param
*  in:  double      -- time
*  in:  char*       -- suffix
*/
void write_cell_collector(parameters* user_param,double time, char* suffix)
{
    char filename[80];
    char scenario_temp[80];
    filename[0] = '\0';
    sprintf(scenario_temp, "paraview_s%d", (*user_param).test_case);
    strcat(filename,scenario_temp);
    strcat(filename,(*user_param).outputfilename);
    strcat(filename,"_cells_");
    strcat(filename,suffix);
    strcat(filename,".vtp");
    (*((*user_param).cell_collector_of))<<"<DataSet timestep=\""<<time<<"\" ";
    (*((*user_param).cell_collector_of))<<"group=\"\" part=\"0\" file=\""<<filename<<"\"/>"<<endl;
}

/* void end_point_collector
*  writes footer of point collector file
*  in:  parameters* -- user_param
*/
void end_point_collector(parameters* user_param)
{
    (*((*user_param).point_collector_of))<<"</Collection>"<<endl;
    (*((*user_param).point_collector_of))<<"</VTKFile>"<<endl;
    (*((*user_param).point_collector_of)).close();
    delete (*user_param).point_collector_of;

}

/* void end_cell_collector
*  writes footer of cell collector file
*  in:  parameters* -- user_param
*/
void end_cell_collector(parameters* user_param)
{
    (*((*user_param).cell_collector_of))<<"</Collection>"<<endl;
    (*((*user_param).cell_collector_of))<<"</VTKFile>"<<endl;
    (*((*user_param).cell_collector_of)).close();
    delete (*user_param).cell_collector_of;

}

