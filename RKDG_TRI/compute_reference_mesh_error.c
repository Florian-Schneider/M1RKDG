#include"compute_reference_mesh_error.h"

/*void compute_relative_error
 *in:uh_tn_sol_vec -- the solution of the current computation
 *in:ref_uh_tn_sol_vec -- the reference solution on a very fine mesh
 *in:mesh_node -- current mesh node data structure
 *in:ref_mesh_node -- mesh data structure that gave reference solution
 *in:mesh_element -- current mesh element data structure
 *in:ref_mesh_element - mesh element structure that gave reference solution
 *in:num_elts -- number of elements in the reference mesh
 *in:num_voxes -- number of voxels of the fine mesh
 *in:s_dim --current dimension 
 *out: l1_error
 *out: linf_error*/
void compute_relative_error(double* uh_tn_sol_vec,
                            double* ref_uh_tn_sol_vec,
                            node* mesh_node,
                            node* ref_mesh_node,
                            element* mesh_element,
                            element* ref_mesh_element,
                            voxel* mesh_voxel, 
                            int nelts,
                            int num_voxels,
                            int s_dim,
                            double* l1_error,
                            double* linf_error)
{
    int E=0;
    int V=0;
    int E_v =0;
    int num_e_in_V=0;
    int k=0;
    int deg=0;
    int Nloc=0;
    double* gpx = NULL;
    double* gpy = NULL;
    double* w = NULL;
    double* inv_BE_T = NULL;
    double  phys_coords [2];
    double  ref_coords [2];
    double* phi_vec = NULL;
    double* grad_phi_mat = NULL;
    double det_BE_T;
    double local_uh_tn_sol[SYSTEM_DIM];
    double ref_local_uh_tn_sol[SYSTEM_DIM];
    int found =0;
    int found_e =0;
    int e_in_v=0;
    int V_found=0;
    int e_found=0;
    double local_error =0.0;
    double sup_error =0.0;
    double new_sup_error=0.0;
    double global_error =0.0;

    gpx = (double*)malloc((ngpts+1)*sizeof(double));
    gpy = (double*)malloc((ngpts+1)*sizeof(double));
    w = (double*)malloc((ngpts+1)*sizeof(double));
    inv_BE_T = (double*)malloc((2*2)*sizeof(double));

    //gauss points on reference element
    gauss_pt(gpx,gpy,w);

    //get polynomial degree
    deg = mesh_element[E].degree;
    Nloc = (deg+1)*(deg+2)/2;
    //allocate memory for basis functions
    phi_vec = (double*)malloc(Nloc*sizeof(double));
    grad_phi_mat = (double*)malloc(2*Nloc*sizeof(double));

    V=1;

    for(E=1;E<nelts+1;E++)
    {
        det_BE_T = get_Element_data(inv_BE_T,mesh_element,mesh_node,E);

        local_error =0.0;

        for(k=1;k<ngpts+1;k++)
        {
            init_zero_d(phys_coords,2);
            // Get physical coordinates of Gauss points on current mesh
            map_to_physical_element(mesh_element,mesh_node,E,phys_coords,gpx[k],gpy[k]);
            //printf("----------------phys_coords --------------%lf %lf \n",phys_coords[0],phys_coords[1]);
            init_zero_d(ref_local_uh_tn_sol,SYSTEM_DIM);
            get_approx_solution(uh_tn_sol_vec,E,mesh_element,mesh_node,gpx[k],gpy[k],local_uh_tn_sol);
            V=1;
            E_v=0;
            found =0;
            while(!(found)&& (V < num_voxels+1))
            {
                found = check_point_in_voxel(phys_coords,mesh_voxel,V);
                if(found)
                {
                    V_found=V;
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
                //printf("e_in_v %d found %d Vfound  %d \n",e_in_v,found,V_found);
                found_e = check_point_in_element(e_in_v,ref_mesh_element,ref_mesh_node,phys_coords);
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

            init_zero_d(ref_local_uh_tn_sol,SYSTEM_DIM);
            // Get triangle coordinates of Gauss points on reference mesh
            map_to_reference_element(ref_mesh_element,ref_mesh_node,e_found,ref_coords,phys_coords[0],phys_coords[1]);
            get_approx_solution(ref_uh_tn_sol_vec,e_found,ref_mesh_element,ref_mesh_node,ref_coords[0],ref_coords[1],ref_local_uh_tn_sol);

            // assumes that the area of the domain is 1 !!!
            // multiply by basis function at gauss point 
            local_error = local_error + det_BE_T*fabs(ref_local_uh_tn_sol[s_dim] - local_uh_tn_sol[s_dim])*w[k];
            new_sup_error = fabs(local_uh_tn_sol[s_dim]-ref_local_uh_tn_sol[s_dim]);

            if(new_sup_error > sup_error)
            {
                sup_error = new_sup_error;
            }

            V_found =0;
            found =0;
            found_e =0;

        }
        global_error = global_error + local_error;
    }//E loop
    *l1_error = global_error;
    *linf_error = sup_error;

    free(grad_phi_mat);
    free(phi_vec);
    free(inv_BE_T);
    free(w);
    free(gpy);
    free(gpx);
}

void compute_cell_error(double* uh_tn_sol_vec,
			    double* ref_uh_tn_sol_vec,
                         int E,
			    node* mesh_node,
			    node* ref_mesh_node,
			    element* mesh_element,
			    element* ref_mesh_element,
			    voxel* mesh_voxel, 
			    int ref_nelts,
			    int num_voxels,
                         double* average_cell_error)
{

    int V=0;
    int E_v =0;
    int num_e_in_V=0;
    int k=0;
    double* gpx = NULL;
    double* gpy = NULL;
    double* w= NULL;
    double  phys_coords [2];
    double  ref_coords [2];
    double local_uh_tn_sol[SYSTEM_DIM];
    double ref_local_uh_tn_sol[SYSTEM_DIM];
    int found =0;
    int found_e =0;
    int e_in_v=0;
    int V_found=0;
    int e_found=0;
    double local_error [3];
    double temp_coords[2];


    gpx = (double*)malloc((ngpts+1)*sizeof(double));
    gpy = (double*)malloc((ngpts+1)*sizeof(double));
    w = (double*)malloc((ngpts+1)*sizeof(double));

 
    //gauss points on reference element
    gauss_pt(gpx,gpy,w);

    for(int i=0;i<SYSTEM_DIM;i++)
    {
        local_error[i]=0;    
    }

    //calculate the error at each point and sum it up
    for(k=1;k<ngpts+1;k++)
    {
        init_zero_d(phys_coords,2);
        // Get physical coordinates of Gauss points on current mesh

        /* old */
        //map_to_physical_element(ref_mesh_element,ref_mesh_node,E,phys_coords,gpx[k],gpy[k]);

        /* new */
        map_to_physical_element(mesh_element,mesh_node,E,phys_coords,gpx[k],gpy[k]);

        init_zero_d(local_uh_tn_sol,SYSTEM_DIM);

        /* old */
        //get_approx_solution(ref_uh_tn_sol_vec,E,ref_mesh_element,ref_mesh_node,gpx[k],gpy[k],ref_local_uh_tn_sol);

        /* new */        
        get_approx_solution(uh_tn_sol_vec,E,mesh_element,mesh_node,gpx[k],gpy[k],local_uh_tn_sol);

        V=1;
        while(!(found)&& (V < num_voxels+1))
        {
            found = check_point_in_voxel(phys_coords,mesh_voxel,V);
            if(found)
            {
                V_found=V;
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
            //printf("e_in_v %d found %d Vfound  %d \n",e_in_v,found,V_found);
            /* changed: mesh_element -> ref_mesh_element, mesh_node-> ref_mesh_node */
            found_e = check_point_in_element(e_in_v,ref_mesh_element,ref_mesh_node,phys_coords);
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

        init_zero_d(ref_local_uh_tn_sol,SYSTEM_DIM);
        // Get triangle coordinates of Gauss points on reference mesh

        /* old */
        //map_to_reference_element(mesh_element,mesh_node,e_found,ref_coords,phys_coords[0],phys_coords[1]);
        //get_approx_solution(uh_tn_sol_vec,e_found,mesh_element,mesh_node,ref_coords[0],ref_coords[1],local_uh_tn_sol);

        /* new */
        map_to_reference_element(ref_mesh_element,ref_mesh_node,e_found,ref_coords,phys_coords[0],phys_coords[1]);
        get_approx_solution(ref_uh_tn_sol_vec,e_found,ref_mesh_element,ref_mesh_node,ref_coords[0],ref_coords[1],ref_local_uh_tn_sol);

        
        for(int i=0;i<SYSTEM_DIM;i++)
        {
            local_error[i] = local_error[i] + fabs(ref_local_uh_tn_sol[i] - local_uh_tn_sol[i])*w[k];
        }
        //std::cout<<"E: "<<E<<"local_error: "<<local_error[0]<<std::endl;
        found =0;
        found_e =0;
        E_v=0;

    }
        //calculate the average of all errors at gauss points
        for(int i=0;i<SYSTEM_DIM;i++)
        {
            average_cell_error[i]=local_error[i]/ngpts;
        }
    
    free(gpy);
    free(gpx);
    free(w);
}

/*load_ref_sol_coefs 
 * out: uh_tn_sol_vec -- current solution
 * in: mesh element data structure
 * in: file of coefficients
 * in: Nloc 
 * in: Nelts
 */
void load_sol_coeffs(double* uh_tn_sol_vec,
			 element* mesh_element,
			 node* mesh_node,
			 FILE* coeff_file,
			 int Nloc,
			 int Nelts)
{
    int idofs=0;
    int s_dim=0;
    int E =0;
    int result =0;
    double sol_coef=0.0;
    int ele_e =0;
    int sys_dim=0;
    double saved_time =0.0;

    for(E=1;E<Nelts+1;E++)
    {
        for(s_dim=0;s_dim < SYSTEM_DIM;s_dim++)
        {
            for(idofs=0;idofs<Nloc;idofs++)
            {
            result = fscanf(coeff_file,"%d  %d %lf",&ele_e,&sys_dim,&sol_coef);

            if(!(result))
            {
                fprintf(stderr,"Error in load coefficients \n");
                exit(1);
            }
            uh_tn_sol_vec[SYSTEM_DIM*(E-1)*Nloc + (Nloc*s_dim)+idofs] = sol_coef;
            }
        }
    }
}

void create_voxel(double alpha,
		          element* mesh_element,
	       	      node* mesh_node,
	              double x0,
                  double y0,
                  double x1,
                  double y1,
		          voxel* mesh_voxel,
		          int* num_voxels)
{
    double xj =0.0;
    double yj =0.0;
    double xjplus1 =0.0;
    double yjplus1 =0.0;

    double max_x =0.0;
    double max_y =0.0;

    int i=0;
    int j=0;

    int N_x =0;
    int N_y =0;
    max_x = fabs(x1-x0);
    max_y = fabs(y1-y0);


    N_x = ceil(max_x/alpha);
    N_y = ceil(max_y/alpha);

    int voxel_count=0;

    for(j=0;j<N_y;j++)
    {
        for(i=0;i<N_x;i++)
        {
            xj = x0 + i*alpha;
            yj = y0 + j*alpha;
            xjplus1 = alpha + i*alpha;
            yjplus1 = alpha + j*alpha;
            voxel_count++;

            //node1
          
 
            mesh_voxel[voxel_count].n1[0] = xj;
            mesh_voxel[voxel_count].n1[1] = yj;
            //node2
            mesh_voxel[voxel_count].n2[0] = xjplus1;
            mesh_voxel[voxel_count].n2[1] = yj;
            //node3
            mesh_voxel[voxel_count].n3[0] = xjplus1;
            mesh_voxel[voxel_count].n3[1] = yjplus1;
            //node4
            mesh_voxel[voxel_count].n4[0] = xj;
            mesh_voxel[voxel_count].n4[1] = yjplus1;


        }
    }

    *num_voxels = voxel_count;
}

void add_elements_to_voxel(int ref_num_elts,
			    int num_voxels,
			    voxel* ref_mesh_voxel,
			    element* ref_mesh_element,
			    node* ref_mesh_node,
			    int* added_elements)
{
    int E=0;
    int V=0;
    int local_element_count=0;
    int global_e_count=0;

    for(V=1;V<num_voxels+1;V++)
    {
        //search for elements that could be partially or fully enclosed in voxel
        local_element_count=0;
        for(E=1;E<ref_num_elts+1;E++)
        {
            if(check_triangle_in_voxel(E,V,ref_mesh_voxel,ref_mesh_node,ref_mesh_element))
            {
                local_element_count++;
                global_e_count++;
                ref_mesh_voxel[V].num_elements = local_element_count;
                ref_mesh_voxel[V].element_list[local_element_count-1] = E;
            }
        }

        if(ref_mesh_voxel[V].num_elements > VOXEL_MAX_NUM_ELTS-1)
        {
            fprintf(stderr,"Number of elements in voxel( %d )  exceeds maximum :::Increase VOXEL_MAX_ELTS in data_structure.h!! \n",local_element_count);
            exit(1);
        }	
    }
    *added_elements = global_e_count;
}


int check_triangle_in_voxel(int E,
			    int V,
			    voxel* mesh_voxel,
			    node* mesh_node,
			    element* mesh_element
			    )
{
    int n1=0;
    int n2=0;
    int n3=0;
    int result=0;
    int i=0;

    double n1_x=0.0;
    double n1_y=0.0;
    double n2_x=0.0;
    double n2_y=0.0;
    double n3_x=0.0;
    double n3_y=0.0;

    int found =0;
    double bb_coords[8];
    double bb_phys_coords[2];
    double voxel_coords[8];
    double vv_phys_coords[2];

    n1 = mesh_element[E].vertex[1];
    n2 = mesh_element[E].vertex[2];
    n3 = mesh_element[E].vertex[3];

    n1_x = mesh_node[n1].coord[0];
    n1_y = mesh_node[n1].coord[1];

    n2_x = mesh_node[n2].coord[0];
    n2_y = mesh_node[n2].coord[1];


    n3_x = mesh_node[n3].coord[0];
    n3_y = mesh_node[n3].coord[1];


    init_zero_d(bb_coords,8);
    get_bounding_box(E,mesh_element,mesh_node,bb_coords);
    i=0;
    found =0;
    while(!(found) && (i <4))
    {
        init_zero_d(bb_phys_coords,2);
        bb_phys_coords[0] = bb_coords[2*i];
        bb_phys_coords[1] = bb_coords[2*i+1];
        found = check_point_in_voxel(bb_phys_coords,mesh_voxel,V);
        i++;
    }
    if(found)
    {
	    result =1;
    }
    else if(found == 0)
    {
        voxel_coords[0] =mesh_voxel[V].n1[0];
        voxel_coords[1] =mesh_voxel[V].n1[1];
        voxel_coords[2] =mesh_voxel[V].n2[0];
        voxel_coords[3] =mesh_voxel[V].n2[1];
        voxel_coords[4] =mesh_voxel[V].n3[0];
        voxel_coords[5] =mesh_voxel[V].n3[1];
        voxel_coords[6] =mesh_voxel[V].n4[0];
        voxel_coords[7] =mesh_voxel[V].n4[1];
        i=0;
        while(!(found) &&(i<4))
        {
            vv_phys_coords[0] = voxel_coords[2*i];
            vv_phys_coords[1] = voxel_coords[2*i+1];
            found = check_point_in_bounding_box(vv_phys_coords,bb_coords);
            i++;
        }
    }
    if(found)
    {
	    result =1;
    }
    else
    {
	    result =0;
    }

    return result;
}

int check_point_in_voxel(double*  p,
			 voxel* mesh_voxel,
			 int V)
{
    int result =0;

    if(((p[0] >= mesh_voxel[V].n1[0]) && (p[0] <= mesh_voxel[V].n2[0])) && ((p[1] >= mesh_voxel[V].n1[1]) && (p[1] <=mesh_voxel[V].n4[1])))
    {
	    result =1;
    }
    else 
    {
	    result =0;
    }

    return result;
}
int check_point_in_bounding_box(double* p,
				double* bb_coords
		)
{
    int result =0;

    if((p[0] >= bb_coords[0]) && (p[0] <= bb_coords[2]) && ((p[1] >= bb_coords[1]) && (p[1] <=bb_coords[7])))
    {
	    result =1;
    }
    else 
    {
	    result =0;
    }
    return result;
}

int check_point_in_element(int E,
			   element* mesh_element,
			   node* mesh_node,
			   double* v)
{
    double v0[2];
    double v1[2];
    double v2[2];

    double p0[2];
    double p1[2];
    double p2[2];

    int n1=0;
    int n2=0;
    int n3=0;

    double a=0.0;
    double b=0.0;
    int result=0;

    n1 = mesh_element[E].vertex[1];
    n2 = mesh_element[E].vertex[2];
    n3 = mesh_element[E].vertex[3];

    p0[0] = mesh_node[n1].coord[0];
    p0[1] = mesh_node[n1].coord[1];
    p1[0] = mesh_node[n2].coord[0];
    p1[1] = mesh_node[n2].coord[1];
    p2[0] = mesh_node[n3].coord[0];
    p2[1] = mesh_node[n3].coord[1];

    //printf("x_coords %lf %lf %lf \n",p0[0],p1[0],p2[0]);
    //printf("y_coords %lf %lf %lf \n",p0[1],p0[1],p2[1]);

    v0[0] = p0[0];
    v0[1] = p0[1];
    v1[0] = p1[0]-p0[0];
    v1[1] = p1[1]-p0[1];
    v2[0] = p2[0]-p0[0];
    v2[1] = p2[1]-p0[1];

    a = (det(v,v2) - det(v0,v2))/det(v1,v2);
    b = -1.0*(det(v,v1) -det(v0,v1))/det(v1,v2);

    // Allow for points to lie on the boundary
    /*if(((a>0) && (b > 0)) && (a+b <1))*/
    if(((a >= 0) && (b >= 0)) && (a+b <= 1))
    {
	    result =1;
    }
    else
    {
	    result =0;
    }

    return result;
}

double det(double* u,
	   double* v)
{
    double d =0;
    d = u[0]*v[1] - u[1]*v[0];


    return d;
}

void get_bounding_box(int E,
		      element* mesh_element,
		      node* mesh_node,
		      double* bb_coords
		)

{
    int n1 =0;
    int n2 =0;
    int n3 =0;

    double p0[2];
    double p1[2];
    double p2[2];

    double xc_min =0.0;
    double xc_max =0.0;
    double yc_min =0.0;
    double yc_max =0.0;


    n1 = mesh_element[E].vertex[1];
    n2 = mesh_element[E].vertex[2];
    n3 = mesh_element[E].vertex[3];


    p0[0] = mesh_node[n1].coord[0];
    p0[1] = mesh_node[n1].coord[1];
    p1[0] = mesh_node[n2].coord[0];
    p1[1] = mesh_node[n2].coord[1];
    p2[0] = mesh_node[n3].coord[0];
    p2[1] = mesh_node[n3].coord[1];


    xc_min = min(p0[0],p1[0]);
    xc_min = min(xc_min,p2[0]);

    xc_max = max(p0[0],p1[0]);
    xc_max = max(xc_max,p2[0]);

    yc_min = min(p0[1],p1[1]);
    yc_min = min(yc_min,p2[1]);

    yc_max = max(p0[1],p1[1]);
    yc_max = max(yc_max,p2[1]);

    //first coord
    bb_coords[0] = xc_min;
    bb_coords[1] = yc_min;
    //second coord
    bb_coords[2] = xc_max;
    bb_coords[3] = yc_min;
    //third coord
    bb_coords[4] = xc_max;
    bb_coords[5] = yc_max;
    //fourth coord
    bb_coords[6] = xc_min;
    bb_coords[7] = yc_max;
}

double max(double a,
	   double b)
{
    double result=0;
    if(a>=b)
    {
	    result =a;
    }
    else if(b>=a)
    {
	    result =b;
    }
    return result;
}

double min(double a,
	   double b)

{
    double result =0;

    if(a<=b)
    {
	    result =a;
    }
    else if(b<=a)
    {
	    result =b;
    }

    return result;
}
