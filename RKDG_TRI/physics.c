/*
Definition of systems and physics
Philipp Monreal
Prince Chidyagwai
*/
#include"physics.h"

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
void assemble_source(int E,
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

    for(k=1;k<ngpts+1;k++) 
    {
        phys_coords[0] = mesh_element[E].el_gpts_x[k];
	    phys_coords[1] = mesh_element[E].el_gpts_y[k];

	    source_function(t_val,phys_coords,user_param,source_fval);
        

        factor = mesh_element[E].el_gpts_w[k]*det_BE_T*source_fval[s_dim_counter];
        for(idofs=0;idofs<Nloc;idofs++) 
        {
	        source_loc_vec[idofs] = source_loc_vec[idofs] + factor*mesh_element[E].el_gpts_basis[(k-1)*Nloc+idofs];
        }//loop over idofs
    }//loop over quadrature nodes
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
void scattering(int E,
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

    for(k=1;k<ngpts+1;k++) 
    {
        //initialize vectors
        init_sigma_vectors(E,uh_tn_sol_vec,mesh_element,mesh_node,user_param,mesh_element[E].el_gpts_ref_x[k],mesh_element[E].el_gpts_ref_y[k],loc_sigma_s_vec,
                loc_sigma_t_vec);

        //evaluate scattering coeff at the gauss points
        sigma_s = eval_sigma_s(mesh_element[E].el_gpts_x[k],mesh_element[E].el_gpts_y[k],t_val,user_param); // Scattering cross section
        sigma_t = eval_sigma_t(mesh_element[E].el_gpts_x[k],mesh_element[E].el_gpts_y[k],t_val,user_param); // Total (scattering + absorption) cross section

        // Henyey-Greenstein model
        factors = (1.0/(4.0*PI))*(sigma_s*det_BE_T*loc_sigma_s_vec[s_dim_counter]*mesh_element[E].el_density[k]);
        factort = -1.0*sigma_t*det_BE_T*loc_sigma_t_vec[s_dim_counter]*mesh_element[E].el_density[k];

        for(idofs=0;idofs<Nloc;idofs++) 
        {
            sigmas_vec[idofs] = sigmas_vec[idofs] + factors*mesh_element[E].el_gpts_w[k]*mesh_element[E].el_gpts_basis[(k-1)*Nloc+idofs];
            sigmat_vec[idofs] = sigmat_vec[idofs] + factort*mesh_element[E].el_gpts_w[k]*mesh_element[E].el_gpts_basis[(k-1)*Nloc+idofs];
        }//loop over idofs
    }//loop over quadrature nodes
}

/* void init_sigma_vectors
 * initialize structure of scattering coefficients at 
 * (x,y) on element E
 * in: element	    		      -- E
 * in: solution vector                -- uh_tn_sol_vec
 * in: mesh element data structure    -- mesh_element
 * in: mesh_node data structure       -- mesh_node
 * in: parameter structure	      -- user_param
 * in: coordinate x 		      -- ref_coordx
 * in: coordinate y 		      -- ref_coordy
 * out: sigma_s coefficients          -- sigma_s_vec
 * out: sigma_t coefficients          -- sigma_t_vec
 */
void init_sigma_vectors(int E,
		      	double* uh_tn_sol_vec,
		     	element* mesh_element,
		     	node* mesh_node,
		      	parameters* user_param,
		      	double ref_coordx,
		      	double ref_coordy,
		      	double* sigma_s_vec,
		      	double* sigma_t_vec)
{
    double local_approx_sol[SYSTEM_DIM];
    int s_dim;

    init_zero_d(local_approx_sol,SYSTEM_DIM);
    get_approx_solution(uh_tn_sol_vec,E,mesh_element,mesh_node,ref_coordx,ref_coordy,local_approx_sol);

    for (s_dim=0; s_dim < SYSTEM_DIM; s_dim++) 
    {
        sigma_t_vec[s_dim] = local_approx_sol[s_dim];
    }
    
    sigma_s_vec[0] = 4*PI*local_approx_sol[0];
   

    if (SYSTEM_DIM > 2)
    {
        sigma_s_vec[1] = 0.0;
        sigma_s_vec[2] = 0.0;
        if (SYSTEM_DIM > 5) 
        {
            sigma_s_vec[3] = 4.0/3.0*PI*local_approx_sol[0];
            sigma_s_vec[4] = 0.0;
            sigma_s_vec[5] = 4.0/3.0*PI*local_approx_sol[0];

            if (SYSTEM_DIM > 9) 
            {
                for (s_dim=6; s_dim < 10; s_dim++)
                {

                    sigma_s_vec[s_dim] = 0.0;
                }

            }
        }
    }
}	

/* void source_function
 * defines the source function
 * ! MULTIPLY BY 4.0*PI FOR ISOTROPIC SOURCES !
 * in: time	    -- t
 * in: coordinates  -- phys_coords
 * in: parameters   -- user_param
 * out: vector-value of source -- source_fval
 */
void source_function(double t,
		     double* phys_coords,
		     parameters* user_param,
		     double* source_fval)
{
    int s_dim;
    double temp;
    double x,y=0.0;
    x = phys_coords[0];
    y = phys_coords[1];

    switch ((*user_param).test_case) {
        case(4): 
            {
                if (SYSTEM_DIM == 3) 
                {
                    source_fval[0] = ((8.0/3.0)*PI*PI*t-1.0/2.0)*exp(-t)*sin(2*PI*(phys_coords[0]+phys_coords[1]));
                    source_fval[1] = 0.0;
                    source_fval[2] = 0.0;
                }
                break;
            }
        case(5): 
            {
                if (phys_coords[0] >= 3.0 && phys_coords[0] <= 4.0 && phys_coords[1] >= 3.0 && phys_coords[1] <= 4.0)
                {
                    source_fval[0] = 4.0*PI;
                    for (s_dim=1;s_dim<SYSTEM_DIM;s_dim++) 
                    {
                        source_fval[s_dim] = 0.0;
                    }
                    if (SYSTEM_DIM > 2)
                    {
                        source_fval[1] = 0.0;
                        source_fval[2] = 0.0;
                        if (SYSTEM_DIM > 5)
                        {
                            source_fval[3] = 4.0/3.0*PI;
                            source_fval[4] = 0.0;
                            source_fval[5] = 4.0/3.0*PI;

                            for (s_dim=3; s_dim < 6; s_dim++) 
                            {
                                source_fval[s_dim] = 0.0;
                            }
                            if (SYSTEM_DIM > 9)
                            {
                                for (s_dim=6; s_dim < 10; s_dim++) 
                                {
                                    source_fval[s_dim] = 0.0;
                                }
                            }
                        }
                    }
                }
                else 
                {
                    for (s_dim=0;s_dim<SYSTEM_DIM;s_dim++) 
                    {
                        source_fval[s_dim] = 0.0; 
                    }
                }
                break;
            }
        case(7): 
            {
                if (SYSTEM_DIM == 3) 
                {
                    temp = 4.0*PI*((phys_coords[0]-0.5)*(phys_coords[0]-0.5)+(phys_coords[1]-0.5)*(phys_coords[1]-0.5));
                    source_fval[0] = 4.0*PI*exp(-t)*4.0/3.0*t*( PI*((1.0-2*phys_coords[0])*(1.0-2*phys_coords[0])  
                                + (1.0-2*phys_coords[1])*(1.0-2*phys_coords[1])) *sin(temp) - cos(temp)) - 4.0*PI*1.0/(8.0*PI)*exp(-t)*sin(temp);
                    source_fval[1] = 0.0;
                    source_fval[2] = 0.0;
                }
                break;
            }
        case(8): 
            {
                if (SYSTEM_DIM == 3) 
                {
                    temp = 4.0*PI*((phys_coords[0]-0.5)*(phys_coords[0]-0.5)+(phys_coords[1]-0.5)*(phys_coords[1]-0.5));
                    source_fval[0] = -exp(-t)*(16.0/3.0*PI*t*(cos(temp) - PI*((2.0*phys_coords[0] - 1.0)*(2.0*phys_coords[0] - 1.0) + 
                                    (2.0*phys_coords[1] - 1.0)*(2.0*phys_coords[1] - 1.0))*sin(temp)) +
                            1.0/2.0*sin(temp)) + 1.0/2.0*(sin(phys_coords[0]-t) + 3.0*cos(phys_coords[0]-t) - 3.0*sin(phys_coords[1]-t) + cos(phys_coords[1]-t));
                    source_fval[1] = 0.0;
                    source_fval[2] = 0.0;
                }
                break;
            }
        case(9): 
            {

                double x_beam, y_beam, x_width, y_width, t_beam, ints, dir_x, dir_y;
                double x = phys_coords[0];
                double y = phys_coords[1];
                init_zero_d(source_fval,SYSTEM_DIM);

                // Beam 1
                x_beam = 0.1;
                y_beam = 0.5;
                x_width = 0.1;
                y_width = 0.1;
                /*y_width = 1.0;*/
                t_beam = 0.1;
                ints = 1.0;
                dir_x = 1.0;
                dir_y = 0.0;
                /*ints = ints*1e+2;*/

                boundary_beam(phys_coords[0],phys_coords[1],t,x_beam,y_beam,x_width,y_width,t_beam,ints,dir_x,dir_y,source_fval);

                // Beam 2
                x_beam = 0.9;
                y_beam = 0.5;
                x_width = 0.1;
                y_width = 0.1;
                /*y_width = 1.0;*/
                t_beam = 0.1;
                ints = 1.0;
                dir_x = -1.0;
                dir_y = 0.0;
                /*ints = ints*1e+2;*/

                boundary_beam(phys_coords[0],phys_coords[1],t,x_beam,y_beam,x_width,y_width,t_beam,ints,dir_x,dir_y,source_fval);

                // Beam 3
                x_beam = 0.5;
                y_beam = 0.1;
                x_width = 0.1;
                y_width = 0.1;
                /*y_width = 1.0;*/
                t_beam = 0.1;
                ints = 1.0;
                dir_x = 0.0;
                dir_y = 1.0;
                /*ints = ints*1e+2;*/

                //boundary_beam(phys_coords[0],phys_coords[1],t,x_beam,y_beam,x_width,y_width,t_beam,ints,dir_x,dir_y,source_fval);
                break;
            }
        case(10): 
            {
                double x_beam, y_beam, x_width, y_width, t_beam, ints, dir_x, dir_y;
                
                init_zero_d(source_fval,SYSTEM_DIM);

                // Beam 1
                x_beam = 0.1;
                y_beam = 0.5;
                x_width = 0.1;
                y_width = 0.1;
                /*y_width = 1.0;*/
                t_beam = 0.1;
                ints = 1.0;
                dir_x = 1.0;
                dir_y = 0.0;
                /*ints = ints*1e+2;*/

                boundary_beam(phys_coords[0],phys_coords[1],t,x_beam,y_beam,x_width,y_width,t_beam,ints,dir_x,dir_y,source_fval);
                break;
            }
        case(12): 
            {
                if (SYSTEM_DIM == 3) 
                {
                    double x = phys_coords[0];
                    double y = phys_coords[1];
                    // modified Maple output
                    source_fval[0] = -exp(-t) * sin(0.2e1 * 0.3141592654e1 * (double) (x + 2 * y)) / 0.2e1 - 0.6e1 * t * exp(-t) * sin(0.2e1 * 0.3141592654e1 * (double) (x + 2 * y)) * 0.3141592654e1 + 0.3e1 / 0.2e1;



		source_fval[1] = -exp(-t) * (0.756e3 * exp(-0.2e1 * t) * 0.314159265358979323846264338328e1 * sqrt(0.2e1) - 0.756e3 * exp(-0.2e1 * t) * 0.314159265358979323846264338328e1 * sqrt(0.2e1) * pow(cos(0.2e1 * 0.314159265358979323846264338328e1 *(x + 2 * y)), 0.2e1) + 0.14e2 * exp(-0.4e1 * t) * 0.314159265358979323846264338328e1 * sqrt(0.2e1) * pow(cos(0.2e1 * 0.314159265358979323846264338328e1 *(x + 2 * y)), 0.4e1) - 0.28e2 * exp(-0.4e1 * t) * 0.314159265358979323846264338328e1 * sqrt(0.2e1) * pow(cos(0.2e1 * 0.314159265358979323846264338328e1 *(x + 2 * y)), 0.2e1) - 0.16e2 * exp(-0.4e1 * t) * 0.314159265358979323846264338328e1 * sqrt((-0.2e1 * exp(-0.2e1 * t) + 0.2e1 * exp(-0.2e1 * t) * pow(cos(0.2e1 * 0.314159265358979323846264338328e1 *(x + 2 * y)), 0.2e1) - 0.12e2 * exp(-t) * sin(0.2e1 * 0.314159265358979323846264338328e1 *(x + 2 * y)) - 0.18e2 + 0.3e1 * t * t * exp(-0.2e1 * t) * pow(cos(0.2e1 * 0.314159265358979323846264338328e1 *(x + 2 * y)), 0.2e1)) / (-exp(-0.2e1 * t) + exp(-0.2e1 * t) * pow(cos(0.2e1 * 0.314159265358979323846264338328e1 *(x + 2 * y)), 0.2e1) - 0.6e1 * exp(-t) * sin(0.2e1 * 0.314159265358979323846264338328e1 *(x + 2 * y)) - 0.9e1)) * pow(cos(0.2e1 * 0.314159265358979323846264338328e1 *(x + 2 * y)), 0.4e1) + 0.21e2 * sqrt(0.2e1) * 0.314159265358979323846264338328e1 * t * t * exp(-0.4e1 * t) + 0.32e2 * exp(-0.4e1 * t) * 0.314159265358979323846264338328e1 * sqrt((-0.2e1 * exp(-0.2e1 * t) + 0.2e1 * exp(-0.2e1 * t) * pow(cos(0.2e1 * 0.314159265358979323846264338328e1 *(x + 2 * y)), 0.2e1) - 0.12e2 * exp(-t) * sin(0.2e1 * 0.314159265358979323846264338328e1 *(x + 2 * y)) - 0.18e2 + 0.3e1 * t * t * exp(-0.2e1 * t) * pow(cos(0.2e1 * 0.314159265358979323846264338328e1 *(x + 2 * y)), 0.2e1)) / (-exp(-0.2e1 * t) + exp(-0.2e1 * t) * pow(cos(0.2e1 * 0.314159265358979323846264338328e1 *(x + 2 * y)), 0.2e1) - 0.6e1 * exp(-t) * sin(0.2e1 * 0.314159265358979323846264338328e1 *(x + 2 * y)) - 0.9e1)) * pow(cos(0.2e1 * 0.314159265358979323846264338328e1 *(x + 2 * y)), 0.2e1) - 0.1728e4 * exp(-t) * 0.314159265358979323846264338328e1 * sqrt((-0.2e1 * exp(-0.2e1 * t) + 0.2e1 * exp(-0.2e1 * t) * pow(cos(0.2e1 * 0.314159265358979323846264338328e1 *(x + 2 * y)), 0.2e1) - 0.12e2 * exp(-t) * sin(0.2e1 * 0.314159265358979323846264338328e1 *(x + 2 * y)) - 0.18e2 + 0.3e1 * t * t * exp(-0.2e1 * t) * pow(cos(0.2e1 * 0.314159265358979323846264338328e1 *(x + 2 * y)), 0.2e1)) / (-exp(-0.2e1 * t) + exp(-0.2e1 * t) * pow(cos(0.2e1 * 0.314159265358979323846264338328e1 *(x + 2 * y)), 0.2e1) - 0.6e1 * exp(-t) * sin(0.2e1 * 0.314159265358979323846264338328e1 *(x + 2 * y)) - 0.9e1)) * sin(0.2e1 * 0.314159265358979323846264338328e1 *(x + 2 * y)) + 0.864e3 * exp(-0.2e1 * t) * 0.314159265358979323846264338328e1 * sqrt((-0.2e1 * exp(-0.2e1 * t) + 0.2e1 * exp(-0.2e1 * t) * pow(cos(0.2e1 * 0.314159265358979323846264338328e1 *(x + 2 * y)), 0.2e1) - 0.12e2 * exp(-t) * sin(0.2e1 * 0.314159265358979323846264338328e1 *(x + 2 * y)) - 0.18e2 + 0.3e1 * t * t * exp(-0.2e1 * t) * pow(cos(0.2e1 * 0.314159265358979323846264338328e1 *(x + 2 * y)), 0.2e1)) / (-exp(-0.2e1 * t) + exp(-0.2e1 * t) * pow(cos(0.2e1 * 0.314159265358979323846264338328e1 *(x + 2 * y)), 0.2e1) - 0.6e1 * exp(-t) * sin(0.2e1 * 0.314159265358979323846264338328e1 *(x + 2 * y)) - 0.9e1)) * pow(cos(0.2e1 * 0.314159265358979323846264338328e1 *(x + 2 * y)), 0.2e1) + 0.567e3 * exp(-0.2e1 * t) * sqrt(0.2e1) * 0.314159265358979323846264338328e1 * t * t - 0.1296e4 * 0.314159265358979323846264338328e1 * sqrt((-0.2e1 * exp(-0.2e1 * t) + 0.2e1 * exp(-0.2e1 * t) * pow(cos(0.2e1 * 0.314159265358979323846264338328e1 *(x + 2 * y)), 0.2e1) - 0.12e2 * exp(-t) * sin(0.2e1 * 0.314159265358979323846264338328e1 *(x + 2 * y)) - 0.18e2 + 0.3e1 * t * t * exp(-0.2e1 * t) * pow(cos(0.2e1 * 0.314159265358979323846264338328e1 *(x + 2 * y)), 0.2e1)) / (-exp(-0.2e1 * t) + exp(-0.2e1 * t) * pow(cos(0.2e1 * 0.314159265358979323846264338328e1 *(x + 2 * y)), 0.2e1) - 0.6e1 * exp(-t) * sin(0.2e1 * 0.314159265358979323846264338328e1 *(x + 2 * y)) - 0.9e1)) + 0.162e3 * exp(-0.2e1 * t) * sqrt((-0.2e1 * exp(-0.2e1 * t) + 0.2e1 * exp(-0.2e1 * t) * pow(cos(0.2e1 * 0.314159265358979323846264338328e1 *(x + 2 * y)), 0.2e1) - 0.12e2 * exp(-t) * sin(0.2e1 * 0.314159265358979323846264338328e1 *(x + 2 * y)) - 0.18e2 + 0.3e1 * t * t * exp(-0.2e1 * t) * pow(cos(0.2e1 * 0.314159265358979323846264338328e1 *(x + 2 * y)), 0.2e1)) / (-exp(-0.2e1 * t) + exp(-0.2e1 * t) * pow(cos(0.2e1 * 0.314159265358979323846264338328e1 *(x + 2 * y)), 0.2e1) - 0.6e1 * exp(-t) * sin(0.2e1 * 0.314159265358979323846264338328e1 *(x + 2 * y)) - 0.9e1)) * pow(cos(0.2e1 * 0.314159265358979323846264338328e1 *(x + 2 * y)), 0.2e1) + 0.14e2 * exp(-0.4e1 * t) * 0.314159265358979323846264338328e1 * sqrt(0.2e1) + 0.567e3 * sqrt(0.2e1) * 0.314159265358979323846264338328e1 * t * t * exp(-t) * sin(0.2e1 * 0.314159265358979323846264338328e1 *(x + 2 * y)) - 0.324e3 * exp(-t) * sqrt((-0.2e1 * exp(-0.2e1 * t) + 0.2e1 * exp(-0.2e1 * t) * pow(cos(0.2e1 * 0.314159265358979323846264338328e1 *(x + 2 * y)), 0.2e1) - 0.12e2 * exp(-t) * sin(0.2e1 * 0.314159265358979323846264338328e1 *(x + 2 * y)) - 0.18e2 + 0.3e1 * t * t * exp(-0.2e1 * t) * pow(cos(0.2e1 * 0.314159265358979323846264338328e1 *(x + 2 * y)), 0.2e1)) / (-exp(-0.2e1 * t) + exp(-0.2e1 * t) * pow(cos(0.2e1 * 0.314159265358979323846264338328e1 *(x + 2 * y)), 0.2e1) - 0.6e1 * exp(-t) * sin(0.2e1 * 0.314159265358979323846264338328e1 *(x + 2 * y)) - 0.9e1)) * sin(0.2e1 * 0.314159265358979323846264338328e1 *(x + 2 * y)) - 0.243e3 * sqrt((-0.2e1 * exp(-0.2e1 * t) + 0.2e1 * exp(-0.2e1 * t) * pow(cos(0.2e1 * 0.314159265358979323846264338328e1 *(x + 2 * y)), 0.2e1) - 0.12e2 * exp(-t) * sin(0.2e1 * 0.314159265358979323846264338328e1 *(x + 2 * y)) - 0.18e2 + 0.3e1 * t * t * exp(-0.2e1 * t) * pow(cos(0.2e1 * 0.314159265358979323846264338328e1 *(x + 2 * y)), 0.2e1)) / (-exp(-0.2e1 * t) + exp(-0.2e1 * t) * pow(cos(0.2e1 * 0.314159265358979323846264338328e1 *(x + 2 * y)), 0.2e1) - 0.6e1 * exp(-t) * sin(0.2e1 * 0.314159265358979323846264338328e1 *(x + 2 * y)) - 0.9e1)) - 0.162e3 * exp(-0.2e1 * t) * sqrt((-0.2e1 * exp(-0.2e1 * t) + 0.2e1 * exp(-0.2e1 * t) * pow(cos(0.2e1 * 0.314159265358979323846264338328e1 *(x + 2 * y)), 0.2e1) - 0.12e2 * exp(-t) * sin(0.2e1 * 0.314159265358979323846264338328e1 *(x + 2 * y)) - 0.18e2 + 0.3e1 * t * t * exp(-0.2e1 * t) * pow(cos(0.2e1 * 0.314159265358979323846264338328e1 *(x + 2 * y)), 0.2e1)) / (-exp(-0.2e1 * t) + exp(-0.2e1 * t) * pow(cos(0.2e1 * 0.314159265358979323846264338328e1 *(x + 2 * y)), 0.2e1) - 0.6e1 * exp(-t) * sin(0.2e1 * 0.314159265358979323846264338328e1 *(x + 2 * y)) - 0.9e1)) - 0.189e3 * sqrt(0.2e1) * 0.314159265358979323846264338328e1 * t * t * exp(-0.3e1 * t) * pow(cos(0.2e1 * 0.314159265358979323846264338328e1 *(x + 2 * y)), 0.2e1) * sin(0.2e1 * 0.314159265358979323846264338328e1 *(x + 2 * y)) - 0.864e3 * exp(-0.2e1 * t) * 0.314159265358979323846264338328e1 * sqrt((-0.2e1 * exp(-0.2e1 * t) + 0.2e1 * exp(-0.2e1 * t) * pow(cos(0.2e1 * 0.314159265358979323846264338328e1 *(x + 2 * y)), 0.2e1) - 0.12e2 * exp(-t) * sin(0.2e1 * 0.314159265358979323846264338328e1 *(x + 2 * y)) - 0.18e2 + 0.3e1 * t * t * exp(-0.2e1 * t) * pow(cos(0.2e1 * 0.314159265358979323846264338328e1 *(x + 2 * y)), 0.2e1)) / (-exp(-0.2e1 * t) + exp(-0.2e1 * t) * pow(cos(0.2e1 * 0.314159265358979323846264338328e1 *(x + 2 * y)), 0.2e1) - 0.6e1 * exp(-t) * sin(0.2e1 * 0.314159265358979323846264338328e1 *(x + 2 * y)) - 0.9e1)) + 0.168e3 * exp(-0.3e1 * t) * 0.314159265358979323846264338328e1 * sqrt(0.2e1) * sin(0.2e1 * 0.314159265358979323846264338328e1 *(x + 2 * y)) + 0.36e2 * exp(-0.3e1 * t) * sqrt((-0.2e1 * exp(-0.2e1 * t) + 0.2e1 * exp(-0.2e1 * t) * pow(cos(0.2e1 * 0.314159265358979323846264338328e1 *(x + 2 * y)), 0.2e1) - 0.12e2 * exp(-t) * sin(0.2e1 * 0.314159265358979323846264338328e1 *(x + 2 * y)) - 0.18e2 + 0.3e1 * t * t * exp(-0.2e1 * t) * pow(cos(0.2e1 * 0.314159265358979323846264338328e1 *(x + 2 * y)), 0.2e1)) / (-exp(-0.2e1 * t) + exp(-0.2e1 * t) * pow(cos(0.2e1 * 0.314159265358979323846264338328e1 *(x + 2 * y)), 0.2e1) - 0.6e1 * exp(-t) * sin(0.2e1 * 0.314159265358979323846264338328e1 *(x + 2 * y)) - 0.9e1)) * pow(cos(0.2e1 * 0.314159265358979323846264338328e1 *(x + 2 * y)), 0.2e1) * sin(0.2e1 * 0.314159265358979323846264338328e1 *(x + 2 * y)) - 0.192e3 * exp(-0.3e1 * t) * 0.314159265358979323846264338328e1 * sqrt((-0.2e1 * exp(-0.2e1 * t) + 0.2e1 * exp(-0.2e1 * t) * pow(cos(0.2e1 * 0.314159265358979323846264338328e1 *(x + 2 * y)), 0.2e1) - 0.12e2 * exp(-t) * sin(0.2e1 * 0.314159265358979323846264338328e1 *(x + 2 * y)) - 0.18e2 + 0.3e1 * t * t * exp(-0.2e1 * t) * pow(cos(0.2e1 * 0.314159265358979323846264338328e1 *(x + 2 * y)), 0.2e1)) / (-exp(-0.2e1 * t) + exp(-0.2e1 * t) * pow(cos(0.2e1 * 0.314159265358979323846264338328e1 *(x + 2 * y)), 0.2e1) - 0.6e1 * exp(-t) * sin(0.2e1 * 0.314159265358979323846264338328e1 *(x + 2 * y)) - 0.9e1)) * sin(0.2e1 * 0.314159265358979323846264338328e1 *(x + 2 * y)) + 0.1134e4 * 0.314159265358979323846264338328e1 * sqrt(0.2e1) - 0.42e2 * sqrt(0.2e1) * 0.314159265358979323846264338328e1 * t * t * exp(-0.4e1 * t) * pow(cos(0.2e1 * 0.314159265358979323846264338328e1 *(x + 2 * y)), 0.2e1) + 0.6e1 * exp(-0.4e1 * t) * sqrt((-0.2e1 * exp(-0.2e1 * t) + 0.2e1 * exp(-0.2e1 * t) * pow(cos(0.2e1 * 0.314159265358979323846264338328e1 *(x + 2 * y)), 0.2e1) - 0.12e2 * exp(-t) * sin(0.2e1 * 0.314159265358979323846264338328e1 *(x + 2 * y)) - 0.18e2 + 0.3e1 * t * t * exp(-0.2e1 * t) * pow(cos(0.2e1 * 0.314159265358979323846264338328e1 *(x + 2 * y)), 0.2e1)) / (-exp(-0.2e1 * t) + exp(-0.2e1 * t) * pow(cos(0.2e1 * 0.314159265358979323846264338328e1 *(x + 2 * y)), 0.2e1) - 0.6e1 * exp(-t) * sin(0.2e1 * 0.314159265358979323846264338328e1 *(x + 2 * y)) - 0.9e1)) * pow(cos(0.2e1 * 0.314159265358979323846264338328e1 *(x + 2 * y)), 0.2e1) - 0.3e1 * exp(-0.4e1 * t) * sqrt((-0.2e1 * exp(-0.2e1 * t) + 0.2e1 * exp(-0.2e1 * t) * pow(cos(0.2e1 * 0.314159265358979323846264338328e1 *(x + 2 * y)), 0.2e1) - 0.12e2 * exp(-t) * sin(0.2e1 * 0.314159265358979323846264338328e1 *(x + 2 * y)) - 0.18e2 + 0.3e1 * t * t * exp(-0.2e1 * t) * pow(cos(0.2e1 * 0.314159265358979323846264338328e1 *(x + 2 * y)), 0.2e1)) / (-exp(-0.2e1 * t) + exp(-0.2e1 * t) * pow(cos(0.2e1 * 0.314159265358979323846264338328e1 *(x + 2 * y)), 0.2e1) - 0.6e1 * exp(-t) * sin(0.2e1 * 0.314159265358979323846264338328e1 *(x + 2 * y)) - 0.9e1)) * pow(cos(0.2e1 * 0.314159265358979323846264338328e1 *(x + 2 * y)), 0.4e1) - 0.16e2 * exp(-0.4e1 * t) * 0.314159265358979323846264338328e1 * sqrt((-0.2e1 * exp(-0.2e1 * t) + 0.2e1 * exp(-0.2e1 * t) * pow(cos(0.2e1 * 0.314159265358979323846264338328e1 *(x + 2 * y)), 0.2e1) - 0.12e2 * exp(-t) * sin(0.2e1 * 0.314159265358979323846264338328e1 *(x + 2 * y)) - 0.18e2 + 0.3e1 * t * t * exp(-0.2e1 * t) * pow(cos(0.2e1 * 0.314159265358979323846264338328e1 *(x + 2 * y)), 0.2e1)) / (-exp(-0.2e1 * t) + exp(-0.2e1 * t) * pow(cos(0.2e1 * 0.314159265358979323846264338328e1 *(x + 2 * y)), 0.2e1) - 0.6e1 * exp(-t) * sin(0.2e1 * 0.314159265358979323846264338328e1 *(x + 2 * y)) - 0.9e1)) - 0.3e1 * exp(-0.4e1 * t) * sqrt((-0.2e1 * exp(-0.2e1 * t) + 0.2e1 * exp(-0.2e1 * t) * pow(cos(0.2e1 * 0.314159265358979323846264338328e1 *(x + 2 * y)), 0.2e1) - 0.12e2 * exp(-t) * sin(0.2e1 * 0.314159265358979323846264338328e1 *(x + 2 * y)) - 0.18e2 + 0.3e1 * t * t * exp(-0.2e1 * t) * pow(cos(0.2e1 * 0.314159265358979323846264338328e1 *(x + 2 * y)), 0.2e1)) / (-exp(-0.2e1 * t) + exp(-0.2e1 * t) * pow(cos(0.2e1 * 0.314159265358979323846264338328e1 *(x + 2 * y)), 0.2e1) - 0.6e1 * exp(-t) * sin(0.2e1 * 0.314159265358979323846264338328e1 *(x + 2 * y)) - 0.9e1)) - 0.36e2 * exp(-0.3e1 * t) * sqrt((-0.2e1 * exp(-0.2e1 * t) + 0.2e1 * exp(-0.2e1 * t) * pow(cos(0.2e1 * 0.314159265358979323846264338328e1 *(x + 2 * y)), 0.2e1) - 0.12e2 * exp(-t) * sin(0.2e1 * 0.314159265358979323846264338328e1 *(x + 2 * y)) - 0.18e2 + 0.3e1 * t * t * exp(-0.2e1 * t) * pow(cos(0.2e1 * 0.314159265358979323846264338328e1 *(x + 2 * y)), 0.2e1)) / (-exp(-0.2e1 * t) + exp(-0.2e1 * t) * pow(cos(0.2e1 * 0.314159265358979323846264338328e1 *(x + 2 * y)), 0.2e1) - 0.6e1 * exp(-t) * sin(0.2e1 * 0.314159265358979323846264338328e1 *(x + 2 * y)) - 0.9e1)) * sin(0.2e1 * 0.314159265358979323846264338328e1 *(x + 2 * y)) + 0.1512e4 * exp(-t) * 0.314159265358979323846264338328e1 * sqrt(0.2e1) * sin(0.2e1 * 0.314159265358979323846264338328e1 *(x + 2 * y)) + 0.21e2 * exp(-0.4e1 * t) * 0.314159265358979323846264338328e1 * sqrt(0.2e1) * pow(cos(0.2e1 * 0.314159265358979323846264338328e1 *(x + 2 * y)), 0.4e1) * t * t - 0.567e3 * sqrt(0.2e1) * 0.314159265358979323846264338328e1 * t * t * exp(-0.2e1 * t) * pow(cos(0.2e1 * 0.314159265358979323846264338328e1 *(x + 2 * y)), 0.2e1) + 0.189e3 * sqrt(0.2e1) * 0.314159265358979323846264338328e1 * t * t * exp(-0.3e1 * t) * sin(0.2e1 * 0.314159265358979323846264338328e1 *(x + 2 * y)) + 0.192e3 * exp(-0.3e1 * t) * 0.314159265358979323846264338328e1 * sqrt((-0.2e1 * exp(-0.2e1 * t) + 0.2e1 * exp(-0.2e1 * t) * pow(cos(0.2e1 * 0.314159265358979323846264338328e1 *(x + 2 * y)), 0.2e1) - 0.12e2 * exp(-t) * sin(0.2e1 * 0.314159265358979323846264338328e1 *(x + 2 * y)) - 0.18e2 + 0.3e1 * t * t * exp(-0.2e1 * t) * pow(cos(0.2e1 * 0.314159265358979323846264338328e1 *(x + 2 * y)), 0.2e1)) / (-exp(-0.2e1 * t) + exp(-0.2e1 * t) * pow(cos(0.2e1 * 0.314159265358979323846264338328e1 *(x + 2 * y)), 0.2e1) - 0.6e1 * exp(-t) * sin(0.2e1 * 0.314159265358979323846264338328e1 *(x + 2 * y)) - 0.9e1)) * pow(cos(0.2e1 * 0.314159265358979323846264338328e1 *(x + 2 * y)), 0.2e1) * sin(0.2e1 * 0.314159265358979323846264338328e1 *(x + 2 * y)) - 0.168e3 * exp(-0.3e1 * t) * 0.314159265358979323846264338328e1 * sqrt(0.2e1) * pow(cos(0.2e1 * 0.314159265358979323846264338328e1 *(x + 2 * y)), 0.2e1) * sin(0.2e1 * 0.314159265358979323846264338328e1 *(x + 2 * y))) * cos(0.2e1 * 0.314159265358979323846264338328e1 *(x + 2 * y)) * pow((-0.2e1 * exp(-0.2e1 * t) + 0.2e1 * exp(-0.2e1 * t) * pow(cos(0.2e1 * 0.314159265358979323846264338328e1 *(x + 2 * y)), 0.2e1) - 0.12e2 * exp(-t) * sin(0.2e1 * 0.314159265358979323846264338328e1 *(x + 2 * y)) - 0.18e2 + 0.3e1 * t * t * exp(-0.2e1 * t) * pow(cos(0.2e1 * 0.314159265358979323846264338328e1 *(x + 2 * y)), 0.2e1)) / (-exp(-0.2e1 * t) + exp(-0.2e1 * t) * pow(cos(0.2e1 * 0.314159265358979323846264338328e1 *(x + 2 * y)), 0.2e1) - 0.6e1 * exp(-t) * sin(0.2e1 * 0.314159265358979323846264338328e1 *(x + 2 * y)) - 0.9e1), -0.1e1 / 0.2e1) / (exp(-0.4e1 * t) - 0.2e1 * pow(cos(0.2e1 * 0.314159265358979323846264338328e1 *(x + 2 * y)), 0.2e1) * exp(-0.4e1 * t) + 0.12e2 * exp(-0.3e1 * t) * sin(0.2e1 * 0.314159265358979323846264338328e1 *(x + 2 * y)) + 0.54e2 * exp(-0.2e1 * t) + pow(cos(0.2e1 * 0.314159265358979323846264338328e1 *(x + 2 * y)), 0.4e1) * exp(-0.4e1 * t) - 0.12e2 * pow(cos(0.2e1 * 0.314159265358979323846264338328e1 *(x + 2 * y)), 0.2e1) * exp(-0.3e1 * t) * sin(0.2e1 * 0.314159265358979323846264338328e1 *(x + 2 * y)) - 0.54e2 * exp(-0.2e1 * t) * pow(cos(0.2e1 * 0.314159265358979323846264338328e1 *(x + 2 * y)), 0.2e1) + 0.108e3 * exp(-t) * sin(0.2e1 * 0.314159265358979323846264338328e1 *(x + 2 * y)) + 0.81e2) / 0.3e1;

		source_fval[2] = -exp(-t) * (-0.6e1 * exp(-0.4e1 * t) * sqrt((-0.2e1 * exp(-0.2e1 * t) + 0.2e1 * exp(-0.2e1 * t) * pow(cos(0.2e1 * 0.314159265358979323846264338328e1 * (double) (x + 2 * y)), 0.2e1) - 0.12e2 * exp(-t) * sin(0.2e1 * 0.314159265358979323846264338328e1 * (double) (x + 2 * y)) - 0.18e2 + 0.3e1 * t * t * exp(-0.2e1 * t) * pow(cos(0.2e1 * 0.314159265358979323846264338328e1 * (double) (x + 2 * y)), 0.2e1)) / (-exp(-0.2e1 * t) + exp(-0.2e1 * t) * pow(cos(0.2e1 * 0.314159265358979323846264338328e1 * (double) (x + 2 * y)), 0.2e1) - 0.6e1 * exp(-t) * sin(0.2e1 * 0.314159265358979323846264338328e1 * (double) (x + 2 * y)) - 0.9e1)) * pow(cos(0.2e1 * 0.314159265358979323846264338328e1 * (double) (x + 2 * y)), 0.2e1) + 0.243e3 * sqrt((-0.2e1 * exp(-0.2e1 * t) + 0.2e1 * exp(-0.2e1 * t) * pow(cos(0.2e1 * 0.314159265358979323846264338328e1 * (double) (x + 2 * y)), 0.2e1) - 0.12e2 * exp(-t) * sin(0.2e1 * 0.314159265358979323846264338328e1 * (double) (x + 2 * y)) - 0.18e2 + 0.3e1 * t * t * exp(-0.2e1 * t) * pow(cos(0.2e1 * 0.314159265358979323846264338328e1 * (double) (x + 2 * y)), 0.2e1)) / (-exp(-0.2e1 * t) + exp(-0.2e1 * t) * pow(cos(0.2e1 * 0.314159265358979323846264338328e1 * (double) (x + 2 * y)), 0.2e1) - 0.6e1 * exp(-t) * sin(0.2e1 * 0.314159265358979323846264338328e1 * (double) (x + 2 * y)) - 0.9e1)) + 0.756e3 * exp(-0.2e1 * t) * 0.314159265358979323846264338328e1 * sqrt((-0.2e1 * exp(-0.2e1 * t) + 0.2e1 * exp(-0.2e1 * t) * pow(cos(0.2e1 * 0.314159265358979323846264338328e1 * (double) (x + 2 * y)), 0.2e1) - 0.12e2 * exp(-t) * sin(0.2e1 * 0.314159265358979323846264338328e1 * (double) (x + 2 * y)) - 0.18e2 + 0.3e1 * t * t * exp(-0.2e1 * t) * pow(cos(0.2e1 * 0.314159265358979323846264338328e1 * (double) (x + 2 * y)), 0.2e1)) / (-exp(-0.2e1 * t) + exp(-0.2e1 * t) * pow(cos(0.2e1 * 0.314159265358979323846264338328e1 * (double) (x + 2 * y)), 0.2e1) - 0.6e1 * exp(-t) * sin(0.2e1 * 0.314159265358979323846264338328e1 * (double) (x + 2 * y)) - 0.9e1)) + 0.162e3 * sqrt((-0.2e1 * exp(-0.2e1 * t) + 0.2e1 * exp(-0.2e1 * t) * pow(cos(0.2e1 * 0.314159265358979323846264338328e1 * (double) (x + 2 * y)), 0.2e1) - 0.12e2 * exp(-t) * sin(0.2e1 * 0.314159265358979323846264338328e1 * (double) (x + 2 * y)) - 0.18e2 + 0.3e1 * t * t * exp(-0.2e1 * t) * pow(cos(0.2e1 * 0.314159265358979323846264338328e1 * (double) (x + 2 * y)), 0.2e1)) / (-exp(-0.2e1 * t) + exp(-0.2e1 * t) * pow(cos(0.2e1 * 0.314159265358979323846264338328e1 * (double) (x + 2 * y)), 0.2e1) - 0.6e1 * exp(-t) * sin(0.2e1 * 0.314159265358979323846264338328e1 * (double) (x + 2 * y)) - 0.9e1)) * exp(-0.2e1 * t) - 0.540e3 * exp(-0.2e1 * t) * 0.314159265358979323846264338328e1 * sqrt(0.2e1) + 0.36e2 * exp(-0.3e1 * t) * sqrt((-0.2e1 * exp(-0.2e1 * t) + 0.2e1 * exp(-0.2e1 * t) * pow(cos(0.2e1 * 0.314159265358979323846264338328e1 * (double) (x + 2 * y)), 0.2e1) - 0.12e2 * exp(-t) * sin(0.2e1 * 0.314159265358979323846264338328e1 * (double) (x + 2 * y)) - 0.18e2 + 0.3e1 * t * t * exp(-0.2e1 * t) * pow(cos(0.2e1 * 0.314159265358979323846264338328e1 * (double) (x + 2 * y)), 0.2e1)) / (-exp(-0.2e1 * t) + exp(-0.2e1 * t) * pow(cos(0.2e1 * 0.314159265358979323846264338328e1 * (double) (x + 2 * y)), 0.2e1) - 0.6e1 * exp(-t) * sin(0.2e1 * 0.314159265358979323846264338328e1 * (double) (x + 2 * y)) - 0.9e1)) * sin(0.2e1 * 0.314159265358979323846264338328e1 * (double) (x + 2 * y)) - 0.162e3 * exp(-0.2e1 * t) * sqrt((-0.2e1 * exp(-0.2e1 * t) + 0.2e1 * exp(-0.2e1 * t) * pow(cos(0.2e1 * 0.314159265358979323846264338328e1 * (double) (x + 2 * y)), 0.2e1) - 0.12e2 * exp(-t) * sin(0.2e1 * 0.314159265358979323846264338328e1 * (double) (x + 2 * y)) - 0.18e2 + 0.3e1 * t * t * exp(-0.2e1 * t) * pow(cos(0.2e1 * 0.314159265358979323846264338328e1 * (double) (x + 2 * y)), 0.2e1)) / (-exp(-0.2e1 * t) + exp(-0.2e1 * t) * pow(cos(0.2e1 * 0.314159265358979323846264338328e1 * (double) (x + 2 * y)), 0.2e1) - 0.6e1 * exp(-t) * sin(0.2e1 * 0.314159265358979323846264338328e1 * (double) (x + 2 * y)) - 0.9e1)) * pow(cos(0.2e1 * 0.314159265358979323846264338328e1 * (double) (x + 2 * y)), 0.2e1) + 0.135e3 * sqrt(0.2e1) * 0.314159265358979323846264338328e1 * t * t * exp(-0.3e1 * t) * pow(cos(0.2e1 * 0.314159265358979323846264338328e1 * (double) (x + 2 * y)), 0.2e1) * sin(0.2e1 * 0.314159265358979323846264338328e1 * (double) (x + 2 * y)) - 0.810e3 * sqrt(0.2e1) * 0.314159265358979323846264338328e1 + 0.1134e4 * 0.314159265358979323846264338328e1 * sqrt((-0.2e1 * exp(-0.2e1 * t) + 0.2e1 * exp(-0.2e1 * t) * pow(cos(0.2e1 * 0.314159265358979323846264338328e1 * (double) (x + 2 * y)), 0.2e1) - 0.12e2 * exp(-t) * sin(0.2e1 * 0.314159265358979323846264338328e1 * (double) (x + 2 * y)) - 0.18e2 + 0.3e1 * t * t * exp(-0.2e1 * t) * pow(cos(0.2e1 * 0.314159265358979323846264338328e1 * (double) (x + 2 * y)), 0.2e1)) / (-exp(-0.2e1 * t) + exp(-0.2e1 * t) * pow(cos(0.2e1 * 0.314159265358979323846264338328e1 * (double) (x + 2 * y)), 0.2e1) - 0.6e1 * exp(-t) * sin(0.2e1 * 0.314159265358979323846264338328e1 * (double) (x + 2 * y)) - 0.9e1)) - 0.405e3 * sqrt(0.2e1) * 0.314159265358979323846264338328e1 * t * t * exp(-t) * sin(0.2e1 * 0.314159265358979323846264338328e1 * (double) (x + 2 * y)) + 0.405e3 * sqrt(0.2e1) * 0.314159265358979323846264338328e1 * t * t * exp(-0.2e1 * t) * pow(cos(0.2e1 * 0.314159265358979323846264338328e1 * (double) (x + 2 * y)), 0.2e1) - 0.1080e4 * exp(-t) * 0.314159265358979323846264338328e1 * sqrt(0.2e1) * sin(0.2e1 * 0.314159265358979323846264338328e1 * (double) (x + 2 * y)) - 0.10e2 * exp(-0.4e1 * t) * 0.314159265358979323846264338328e1 * sqrt(0.2e1) + 0.324e3 * exp(-t) * sqrt((-0.2e1 * exp(-0.2e1 * t) + 0.2e1 * exp(-0.2e1 * t) * pow(cos(0.2e1 * 0.314159265358979323846264338328e1 * (double) (x + 2 * y)), 0.2e1) - 0.12e2 * exp(-t) * sin(0.2e1 * 0.314159265358979323846264338328e1 * (double) (x + 2 * y)) - 0.18e2 + 0.3e1 * t * t * exp(-0.2e1 * t) * pow(cos(0.2e1 * 0.314159265358979323846264338328e1 * (double) (x + 2 * y)), 0.2e1)) / (-exp(-0.2e1 * t) + exp(-0.2e1 * t) * pow(cos(0.2e1 * 0.314159265358979323846264338328e1 * (double) (x + 2 * y)), 0.2e1) - 0.6e1 * exp(-t) * sin(0.2e1 * 0.314159265358979323846264338328e1 * (double) (x + 2 * y)) - 0.9e1)) * sin(0.2e1 * 0.314159265358979323846264338328e1 * (double) (x + 2 * y)) + 0.14e2 * exp(-0.4e1 * t) * 0.314159265358979323846264338328e1 * sqrt((-0.2e1 * exp(-0.2e1 * t) + 0.2e1 * exp(-0.2e1 * t) * pow(cos(0.2e1 * 0.314159265358979323846264338328e1 * (double) (x + 2 * y)), 0.2e1) - 0.12e2 * exp(-t) * sin(0.2e1 * 0.314159265358979323846264338328e1 * (double) (x + 2 * y)) - 0.18e2 + 0.3e1 * t * t * exp(-0.2e1 * t) * pow(cos(0.2e1 * 0.314159265358979323846264338328e1 * (double) (x + 2 * y)), 0.2e1)) / (-exp(-0.2e1 * t) + exp(-0.2e1 * t) * pow(cos(0.2e1 * 0.314159265358979323846264338328e1 * (double) (x + 2 * y)), 0.2e1) - 0.6e1 * exp(-t) * sin(0.2e1 * 0.314159265358979323846264338328e1 * (double) (x + 2 * y)) - 0.9e1)) + 0.3e1 * exp(-0.4e1 * t) * sqrt((-0.2e1 * exp(-0.2e1 * t) + 0.2e1 * exp(-0.2e1 * t) * pow(cos(0.2e1 * 0.314159265358979323846264338328e1 * (double) (x + 2 * y)), 0.2e1) - 0.12e2 * exp(-t) * sin(0.2e1 * 0.314159265358979323846264338328e1 * (double) (x + 2 * y)) - 0.18e2 + 0.3e1 * t * t * exp(-0.2e1 * t) * pow(cos(0.2e1 * 0.314159265358979323846264338328e1 * (double) (x + 2 * y)), 0.2e1)) / (-exp(-0.2e1 * t) + exp(-0.2e1 * t) * pow(cos(0.2e1 * 0.314159265358979323846264338328e1 * (double) (x + 2 * y)), 0.2e1) - 0.6e1 * exp(-t) * sin(0.2e1 * 0.314159265358979323846264338328e1 * (double) (x + 2 * y)) - 0.9e1)) * pow(cos(0.2e1 * 0.314159265358979323846264338328e1 * (double) (x + 2 * y)), 0.4e1) + 0.3e1 * exp(-0.4e1 * t) * sqrt((-0.2e1 * exp(-0.2e1 * t) + 0.2e1 * exp(-0.2e1 * t) * pow(cos(0.2e1 * 0.314159265358979323846264338328e1 * (double) (x + 2 * y)), 0.2e1) - 0.12e2 * exp(-t) * sin(0.2e1 * 0.314159265358979323846264338328e1 * (double) (x + 2 * y)) - 0.18e2 + 0.3e1 * t * t * exp(-0.2e1 * t) * pow(cos(0.2e1 * 0.314159265358979323846264338328e1 * (double) (x + 2 * y)), 0.2e1)) / (-exp(-0.2e1 * t) + exp(-0.2e1 * t) * pow(cos(0.2e1 * 0.314159265358979323846264338328e1 * (double) (x + 2 * y)), 0.2e1) - 0.6e1 * exp(-t) * sin(0.2e1 * 0.314159265358979323846264338328e1 * (double) (x + 2 * y)) - 0.9e1)) - 0.10e2 * exp(-0.4e1 * t) * 0.314159265358979323846264338328e1 * sqrt(0.2e1) * pow(cos(0.2e1 * 0.314159265358979323846264338328e1 * (double) (x + 2 * y)), 0.4e1) - 0.28e2 * exp(-0.4e1 * t) * 0.314159265358979323846264338328e1 * sqrt((-0.2e1 * exp(-0.2e1 * t) + 0.2e1 * exp(-0.2e1 * t) * pow(cos(0.2e1 * 0.314159265358979323846264338328e1 * (double) (x + 2 * y)), 0.2e1) - 0.12e2 * exp(-t) * sin(0.2e1 * 0.314159265358979323846264338328e1 * (double) (x + 2 * y)) - 0.18e2 + 0.3e1 * t * t * exp(-0.2e1 * t) * pow(cos(0.2e1 * 0.314159265358979323846264338328e1 * (double) (x + 2 * y)), 0.2e1)) / (-exp(-0.2e1 * t) + exp(-0.2e1 * t) * pow(cos(0.2e1 * 0.314159265358979323846264338328e1 * (double) (x + 2 * y)), 0.2e1) - 0.6e1 * exp(-t) * sin(0.2e1 * 0.314159265358979323846264338328e1 * (double) (x + 2 * y)) - 0.9e1)) * pow(cos(0.2e1 * 0.314159265358979323846264338328e1 * (double) (x + 2 * y)), 0.2e1) + 0.14e2 * exp(-0.4e1 * t) * 0.314159265358979323846264338328e1 * sqrt((-0.2e1 * exp(-0.2e1 * t) + 0.2e1 * exp(-0.2e1 * t) * pow(cos(0.2e1 * 0.314159265358979323846264338328e1 * (double) (x + 2 * y)), 0.2e1) - 0.12e2 * exp(-t) * sin(0.2e1 * 0.314159265358979323846264338328e1 * (double) (x + 2 * y)) - 0.18e2 + 0.3e1 * t * t * exp(-0.2e1 * t) * pow(cos(0.2e1 * 0.314159265358979323846264338328e1 * (double) (x + 2 * y)), 0.2e1)) / (-exp(-0.2e1 * t) + exp(-0.2e1 * t) * pow(cos(0.2e1 * 0.314159265358979323846264338328e1 * (double) (x + 2 * y)), 0.2e1) - 0.6e1 * exp(-t) * sin(0.2e1 * 0.314159265358979323846264338328e1 * (double) (x + 2 * y)) - 0.9e1)) * pow(cos(0.2e1 * 0.314159265358979323846264338328e1 * (double) (x + 2 * y)), 0.4e1) + 0.540e3 * exp(-0.2e1 * t) * 0.314159265358979323846264338328e1 * sqrt(0.2e1) * pow(cos(0.2e1 * 0.314159265358979323846264338328e1 * (double) (x + 2 * y)), 0.2e1) - 0.405e3 * sqrt(0.2e1) * 0.314159265358979323846264338328e1 * t * t * exp(-0.2e1 * t) + 0.1512e4 * exp(-t) * 0.314159265358979323846264338328e1 * sqrt((-0.2e1 * exp(-0.2e1 * t) + 0.2e1 * exp(-0.2e1 * t) * pow(cos(0.2e1 * 0.314159265358979323846264338328e1 * (double) (x + 2 * y)), 0.2e1) - 0.12e2 * exp(-t) * sin(0.2e1 * 0.314159265358979323846264338328e1 * (double) (x + 2 * y)) - 0.18e2 + 0.3e1 * t * t * exp(-0.2e1 * t) * pow(cos(0.2e1 * 0.314159265358979323846264338328e1 * (double) (x + 2 * y)), 0.2e1)) / (-exp(-0.2e1 * t) + exp(-0.2e1 * t) * pow(cos(0.2e1 * 0.314159265358979323846264338328e1 * (double) (x + 2 * y)), 0.2e1) - 0.6e1 * exp(-t) * sin(0.2e1 * 0.314159265358979323846264338328e1 * (double) (x + 2 * y)) - 0.9e1)) * sin(0.2e1 * 0.314159265358979323846264338328e1 * (double) (x + 2 * y)) - 0.756e3 * exp(-0.2e1 * t) * 0.314159265358979323846264338328e1 * sqrt((-0.2e1 * exp(-0.2e1 * t) + 0.2e1 * exp(-0.2e1 * t) * pow(cos(0.2e1 * 0.314159265358979323846264338328e1 * (double) (x + 2 * y)), 0.2e1) - 0.12e2 * exp(-t) * sin(0.2e1 * 0.314159265358979323846264338328e1 * (double) (x + 2 * y)) - 0.18e2 + 0.3e1 * t * t * exp(-0.2e1 * t) * pow(cos(0.2e1 * 0.314159265358979323846264338328e1 * (double) (x + 2 * y)), 0.2e1)) / (-exp(-0.2e1 * t) + exp(-0.2e1 * t) * pow(cos(0.2e1 * 0.314159265358979323846264338328e1 * (double) (x + 2 * y)), 0.2e1) - 0.6e1 * exp(-t) * sin(0.2e1 * 0.314159265358979323846264338328e1 * (double) (x + 2 * y)) - 0.9e1)) * pow(cos(0.2e1 * 0.314159265358979323846264338328e1 * (double) (x + 2 * y)), 0.2e1) - 0.120e3 * exp(-0.3e1 * t) * 0.314159265358979323846264338328e1 * sqrt(0.2e1) * sin(0.2e1 * 0.314159265358979323846264338328e1 * (double) (x + 2 * y)) - 0.36e2 * exp(-0.3e1 * t) * sqrt((-0.2e1 * exp(-0.2e1 * t) + 0.2e1 * exp(-0.2e1 * t) * pow(cos(0.2e1 * 0.314159265358979323846264338328e1 * (double) (x + 2 * y)), 0.2e1) - 0.12e2 * exp(-t) * sin(0.2e1 * 0.314159265358979323846264338328e1 * (double) (x + 2 * y)) - 0.18e2 + 0.3e1 * t * t * exp(-0.2e1 * t) * pow(cos(0.2e1 * 0.314159265358979323846264338328e1 * (double) (x + 2 * y)), 0.2e1)) / (-exp(-0.2e1 * t) + exp(-0.2e1 * t) * pow(cos(0.2e1 * 0.314159265358979323846264338328e1 * (double) (x + 2 * y)), 0.2e1) - 0.6e1 * exp(-t) * sin(0.2e1 * 0.314159265358979323846264338328e1 * (double) (x + 2 * y)) - 0.9e1)) * pow(cos(0.2e1 * 0.314159265358979323846264338328e1 * (double) (x + 2 * y)), 0.2e1) * sin(0.2e1 * 0.314159265358979323846264338328e1 * (double) (x + 2 * y)) + 0.168e3 * exp(-0.3e1 * t) * 0.314159265358979323846264338328e1 * sqrt((-0.2e1 * exp(-0.2e1 * t) + 0.2e1 * exp(-0.2e1 * t) * pow(cos(0.2e1 * 0.314159265358979323846264338328e1 * (double) (x + 2 * y)), 0.2e1) - 0.12e2 * exp(-t) * sin(0.2e1 * 0.314159265358979323846264338328e1 * (double) (x + 2 * y)) - 0.18e2 + 0.3e1 * t * t * exp(-0.2e1 * t) * pow(cos(0.2e1 * 0.314159265358979323846264338328e1 * (double) (x + 2 * y)), 0.2e1)) / (-exp(-0.2e1 * t) + exp(-0.2e1 * t) * pow(cos(0.2e1 * 0.314159265358979323846264338328e1 * (double) (x + 2 * y)), 0.2e1) - 0.6e1 * exp(-t) * sin(0.2e1 * 0.314159265358979323846264338328e1 * (double) (x + 2 * y)) - 0.9e1)) * sin(0.2e1 * 0.314159265358979323846264338328e1 * (double) (x + 2 * y)) + 0.20e2 * exp(-0.4e1 * t) * 0.314159265358979323846264338328e1 * sqrt(0.2e1) * pow(cos(0.2e1 * 0.314159265358979323846264338328e1 * (double) (x + 2 * y)), 0.2e1) - 0.15e2 * sqrt(0.2e1) * 0.314159265358979323846264338328e1 * t * t * exp(-0.4e1 * t) - 0.135e3 * sqrt(0.2e1) * 0.314159265358979323846264338328e1 * t * t * exp(-0.3e1 * t) * sin(0.2e1 * 0.314159265358979323846264338328e1 * (double) (x + 2 * y)) + 0.120e3 * exp(-0.3e1 * t) * 0.314159265358979323846264338328e1 * sqrt(0.2e1) * pow(cos(0.2e1 * 0.314159265358979323846264338328e1 * (double) (x + 2 * y)), 0.2e1) * sin(0.2e1 * 0.314159265358979323846264338328e1 * (double) (x + 2 * y)) - 0.168e3 * exp(-0.3e1 * t) * 0.314159265358979323846264338328e1 * sqrt((-0.2e1 * exp(-0.2e1 * t) + 0.2e1 * exp(-0.2e1 * t) * pow(cos(0.2e1 * 0.314159265358979323846264338328e1 * (double) (x + 2 * y)), 0.2e1) - 0.12e2 * exp(-t) * sin(0.2e1 * 0.314159265358979323846264338328e1 * (double) (x + 2 * y)) - 0.18e2 + 0.3e1 * t * t * exp(-0.2e1 * t) * pow(cos(0.2e1 * 0.314159265358979323846264338328e1 * (double) (x + 2 * y)), 0.2e1)) / (-exp(-0.2e1 * t) + exp(-0.2e1 * t) * pow(cos(0.2e1 * 0.314159265358979323846264338328e1 * (double) (x + 2 * y)), 0.2e1) - 0.6e1 * exp(-t) * sin(0.2e1 * 0.314159265358979323846264338328e1 * (double) (x + 2 * y)) - 0.9e1)) * pow(cos(0.2e1 * 0.314159265358979323846264338328e1 * (double) (x + 2 * y)), 0.2e1) * sin(0.2e1 * 0.314159265358979323846264338328e1 * (double) (x + 2 * y)) - 0.15e2 * exp(-0.4e1 * t) * 0.314159265358979323846264338328e1 * sqrt(0.2e1) * pow(cos(0.2e1 * 0.314159265358979323846264338328e1 * (double) (x + 2 * y)), 0.4e1) * t * t + 0.30e2 * sqrt(0.2e1) * 0.314159265358979323846264338328e1 * t * t * exp(-0.4e1 * t) * pow(cos(0.2e1 * 0.314159265358979323846264338328e1 * (double) (x + 2 * y)), 0.2e1)) * cos(0.2e1 * 0.314159265358979323846264338328e1 * (double) (x + 2 * y)) * pow((-0.2e1 * exp(-0.2e1 * t) + 0.2e1 * exp(-0.2e1 * t) * pow(cos(0.2e1 * 0.314159265358979323846264338328e1 * (double) (x + 2 * y)), 0.2e1) - 0.12e2 * exp(-t) * sin(0.2e1 * 0.314159265358979323846264338328e1 * (double) (x + 2 * y)) - 0.18e2 + 0.3e1 * t * t * exp(-0.2e1 * t) * pow(cos(0.2e1 * 0.314159265358979323846264338328e1 * (double) (x + 2 * y)), 0.2e1)) / (-exp(-0.2e1 * t) + exp(-0.2e1 * t) * pow(cos(0.2e1 * 0.314159265358979323846264338328e1 * (double) (x + 2 * y)), 0.2e1) - 0.6e1 * exp(-t) * sin(0.2e1 * 0.314159265358979323846264338328e1 * (double) (x + 2 * y)) - 0.9e1), -0.1e1 / 0.2e1) / (-exp(-0.4e1 * t) + 0.2e1 * pow(cos(0.2e1 * 0.314159265358979323846264338328e1 * (double) (x + 2 * y)), 0.2e1) * exp(-0.4e1 * t) - 0.12e2 * exp(-0.3e1 * t) * sin(0.2e1 * 0.314159265358979323846264338328e1 * (double) (x + 2 * y)) - 0.54e2 * exp(-0.2e1 * t) - pow(cos(0.2e1 * 0.314159265358979323846264338328e1 * (double) (x + 2 * y)), 0.4e1) * exp(-0.4e1 * t) + 0.12e2 * pow(cos(0.2e1 * 0.314159265358979323846264338328e1 * (double) (x + 2 * y)), 0.2e1) * exp(-0.3e1 * t) * sin(0.2e1 * 0.314159265358979323846264338328e1 * (double) (x + 2 * y)) + 0.54e2 * exp(-0.2e1 * t) * pow(cos(0.2e1 * 0.314159265358979323846264338328e1 * (double) (x + 2 * y)), 0.2e1) - 0.108e3 * exp(-t) * sin(0.2e1 * 0.314159265358979323846264338328e1 * (double) (x + 2 * y)) - 0.81e2) / 0.3e1;
	    }
	    break;
            }
        case(14): 
            {
                if (SYSTEM_DIM == 3) 
                {
                    double x = phys_coords[0];
                    double y = phys_coords[1];
                    source_fval[0] = ((20.0/3.0)*PI*PI*t-1.0/2.0)*exp(-t)*sin(2*PI*(phys_coords[0]+2.*phys_coords[1]));

		            source_fval[1] = 0.0;
		            source_fval[2] = 0.0;
                }
                break;
            }
        case(16): 
            {
	    if (SYSTEM_DIM == 3) {
	        double x = phys_coords[0];
	        double y = phys_coords[1];

		source_fval[0] = -exp(-t) * sin(0.2e1 * 0.314159265358979323846264338328e1 * (double) (x + 2 * y)) / 0.2e1 - 0.6e1 * t * exp(-t) * sin(0.2e1 * 0.314159265358979323846264338328e1 * (double) (x + 2 * y)) * 0.314159265358979323846264338328e1 + 0.3e1 / 0.2e1;

		/*source_fval[0] = -(0.48e2 * t * exp(-t) * sin(0.2e1 * 0.314159265358979323846264338328e1 * (double) (x + 2 * y)) * 0.314159265358979323846264338328e1 * PI - 0.24e2 * PI + exp(-t) * sin(0.2e1 * 0.314159265358979323846264338328e1 * (double) (x + 2 * y)) + 0.3e1) / PI / 0.8e1;*/

		double t1 = exp(-t);
		double t5 = 0.2e1 * 0.314159265358979323846264338328e1 * (double) (x + 2 * y);
		double t6 = cos(t5);
		double t8 = t * t;
		double t9 = t1 * t1;
		double t10 = t8 * t9;
		double t11 = t10 * t6;
		double t12 = sin(t5);
		double t14 = t1 * t12 + 0.3e1;
		double t15 = t14 * t14;
		double t16 = 0.1e1 / t15;
		double t21 = t8 * t9 * t1;
		double t22 = t6 * t6;
		double t23 = t22 * t6;
		source_fval[1] = t1 * t6 + 0.8e1 / 0.3e1 * (-t11 * t16 * t12 * 0.314159265358979323846264338328e1 - t21 * t23 / t15 / t14 * 0.314159265358979323846264338328e1) * t14 + 0.2e1 * (0.1e1 / 0.3e1 + 0.2e1 / 0.3e1 * t10 * t22 * t16) * t1 * t6 * 0.314159265358979323846264338328e1 - 0.16e2 / 0.3e1 * t11 / t14 * t12 * 0.314159265358979323846264338328e1 - 0.8e1 / 0.3e1 * t21 * t23 * t16 * 0.314159265358979323846264338328e1;

		double r1 = exp(-t);
		double r5 = 0.2e1 * 0.314159265358979323846264338328e1 * (double) (x + 2 * y);
		double r6 = cos(r5);
		double r8 = t * t;
		double r9 = r1 * r1;
		double r10 = r8 * r9;
		double r11 = r10 * r6;
		double r12 = sin(r5);
		double r14 = r1 * r12 + 0.3e1;
		double r21 = r8 * r9 * r1;
		double r22 = r6 * r6;
		double r23 = r22 * r6;
		double r24 = r14 * r14;
		double r25 = 0.1e1 / r24;
		source_fval[2] = r1 * r6 - 0.8e1 / 0.3e1 * r11 / r14 * r12 * 0.314159265358979323846264338328e1 - 0.4e1 / 0.3e1 * r21 * r23 * r25 * 0.314159265358979323846264338328e1 + 0.16e2 / 0.3e1 * (-r11 * r25 * r12 * 0.314159265358979323846264338328e1 - r21 * r23 / r24 / r14 * 0.314159265358979323846264338328e1) * r14 + 0.4e1 * (0.1e1 / 0.3e1 + 0.2e1 / 0.3e1 * r10 * r22 * r25) * r1 * r6 * 0.314159265358979323846264338328e1;
	    }
	    break;
	}

        case(18): 
	{
	    if (SYSTEM_DIM == 3) {
		double x = phys_coords[0];
		double y = phys_coords[1];

		if ( (((x > 1.0 && x < 2.0) || (x > 5.0 && x < 6.0)) &&
		    ((y > 1.0 && y < 2.0) || (y > 3.0 && y < 4.0) || (y > 5.0 && y < 6.0))) || 
		     (((x > 2.0 && x < 3.0) || (x > 4.0 && x < 5.0)) &&
		    ((y > 2.0 && y < 3.0) || (y > 4.0 && y < 5.0))) || 
		    ((x > 3.0 && x < 4.0) && (y > 1.0 && y < 2.0)) ) {

		    source_fval[0] = 0.9e1 * exp(-t) * sin(0.2e1 * 0.314159265358979323846264338328e1 * (double) (x + 2 * y)) - 0.6e1 * t * exp(-t) * sin(0.2e1 * 0.314159265358979323846264338328e1 * (double) (x + 2 * y)) * 0.314159265358979323846264338328e1 + 0.30e2;

		    double t1 = exp(-t);
		    double t6 = cos(0.2e1 * 0.314159265358979323846264338328e1 * (double) (x + 2 * y));
		    double t7 = t1 * t6;
		    source_fval[1] = t7 + 0.9e1 * t * t1 * t6 + 0.2e1 / 0.3e1 * t7 * 0.314159265358979323846264338328e1;

		    double r1 = exp(-t);
		    double r6 = cos(0.2e1 * 0.314159265358979323846264338328e1 * (double) (x + 2 * y));
		    double r7 = r1 * r6;
		    source_fval[2] = r7 + 0.9e1 * t * r1 * r6 + 0.4e1 / 0.3e1 * r7 * 0.314159265358979323846264338328e1;
		}
		else {
		    source_fval[0] = -exp(-t) * sin(0.2e1 * 0.314159265358979323846264338328e1 * (double) (x + 2 * y)) - 0.6e1 * t * exp(-t) * sin(0.2e1 * 0.314159265358979323846264338328e1 * (double) (x + 2 * y)) * 0.314159265358979323846264338328e1;

		    double t1 = exp(-t);
		    double t6 = cos(0.2e1 * 0.314159265358979323846264338328e1 * (double) (x + 2 * y));
		    double t7 = t1 * t6;
		    source_fval[1] = t7 + 0.2e1 / 0.3e1 * t7 * 0.314159265358979323846264338328e1;

		    double r1 = exp(-t);
		    double r6 = cos(0.2e1 * 0.314159265358979323846264338328e1 * (double) (x + 2 * y));
		    double r7 = r1 * r6;
		    source_fval[2] = r7 + 0.4e1 / 0.3e1 * r7 * 0.314159265358979323846264338328e1;
		}
	    }
	    break;
	}

        case(19): 
	{
	    if (SYSTEM_DIM == 3) {
		double x = phys_coords[0];
		double y = phys_coords[1];

		if ((x > 0.25 && x < 0.5) && (y > 0.25 && y < 0.5)) {

		    source_fval[0] = 0.9e1 * exp(-t) * sin(0.2e1 * 0.314159265358979323846264338328e1 * (double) (x + 2 * y)) - 0.6e1 * t * exp(-t) * sin(0.2e1 * 0.314159265358979323846264338328e1 * (double) (x + 2 * y)) * 0.314159265358979323846264338328e1 + 0.30e2;

		    double t1 = exp(-t);
		    double t6 = cos(0.2e1 * 0.314159265358979323846264338328e1 * (double) (x + 2 * y));
		    double t7 = t1 * t6;
		    source_fval[1] = t7 + 0.9e1 * t * t1 * t6 + 0.2e1 / 0.3e1 * t7 * 0.314159265358979323846264338328e1;

		    double r1 = exp(-t);
		    double r6 = cos(0.2e1 * 0.314159265358979323846264338328e1 * (double) (x + 2 * y));
		    double r7 = r1 * r6;
		    source_fval[2] = r7 + 0.9e1 * t * r1 * r6 + 0.4e1 / 0.3e1 * r7 * 0.314159265358979323846264338328e1;
		}
		else {
		    source_fval[0] = -exp(-t) * sin(0.2e1 * 0.314159265358979323846264338328e1 * (double) (x + 2 * y)) - 0.6e1 * t * exp(-t) * sin(0.2e1 * 0.314159265358979323846264338328e1 * (double) (x + 2 * y)) * 0.314159265358979323846264338328e1;

		    double t1 = exp(-t);
		    double t6 = cos(0.2e1 * 0.314159265358979323846264338328e1 * (double) (x + 2 * y));
		    double t7 = t1 * t6;
		    source_fval[1] = t7 + 0.2e1 / 0.3e1 * t7 * 0.314159265358979323846264338328e1;

		    double r1 = exp(-t);
		    double r6 = cos(0.2e1 * 0.314159265358979323846264338328e1 * (double) (x + 2 * y));
		    double r7 = r1 * r6;
		    source_fval[2] = r7 + 0.4e1 / 0.3e1 * r7 * 0.314159265358979323846264338328e1;
		}
	    }
	    break;
	}

	case(20): {
	    if (phys_coords[0] >= 3.0 && phys_coords[0] <= 4.0 && phys_coords[1] >= 3.0 && phys_coords[1] <= 4.0)
	    {
		source_fval[0] = 4.*PI*pow(sin(PI*(phys_coords[0]-3.)),4.)*pow(sin(PI*(phys_coords[1]-3.)),4.);
	    }
	    else
	    {
		source_fval[0] = 0.;
	    }
		    
	    for (s_dim=1;s_dim<SYSTEM_DIM;s_dim++) 
		source_fval[s_dim] = 0.0; 

	    break;
	}

        case(21): {
	    if (SYSTEM_DIM == 3) {
		double x = phys_coords[0];
		source_fval[1] = 0.0;
		source_fval[2] = 0.0;

		source_fval[0] = 0.0;

		if ((x >= 0. && x < 2./16.)) {
		    source_fval[0] = 0.0;
		}
		if ((x >= 2./16. && x < 3./16.)) {
		    source_fval[0] = 1.0;
		}
		if ((x >= 3./16. && x < 5./16.)) {
		    source_fval[0] = 0.0;
		}
		if ((x >= 5./16. && x < 6./16.)) {
		    source_fval[0] = 0.0;
		}
		if ((x >= 6./16. && x < 8./16.)) {
		    source_fval[0] = 50.0;
		}
		if ((x >= 8./16. && x <= 10./16.)) {
		    source_fval[0] = 50.0;
		}
		if ((x > 10./16. && x <= 11./16.)) {
		    source_fval[0] = 0.0;
		}
		if ((x > 11./16. && x <= 13./16.)) {
		    source_fval[0] = 0.0;
		}
		if ((x > 13./16. && x <= 14./16.)) {
		    source_fval[0] = 1.0;
		}
		if ((x > 14./16. && x <= 1.)) {
		    source_fval[0] = 0.0;
		}
		source_fval[0] = source_fval[0]*2.0;
	    }
	    break;
	}

        case(22): {
	    double x_beam, y_beam, x_width, y_width, t_beam, ints, dir_x, dir_y;
	    double x = phys_coords[0];
	    double y = phys_coords[1];

	    if (SYSTEM_DIM == 3) {
		source_fval[0] = 0.0;
		source_fval[1] = 0.0;
		source_fval[2] = 0.0;

		// Beam 1
		x_beam = 0.1;
		y_beam = 0.0;
		x_width = 0.05;
		y_width = 0.01;
		t_beam = t;
		ints = 1.0;
		dir_x = 0.9;
		dir_y = 1.0;

		boundary_beam(x,y,t,x_beam,y_beam,x_width,y_width,t_beam,ints,dir_x,dir_y,source_fval);

		// Beam 2
		x_beam = 0.0;
		y_beam = 0.1;
		x_width = 0.01;
		y_width = 0.05;
		t_beam = t;
		ints = 1.0;
		dir_x = 1.0;
		dir_y = 0.9;

		boundary_beam(x,y,t,x_beam,y_beam,x_width,y_width,t_beam,ints,dir_x,dir_y,source_fval);

		// Beam 3
		x_beam = 1.0;
		y_beam = 0.9;
		x_width = 0.01;
		y_width = 0.05;
		t_beam = t;
		ints = 1.0;
		dir_x = -1.0;
		dir_y = -0.9;

		boundary_beam(x,y,t,x_beam,y_beam,x_width,y_width,t_beam,ints,dir_x,dir_y,source_fval);

		// Beam 4
		x_beam = 0.9;
		y_beam = 1.0;
		x_width = 0.05;
		y_width = 0.01;
		t_beam = t;
		ints = 1.0;
		dir_x = -0.9;
		dir_y = -1.0;

		boundary_beam(x,y,t,x_beam,y_beam,x_width,y_width,t_beam,ints,dir_x,dir_y,source_fval);
	    }
	    break;
	}
        case(24):
                  {
                      //if((pow(x,2.0) + pow(y,2.0))<1.0)
                      //if(((x*x>0.0) && (x*x<2.0)) && ((y*y>0.0)&&(y*y<2.0)))
                          
                     
                      if((pow(x,2.0) + pow(y,2.0)) <=1.0)
                      {
                          source_fval[0] = 1.0;
                          source_fval[1] = 0.0;
                          source_fval[2] = 0.0;
                       
                      }
                      else
                      {
                          source_fval[0] = 0.0;
                          source_fval[1] = 0.0;
                          source_fval[2] = 0.0;
                      }
                      break;
                  }
        case(25):
                  {
                      source_fval[0] = 0.0;
                      source_fval[1] = 0.0;
                      source_fval[2] = 0.0;
                      break;
                  }
        case(26):
                  {
                      source_fval[0] = 0.0;
                      source_fval[1] = 0.0;
                      source_fval[2] = 0.0;
                      break;
                  }
        default: {
	    for (s_dim=0;s_dim<SYSTEM_DIM;s_dim++) source_fval[s_dim] = 0.0;
	}
    }
}

/* double eval_sigma_s
 * evaluates scattering cross section
 * in: x coordinate -- x
 * in: y coordinate -- y
 * in: time	    -- t
 * in: parameters   -- user_param
 * out: sigma_s at x,y,t
 */

double eval_sigma_s(double x,
		    double y,
		    double t,
		    parameters* user_param)
{
    double sigma_s;

    switch ((*user_param).test_case) 
    {
        case(4):
            {
                sigma_s = 0.5;
                break;
            }
        case(5): 
            {
                // Checkerboard test-case
	            sigma_s = 1.0; 

                // First and fifth column 
                if (((x >= 1.0 && x <= 2.0) || (x >= 5.0 && x <= 6.0)) &&
                    ((y >= 1.0 && y <= 2.0) || (y >= 3.0 && y <= 4.0) || (y >= 5.0 && y <= 6.0))) 
                {
                     sigma_s = 0.0;
                }
                // Second and fourth column 
                if (((x >= 2.0 && x <= 3.0) || (x >= 4.0 && x <= 5.0)) &&
                    ((y >= 2.0 && y <=3.0) || (y >= 4.0 && y <= 5.0))) 
                {
                    sigma_s = 0.0;
                }
                // Third column 
                if ((x >= 3.0 && x <= 4.0) && (y >= 1.0 && y <= 2.0)) 
                {
                    sigma_s = 0.0;
	            }
                break;
            }
        case(6):
            {
                sigma_s = 0.0;
                break;
            }
        case(7):
            {
                sigma_s = 0.5;
                break;
            }
        case(8):
            {
                sigma_s = 0.5;
	            break;
            }
        case(9):
            {
                sigma_s = 1.0/40.0*1.0/10.0;
                break;
            }
        case(10):
            {
                sigma_s = 0.1;
                break;
            }
        case(11):
            {
                sigma_s = 0.0;
                break;
            }
        case(12):
            {
                sigma_s = 0.5;
                break;
            }
        case(13):
            {
                sigma_s = 0.5;
                break;
            }
        case(14):
            {
                sigma_s = 0.5;
                break;
            }
        case(15):
            {
                sigma_s = 0.0;
                break;
            }
        case(16):
            {
                sigma_s = 0.5;
                break;
            }
        case(17):
            {
                sigma_s = 0.0;
                break;
            }
        case(18):
            {
                sigma_s = 1.0; 

                // First and fifth column 
                if (((x > 1.0 && x < 2.0) || (x > 5.0 && x < 6.0)) &&
                    ((y > 1.0 && y < 2.0) || (y > 3.0 && y < 4.0) || (y > 5.0 && y < 6.0))) 
                {
                    sigma_s = 0.0;
                }
                // Second and fourth column 
                if (((x > 2.0 && x < 3.0) || (x > 4.0 && x < 5.0)) &&
                    ((y > 2.0 && y < 3.0) || (y > 4.0 && y < 5.0))) 
                {
                    sigma_s = 0.0;
                }
                // Third column 
                if ((x > 3.0 && x < 4.0) && (y > 1.0 && y < 2.0)) 
                {
                    sigma_s = 0.0;
                }
                break;
            }
        case(19):
            {
                sigma_s = 1.0; 
                if ((x > 0.25 && x < 0.5) && (y > 0.25 && y < 0.5)) 
                {
                    sigma_s = 0.0;
                }
                break;
            }
        case(20): 
            {
                // Checkerboard test-case
                sigma_s = 1.0; 

                // First and fifth column 
                if (((x > 1.0 && x < 2.0) || (x > 5.0 && x < 6.0)) &&
                    ((y > 1.0 && y < 2.0) || (y > 3.0 && y < 4.0) || (y > 5.0 && y < 6.0))) 
                {
                    sigma_s = 0.0;
                }
                // Second and fourth column 
                if (((x > 2.0 && x < 3.0) || (x > 4.0 && x < 5.0)) &&
                    ((y > 2.0 && y < 3.0) || (y > 4.0 && y < 5.0))) 
                {
                    sigma_s = 0.0;
                }
                // Third column 
                if ((x > 3.0 && x < 4.0) && (y > 1.0 && y < 2.0)) 
                {
                    sigma_s = 0.0;
                }
                break;
            }
        case(21): 
            {
                sigma_s = 0.0; 

                if ((x >= 0. && x < 2./16.)) {
                sigma_s = 0.9;
                }
                if ((x >= 2./16. && x < 3./16.)) {
                sigma_s = 0.1;
                }
                if ((x >= 3./16. && x < 5./16.)) {
                sigma_s = 0.0;
                }
                if ((x >= 5./16. && x < 6./16.)) {
                sigma_s = 0.0;
                }
                if ((x >= 6./16. && x < 8./16.)) {
                sigma_s = 0.0;
                }
                if ((x >= 8./16. && x <= 10./16.)) {
                sigma_s = 0.0;
                }
                if ((x > 10./16. && x <= 11./16.)) {
                sigma_s = 0.0;
                }
                if ((x > 11./16. && x <= 13./16.)) {
                sigma_s = 0.0;
                }
                if ((x > 13./16. && x <= 14./16.)) {
                sigma_s = 0.1;
                }
                if ((x > 14./16. && x <= 1.)) {
                sigma_s = 0.9;
                }
                /*sigma_s = sigma_s*4.0*PI;*/
                break;
            }
        case(22): 
            {
                sigma_s = 1.0; 
                double dist = sqrt( (x - 0.5)*(x - 0.5) + (y - 0.5)*(y - 0.5));
                if ( dist < 0.3 ) 
                {
                    if ( dist > 0.2 ) 
                    {
                        sigma_s = 0.2;
                    }
                }
                break;
            }
        case(24):
            {
                if((pow(x,2.0) + pow(y,2.0)) <= 1.0)
                {
                    sigma_s = 0.0;
                }
                else
                {
                    sigma_s = 0.0;
                }
                break;
            }
        case(25):
            {
                sigma_s =0.0;
                break;
            }
        case(26):
            {
                sigma_s =0.0;
                break;
            }
        case(27):
            {
                sigma_s = 0.01;
                break;
            }
        default: 
            {
                sigma_s = 1.0;
                break;
            }
    }
    return (sigma_s);
}

/* double eval_sigma_t
 * evaluates total cross section
 * in: x coordinate -- x
 * in: y coordinate -- y
 * in: time	    -- t
 * in: parameters   -- user_param
 * out: sigma_t at x,y,t
 */
double eval_sigma_t(double x,
		    double y,
		    double t,
		    parameters* user_param)
{
    double sigma_t;

    switch ((*user_param).test_case) 
    {
        case(5): 
            {
                // Checkerboard test-case
                sigma_t = 1.0; 
                // First and fifth column 
                if (((x > 1.0 && x < 2.0) || (x > 5.0 && x < 6.0)) &&
                ((y > 1.0 && y < 2.0) || (y > 3.0 && y < 4.0) || (y > 5.0 && y < 6.0))) 
                {
                    sigma_t = 10.0;
                }
                // Second and fourth column
                if (((x > 2.0 && x < 3.0) || (x > 4.0 && x < 5.0)) &&
                        ((y > 2.0 && y < 3.0) || (y > 4.0 && y < 5.0))) 
                {
                    sigma_t = 10.0;
                }
                // Third column 
	            if ((x > 3.0 && x < 4.0) && (y > 1.0 && y < 2.0)) 
                {
                    sigma_t = 10.0;
                }
                break;
            }
        case(6):
            {
                sigma_t = 0.0;
                break;
            }
        case(9):
            {
                sigma_t = 1.0/40.0; // sigma_t * length(domain) = 1/2
	            break;
            }
        case(10):
            {
                sigma_t = 0.2;
                break;
            }
        case(11):
            {
                sigma_t = 0.0;
                break;
            }
        case(12):
            {
                sigma_t = 1.0;
                break;
            }
        case(13):
            {
                sigma_t = 0.5;
                break;
            }
        case(14):
            {
                sigma_t = 1.0;
                break;
            }
        case(15):
            {
                if (x >= 0.5) 
                {
                    sigma_t = 1.0;
                }
                else 
                {
                    sigma_t = 0.0;
                }
                break;
            }
        case(16):
            {
                sigma_t = 1.0;
                break;
            }
        case(17):
            {
                sigma_t = 0.0;
                break;
            }
        case(18): 
            {
                sigma_t = 1.0; 

                // First and fifth column 
                if (((x > 1.0 && x < 2.0) || (x > 5.0 && x < 6.0)) &&
                    ((y > 1.0 && y < 2.0) || (y > 3.0 && y < 4.0) || (y > 5.0 && y < 6.0))) 
                {
                    sigma_t = 10.0;
                }
                // Second and fourth column 
                if (((x > 2.0 && x < 3.0) || (x > 4.0 && x < 5.0)) &&
                    ((y > 2.0 && y < 3.0) || (y > 4.0 && y < 5.0))) 
                {
                    sigma_t = 10.0;
                }
                // Third column 
                if ((x > 3.0 && x < 4.0) && (y > 1.0 && y < 2.0)) 
                {
                    sigma_t = 10.0;
                }
                break;
            }
        case(19): 
            {
                sigma_t = 1.0; 
                if ((x > 0.25 && x < 0.5) && (y > 0.25 && y < 0.5)) 
                {
                    sigma_t = 10.0;
                }
                break;
            }
        case(20): 
            {
                // Checkerboard test-case
                sigma_t = 1.0; 

                // First and fifth column 
                if (((x > 1.0 && x < 2.0) || (x > 5.0 && x < 6.0)) &&
                    ((y > 1.0 && y < 2.0) || (y > 3.0 && y < 4.0) || (y > 5.0 && y < 6.0)))
                {
                    sigma_t = 10.0;
                }
                // Second and fourth column 
                if (((x > 2.0 && x < 3.0) || (x > 4.0 && x < 5.0)) &&
                    ((y > 2.0 && y < 3.0) || (y > 4.0 && y < 5.0))) 
                {
                    sigma_t = 10.0;
                }
                // Third column 
                if ((x > 3.0 && x < 4.0) && (y > 1.0 && y < 2.0)) 
                {
                    sigma_t = 10.0;
                }
                break;
            }
        case(21): 
            {
                sigma_t = 0.0; 

                if ((x >= 0. && x < 2./16.)) 
                {
                    sigma_t = 1.0;
                }
                if ((x >= 2./16. && x < 3./16.)) 
                {
                    sigma_t = 0.2;
                }
                if ((x >= 3./16. && x < 5./16.)) 
                {
                    sigma_t = 0.0;
                }
                if ((x >= 5./16. && x < 6./16.)) 
                {
                    sigma_t = 5.0;
                }
                if ((x >= 6./16. && x < 8./16.)) 
                {
                    sigma_t = 50.0;
                }
                if ((x >= 8./16. && x <= 10./16.)) 
                {
                    sigma_t = 50.0;
                }
                if ((x > 10./16. && x <= 11./16.)) 
                {
                    sigma_t = 5.0;
                }
                if ((x > 11./16. && x <= 13./16.)) 
                {
                    sigma_t = 0.0;
                }
                if ((x > 13./16. && x <= 14./16.)) 
                {
                    sigma_t = 0.2;
                }
                if ((x > 14./16. && x <= 1.)) 
                {
                    sigma_t = 1.0;
                }
                break;
            }
        case(22):
            {
                sigma_t = 1.0; 
                double dist = sqrt( (x - 0.5)*(x - 0.5) + (y - 0.5)*(y - 0.5));
                if ( dist < 0.3 ) 
                {
                    if ( dist > 0.2 ) 
                    {
                        sigma_t = 0.2;
                    }
                }
                break;
            }
        case(24):
            {
                if((pow(x,2.0) + pow(y,2.0)) <=1.0)
                {
                    sigma_t = 10.0;
                  
                }
                else
                {
                    sigma_t =0.0;
                }
                break;
            }
        case(25):
            {
                sigma_t =0.0;
                break;
            }
        case(26):
            {
                if(((x>=2.0)&&(x<=3.0))&&((y>=0.0)&&(y<=2.0)))
                {
                    sigma_t =50.0;
                    break;
                }
            }
        case(27):
            {
                sigma_t = 0.01;
                break;
            }
        default: 
            {
                sigma_t = 1.0;
                break;
            }
    }

    return (sigma_t);
}

