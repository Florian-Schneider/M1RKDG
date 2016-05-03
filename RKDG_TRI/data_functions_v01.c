/* data_functions.c
defines initial and boundary conditions and the flux function
Prince Chidyagwai
Philipp Monreal

Test cases:
  1 : Constant initial and boundary conditions
  2 : Initial Gaussian 
  3 : Initial plane source
  4 : Manufactured solutions: constant along along x+y = c
  5 : Checkerboard
  6 : Linesource
  7 : Manufactured solutions: P1 closure; radially symmetric
  8 : Manufactured solutions: P1 closure; radially symmetric with time- and space-dependent flux
  9 : 3 beamS
  10 : 1 beam
  11 : Quasi 1D Riemann-problem
  12 : Manufactured solutions: M1 closure; constant along along x+2y = c
  13 : Two initial Gaussians
  14 : Manufactured solutions: P1 closure; constant along along x+2y = c
  15 : Quasi 1D test-case with discontinuous coefficients
  16 : Manufactured solutions: K1 closure; constant along along x+2y = c
  17 : Reflective arc
  18 : Manufactured solutions: P1 closure on checkerboard
  19 : Manufactured solutions: P1 closure; discont. parameters in a box
  20 : Smooth checkerboard
  21 : Reed's problem
  22 : Beams enter vacuum ring
  23 : Initial atan

Closures:
    1: P1
    2: K1
    3: M1
    4: P2
    5: K2
    6: P3
    7: Simple nonlinear flux function
*/

#include"data_functions.h"
#include"mesh_utils.h"
#include"hyperbolic_solver.h"

/* double u_zero_val
 * defines initial condition
 * in: x coordinate -- x
 * in: y coordinate -- y
 * in: parameters   -- user_param
 * out: value of u0 at x,y
 */
void u_zero_val(double x,
                        double y,
                        parameters* user_param,
                        double* u_zero_val)
{
    double sigma_x;
    double sigma_y;
    double x_0;
    double y_0;

    int s_dim;

    //double x_0=0.25+((1.0)*(0.5));
    //double y_0=0.25+(sqrt(2.0)*(0.5));
    
    switch ((*user_param).test_case) 
    {
        case(1): 
            {
                if(SYSTEM_DIM ==1) 
                {
                    u_zero_val[0] =1.0; 
	            }
                else if(SYSTEM_DIM ==2) 
                {
                    u_zero_val[0] =1.0;
		            u_zero_val[1] =1.0;
	            }
                else if(SYSTEM_DIM ==3) 
                {
                    u_zero_val[0] =1.0;
		            u_zero_val[1] =1.0;
		            u_zero_val[2] =1.0;
                }
                break;
            }
        case(2):
            {
                sigma_x = 0.12;
                sigma_y = 0.12;
                sigma_x = 0.2;
                sigma_y = 0.2;
                x_0 = 0.5;
                y_0 = 0.5;
                /*x_0 = 0.45;*/
                /*y_0 = 0.35;*/

                u_zero_val[0] = 1.0*exp(-1.0*(((pow((x-x_0),2.0))/(2.0*pow(sigma_x,2.0))) + ((pow((y-y_0),2.0))/(2.0*pow(sigma_y,2.0)))));
                for (s_dim=1;s_dim<SYSTEM_DIM;s_dim++) 
                {
                    u_zero_val[s_dim] = 0.0; 
                }
                break;
            }
        case(3):
            {
                double val = 1.0;
                if(((x<=0.625) && (x>=0.375)) &&((y>=0.375) &&(y<=0.625))) 
                {
                    u_zero_val[0] = val;

                    if (SYSTEM_DIM > 2)
                    {
                        u_zero_val[1] = 0.0;
                        u_zero_val[2] = 0.0;
                        if (SYSTEM_DIM > 5) 
                        {
                            u_zero_val[3] = 1.0/3.0*val;
                            //u_zero_val[4] = 1.0/3.0*val;
                            u_zero_val[4] = 0.0;
                            u_zero_val[5] = 1.0/3.0*val;

                            if (SYSTEM_DIM > 9) 
                            {
                                for (s_dim=6; s_dim < 10; s_dim++) u_zero_val[s_dim] = 0.0;

                            }
                        }
                    }
                }
                else 
                {
                    for (s_dim=0;s_dim<SYSTEM_DIM;s_dim++) 
                    {
                        u_zero_val[s_dim] = 0.0; 
                    }
                }
                break;
            }
        case(4):
            {
                if(SYSTEM_DIM ==3)
                {
                    u_zero_val[0] = sin(2*PI*(x+y));
		            u_zero_val[1] = 0.0;
		            u_zero_val[2] = 0.0;
                }
                break;
            }
        case(6):
            {
                /*sigma_x = 0.06; */
                /*sigma_y = 0.06;*/
                //sigma_x = 0.28;
                //sigma_y = 0.28;
                sigma_x = 0.02;
                sigma_y = 0.02;
                x_0=0.00;
                y_0=0.00;
                /*y_0=0.1;*/

                double d=0.25;

                /*u_zero_val[0] = 1e-2;*/
                /*u_zero_val[1] = 1e-2;*/
                /*u_zero_val[2] = 1e-2;*/

                //u_zero_val[0] = u_zero_val[0] + 1000.0*exp(-10.0*(((pow((x-x_0),2.0))/(2.0*pow(sigma_x,2.0))) + ((pow((y-y_0),2.0))/(2.0*pow(sigma_y,2.0)))));
                u_zero_val[0] = exp(-10.0*(((pow((x-x_0),2.0))/(2.0*pow(sigma_x,2.0))) + ((pow((y-y_0),2.0))/(2.0*pow(sigma_y,2.0)))));
                
                if(u_zero_val[0] < 1.0e-04)
                {
                    u_zero_val[0] = 1e-04;
                }


               //u_zero_val[0] = 1.0*exp(-10.0*(((pow((y-y_0),2.0))/(2.0*pow(sigma_y,2.0)))));// + ((pow((y-y_0),2.0))/(2.0*pow(sigma_y,2.0)))));
                /*if((fabs(x) <= 0.25))
                    u_zero_val[0] = 1.0;
                else
                {

                    if(x>d)
                        u_zero_val[0] = 1.0*exp(-0.5*(pow((x-d),2.0)));// + ((pow((y-y_0),2.0))/(2.0*pow(sigma_y,2.0)))));
                    else if(x<-d)
                        u_zero_val[0] = 1.0*exp(-0.5*(pow((x+d),2.0)));// + ((pow((y-y_0),2.0))/(2.0*pow(sigma_y,2.0)))));
                }*/
                    

                if (SYSTEM_DIM > 2)
                {
                    u_zero_val[1] = 0.0;
                    u_zero_val[2] = 0.0;
                    if (SYSTEM_DIM > 5) 
                    {
                        u_zero_val[3] = u_zero_val[0]/3.0;
                        u_zero_val[4] = 0.0;
                        u_zero_val[5] = u_zero_val[0]/3.0;

                        /*u_zero_val[3] = u_zero_val[0]*4.0/3.0*PI;*/
                        /*u_zero_val[4] = 0.0;*/
                        /*u_zero_val[5] = u_zero_val[0]*4.0/3.0*PI;*/
                        if (SYSTEM_DIM > 9)
                        {
                            for (s_dim=6; s_dim < 10; s_dim++) 
                            {
                                u_zero_val[s_dim] = 0.0;
                            }
                        }
                    }
                }
                break;
            }
        case(7): 
            {
                if(SYSTEM_DIM ==3)
                {
                    u_zero_val[0] = sin(4*PI*((x-0.5)*(x-0.5)+(y-0.5)*(y-0.5)));
                    u_zero_val[1] = 0.0;
                    u_zero_val[2] = 0.0;
                }
                break;
            }
        case(8): 
            {
                if(SYSTEM_DIM ==3)
                {
                    u_zero_val[0] = sin(2*PI*((x-0.5)*(x-0.5)+(y-0.5)*(y-0.5)));
                    u_zero_val[1] = 1.0/2.0*(sin(x) - cos(x));
                    u_zero_val[2] = 1.0/2.0*(sin(y) + cos(y));
                }
                break;
            }
        case(11): 
            {
                if(SYSTEM_DIM ==3)
                {   
                    /*if ((x <= 0.5) && (x>= 0.45))*/
                    if (x <= 0.5) 
                    {
                        u_zero_val[0] = 1.0;
                        u_zero_val[1] = 0.99;
                        u_zero_val[2] = 0.0;
                    }
                    else 
                    {
                        u_zero_val[0] = 0.0;//0.5;
		                u_zero_val[1] = 0.0;
		                u_zero_val[2] = 0.0;
                    }
                }
                break;
            }
        case(12): 
            {
                if(SYSTEM_DIM ==3)
                {
                    u_zero_val[0] = sin(2.*PI*(x+2.*y)) + 3.;
                    u_zero_val[1] = 0.0;
                    u_zero_val[2] = 0.0;
                }
                break;
            }
        case(13): 
            {
                sigma_x = 0.12;
                sigma_y = 0.12;
                x_0 = 0.25;
                y_0 = 0.25;

                u_zero_val[0] = 1.0*exp(-1.0*(((pow((x-x_0),2.0))/(2.0*pow(sigma_x,2.0))) + ((pow((y-y_0),2.0))/(2.0*pow(sigma_y,2.0)))));

                sigma_x = 0.12;
                sigma_y = 0.12;
                x_0 = 0.6;
                y_0 = 0.6;
                
                u_zero_val[0] = u_zero_val[0] + 1.0*exp(-1.0*(((pow((x-x_0),2.0))/(2.0*pow(sigma_x,2.0))) + ((pow((y-y_0),2.0))/(2.0*pow(sigma_y,2.0)))));
                
                for (s_dim=1;s_dim<SYSTEM_DIM;s_dim++) 
                {
                    u_zero_val[s_dim] = 0.0; 
                }
                break;
            }
        case(14): 
            {
                if(SYSTEM_DIM ==3)
                {
                    u_zero_val[0] = sin(2*PI*(x+2.*y));
                    u_zero_val[1] = 0.0;
                    u_zero_val[2] = 0.0;
                }
                break;
            }
        case(15): 
            {
                if(SYSTEM_DIM ==3)
                {
                    if (x >= 0.25 && x <= 0.5) 
                    {
                        u_zero_val[0] = pow(sin(4*PI*(x-0.5)),8.0);
                    }
                    else 
                    {
                        u_zero_val[0] = 0;
                    }
                    u_zero_val[1] = 0.0;
                    u_zero_val[2] = 0.0;
                }
                break;
            }
        case(16): 
            {
                if(SYSTEM_DIM ==3)
                {
                    u_zero_val[0] = sin(2.*PI*(x+2.*y)) + 3.;
                    u_zero_val[1] = 0.0;
                    u_zero_val[2] = 0.0;
                }
                break;
            }
        case(17): 
            {
                if(SYSTEM_DIM ==3)
                {
                    if (x >= 10.9) 
                    {
                        if ( (y >= -5.55) && (y <= -5.35) ) 
                        {
                            u_zero_val[0] = 1.0;
                            u_zero_val[1] = -0.99;
                        }
                    }
                    else 
                    {
                        u_zero_val[0] = 0;
		                u_zero_val[1] = 0.0;
		            }
                    u_zero_val[2] = 0.0;
                }
                break;
            }
        case(18): 
            {
                if(SYSTEM_DIM ==3)
                {
                    u_zero_val[0] = sin(2.*PI*(x+2.*y)) + 3.;
                    u_zero_val[1] = 0.0;
                    u_zero_val[2] = 0.0;
	            }
                break;
            }
        case(19): 
            {
                if(SYSTEM_DIM ==3)
                {
                    u_zero_val[0] = sin(2.*PI*(x+2.*y)) + 3.;
		            u_zero_val[1] = 0.0;
		            u_zero_val[2] = 0.0;
	            }
                break;
            }
        case(23): 
            {
                if(SYSTEM_DIM ==3)
                {
                    u_zero_val[0] = atan(200.*(y-0.51)) + 1.;
		            u_zero_val[1] = 0.0;
		            u_zero_val[2] = 0.0;
	            }
                break;
            }
    case(24):{
                 u_zero_val[0] = 0.0;
                 u_zero_val[1] = 0.0;
                 u_zero_val[2] = 0.0;
                 break;
             }
    case(25):
        {
            if((pow(x,2.0) +pow(y,2.0)) <=0.25)
             {
                 u_zero_val[0] = 1.0;
                 u_zero_val[1] = 0.9;
                 u_zero_val[2] = 0.0;
             }
            else
            {
                u_zero_val[0] = 0.0;
                u_zero_val[1] = 0.0;
                u_zero_val[2] = 0.0;
            }
            break;
        }
    case(26):
        {
             u_zero_val[0] = 0.0;
             u_zero_val[1] = 0.0;
             u_zero_val[2] = 0.0;
             break;
        }
    case(50):
        {
            u_zero_val[0]= 1.0e-10/4/PI;
            u_zero_val[1]= 1.0e-10;
            u_zero_val[2] = 0.0;
            break;
        }
	default: 
        {
            for (s_dim=0;s_dim<SYSTEM_DIM;s_dim++)
            {
                u_zero_val[s_dim] = 0.0;
            }
            break;
        }
    }
}

/* double bdry_function_val
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
void  bdry_function_val(double x,
                        double y,
                        double x_ref,
                        double y_ref,
                        double t,
                        int E,
                        int iedge,
                        element* mesh_element,
                        node* mesh_node,
                        edge* mesh_edge,
                        double* uh_tn_sol_vec,
                        parameters* user_param,
                        double* bdry_f_val)
{
    int s_dim;
    double temp;
    int node_a;
    int node_b;
    double na_xc,na_yc;
    double nb_xc,nb_yc;

    node_a = mesh_edge[iedge].vertex[1];
    node_b = mesh_edge[iedge].vertex[2];

    na_xc = mesh_node[node_a].coord[0];
    na_yc = mesh_node[node_a].coord[1];

    nb_xc = mesh_node[node_b].coord[0];
    nb_yc = mesh_node[node_b].coord[1];



    switch ((*user_param).test_case) 
    {
        case(1): 
            {
                for (s_dim=0;s_dim<SYSTEM_DIM;s_dim++) 
                {
                    bdry_f_val[s_dim] = 1.0; 
                }
                break;
	}
        case(4): 
            {
                if(SYSTEM_DIM ==3) 
                {
                    bdry_f_val[0] = exp(-t)*sin(2*PI*(x+y));
                    bdry_f_val[1] = -1.0*exp(-t)*t*(2.0/3.0)*PI*cos(2*PI*(x+y));
                    bdry_f_val[2] = -1.0*exp(-t)*t*(2.0/3.0)*PI*cos(2*PI*(x+y));
                }
                break;
            }
        case(6):
            {
                bdry_f_val[0] = 1.0e-04;
                bdry_f_val[1] = 0.0;
                bdry_f_val[2] = 0.0;
                break;
            }
        case(7): 
            {
                if(SYSTEM_DIM ==3)
                {
                    bdry_f_val[0] = exp(-t)*sin(4*PI*((x-0.5)*(x-0.5)+(y-0.5)*(y-0.5)));
                    bdry_f_val[1] = -4./3.*PI*t*exp(-t)*(2.0*x-1.0)*cos(4*PI*((x-0.5)*(x-0.5)+(y-0.5)*(y-0.5)));
                    bdry_f_val[2] = -4./3.*PI*t*exp(-t)*(2.0*y-1.0)*cos(4*PI*((x-0.5)*(x-0.5)+(y-0.5)*(y-0.5)));
                }
                break;
            }
        case(8): 
            {
                if(SYSTEM_DIM ==3)
                {
                    temp = -4./3.*PI*t*exp(-t)*cos(4*PI*((x-0.5)*(x-0.5)+(y-0.5)*(y-0.5)));
                    bdry_f_val[0] = exp(-t)*sin(4*PI*((x-0.5)*(x-0.5)+(y-0.5)*(y-0.5)));
                    bdry_f_val[1] = (2.0*x-1.0)*temp + 1.0/2.0*(sin(x-t) - cos(x-t));
                    bdry_f_val[2] = (2.0*y-1.0)*temp + 1.0/2.0*(sin(y-t) + cos(y-t));
                }
                break;
            }
        case(11): 
            {
                if(SYSTEM_DIM ==3)
                {
                    if (x <= 0.00) 
                    {
                        bdry_f_val[0] = 1.0;
                        bdry_f_val[1] = 0.99;
                        bdry_f_val[2] = 0.0;
                    }
                    if (x >= 1.0) 
                    {
                        bdry_f_val[0] = 0.5;
                        bdry_f_val[1] = 0.0;
                        bdry_f_val[2] = 0.0;
                    }
                }
                break;
            }
        case(12): 
            {
                if(SYSTEM_DIM ==3) 
                {
                    bdry_f_val[0] = exp(-t)*sin(2.*PI*(x+2.*y)) + 3.;
                    bdry_f_val[1] = exp(-t)*t*cos(2.*PI*(x+2.*y));
                    bdry_f_val[2] = exp(-t)*t*cos(2.*PI*(x+2.*y));
                }
                break;
            }
        case(14): 
            {
                if(SYSTEM_DIM ==3) 
                {
                    bdry_f_val[0] = exp(-t)*sin(2*PI*(x+2.*y));
                    bdry_f_val[1] = -1.0*exp(-t)*t*(2.0/3.0)*PI*cos(2*PI*(x+2.*y));
                    bdry_f_val[2] = -1.0*exp(-t)*t*(4.0/3.0)*PI*cos(2*PI*(x+2.*y));
	            }
                break;
	        }
        case(16): 
            {
                if(SYSTEM_DIM ==3) 
                {
                    bdry_f_val[0] = exp(-t)*sin(2.*PI*(x+2.*y)) + 3.;
                    bdry_f_val[1] = exp(-t)*t*cos(2.*PI*(x+2.*y));
                    bdry_f_val[2] = exp(-t)*t*cos(2.*PI*(x+2.*y));
                }
                break;
            }
        case(18): 
            {
                if(SYSTEM_DIM ==3) 
                {
		            bdry_f_val[0] = exp(-t)*sin(2.*PI*(x+2.*y)) + 3.;
		            bdry_f_val[1] = exp(-t)*t*cos(2.*PI*(x+2.*y));
		            bdry_f_val[2] = exp(-t)*t*cos(2.*PI*(x+2.*y));
                }
                break;
            }
        case(19): 
            {
                if(SYSTEM_DIM ==3) 
                {
                    bdry_f_val[0] = exp(-t)*sin(2.*PI*(x+2.*y)) + 3.;
		            bdry_f_val[1] = exp(-t)*t*cos(2.*PI*(x+2.*y));
		            bdry_f_val[2] = exp(-t)*t*cos(2.*PI*(x+2.*y));
                }
                break;
            }
        case(23): 
            {
                if(SYSTEM_DIM ==3)
                {
                    bdry_f_val[0] = atan(200.*(y-0.51)) + 1;;
                    bdry_f_val[1] = 0.0;
                    bdry_f_val[2] = 0.0;
                }
                break;
            }
        case(24):
            {
                 bdry_f_val[0] = 0.0;
                 bdry_f_val[1] = 0.0;
                 bdry_f_val[2] = 0.0;
                 break;
            }
        case(25):
            {
                 bdry_f_val[0] = 0.0;
                 bdry_f_val[1] = 0.0;
                 bdry_f_val[2] = 0.0;
                 break;
            }
        case(26):
            {
                if((na_xc < EPSILON) && (nb_xc < EPSILON))
                {
                    bdry_f_val[0] = 1.0; 
                    bdry_f_val[1] = 0.99;
                    bdry_f_val[2] = 0.0;
                }
                else
                {
                     bdry_f_val[0] = 0.0;
                     bdry_f_val[1] = 0.0;
                     bdry_f_val[2] = 0.0;
                }
                break;
             }
        case(27):
            {
                //Spot 1
                if((y >=3.0) && (y <=4.0) && (fabs(x) < EPSILON))
                {
                    bdry_f_val[0] = 62.733495998904842849697161000222;
                    bdry_f_val[1] = 56.451755972807482919506583129987;
                    bdry_f_val[2] = 0.000000010684557813062850177457212275112;
                 }//spot 2
                 else if((x >=3.0) && (x <=4.0) && (fabs(y) < EPSILON))
                 {
                     bdry_f_val[0] = 62.733495998904842849697161000222;
                     bdry_f_val[1] = 0.000000010684557813062850177457212275112;
                     bdry_f_val[2] = 56.451755972807482919506583129987;
                 }
                 else //spot 3
                 {
                     bdry_f_val[0] = 1e-10;
                     bdry_f_val[1] = 0.0;
                     bdry_f_val[2] = 0.0;
                 }
                 break;
             }
        case(50):
            {
                //if((fabs(x-0.5 < EPSILON)) && ((y>=-0.2143) && (y<=0.2143)))
              
                if((fabs(x-0.5 < EPSILON)) && ((y>=(-0.05 +0.0)) && (y<=(0.05+0.0))))
                { 
                    bdry_f_val[0] = 100.00;
                    bdry_f_val[1] = -99.999;
                    bdry_f_val[2] = 0.0;
                }
                else
                {
                    bdry_f_val[0] = 1.0e-10/4/PI;
                    bdry_f_val[1] = 1.0e-10;
                    bdry_f_val[2] = 0.0;
                }
                break;
            }
        default: 
            {
                for (s_dim=0;s_dim<SYSTEM_DIM;s_dim++) 
                {
                    bdry_f_val[s_dim] = 0.0;
                }
                break;
            }
    }
}

/* void flux_function
 * defines the flux in x and y direction
 * in: solution vector -- uh_tn_sol
 * in: density	       -- density at current position
 * in: parameters      -- user_param
 * out: flux funcion   -- flux_val
 */
void flux_function(double* uh_tn_sol,
                   double x,
                   double y,
                   double t,
                   double density,
                   parameters* user_param,
                   int s_dim_count,
                   double* flux_val)
{
    // m1 Eddington factor
    double D[2];
    double n[2];
    double norm, f, fac1, fac2, Chi;

    if ((*user_param).test_case == 4 || (*user_param).test_case == 7 || (*user_param).test_case == 14 || (*user_param).test_case == 15 || (*user_param).test_case == 18 || (*user_param).test_case == 19) 
    {
        // Enforce that P1 Closure is used
        (*user_param).closure = 1;
    }
    else if ((*user_param).test_case == 8)
    {
        // modified nonlinear P1 Closure 
	    if(s_dim_count == 0)
        {
	        flux_val[0] = uh_tn_sol[1] + sin(x-t); 
	        flux_val[1] = uh_tn_sol[2] + cos(y-t);
	    }
        else if(s_dim_count == 1)
        {
            flux_val[0] = uh_tn_sol[0]/3.0 + sin(x-t);
            flux_val[1] = 0.0;
        }
        else if(s_dim_count == 2)
        {
            flux_val[0] = 0.0;
            flux_val[1] = uh_tn_sol[0]/3.0 + cos(y-t);
        }
        return;
    }
    else if ((*user_param).test_case == 12)
    {
        // Enforce that M1 Closure is used
	    (*user_param).closure = 3;
    }
    else if ((*user_param).test_case == 16)
    {
	    // Enforce that K1 Closure is used
	    (*user_param).closure = 2;
    }
    if(SYSTEM_DIM ==3)
    {
        switch ((*user_param).closure) 
        {
            // P1 Closure
            case (1): 
                {
                    /*ensure_realizability_ortho_proj(uh_tn_sol);*/
		            /*ensure_realizability_vert_proj(1e-8,uh_tn_sol);*/

                    if(s_dim_count == 0)
                    {
                        flux_val[0] = uh_tn_sol[1];
                        flux_val[1] = uh_tn_sol[2];
                    }
                    else if(s_dim_count == 1)
                    {
                        flux_val[0] = uh_tn_sol[0]/3.;
                        flux_val[1] = 0.0;
                    }
                    else if(s_dim_count == 2)
                    {
                        flux_val[0] = 0.0;
                        flux_val[1] = uh_tn_sol[0]/3.;
                    }
                    break;
                }// P1 closure
            case (2): 
                {
                    ensure_realizability_vert_proj(1e-8,uh_tn_sol);
                    double u0tol = 1e-7;
                    if(s_dim_count == 0)
                    {
                        flux_val[0] = uh_tn_sol[1];
                        flux_val[1] = uh_tn_sol[2];
                    }
                    else if(s_dim_count == 1)
                    {
                        flux_val[0] = uh_tn_sol[0]/3.0;
                        flux_val[1] = 0.0;

                        if (uh_tn_sol[0] > u0tol) 
                        {
                            flux_val[0] += 2.0/3.0*uh_tn_sol[1]*uh_tn_sol[1] / uh_tn_sol[0];
                            flux_val[1] += 2.0/3.0*uh_tn_sol[1]*uh_tn_sol[2] / uh_tn_sol[0];
                        }
                        if (flux_val[0] != flux_val[0]) 
                        {
                            printf("u0 = %lf\t, u1 = %lf\t, u2 = %lf\n",uh_tn_sol[0],uh_tn_sol[1],uh_tn_sol[2]);
                            printf("f0 = %lf\t, f1 = %lf\t",flux_val[0],flux_val[1]);
                            getchar();
                        }
                    }
                    else if(s_dim_count == 2)
                    {
                        flux_val[0] = 0.0;
                        flux_val[1] = uh_tn_sol[0]/3.0;

                        if (uh_tn_sol[0] > u0tol) 
                        {
                            flux_val[0] += 2.0/3.0*uh_tn_sol[1]*uh_tn_sol[2] / uh_tn_sol[0];
                            flux_val[1] += 2.0/3.0*uh_tn_sol[2]*uh_tn_sol[2] / uh_tn_sol[0];
                        }
                        if (flux_val[0] != flux_val[0]) 
                        {
                            printf("u0 = %lf\t, u1 = %lf\t, u2 = %lf\n",uh_tn_sol[0],uh_tn_sol[1],uh_tn_sol[2]);
                            printf("f0 = %lf\t, f1 = %lf\t",flux_val[0],flux_val[1]);
                            getchar();
                        }
                    }
                }// K1 closure
            case (3): 
             {
                 double u0tol = 1e-12;
                 double normtol = 1e-12;

                //double u0tol = 1e-06;
                //double normtol = 1e-06;


                if(!(*user_param).flag_ensure_positivity)
                {
                    realizability_hack(1e-05,uh_tn_sol);
                }
                        
                int rel_flag=0;
                rel_flag = is_realizable(uh_tn_sol);
                if(rel_flag ==0)
                {
                    //printf("rel_flag in M1 closure  %d ============>solution %10.16e %10.16e %10.16e  \n",rel_flag,uh_tn_sol[0],uh_tn_sol[1],uh_tn_sol[2]);
                    //exit(1);
                }
                if(s_dim_count == 0) 
                {
                    flux_val[0] = uh_tn_sol[1];
                    flux_val[1] = uh_tn_sol[2];
                }
                else if(s_dim_count == 1) 
                {
                    norm = sqrt(uh_tn_sol[1]*uh_tn_sol[1] + uh_tn_sol[2]*uh_tn_sol[2]);
                    if ( (uh_tn_sol[0] <= u0tol ) || ( norm <= normtol ) ) 
                    {
                        flux_val[0] = uh_tn_sol[0]/3.;
                        flux_val[1] = 0.;
                    }
                    else 
                    {
                        f = norm / uh_tn_sol[0];
                        n[0] = uh_tn_sol[1] / norm;
                        n[1] = uh_tn_sol[2] / norm;

                        if (f > 1.1547) 
                        {
                            //printf(">>>>>>>>>>>>>>>>>>>alert f %lf norm %10.16e  uhtnsol %10.16e >>>>>>>>>>>>>>>\n",f,norm,uh_tn_sol[0]);
                            f = 1.1547;
                        }

                        /*Chi = (5. - 2*sqrt(4. - 3.*f*f))/3.;*/
                        Chi = (3. + 4.*f*f) / (5. + 2.*sqrt(4. - 3.*f*f));
                        fac1 = 1. - Chi;
                        fac2 = 3. * Chi - 1.;

                        D[0] = fac2 * n[0]*n[0] + fac1;  
                        D[1] = fac2 * n[0]*n[1];
                   

                        flux_val[0] = D[0]/2. * uh_tn_sol[0];
                        flux_val[1] = D[1]/2. * uh_tn_sol[0];
                    }
                }
                else if(s_dim_count == 2) 
                {
                    norm = sqrt(uh_tn_sol[1]*uh_tn_sol[1] + uh_tn_sol[2]*uh_tn_sol[2]);
                    if ( (uh_tn_sol[0] <= u0tol) || ( norm <= normtol ) ) 
                    {
                        flux_val[0] = 0.;
			            flux_val[1] = uh_tn_sol[0]/3.;
                    }
                    else 
                    {
                        f = norm / uh_tn_sol[0];
                        n[0] = uh_tn_sol[1] / norm;
                        n[1] = uh_tn_sol[2] / norm;

                        if (f > 1.1547) 
                        {
                            //printf(">>>>>>>>>>>>>>>>>>>alert f %lf norm %10.16e  uhtnsol %10.16e >>>>>>>>>>>>>>>\n",f,norm,uh_tn_sol[0]);
                            f = 1.1547;
                        }

                        /*Chi = (5. - 2*sqrt(4. - 3.*f*f))/3.;*/
                        Chi = (3. + 4.*f*f) / (5. + 2.*sqrt(4. - 3.*f*f));
                        fac1 = 1. - Chi;
                        fac2 = 3. * Chi - 1.;

                        D[0] = fac2 * n[0]*n[1];  
                        D[1] = fac2 * n[1]*n[1] + fac1;


                        flux_val[0] = D[0]/2. * uh_tn_sol[0];
                        flux_val[1] = D[1]/2. * uh_tn_sol[0];
                    }
                }
                break;
             }// M1 closure
            case (7): 
             {
                 /*ensure_realizability_ortho_proj(uh_tn_sol);*/

                if(s_dim_count == 0)
                {
                    flux_val[0] = uh_tn_sol[1];
                    flux_val[1] = uh_tn_sol[2];
                }
                else if(s_dim_count == 1)
                {
                    flux_val[0] = pow(uh_tn_sol[0],3.0)/3.;
                    flux_val[1] = 0.0;
                }
                else if(s_dim_count == 2)
                {
                    flux_val[0] = 0.0;
                    flux_val[1] = pow(uh_tn_sol[0],3.0)/3.;
                }
                break;
             }// Simple nonelinear closure
        }//case closure
    }
    else if(SYSTEM_DIM == 6)
    {
	switch ((*user_param).closure) {

	    case (4): {
		// P2 Closure
		switch (s_dim_count) {

		    case (0): {
			flux_val[0] = uh_tn_sol[1];
			flux_val[1] = uh_tn_sol[2];
			break;
		    }

		    case (1): {
			flux_val[0] = uh_tn_sol[3];
			flux_val[1] = uh_tn_sol[4];
			break;
		    }

		    case (2): {
			flux_val[0] = uh_tn_sol[4];
			flux_val[1] = uh_tn_sol[5];
			break;
		    }

		    case (3): {
			flux_val[0] = 3.0/5.0*uh_tn_sol[1];
			flux_val[1] = 1.0/5.0*uh_tn_sol[2];
			break;
		    }

		    case (4): {
			flux_val[0] = 1.0/5.0*uh_tn_sol[2];
			flux_val[1] = 1.0/5.0*uh_tn_sol[1];
			break;
		    }

		    case (5): {
			flux_val[0] = 1.0/5.0*uh_tn_sol[1];
			flux_val[1] = 3.0/5.0*uh_tn_sol[2];
			break;
		    }
		    /*case (3): {*/
			/*flux_val[0] = uh_tn_sol[0]*3.0/5.0*uh_tn_sol[1];*/
			/*flux_val[1] = uh_tn_sol[0]*1.0/5.0*uh_tn_sol[2];*/
			/*break;*/
		    /*}*/

		    /*case (4): {*/
			/*flux_val[0] = uh_tn_sol[0]*1.0/5.0*uh_tn_sol[2];*/
			/*flux_val[1] = uh_tn_sol[0]*1.0/5.0*uh_tn_sol[1];*/
			/*break;*/
		    /*}*/

		    /*case (5): {*/
			/*flux_val[0] = uh_tn_sol[0]*1.0/5.0*uh_tn_sol[1];*/
			/*flux_val[1] = uh_tn_sol[0]*3.0/5.0*uh_tn_sol[2];*/
			/*break;*/
		    /*}*/
		}
		return;
	    }// P2 closure

	    case (5): {
		// K2 Closure

		double u0tol = 1e-12;

		switch (s_dim_count) {

		    case (0): {
			flux_val[0] = uh_tn_sol[1];
			flux_val[1] = uh_tn_sol[2];
			break;
		    }

		    case (1): {
			flux_val[0] = uh_tn_sol[3];
			flux_val[1] = uh_tn_sol[4];
			break;
		    }

		    case (2): {
			flux_val[0] = uh_tn_sol[4];
			flux_val[1] = uh_tn_sol[5];
			break;
		    }

		    case (3): {

			if ( uh_tn_sol[0] <= u0tol ) {
			    flux_val[0] = 0.;
			    flux_val[1] = 0.;
			}
			else {
			    flux_val[0] = uh_tn_sol[1]*uh_tn_sol[3]/uh_tn_sol[0];
			    flux_val[1] = (2.0/3.0*uh_tn_sol[1]*uh_tn_sol[4] + 1.0/3.0*uh_tn_sol[2]*uh_tn_sol[3])/uh_tn_sol[0];
			}
			break;
		    }

		    case (4): {
			if ( uh_tn_sol[0] <= u0tol ) {
			    flux_val[0] = 0.;
			    flux_val[1] = 0.;
			}
			else {
			    flux_val[0] = (2.0/3.0*uh_tn_sol[1]*uh_tn_sol[4] + 1.0/3.0*uh_tn_sol[2]*uh_tn_sol[3])/uh_tn_sol[0];
			    flux_val[1] = (1.0/3.0*uh_tn_sol[1]*uh_tn_sol[4] + 2.0/3.0*uh_tn_sol[2]*uh_tn_sol[3])/uh_tn_sol[0];
			}
			break;
		    }

		    case (5): {
			if ( uh_tn_sol[0] <= u0tol ) {
			    flux_val[0] = 0.;
			    flux_val[1] = 0.;
			}
			else {
			    flux_val[0] = (1.0/3.0*uh_tn_sol[1]*uh_tn_sol[4] + 2.0/3.0*uh_tn_sol[2]*uh_tn_sol[3])/uh_tn_sol[0];
			    flux_val[1] = uh_tn_sol[2]*uh_tn_sol[5]/uh_tn_sol[0];
			}
			break;
		    }
		}
	    }// K2 closure
	    return;
	}
    }
    else if(SYSTEM_DIM == 10)
    {
	switch ((*user_param).closure) {

	    case (6): {
		// P3 Closure
		switch (s_dim_count) {

		    case (0): {
			flux_val[0] = uh_tn_sol[1];
			flux_val[1] = uh_tn_sol[2];
			break;
		    }

		    case (1): {
			flux_val[0] = uh_tn_sol[3];
			flux_val[1] = uh_tn_sol[4];
			break;
		    }

		    case (2): {
			flux_val[0] = uh_tn_sol[4];
			flux_val[1] = uh_tn_sol[5];
			break;
		    }

		    case (3): {
			flux_val[0] = uh_tn_sol[6];
			flux_val[1] = uh_tn_sol[7];
			break;
		    }

		    case (4): {
			flux_val[0] = uh_tn_sol[7];
			flux_val[1] = uh_tn_sol[8];
			break;
		    }

		    case (5): {
			flux_val[0] = uh_tn_sol[8];
			flux_val[1] = uh_tn_sol[9];
			break;
		    }

		    case (6): {
			flux_val[0] = -3.0/35.0*uh_tn_sol[0] + 6.0/7.0*uh_tn_sol[3];
			flux_val[1] = 3.0/7.0*uh_tn_sol[4];
			break;
		    }

		    case (7): {
			flux_val[0] = 3.0/7.0*uh_tn_sol[4];
			flux_val[1] = -1.0/35.0*uh_tn_sol[0] + 1.0/7.0*uh_tn_sol[3] + 1.0/7.0*uh_tn_sol[5];
			break;
		    }

		    case (8): {
			flux_val[0] = -1.0/35.0*uh_tn_sol[0] + 1.0/7.0*uh_tn_sol[3] + 1.0/7.0*uh_tn_sol[5];
			flux_val[1] = 3.0/7.0*uh_tn_sol[4];
			break;
		    }

		    case (9): {
			flux_val[0] = 3.0/7.0*uh_tn_sol[4];
			flux_val[1] = -3.0/35.0*uh_tn_sol[0] + 6.0/7.0*uh_tn_sol[5];
			break;
		    }
		}
	    }// P3 closure
	    break;
	}
    }

    // Flux after transformation 
    // do not transform when using Henyey-Greenstein model!
    /*flux_val[0] = flux_val[0] / density;*/
    /*flux_val[1] = flux_val[1] / density;*/
}

/* This function seems to be obsolete, as long as upwind_flux 
 * in hyperbolic_solver is not used!
 */
void advection_v(double t_val,
                  double* phys_coords,
                  parameters* user_param,
		          int s_dim_counter,
                  double* a_vector)
{
    if ((*user_param).test_case ==1)
    {
        a_vector[0] = 1.0;
        a_vector[1] = 1.0;
    }
    else if ((*user_param).test_case == 2)
    {

        a_vector[0] = cos((PI*t_val)/T)*(pow(sin(PI*phys_coords[0]),2.0)*sin(2*PI*phys_coords[1]));
        a_vector[1] = cos((PI*t_val)/T)*(-1.0*sin(2*PI*phys_coords[0])*pow(sin(PI*phys_coords[1]),2.0));
    }
    else if ((*user_param).test_case ==3)
    {
	a_vector[0] = 0.0;
	a_vector[1] = 0.0;
    }
    else if ( (*user_param).test_case == 4 || (*user_param).test_case == 7 || (*user_param).test_case == 8)
    {
        if(s_dim_counter==0)
        {
           a_vector[0] = 1.0;
           a_vector[1] = 1.0;
        }
        else if(s_dim_counter==1)
        {
            a_vector[0] = 1.0/3.0;
            a_vector[1] = 0.0;
        }
        else if(s_dim_counter==2)
        {
            a_vector[0] = 0.0;
            a_vector[1] = 1.0/3.0;
        }
    }

}

/* void transf_to_char_var
 * transforms vector into the characteristic field,
 * assuming n has length 1 and solution is evaluated at barycenter
 * in: solution vector at the barycenter -- uh_sol
 * in: vector of length SYSTEM_DIM -- vec
 * in: unit vector from barycenter to edge midpoint -- n
 * in: parameters      -- user_param
 * out: transformed vector -- char_vec
 */
void transf_to_char_var(double* uh_sol,
		            	double* vec,
                        double* n,
                        parameters* user_param,
                        double* char_vec)
{
    // Transformation matrix
    double inv_R[3][3];

    int verbose=0;

    int i;

    switch ((*user_param).closure) {

	case (1): {
	    if (SYSTEM_DIM == 3){
	        inv_R[0][0] = -sqrt(3.0)/6.0; inv_R[0][1] = n[0]/2.0;  inv_R[0][2] = n[1]/2.0;
	        inv_R[1][0] = sqrt(3.0)/6.0;  inv_R[1][1] = n[0]/2.0;  inv_R[1][2] = n[1]/2.0;
	        inv_R[2][0] = 0.0;            inv_R[2][1] = -n[1];     inv_R[2][2] = n[0];

		for (i=0;i<3;i++) char_vec[i] = vec[0]*inv_R[i][0] + vec[1]*inv_R[i][1] + vec[2]*inv_R[i][2];
	    }
	    break;
	}

	// K1
	case (2): {
	    if (SYSTEM_DIM == 3){
	        double N1x, N1y;

		/*ensure_realizability_ortho_proj(uh_sol);*/
		//ensure_realizability_vert_proj(1e-6,uh_sol);

		if (uh_sol[0] > 1e-12){
		    N1x = uh_sol[1]/uh_sol[0];
		    N1y = uh_sol[2]/uh_sol[0];
		}
		else {
		    N1x = 0.0;
		    N1y = 0.0;
		}

		double nnx = n[0]*N1x + n[1]*N1y;
		double nny = n[0]*N1y - n[1]*N1x;

		if (nnx > 1.00) {
		    printf("Cutting off! This shouldn't happen! ensure_realizability_ortho_proj working?\n");
		    nnx = 1.00;
		}

		double root = sqrt(3.0 - 2.0*nnx*nnx);

		double l1 = (2.0 * nnx + root)/3.0;
		double l2 = (2.0 * nnx - root)/3.0;

		double a1 = 2.0/3.0*nny*(1 - nnx/root);
		double a2 = 2.0/3.0*nny*(1 + nnx/root);

		double nl = 2.0/3.0*root;

	        inv_R[0][0] = -l2/nl;             inv_R[0][1] =  n[0]/nl; 		   inv_R[0][2] =  n[1]/nl;
	        inv_R[1][0] =  l1/nl;             inv_R[1][1] = -n[0]/nl;   		   inv_R[1][2] = -n[1]/nl;
	        inv_R[2][0] = (a1*l2-a2*l1)/nl;   inv_R[2][1] =  n[0]*(a2-a1)/nl - n[1];   inv_R[2][2] =  n[1]*(a2-a1)/nl + n[0];

		for (i=0;i<3;i++) {
		    char_vec[i] = vec[0]*inv_R[i][0] + vec[1]*inv_R[i][1] + vec[2]*inv_R[i][2];
		}
	    }
	    break;
	}

	// M1
	case (3): {
	    if (SYSTEM_DIM == 3){
		/*ensure_realizability_ortho_proj(uh_sol);*/
	   //ensure_realizability_vert_proj(1e-4,uh_sol);
        int rel_flag=0;
        //rel_flag = is_realizable(uh_sol);
        if(rel_flag ==0)
        {
            //printf("rel_flag %d ============>solution %10.16e %10.16e %10.16e  \n",rel_flag,uh_sol[0],uh_sol[1],uh_sol[2]);
        }

		double N1x, N1y;

		N1x = uh_sol[1]/uh_sol[0];
		N1y = uh_sol[2]/uh_sol[0];

        //printf("N1_x %lf N1_y %lf  %lf \n",N1x,N1y,uh_sol[0]);

		double nnx = n[0]*N1x + n[1]*N1y;
		double nny = n[0]*N1y - n[1]*N1x;

		int flag_use_table = 0;

		if (flag_use_table == 1) {
		    double f, phi;
		    f = sqrt(uh_sol[1]*uh_sol[1] + uh_sol[2]*uh_sol[2]) / uh_sol[0];
		    phi = atan2(nny,nnx) + PI;

		    int i,j;
		    i = round(f/(*user_param).fmax*((*user_param).nf-1));
		    j = round(phi/(2*PI)*((*user_param).nphi-1));

	            int r,c;
		    for (r=0;r<3;r++) {
			for (c=0;c<3;c++) {
			    inv_R[r][c] = (*user_param).inv_R[i][j][r][c];
			}
		    }

// invR = [S11(i,j) S12(i,j) S13(i,j); S21(i,j) S22(i,j) S23(i,j); S31(i,j) S32(i,j) S33(i,j)] * [1 0 0; 0 nx ny; 0 -ny nx];

		    double trafo[3][3];

		    for (c=0;c<3;c++) {
			trafo[c][0] = inv_R[c][0]; 
			trafo[c][1] = n[0]*inv_R[c][1] - n[1]*inv_R[c][2];
			trafo[c][2] = n[0]*inv_R[c][2] + n[1]*inv_R[c][1];
		    }

		    for (i=0;i<3;i++) char_vec[i] = vec[0]*trafo[i][0] + vec[1]*trafo[i][1] + vec[2]*trafo[i][2];
		    /*for (i=0;i<3;i++) char_vec[i] = vec[0]*inv_R[i][0] + vec[1]*inv_R[i][1] + vec[2]*inv_R[i][2];*/

		    //if (char_vec[0] != char_vec[0]) 
            if(verbose){
			printf("TO!\n");
			printf("f = %lf\t phi = %lf\n",f, phi);
			printf("inv_R = [%lf\t%lf\t%lf]\n",inv_R[0][0],inv_R[0][1],inv_R[0][2]);
			printf("inv_R = [%lf\t%lf\t%lf]\n",inv_R[1][0],inv_R[1][1],inv_R[1][2]);
			printf("inv_R = [%lf\t%lf\t%lf]\n",inv_R[2][0],inv_R[2][1],inv_R[2][2]);
			printf("--------------------------------------------------------------\n");

			printf("inv_R_trafo = [%lf\t%lf\t%lf]\n",trafo[0][0],trafo[0][1],trafo[0][2]);
			printf("inv_R_trafo = [%lf\t%lf\t%lf]\n",trafo[1][0],trafo[1][1],trafo[1][2]);
			printf("inv_R_trafo = [%lf\t%lf\t%lf]\n",trafo[2][0],trafo[2][1],trafo[2][2]);
			printf("--------------------------------------------------------------\n");

			printf("vec = [%lf\t%lf\t%lf]\n",vec[0],vec[1],vec[2]);
			printf("char_vec = [%lf\t%lf\t%lf]\n\n",char_vec[0],char_vec[1],char_vec[2]);
			getchar();
		    }
		}
		else {
		    if (nnx > 1.00) {
                printf("nnx %lf N1_x %lf N1_y %lf  %lf \n",nnx,N1x,N1y,uh_sol[0]);
			printf("Cutting off! This shouldn't happen! ensure_realizability_ortho_proj working?\n");
			nnx = 1.00;
		    }

		    double lambda1 = nnx*nnx + nny*nny;
		    double lambda2 = nny*nny - nnx*nnx;

		    inv_R[0][0] = 0.0;             
		    inv_R[0][1] = -(((pow(lambda1,2) + lambda2)*(nny*n[0] + nnx*n[1]))/(nnx*(pow(lambda1,2) + lambda2 - 2*pow(nny,2)))); 		   
		    inv_R[0][2] = ((pow(lambda1,2) + lambda2)*(nnx*n[0] - nny*n[1]))/(nnx*(pow(lambda1,2) + lambda2 - 2*pow(nny,2)));

		    inv_R[1][0] = -((pow(lambda1,6) - 4*pow(lambda2,2))*nny)/(4.*(pow(lambda1,5) - 2*sqrt(lambda1)*pow(lambda2,2)));             
		    inv_R[1][1] = -((pow(lambda1,3) - 2*lambda2)*nny*((-(pow(lambda1,1.5)*lambda2) - pow(lambda1,2)*lambda2 - pow(lambda2,2) + 2*sqrt(lambda1)*lambda2*pow(nny,2) + pow(lambda1,3.5)*(-1 + pow(nny,2)))*n[0] + (-2*pow(lambda1,1.5) + pow(lambda1,3.5) - 2*lambda2 + 2*sqrt(lambda1)*lambda2)*nnx*nny*n[1]))/ (2.*(pow(lambda1,4.5) - 2*pow(lambda2,2))*nnx*(pow(lambda1,2) + lambda2 - 2*pow(nny,2)));   		   
		    inv_R[1][2] = ((pow(lambda1,3) - 2*lambda2)*nny*((-2*pow(lambda1,1.5) + pow(lambda1,3.5) - 2*lambda2 + 2*sqrt(lambda1)*lambda2)*nnx*nny*n[0] + (pow(lambda1,1.5)*lambda2 + pow(lambda1,2)*lambda2 + pow(lambda2,2) - 2*sqrt(lambda1)*lambda2*pow(nny,2) - pow(lambda1,3.5)*(-1 + pow(nny,2)))*n[1]))/ (2.*(pow(lambda1,4.5) - 2*pow(lambda2,2))*nnx*(pow(lambda1,2) + lambda2 - 2*pow(nny,2))); 

		    inv_R[2][0] = ((pow(lambda1,6) - 4*pow(lambda2,2))*nny)/(4.*(pow(lambda1,5) - 2*sqrt(lambda1)*pow(lambda2,2)));   
		    inv_R[2][1] = ((pow(lambda1,3) + 2*lambda2)*nny*((pow(lambda1,1.5)*lambda2 - pow(lambda1,2)*lambda2 - pow(lambda2,2) - 2*sqrt(lambda1)*lambda2*pow(nny,2) + pow(lambda1,3.5)*(1 + pow(nny,2)))*n[0] + (2*pow(lambda1,1.5) + pow(lambda1,3.5) - 2*lambda2 - 2*sqrt(lambda1)*lambda2)*nnx*nny*n[1]))/ (2.*(pow(lambda1,4.5) - 2*pow(lambda2,2))*nnx*(pow(lambda1,2) + lambda2 - 2*pow(nny,2)));   
		    inv_R[2][2] = ((pow(lambda1,3) + 2*lambda2)*nny*((-2*pow(lambda1,1.5) - pow(lambda1,3.5) + 2*lambda2 + 2*sqrt(lambda1)*lambda2)*nnx*nny*n[0] + (pow(lambda1,1.5)*lambda2 - pow(lambda1,2)*lambda2 - pow(lambda2,2) - 2*sqrt(lambda1)*lambda2*pow(nny,2) + pow(lambda1,3.5)*(1 + pow(nny,2)))*n[1]))/ (2.*(pow(lambda1,4.5) - 2*pow(lambda2,2))*nnx*(pow(lambda1,2) + lambda2 - 2*pow(nny,2))); 

		    for (i=0;i<3;i++) {
			char_vec[i] = vec[0]*inv_R[i][0] + vec[1]*inv_R[i][1] + vec[2]*inv_R[i][2];
		    }

		    /*if (vec[0] > 1e-12) {*/
		    if (char_vec[0] != char_vec[0]) {
			printf("TO!\n");
			printf("n = [%lf \t %lf]\tnorm = %lf\n",n[0],n[1],sqrt(n[0]*n[0]+n[1]*n[1]));
			printf("nnx = %lf \t nny = %lf\n",nnx,nny);
			printf("lambda1 = %lf \t lambda2 = %lf\n",lambda1,lambda2);
			printf("vec = [%lf\t%lf\t%lf]\n",vec[0],vec[1],vec[2]);
			printf("char_vec = [%lf\t%lf\t%lf]\n\n",char_vec[0],char_vec[1],char_vec[2]);
		    }
		}
	    }
	    break;
	}

	// No transformation by default
	default: {
	    for (i=0;i<SYSTEM_DIM;i++) char_vec[i] = vec[i];
	    break;
	}
    }
}

/* void transf_from_char_var
 * transforms vector from the characteristic field
 * to the original coordinates, assuming n has length 1
 * and solution is evaluated at barycenter
 * in: solution vector at the barycenter -- uh_sol
 * in: vector -- char_vec
 * in: unit vector from barycenter to edge midpoint -- n
 * in: parameters      -- user_param
 * out: transformed vector -- vec
 */
void transf_from_char_var(double* uh_sol,
                          double* char_vec,
                          double* n,
                          parameters* user_param,
                          double* vec)
{
    // Transformation matrix
    double R[3][3];

    int i;
    double tol = 1.0e-05;

    switch ((*user_param).closure) {

	// P1
	case (1): {
	    if (SYSTEM_DIM == 3){
	        R[0][0] = -sqrt(3.0);  R[0][1] = sqrt(3.0); R[0][2] = 0.0;
	        R[1][0] = n[0];        R[1][1] = n[0];      R[1][2] = -n[1];
	        R[2][0] = n[1];        R[2][1] = n[1];      R[2][2] = n[0];

		for (i=0;i<3;i++) vec[i] = char_vec[0]*R[i][0] + char_vec[1]*R[i][1] + char_vec[2]*R[i][2];
	    }
	    break;
	}

	// K1
	case (2): {
	    if (SYSTEM_DIM == 3){
	        double N1x, N1y;

		/*ensure_realizability_ortho_proj(uh_sol);*/
		//ensure_realizability_vert_proj(1e-4,uh_sol);

		N1x = uh_sol[1]/uh_sol[0];
		N1y = uh_sol[2]/uh_sol[0];
		

		double nnx = n[0]*N1x + n[1]*N1y;
		double nny = n[0]*N1y - n[1]*N1x;

		double root = sqrt(3.0 - 2.0*nnx*nnx);

		double l1 = (2.0 * nnx + root)/3.0;
		double l2 = (2.0 * nnx - root)/3.0;

		double a1 = 2.0/3.0*nny*(1 - nnx/root);
		double a2 = 2.0/3.0*nny*(1 + nnx/root);

	        R[0][0] = 1.0 ;              R[0][1] = 1.0; 		  R[0][2] = 0.0;
	        R[1][0] = l1*n[0]-a1*n[1];   R[1][1] = l2*n[0]-a2*n[1];   R[1][2] = -n[1];
	        R[2][0] = l1*n[1]+a1*n[0];   R[2][1] = l2*n[1]+a2*n[0];   R[2][2] = n[0];

		for (i=0;i<3;i++) {
		    vec[i] = char_vec[0]*R[i][0] + char_vec[1]*R[i][1] + char_vec[2]*R[i][2];
		}
	    }
	    break;
	}

	// M1
	case (3): {
	    if (SYSTEM_DIM == 3){

		/*ensure_realizability_ortho_proj(uh_sol);*/
		//ensure_realizability_vert_proj(1e-6,uh_sol);

        int rel_flag=0;
        //rel_flag = is_realizable(uh_sol);
        
        
       if(!((*user_param).flag_ensure_positivity))
       {
           realizability_hack(tol,uh_sol);
       }

        
        if(rel_flag ==0)
        {
            //printf("rel_flag in transformation %d ============>solution %10.16e %10.16e %10.16e  \n",rel_flag,uh_sol[0],uh_sol[1],uh_sol[2]);
        }

		double N1x, N1y;

		if (uh_sol[0] > 1e-14){
		    N1x = uh_sol[1]/uh_sol[0];
		    N1y = uh_sol[2]/uh_sol[0];
		}
		else {
		    N1x = 0.0;
		    N1y = 0.0;
		}

		double nnx = n[0]*N1x + n[1]*N1y;
		double nny = n[0]*N1y - n[1]*N1x;

		int flag_use_table = 0;

		if (flag_use_table == 1) {
		    double f, phi;
		    f = sqrt(uh_sol[1]*uh_sol[1] + uh_sol[2]*uh_sol[2]) / uh_sol[0];
		    //phi = atan2(uh_sol[2],uh_sol[1]) + PI;
		    phi = atan2(nny,nnx) + PI;

		    int i,j;
		    i = round(f/(*user_param).fmax*((*user_param).nf-1));
		    j = round(phi/(2*PI)*((*user_param).nphi-1));
		
		    int r,c;
		    for (r=0;r<3;r++) {
			for (c=0;c<3;c++) {
			    R[r][c] = (*user_param).R[i][j][r][c];
			}
		    }

//R = [1 0 0; 0 nx -ny; 0 ny nx] * [R11(i,j) R12(i,j) R13(i,j); R21(i,j) R22(i,j) R23(i,j); R31(i,j) R32(i,j) R33(i,j)];

		    double trafo[3][3];

		    for (c=0;c<3;c++) {
			trafo[0][c] = R[0][c]; 
			trafo[1][c] = n[0]*R[1][c] - n[1]*R[2][c];
			trafo[2][c] = n[1]*R[1][c] + n[0]*R[2][c];
		    }

		    for (i=0;i<3;i++) vec[i] = char_vec[0]*trafo[i][0] + char_vec[1]*trafo[i][1] + char_vec[2]*trafo[i][2];
		    /*for (i=0;i<3;i++) vec[i] = char_vec[0]*R[i][0] + char_vec[1]*R[i][1] + char_vec[2]*R[i][2];*/

		    if (char_vec[0] != char_vec[0]) {
			printf("BACK!\n");
			printf("f = %lf\t phi = %lf\n",f, phi);
			printf("R = [%lf\t%lf\t%lf]\n",R[0][0],R[0][1],R[0][2]);
			printf("R = [%lf\t%lf\t%lf]\n",R[1][0],R[1][1],R[1][2]);
			printf("R = [%lf\t%lf\t%lf]\n",R[2][0],R[2][1],R[2][2]);
			printf("--------------------------------------------------------------\n");

			printf("R_trafo = [%lf\t%lf\t%lf]\n",trafo[0][0],trafo[0][1],trafo[0][2]);
			printf("R_trafo = [%lf\t%lf\t%lf]\n",trafo[1][0],trafo[1][1],trafo[1][2]);
			printf("R_trafo = [%lf\t%lf\t%lf]\n",trafo[2][0],trafo[2][1],trafo[2][2]);
			printf("--------------------------------------------------------------\n");

			printf("char_vec = [%lf\t%lf\t%lf]\n\n",char_vec[0],char_vec[1],char_vec[2]);
			printf("vec = [%lf\t%lf\t%lf]\n",vec[0],vec[1],vec[2]);
			getchar();
		    }
		}
		else {
		    double lambda1 = nnx*nnx + nny*nny;
		    double lambda2 = nny*nny - nnx*nnx;

		    R[0][0] = (2*lambda1*nny)/(pow(lambda1,2) + lambda2);              
		    R[0][1] = (-2*(pow(lambda1,2) - sqrt(lambda1)*lambda2))/((pow(lambda1,3) - 2*lambda2)*nny);
		    R[0][2] = (2*(pow(lambda1,2) + sqrt(lambda1)*lambda2))/((pow(lambda1,3) + 2*lambda2)*nny);

		    R[1][0] = (2*nnx*nny*n[0])/(pow(lambda1,2) + lambda2) - n[1];   
		    R[1][1] = (nnx*n[0])/nny - n[1];   
		    R[1][2] = (nnx*n[0])/nny - n[1];

		    R[2][0] = n[0] + (2*nnx*nny*n[1])/(pow(lambda1,2) + lambda2);   
		    R[2][1] = n[0] + (nnx*n[1])/nny;   
		    R[2][2] = n[0] + (nnx*n[1])/nny;

		    for (i=0;i<3;i++) {
			vec[i] = char_vec[0]*R[i][0] + char_vec[1]*R[i][1] + char_vec[2]*R[i][2];
		    }

		    /*if (vec[0] > 1e-4) {*/
		    if (vec[0] != vec[0]) {
			printf("BACK!\n");
			printf("n = [%lf \t %lf]\tnorm = %lf\n",n[0],n[1],sqrt(n[0]*n[0]+n[1]*n[1]));
			printf("nnx = %lf \t nny = %lf\n",nnx,nny);
			printf("lambda1 = %lf \t lambda2 = %lf\n",lambda1,lambda2);
			printf("char_vec = [%lf\t%lf\t%lf]\n",char_vec[0],char_vec[1],char_vec[2]);
			printf("vec = [%lf\t%lf\t%lf]\n\n",vec[0],vec[1],vec[2]);
		    }
		}
	    }
	    break;
	}

	// No transformation by default
	default: {
	    for (i=0;i<SYSTEM_DIM;i++) vec[i] = char_vec[i];
	    break;
	}
    }
}

/* void ensure_realizability_ortho_proj
 * ensure realizability for up two first order moments
 * inout: vector containing moments -- uh_sol
 */
void ensure_realizability_ortho_proj(double epsilon, double delta, double* uh_sol)
{
    /*double epsilon = 1e-9;*/
    /*double delta = 1e-9;*/
    double n1_norm,norm_new;
    double phi;

    if (SYSTEM_DIM > 1) {
	n1_norm = sqrt(uh_sol[1]*uh_sol[1] + uh_sol[2]*uh_sol[2]);

	if (uh_sol[0] < epsilon) {
	    uh_sol[0] = epsilon;
	    uh_sol[1] = 0.0;
	    uh_sol[2] = 0.0;
	    return;
	}
	else if (n1_norm > uh_sol[0]) {
	    uh_sol[0] = (uh_sol[0] + n1_norm)/2.0 + delta;

	    /*phi = atan2(uh_sol[2],uh_sol[1]) + PI;*/
	    /*uh_sol[1] = uh_sol[1] * cos(phi);*/
	    /*uh_sol[2] = uh_sol[2] * sin(phi);*/
	    /*double varphi = atan2(uh_sol[2],uh_sol[1]) + PI;*/

	    phi = (uh_sol[0] - 2*delta)/n1_norm;
	    uh_sol[1] = uh_sol[1] * phi;
	    uh_sol[2] = uh_sol[2] * phi;

	    /*norm_new = sqrt(uh_sol[1]*uh_sol[1] + uh_sol[2]*uh_sol[2]);*/
	    /*if (norm_new > uh_sol[0]) {*/
		/*printf("SHOULDN'T HAPPEN! ERROR IN ENSURE_REALIZABILITY\n");*/
		/*printf("uh_sol[0] = %lf\tuh_sol[1] = %lf\tuh_sol[2] = %lf\n",uh_sol[0],uh_sol[1],uh_sol[2]);*/
		/*printf("n1_norm = %lf\tn1_norm_new = %lf\tphi = %lf\n",n1_norm,norm_new,phi);*/
		/*getchar();*/
	    /*}*/

	    /*double varphi2 = atan2(uh_sol[2],uh_sol[1]) + PI;*/
	    /*if ( fabs(varphi2- varphi) > 1e-12) {*/
		/*printf("SHOULDN'T HAPPEN! ERROR IN ENSURE_REALIZABILITY\n");*/
		/*printf("uh_sol[0] = %lf\tuh_sol[1] = %lf\tuh_sol[2] = %lf\n",uh_sol[0],uh_sol[1],uh_sol[2]);*/
		/*printf("n1_norm = %lf\tn1_norm_new = %lf\tphi = %lf\n",n1_norm,norm_new,phi);*/
		/*printf("varphi_pre = %lf\tvarphi_post = %lf\n",varphi, varphi2);*/
		/*getchar();*/
	    /*}*/
	    return;
	}
    }
}


void realizability_hack(double tol,double* uh_sol)
{
    double norm_psi_one =0.0;
    double normalized_first_mom =0.0;


    if(uh_sol[0] < EPSILON)
    {
        uh_sol[0] = EPSILON;
    }
    norm_psi_one = sqrt(pow(uh_sol[1],2.0) + pow(uh_sol[2],2.0));
    normalized_first_mom =  norm_psi_one/uh_sol[0];

    if(normalized_first_mom > (1.0-tol))
    {

        //printf("before::: uh_sol %lf %lf %lf \n",uh_sol[0],uh_sol[1],uh_sol[2]);
        uh_sol[1] = (uh_sol[1]/normalized_first_mom)*(1.0-tol);
        uh_sol[2] = (uh_sol[2]/normalized_first_mom)*(1.0-tol);
        norm_psi_one = sqrt(pow(uh_sol[1],2.0) + pow(uh_sol[2],2.0));
    }
    else if(normalized_first_mom < EPSILON)
    {
        uh_sol[1] = 0.0;
        uh_sol[2] = 0.0;
    }
 
}




/* void ensure_realizability_vert_proj
 * ensure realizability for up two first order moments
 * in: 	  safety distance _in horizontal direction_ from boundary -- delta
 * inout: vector containing moments -- uh_sol
 */
void ensure_realizability_vert_proj(double delta, double* uh_sol)
{
    if (SYSTEM_DIM > 1) {
	if (uh_sol[0] < delta*sqrt(2)) {
	    uh_sol[0] = delta*sqrt(2);
	    uh_sol[1] = 0.0;
	    uh_sol[2] = 0.0;
	    return;
	}
	else{
	    double n1_norm = sqrt(uh_sol[1]*uh_sol[1] + uh_sol[2]*uh_sol[2]);

	    if (n1_norm > uh_sol[0]) {
		double f_old = n1_norm/uh_sol[0];
		double phi = (1-delta) / f_old;
		uh_sol[1] = uh_sol[1] * phi;
		uh_sol[2] = uh_sol[2] * phi;

		/*double norm_new = sqrt(uh_sol[1]*uh_sol[1] + uh_sol[2]*uh_sol[2]);*/
		/*if (norm_new > uh_sol[0]) {*/
		    /*printf("SHOULDN'T HAPPEN! ERROR IN ENSURE_REALIZABILITY\n");*/
		    /*printf("uh_sol[0] = %lf\tuh_sol[1] = %lf\tuh_sol[2] = %lf\n",uh_sol[0],uh_sol[1],uh_sol[2]);*/
		    /*printf("n1_norm = %lf\tn1_norm_new = %lf\tphi = %lf\n",n1_norm,norm_new,phi);*/
		    /*getchar();*/
		/*}*/
	    }
	}
    }
    return;
}

/* void ensure_positivity_sol
 * ensure positivity of the soluion
 * in: 	  current element -- E 
 * in: 	  minimum value for energy density -- delta
 * in: 	  number of basis function -- Nloc
 * inout: dg solution vector -- uh_sol
 */
void ensure_positivity_sol(double delta, 
			   int E, 
			   element* mesh_element,
			   node* mesh_node,
			   int Nloc,
			   double* uh_sol_vec)
{
    if (SYSTEM_DIM > 1) {
	double uh_E[SYSTEM_DIM];

	cell_average(uh_sol_vec,E,mesh_element,mesh_node,uh_E);

	if ( uh_E[0] < delta ) {
	
	    int s_dim, idofs;
	    for(s_dim=0;s_dim < SYSTEM_DIM;s_dim++) {	
		for(idofs=0;idofs<Nloc;idofs++)
		{
		    uh_sol_vec[SYSTEM_DIM*(E-1)*Nloc+(s_dim*Nloc+idofs)] = 0.0;
		}
	    }

	    uh_sol_vec[SYSTEM_DIM*(E-1)*Nloc] = delta;

	    cell_average(uh_sol_vec,E,mesh_element,mesh_node,uh_E);
	}
    }
    return;
}

/* void boundary_beam
 * Defines a beam of radiation
 */
void boundary_beam(double x,
		   double y,
		   double t,
		   double x_beam,
		   double y_beam,
		   double x_width,
		   double y_width,
		   double t_beam,
		   double ints,
		   double dir_x,
		   double dir_y,
		   double *f_val)
{
    double rho[4];

    /*rho[0] = 0.088622692545276; // int_-1^1 1  *exp(-100*(1-x)^2) */
    /*rho[1] = 0.083622692545276; // int_-1^1 x  *exp(-100*(1-x)^2) */
    /*rho[2] = 0.079065806008002; // int_-1^1 x^2*exp(-100*(1-x)^2) */

    int s_dim;

    rho[0] = 1.0;
    rho[1] = 0.99;
    rho[2] = 0.99;
    rho[3] = 0.99;

    double n[2];
    double norm = sqrt(dir_x*dir_x + dir_y*dir_y);
    n[0] = dir_x / norm;
    n[1] = dir_y / norm;
    

    if ( (x >= x_beam - x_width/2.0) && (x <= x_beam + x_width/2.0) 
       &&(y >= y_beam - y_width/2.0) && (y <= y_beam + y_width/2.0) ) {
	if (SYSTEM_DIM > 2) {
	    f_val[0] = f_val[0] + rho[0]*ints*exp(-200*(t-t_beam)*(t-t_beam));
	    f_val[1] = f_val[1] + n[0]*rho[1]*ints*exp(-200*(t-t_beam)*(t-t_beam));
	    f_val[2] = f_val[2] + n[1]*rho[1]*ints*exp(-200*(t-t_beam)*(t-t_beam));

	    if (SYSTEM_DIM > 5) {
		for (s_dim=3; s_dim < 6; s_dim++) f_val[s_dim] = 0.0;
		/*f_val[3] = f_val[3] + n[0]*n[0]*rho[2]*ints*exp(-200*(t-t_beam)*(t-t_beam));*/
		/*f_val[4] = f_val[4] + n[0]*n[1]*rho[2]*ints*exp(-200*(t-t_beam)*(t-t_beam));*/
		/*f_val[5] = f_val[5] + n[1]*n[1]*rho[2]*ints*exp(-200*(t-t_beam)*(t-t_beam));*/

		if (SYSTEM_DIM > 9) {
		    for (s_dim=6; s_dim < 10; s_dim++) f_val[s_dim] = 0.0;
		    /*f_val[6] = f_val[6] + n[0]*n[0]*n[0]*rho[3]*ints*exp(-200*(t-t_beam)*(t-t_beam));*/
		    /*f_val[7] = f_val[7] + n[0]*n[0]*n[1]*rho[3]*ints*exp(-200*(t-t_beam)*(t-t_beam));*/
		    /*f_val[8] = f_val[8] + n[0]*n[1]*n[1]*rho[3]*ints*exp(-200*(t-t_beam)*(t-t_beam));*/
		    /*f_val[9] = f_val[9] + n[1]*n[1]*n[1]*rho[3]*ints*exp(-200*(t-t_beam)*(t-t_beam));*/
		}
	    }
	}
    }
}

/* double sign
 * returnes the sign of a double
 */
double sign(double x)
{
    if (x > 0) return 1;
    if (x < 0) return -1;
    return 0;
}

/* void estimate_largest_eval
 * in: solution vector 	     -- uh_sol
 * in: normal vector on edge -- n
 * in: parameters      	     -- user_param
 * out: estimate of the largest absolute value of the eigenvectors of the Jacobian -- vec
 */
double estimate_largest_eval(double* uh_sol,
                             double* n,
			     parameters* user_param)
{

    switch ((*user_param).closure) {

	// P1
	case (1): {
	    if (SYSTEM_DIM == 3){
	        return(1./sqrt(3.));
	    }
	}

	// K1
	case (2): {
	    if (SYSTEM_DIM == 3){
	        double N1x, N1y;

		ensure_realizability_ortho_proj(1e-9,1e-6,uh_sol);

		if (uh_sol[0] > 1e-16){
		    N1x = uh_sol[1]/uh_sol[0];
		    N1y = uh_sol[2]/uh_sol[0];
		}
		else {
		    N1x = 0.0;
		    N1y = 0.0;
		}

		double root = sqrt(-0.2e1 * n[0] * n[0] * N1x * N1x + 0.2e1 * n[0] * N1x * n[1] * N1y - 0.2e1 * n[1] * n[1] * N1y * N1y - 0.3e1 * n[1] * n[1] * N1x * N1x + 0.3e1 - 0.3e1 * n[0] * n[0] * N1y * N1y) / 0.3e1;

		double t = 0.2e1 / 0.3e1 * n[0] * N1x + 0.2e1 / 0.3e1 * n[1] * N1y;

		return(max(t+root,t-root));
	    }
	    break;
	}

	// M1
	case (3): {
	    if (SYSTEM_DIM == 3){
		/*double N1x, N1y;*/

		/*ensure_realizability_ortho_proj(uh_sol);*/

		/*if (uh_sol[0] > 1e-12){*/
		    /*N1x = uh_sol[1]/uh_sol[0];*/
		    /*N1y = uh_sol[2]/uh_sol[0];*/
		/*}*/
		/*else {*/
		    /*N1x = 0.0;*/
		    /*N1y = 0.0;*/
		/*}*/
		return(1.);
	    }
	    break;
	}

	// No transformation by default
	default: {
	    return(1.);
	}
    }
    return(1.);
}

void transform_to_char_var(double* uh_sol,
		            	double* vec,
                        double* normal,
                        parameters* user_param,
                        double* transformed_vec,
                        int flag)
{
    double f =0.0;      //norm f
    double chi_f=0.0;
    double xi =0.0;
    double F[2];        //flux
    double E =0.0;      //first moment
    double chi_f_prime =0.0;
    double fx,fy=0.0;
    int s=0;


    double JX[3][3];
    double JY[3][3];
    double RMAT[3][3];
    double J[3][3];
    double R[3][3];
    double R_inv[3][3];
    double inv_R[9];
    int i,j=0;
    double theta=0.0;

    double J_2D[2][2];
    double R_2D[2][2];
    double R_2D_inv[2][2];
    double inv_R_2D[4];

    int flag_2d =0;
    int flag_3d =0;
    int fx_zero_flag =0;
    int fy_zero_flag =0;

    int rel_flag=0;
    double tol = 1.0e-06;
    double nfirst_mom =0.0;

    //rel_flag = is_realizable(uh_sol);
    
    
                
    E = uh_sol[0];
    F[0] = uh_sol[1];
    F[1] = uh_sol[2];


    fx = F[0]/E;
    fy = F[1]/E;



    f = sqrt(pow(fx,2.0) + pow(fy,2.0));
  
    if((E > CHAR_EPSILON) && (fabs(fx) > CHAR_EPSILON) && (fabs(fy) > CHAR_EPSILON) && (f > CHAR_EPSILON))
    {
        if(!((*user_param).flag_ensure_positivity))
        {
            if(uh_sol[0]  <tol)
            {
                E = tol;
            }
            nfirst_mom = f;
            if(nfirst_mom > (1.0-tol))
            {
                //F[0] = (F[0]/f)*(1.0-tol);
                //F[1] = (F[1]/f)*(1.0-tol);

                fx = (fx/f)*(1-tol);
                fy = (fy/f)*(1-tol);

                //fx = F[0]/E;
                //fy = F[1]/E;
                f = sqrt(pow(fx,2.0) + pow(fy,2.0));
            }
        }
        nfirst_mom = f/E;

        //printf("nfirst_mom %lf \n",nfirst_mom);
        //assert(nfirst_mom <1.0);
   

      
        rel_flag = is_realizable(uh_sol);

        if(rel_flag ==0)
        {
            //printf(" REL_FLAG 0 E %10.16e fx %10.16e fy %10.16e \n",E,fx,fy);
            //exit(1);
        }

        //assert(rel_flag);
        //printf(" CASE 1:::E %10.16e fx %10.16e fy %10.16e rel_flag %d \n",E,fx,fy,rel_flag);

        assert(f>0);
      
 
        if((4.0-3.0*pow(f,2.0)<CHAR_EPSILON))
        {
            xi = CHAR_EPSILON;
            printf("cut off \n");
            //exit(1);   
        }
        else
        {
            xi = sqrt(4.0 - 3.0*pow(f,2.0));
        }
        //printf("xi %lf f %lf \n",xi,f);
        chi_f = (3.0 + 4.0*pow(f,2.0))/(5.0+2.0*xi);
        
        chi_f_prime = (2.0*f)/xi;
        //chi_f_prime = ((5.0+2.0*xi)*(8.0*f) + (3.0+4.0*pow(f,2.0))*((6.0*f)/xi))/(pow((5.0+2.0*xi),2.0));
        
        theta = 2.0+3.0*f*chi_f_prime - 6.0*chi_f;
       

        double n = fx*fx+fy*fy; //norm(f)^2
        //double xi = sqrt(4.0-3.0*n);
        double chi = (3.0+4.0*n)/(5.0+2.0*xi);
        double dchi = 2.0*f/xi; //dchi/df
        double nx = fx/f;
        double ny = fy/f;
        double a = chi-f*dchi;
        double b = dchi;
       
            
            
        double dPxxdE = (1.0-a)/2.0+(3.0*a-1.0)/2.0*nx*nx;
        double dPxydE = (3.0*a-1.0)/2.0*nx*ny;
        double dPyydE = (1.0-a)/2.0+(3.0*a-1.0)/2.0*ny*ny;
        double dPxxdFx = (-b/2.0+3.0*b/2.0*nx*nx+(3.0*a+3.0*f*b-1.0)/f*ny*ny)*nx;
        double dPxydFx = (-(3.0*a-1.0)/2.0/f*nx*nx+(3.0*a-1.0+3.0*f*b)/2.0/f*ny*ny)*ny;
        double dPyydFx = (-b/2.0-(6.0*a+3.0*f*b-2.0)/2.0/f*ny*ny)*nx;
        double dPxxdFy = (-b/2.0-(6.0*a+3.0*f*b-2.0)/2.0/f*nx*nx)*ny;
        double dPxydFy = ((3.0*a-1.0+3.0*f*b)/2.0/f*nx*nx-(3.0*a-1.0)/2.0/f*ny*ny)*nx;
        double dPyydFy = (-b/2.0+3.0*b/2.0*ny*ny+(3.0*a+3.0*f*b-1.0)/f*nx*nx)*ny;
        

        JX[0][0] = 0.0;
        JX[0][1] = 1.0;
        JX[0][2] = 0.0;

        JX[1][0] = dPxxdE;
        JX[1][1] = dPxxdFx;
        JX[1][2] = dPxxdFy;

        JX[2][0] = dPxydE;
        JX[2][1] = dPxydFx;
        JX[2][2] = dPxydFy;





        
        /*JX[1][0] = ((2.0*pow(f,2.0) - 3.0*pow(fy,2.0))*(chi_f - f*chi_f_prime) + pow(fy,2.0))/(2.0*pow(f,2.0));
        JX[1][1] = -1.0*((fx*pow(fy,2.0))/(2.0*pow(f,4.0)))*theta + (fx*chi_f_prime)/f;
        JX[1][2] = ((pow(fx,2.0)*fy)/(2.0*pow(f,4.0)))*theta - (pow(fy,2.0)*chi_f_prime)/(2.0*f);

        
        JX[2][0] = ((fx*fy)/(2.0*pow(f,2.0)))*(3.0*chi_f - 3.0*f*chi_f_prime-1.0);
        JX[2][1] = ((pow(fx,2.0)*fy)/(2.0*pow(f,4.0)))*theta + (fy/(2.0*pow(f,2.0)))*(3.0*chi_f-1.0);
        JX[2][2] = ((fx*pow(fy,2.0))/(4.0*pow(f,4.0)))*theta + (fx/(2.0*pow(f,2.0)))*(3.0*chi_f-1.0);
        */



        JY[0][0] = 0.0;
        JY[0][1] = 0.0;
        JY[0][2] = 1.0;

        JY[1][0] = dPxydE;
        JY[1][1] = dPxydFx;
        JY[1][2] = dPxydFy;

        JY[2][0] = dPyydE;
        JY[2][1] = dPyydFx;
        JY[2][2] = dPyydFy;





        /*JY[1][0] = ((fx*fy)/(2.0*pow(f,2.0)))*(3.0*chi_f-3*f*chi_f_prime -1.0);
        JY[1][1] = ((fy*(3.0*chi_f-1.0))/(2.0*pow(f,2.0))) + ((pow(fx,2.0)*fy)/(4.0*pow(f,4.0)))*theta;
        JY[1][2] = (fx*(3.0*chi_f-1.0))/(2.0*pow(f,2.0)) +  ((fx*pow(fy,2.0))/(2.0*pow(f,4.0)))*theta;

        JY[2][0] = ((2.0*pow(f,2.0) -3.0*pow(fx,2.0))*(chi_f- f*chi_f_prime) + pow(fx,2.0))/(pow(f,2.0));
        JY[2][1] = -1.0*((fx*chi_f_prime)/(2.0*f)) + ((fx*pow(fy,2.0))/(2.0*pow(f,4.0)))*theta;
        JY[2][2] = (fy*chi_f_prime)/f - ((pow(fx,2.0)*fy)/(2.0*pow(f,4.0)))*theta;*/
        
        for(i=0;i<3;i++)
        {
            for(j=0;j<3;j++)
            {
                J[i][j] = normal[0]*JX[i][j] + normal[1]*JY[i][j];
            }
        }

        
        flag_3d =1;
    }
   /*else if((E > CHAR_EPSILON) && (((fabs(fx) < CHAR_EPSILON) && (fabs(fy) > CHAR_EPSILON)) || ((fabs(fy) < CHAR_EPSILON) && (fabs(fx) > CHAR_EPSILON)))) 
   {

        rel_flag = is_realizable(uh_sol);
        //assert(rel_flag);

       if(rel_flag ==0)
       {
           printf(" WARNING :: CASE 2:::::::E %10.16e fx %10.16e fy %10.16e rel_flag %d \n",E,fx,fy,rel_flag);
       }


       if(fabs(fx) < CHAR_EPSILON)
       {
           f = fabs(fy);
           fx_zero_flag =1;
       }
       else
       {
           f=fabs(fx);
           fy_zero_flag =1;
       }

      

       xi = sqrt(4.0 - 3.0*pow(f,2.0));
       if((4.0-3.0*pow(f,2.0)<0))
       {
           xi =CHAR_EPSILON;
       }
       else
       {
            xi = sqrt(4.0 - 3.0*pow(f,2.0));
       }
       
       chi_f = (3.0+4.0*pow(f,2.0))/(5.0+2.0*xi);
     

     
 
       //printf("f %lf xi %lf \n",f,xi);

       assert(xi>0);
       chi_f_prime = (2.0*f)/xi;

       if(fx_zero_flag)
       {
           J_2D[0][0] = normal[1]*0.0;
           J_2D[0][1] = normal[1]*1.0;

           J_2D[1][0] = normal[1]*chi_f*chi_f_prime;
           J_2D[1][1] = normal[1]*chi_f_prime;
           flag_2d =1;
       }
       else
       {
           J_2D[0][0] = normal[0]*0.0;
           J_2D[0][1] = normal[0]*1.0;

           J_2D[1][0] = normal[0]*chi_f*chi_f_prime;
           J_2D[1][1] = normal[0]*chi_f_prime;
           flag_2d =1;
       }


   }*/
   else
   {

        rel_flag = is_realizable(uh_sol);
        /*if(rel_flag ==0)
        {
            printf(" CASE 3:::::::E %10.16e fx %10.16e fy %10.16e rel_flag %d \n",E,fx,fy,rel_flag);
        }*/
      
       //printf(" CASE 3:::::::E %10.16e fx %10.16e fy %10.16e rel_flag %d \n",E,fx,fy,rel_flag);
 
        //assert(rel_flag);
        JX[0][0] = 0.0;
        JX[0][1] = 1.0;
        JX[0][2] = 0.0;

        JX[1][0] = 1.0/3.0;
        JX[1][1] = 0.0;
        JX[1][2] = 0.0;

        JX[2][0] = 0.0;
        JX[2][1] = 0.0;
        JX[2][2] = 0.0;

        JY[0][0] = 0.0;
        JY[0][1] = 0.0;
        JY[0][2] = 1.0;

        JY[1][0] = 0.0;
        JY[1][1] = 0.0;
        JY[1][2] = 0.0;

        JY[2][0] = 1.0/3.0;
        JY[2][1] = 0.0;
        JY[2][2] = 0.0;
        flag_3d =1;


        for(i=0;i<3;i++)
        {
            for(j=0;j<3;j++)
            {
                J[i][j] = normal[0]*JX[i][j] + normal[1]*JY[i][j];
            }
        }

        flag_3d =1;
    }
   if(flag_3d)
   {
       
 

        double J_vec[]={J[0][0],J[0][1],J[0][2],
           J[1][0],J[1][1],J[1][2],
           J[2][0],J[2][1],J[2][2]};

        gsl_matrix_view JAC = gsl_matrix_view_array(J_vec,3,3);

        gsl_vector_complex *eval = gsl_vector_complex_alloc(3);
        gsl_matrix_complex *evec = gsl_matrix_complex_alloc(3,3);


        gsl_eigen_nonsymmv_workspace *w = gsl_eigen_nonsymmv_alloc(3);
        gsl_eigen_nonsymmv(&JAC.matrix, eval, evec,w);
        gsl_eigen_nonsymmv_free(w);

        /*printf("normal %lf %lf \n",n[0],n[1]);
        printf("J = \n");
        for(i=0;i<3;i++)
        {
            for(j=0;j<3;j++)
            {
                printf("%lf ",J[i][j]);
            }
            printf("\n");
        }
        printf("] \n");*/
        

        
        //printf("R =[");
        gsl_eigen_nonsymmv_sort(eval,evec,GSL_EIGEN_SORT_ABS_DESC);
        {
            for(i=0;i<3;i++)
            {
                gsl_complex eval_i = gsl_vector_complex_get(eval,i);
                gsl_vector_complex_view evec_i = gsl_matrix_complex_column(evec,i);
                for(j=0;j<3;j++)
                {
                    gsl_complex z = gsl_vector_complex_get(&evec_i.vector,j);
                    R[j][i] = GSL_REAL(z); 
                   // printf("%lf ",R[j][i]);
                }
               //printf("\n");

            }
            //printf(" ] \n");
        }

        double R_data[] = {R[0][0],R[0][1],R[0][2],
                           R[1][0],R[1][1],R[1][2],
                           R[2][0],R[2][1],R[2][2]};

        gsl_matrix_view rmat = gsl_matrix_view_array(R_data,3,3);
        gsl_matrix_view inv_rmat = gsl_matrix_view_array(inv_R,3,3);
        gsl_permutation *p = gsl_permutation_alloc(3);

        gsl_linalg_LU_decomp(&rmat.matrix,p,&s);
        gsl_linalg_LU_invert(&rmat.matrix,p,&inv_rmat.matrix);

        for(i=0;i<3;i++)
        {
            for(j=0;j<3;j++)
            {
                R_inv[i][j] = gsl_matrix_get(&inv_rmat.matrix,i,j);
            }
        }

        if(flag == TRANS_FROM_CHAR)
        {
            for (i=0;i<3;i++) 
            {
                transformed_vec[i] = vec[0]*R[i][0] +vec[1]*R[i][1] + vec[2]*R[i][2];
                if(isnan(transformed_vec[i]))
                {
                    printf("tvec %10.16e :::: vec %10.16e \n",transformed_vec[i],vec[i]);
                    getchar();
                }
            }
        }
        else if(flag == TRANS_TO_CHAR)
        {
            for (i=0;i<3;i++) 
            {
                transformed_vec[i] = vec[0]*R_inv[i][0] + vec[1]*R_inv[i][1] + vec[2]*R_inv[i][2];
               if(isnan(transformed_vec[i]))
                {
                    printf("tvec %10.16e :::: vec %10.16e \n",transformed_vec[i],vec[i]);
                    getchar();
                }
            }
        }
        gsl_vector_complex_free(eval);
        gsl_matrix_complex_free(evec);
        gsl_permutation_free(p);
   }//3d_flag
   else if(flag_2d)
   {
       double J2D_vec[]={J_2D[0][0],J_2D[0][1],
                         J_2D[1][0],J_2D[1][1]};

        gsl_matrix_view JAC_2D = gsl_matrix_view_array(J2D_vec,2,2);

        gsl_vector_complex *eval_2d = gsl_vector_complex_alloc(2);
        gsl_matrix_complex *evec_2d = gsl_matrix_complex_alloc(2,2);

        gsl_eigen_nonsymmv_workspace *w_2d = gsl_eigen_nonsymmv_alloc(2);
        gsl_eigen_nonsymmv(&JAC_2D.matrix, eval_2d, evec_2d,w_2d);
        gsl_eigen_nonsymmv_free(w_2d);

        gsl_eigen_nonsymmv_sort(eval_2d,evec_2d,GSL_EIGEN_SORT_ABS_DESC);
        {
          
          // printf("R2d =[");
 
            for(i=0;i<2;i++)
            {
                gsl_complex eval_2d_i = gsl_vector_complex_get(eval_2d,i);
                gsl_vector_complex_view evec_2d_i = gsl_matrix_complex_column(evec_2d,i);
                for(j=0;j<2;j++)
                {
                    gsl_complex z_2d = gsl_vector_complex_get(&evec_2d_i.vector,j);
                  
                    R_2D[j][i] = GSL_REAL(z_2d); 
                //    printf("%lf ",R_2D[j][i]);
                    
                    
                }
              //  printf("\n");

            }
            //printf("]\n");
        }
        
        double R_data [] = {R_2D[0][0],R_2D[0][1],
                            R_2D[1][0],R_2D[1][1]};
        
        gsl_matrix_view rmat = gsl_matrix_view_array(R_data,2,2);
        gsl_matrix_view inv_rmat = gsl_matrix_view_array(inv_R_2D,2,2);
        gsl_permutation *p = gsl_permutation_alloc(2);

        gsl_linalg_LU_decomp(&rmat.matrix,p,&s);
        gsl_linalg_LU_invert(&rmat.matrix,p,&inv_rmat.matrix);
        
        
        for(i=0;i<2;i++)
        {
            for(j=0;j<2;j++)
            {
                R_2D_inv[i][j] = gsl_matrix_get(&inv_rmat.matrix,i,j);
            }
        }
        //printf("Rinv computed \n'");
         
        if(flag == TRANS_FROM_CHAR)
        {
            if(fx_zero_flag)
            {
                transformed_vec[0] = vec[0]*R_2D[0][0] + vec[2]*R_2D[0][1];
                transformed_vec[1] = vec[1];
                transformed_vec[2] = vec[0]*R_2D[1][0] + vec[2]*R_2D[1][1];
            }
            else if(fy_zero_flag)
            {
                transformed_vec[0] = vec[0]*R_2D[0][0] + vec[1]*R_2D[0][1];
                transformed_vec[1] = vec[0]*R_2D[1][0] + vec[1]*R_2D[1][1];
                transformed_vec[2] = vec[2];
            }
            
            if(isnan(transformed_vec[i]))
            {
                printf("tvec %10.16e :::: vec %10.16e \n",transformed_vec[i],vec[i]);
                //getchar();
            }
            
        }
        else if(flag == TRANS_TO_CHAR)
        {
            if(fx_zero_flag)
            {
                transformed_vec[0] = vec[0]*R_2D_inv[0][0] + vec[2]*R_2D_inv[0][1];
                transformed_vec[1] = vec[1];
                transformed_vec[2] = vec[0]*R_2D_inv[1][0] + vec[2]*R_2D_inv[1][1];
            }
            else if(fy_zero_flag)
            {
                transformed_vec[0] = vec[0]*R_2D_inv[0][0] + vec[1]*R_2D_inv[0][1];
                transformed_vec[1] = vec[0]*R_2D_inv[1][0] + vec[1]*R_2D_inv[1][1];
                transformed_vec[2] = vec[2];
            }                
                   /* printf("tvec %10.16e :::: vec %10.16e \n",transformed_vec[i],vec[i]);
                     printf(" E %10.16e fx %10.16e fy %10.16e \n",E,fx,fy);
                     getchar();*/
        }

    gsl_vector_complex_free(eval_2d);
    gsl_matrix_complex_free(evec_2d);
    gsl_permutation_free(p);
   }    
}

