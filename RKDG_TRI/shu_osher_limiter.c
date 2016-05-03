#include"shu_osher_limiter.h"
void min_mod_limiter(int K0,
                     element* mesh_element,
                     node* mesh_node,
                     edge* mesh_edge,
                     parameters* user_param,
                     int Nelts,
                     int Nloc,
                     double* uh_linear_unlimited,
                     double* uh_limit)
{
    //barycenter of triangles
    double b0[2];
    double b1[2];
    double b2[2];
    double b3[2];
    
    double ref_b0[2];
    double ref_b1[2];
    double ref_b2[2];
    double ref_b3[2];

    //mid points
    double m1[2];
    double m2[2];
    double m3[2];
    //ref coords
    double ref_m1[2];
    double ref_m2[2];
    double ref_m3[2];

    //neighbouring elements
    int K1,K2,K3; 
    double uh_K0_bar[SYSTEM_DIM];
    double uh_K1_bar[SYSTEM_DIM];
    double uh_K2_bar[SYSTEM_DIM];
    double uh_K3_bar[SYSTEM_DIM];
    
    double uhm1_K0[SYSTEM_DIM]; //linear solution at midpoint edge1
    double uhm2_K0[SYSTEM_DIM]; //linear solution at midpoint edge2
    double uhm3_K0[SYSTEM_DIM]; //linear solution at midpoint edge3

    
    int edge1,edge2,edge3;
    
    double alphas_edge1[2];
    double alphas_edge2[2];
    double alphas_edge3[2];

    int alpha1_elements[2];
    int alpha2_elements[2];
    int alpha3_elements[2];

    double tilde_uh_m1_K0[SYSTEM_DIM];
    double tilde_uh_m2_K0[SYSTEM_DIM];
    double tilde_uh_m3_K0[SYSTEM_DIM];

    double transf_tilde_uh_m1_K0[SYSTEM_DIM];
    double transf_tilde_uh_m2_K0[SYSTEM_DIM];
    double transf_tilde_uh_m3_K0[SYSTEM_DIM];


    double delta_u_bar_m1_K0[SYSTEM_DIM];
    double delta_u_bar_m2_K0[SYSTEM_DIM];
    double delta_u_bar_m3_K0[SYSTEM_DIM];

    double transf_delta_u_bar_m1_K0[SYSTEM_DIM];
    double transf_delta_u_bar_m2_K0[SYSTEM_DIM];
    double transf_delta_u_bar_m3_K0[SYSTEM_DIM];


    double delta_1[SYSTEM_DIM];
    double delta_2[SYSTEM_DIM];
    double delta_3[SYSTEM_DIM];

    double transf_delta_1[SYSTEM_DIM];
    double transf_delta_2[SYSTEM_DIM];
    double transf_delta_3[SYSTEM_DIM];


    double delta_1_hat=0.0;
    double delta_2_hat=0.0;
    double delta_3_hat=0.0;
   
    double tdelta_1_hat[SYSTEM_DIM];
    double tdelta_2_hat[SYSTEM_DIM];
    double tdelta_3_hat[SYSTEM_DIM];

    double test_debug=0.0;

    double vec_m1minusb0[2];
    double vec_m2minusb0[2];
    double vec_m3minusb0[2];


    

    double* limited_dg_coefs = NULL;

    double* limited_dg_coefs1 = NULL;

    Nloc = (mesh_element[K0].degree+1)*(mesh_element[K0].degree+2)/2;
    //printf("Nloc %d \n",Nloc);getchar();
    limited_dg_coefs1 = (double*)malloc(Nloc*sizeof(double));

    limited_dg_coefs = (double*)malloc(Nloc*sizeof(double));


    int idofs=0;

    int s_dim=0;
    double pos=0.0;
    double neg=0.0;

    double theta_plus=0.0;
    double theta_minus=0.0;

    double grad_x=0.0;
    double grad_y=0.0;

    double grad=0.0;
    double limit_grad=0.0;

    double limit_grad_x =0.0;
    double limit_grad_y =0.0;
    double diff;

    edge1 = mesh_element[K0].edge[1];
    edge2 = mesh_element[K0].edge[2];
    edge3 = mesh_element[K0].edge[3];


    //mid points
    get_edge_mid_point(mesh_edge,edge1,mesh_node,m1);
    get_edge_mid_point(mesh_edge,edge2,mesh_node,m2);
    get_edge_mid_point(mesh_edge,edge3,mesh_node,m3);

    


    //barycenter of element E
    get_barycenter(K0,mesh_element,mesh_node,b0);
    
    //solution K0 @ b0
    map_to_reference_element(mesh_element,mesh_node,K0,ref_b0,b0[0],b0[1]);
    get_approx_sol_linear_basis(K0,mesh_element,mesh_node,uh_linear_unlimited,ref_b0[0],ref_b0[1],uh_K0_bar);


    //solution K0 @ mid points
    map_to_reference_element(mesh_element,mesh_node,K0,ref_m1,m1[0],m1[1]);
    map_to_reference_element(mesh_element,mesh_node,K0,ref_m2,m2[0],m2[1]);
    map_to_reference_element(mesh_element,mesh_node,K0,ref_m3,m3[0],m3[1]);
    
    get_approx_sol_linear_basis(K0,mesh_element,mesh_node,uh_linear_unlimited,ref_m1[0],ref_m1[1],uhm1_K0);
    get_approx_sol_linear_basis(K0,mesh_element,mesh_node,uh_linear_unlimited,ref_m2[0],ref_m2[1],uhm2_K0);
    get_approx_sol_linear_basis(K0,mesh_element,mesh_node,uh_linear_unlimited,ref_m3[0],ref_m3[1],uhm3_K0);

    double dummy[SYSTEM_DIM];
    get_approx_solution(uh_limit,K0,mesh_element,mesh_node,ref_m1[0],ref_m1[1],dummy);
    //assert(fabs(dummy[0] - uhm1_K0[0]) < 1e-10);

    int i;

    double dg_sol_av1[SYSTEM_DIM];
    dg_quadratic_basis_cell_average(K0,mesh_element,mesh_node,uh_limit,dg_sol_av1);


    //edge1 
    if(mesh_edge[edge1].edge_type == INTERIOR)
    {
        if(mesh_edge[edge1].neighbour[1] ==K0)
        {
            K1 = mesh_edge[edge1].neighbour[2];
        }
        else
        {
            K1 = mesh_edge[edge1].neighbour[1];
        }

        get_barycenter(K1,mesh_element,mesh_node,b1);
        map_to_reference_element(mesh_element,mesh_node,K1,ref_b1,b1[0],b1[1]);
        get_approx_sol_linear_basis(K1,mesh_element,mesh_node,uh_linear_unlimited,ref_b1[0],ref_b1[1],uh_K1_bar);


     //  printf("iedge @ interior K1 %d  b1 (%lf %lf ) m1 (%lf %lf) \n",K1,b1[0],b1[1],m1[0],m1[1]); getchar();

    }
    else if(mesh_edge[edge1].edge_type == EXTERIOR)
    {
        if((*user_param).periodic)
        {
            K1 = mesh_edge[edge1].neighbour[2];
            //get barycenter for K1 (physical)
            get_barycenter(K1,mesh_element,mesh_node,b1);

            //printf("iedge @ exterior K1 %d  b1 (%lf %lf ) m1 (%lf %lf) \n",K1,b1[0],b1[1],m1[0],m1[1]); //getchar();
            
            //solution K1 @ b1
            map_to_reference_element(mesh_element,mesh_node,K1,ref_b1,b1[0],b1[1]);
            get_approx_sol_linear_basis(K1,mesh_element,mesh_node,uh_linear_unlimited,ref_b1[0],ref_b1[1],uh_K1_bar);

            //modify K1 so it is ghost neighbour of K0 to compute alphas
            if(mesh_edge[edge1].slope == VERTICAL)
            {
         //       printf("iedge vertical \n");
                //ghost element is left => translate center right
                if(mesh_edge[edge1].boundary_side== RIGHT)
                {
                    b1[0] = b1[0]+ (*user_param).periodic_domain_size_x;
                }
                else if(mesh_edge[edge1].boundary_side == LEFT)
                {
                    b1[0] = b1[0]- (*user_param).periodic_domain_size_x;
                }
            }
            else if(mesh_edge[edge1].slope == HORIZONTAL)
            {
                if(mesh_edge[edge1].boundary_side == BOTTOM)
                {
                    b1[1] = b1[1] - (*user_param).periodic_domain_size_y;
                }
                else if(mesh_edge[edge1].boundary_side == TOP)
                {
                    b1[1] = b1[1] + (*user_param).periodic_domain_size_y;
                }
            }         
        //printf("b1 (%lf %lf ) \n",b1[0],b1[1]); getchar();
        }//periodic
    }
    //edge2
    if(mesh_edge[edge2].edge_type == INTERIOR)
    {
        if(mesh_edge[edge2].neighbour[1] ==K0)
        {
            K2 = mesh_edge[edge2].neighbour[2];
        } 
        else
        {
            K2 = mesh_edge[edge2].neighbour[1];
        }

        get_barycenter(K2,mesh_element,mesh_node,b2);
        map_to_reference_element(mesh_element,mesh_node,K2,ref_b2,b2[0],b2[1]);
        get_approx_sol_linear_basis(K2,mesh_element,mesh_node,uh_linear_unlimited,ref_b2[0],ref_b2[1],uh_K2_bar);

        //printf("K2 %d \n",K2); getchar();
        //printf("iedge @ interior K2 %d  b2 (%lf %lf ) m2 (%lf %lf ) \n",K2,b2[0],b2[1],m2[0],m2[1]); getchar();
    }
    else if(mesh_edge[edge2].edge_type == EXTERIOR)
    {
        if((*user_param).periodic)
        {
            K2 = mesh_edge[edge2].neighbour[2];
            //get barycenter for K2 (physical)
            get_barycenter(K2,mesh_element,mesh_node,b2);

           //printf("iedge @ exterior K2 %d  b2 (%lf %lf ) m2 (%lf %lf) \n",K2,b2[0],b2[1],m2[0],m2[1]); //getchar();
            //solution K2 @ b2
            map_to_reference_element(mesh_element,mesh_node,K2,ref_b2,b2[0],b2[1]);
            get_approx_sol_linear_basis(K2,mesh_element,mesh_node,uh_linear_unlimited,ref_b2[0],ref_b2[1],uh_K2_bar);

            //modify K1 so it is ghost neighbour of K0 to compute alphas
            if(mesh_edge[edge2].slope == VERTICAL)
            {
                //ghost element is left => translate center right
                if(mesh_edge[edge2].boundary_side== RIGHT)
                {
                    b2[0] = b2[0]+ (*user_param).periodic_domain_size_x;
                }
                else if(mesh_edge[edge2].boundary_side == LEFT)
                {
                    b2[0] = b2[0]- (*user_param).periodic_domain_size_x;
                }
            }
            else if(mesh_edge[edge2].slope == HORIZONTAL)
            {
                if(mesh_edge[edge2].boundary_side == BOTTOM)
                {
                    b2[1] = b2[1] - (*user_param).periodic_domain_size_y;
                }
                else if(mesh_edge[edge2].boundary_side == TOP)
                {
                    b2[1] = b2[1] + (*user_param).periodic_domain_size_y;
                }
            }         
        }//periodic
    }
    //edge3
    if(mesh_edge[edge3].edge_type == INTERIOR)
    {

        if(mesh_edge[edge3].neighbour[1] ==K0)
        {
            K3 = mesh_edge[edge3].neighbour[2];
        }
        else
        {
            K3 = mesh_edge[edge3].neighbour[1];
        }

        get_barycenter(K3,mesh_element,mesh_node,b3);
        //solution K3 @ b3
        map_to_reference_element(mesh_element,mesh_node,K3,ref_b3,b3[0],b3[1]);
        get_approx_sol_linear_basis(K3,mesh_element,mesh_node,uh_linear_unlimited,ref_b3[0],ref_b3[1],uh_K3_bar);


        //printf("iedge @ interior K3 %d  b3 (%lf %lf ) m3 (%lf %lf ) \n",K3,b3[0],b3[1],m3[0],m3[1]); getchar();
    }
    else if(mesh_edge[edge3].edge_type == EXTERIOR)
    {
        if((*user_param).periodic)
        {
            K3 = mesh_edge[edge3].neighbour[2];
            //get barycenter for K3 (physical)
            get_barycenter(K3,mesh_element,mesh_node,b3);

            //solution K3 @ b3
            map_to_reference_element(mesh_element,mesh_node,K3,ref_b3,b3[0],b3[1]);
            get_approx_sol_linear_basis(K3,mesh_element,mesh_node,uh_linear_unlimited,ref_b3[0],ref_b3[1],uh_K3_bar);

            //modify K1 so it is ghost neighbour of K0 to compute alphas
            if(mesh_edge[edge3].slope == VERTICAL)
            {
                //ghost element is left => translate center right
                if(mesh_edge[edge3].boundary_side== RIGHT)
                {
                    b3[0] = b3[0]+ (*user_param).periodic_domain_size_x;
                }
                else if(mesh_edge[edge3].boundary_side == LEFT)
                {
                    b3[0] = b3[0]- (*user_param).periodic_domain_size_x;
                }
            }
            else if(mesh_edge[edge3].slope == HORIZONTAL)
            {
                if(mesh_edge[edge3].boundary_side == BOTTOM)
                {
                    b3[1] = b3[1] - (*user_param).periodic_domain_size_y;
                }
                else if(mesh_edge[edge3].boundary_side == TOP)
                {
                    b3[1] = b3[1] + (*user_param).periodic_domain_size_y;
                }
            }
        }

           // printf("iedge @ exterior K3 %d  b3 (%lf %lf ) m3 (%lf %lf) \n",K3,b3[0],b3[1],m3[0],m3[1]); getchar();
    }//exterior

    //printf(" K0 %d K1  %d K2 %d K3 %d  \n",K0,K1,K2,K3);
    //printf(" K0_bar %lf K1_bar %lf K2_bar %lf K3_bar %lf \n",uh_K0_bar[0],uh_K1_bar[0],uh_K2_bar[0],uh_K3_bar[0]);
    //printf("uh_m1_K0 %lf uh_m2_K0 %lf uh_m3_K0 %lf \n",uhm1_K0[0],uhm2_K0[0],uhm3_K0[0]);//getchar();
    ////getchar();
    //compute alpha1 and alpha2 edge1
    compute_edge_alphas(m1,b0,b1,b2,b3,1,alphas_edge1,alpha1_elements);
    if(alpha1_elements[1] ==2)
    {
        
        test_debug = ((m1[0]-b0[0]) - (alphas_edge1[0]*(b1[0]-b0[0]) +alphas_edge1[1]*(b2[0]-b0[0])));
        //printf("test_bebug x  %10.16e alpha1 %10.16e alpha2 %10.16e \n",test_debug,alphas_edge1[0],alphas_edge1[1]);
        assert(fabs(test_debug) < 1e-10);
        test_debug = ((m1[1]-b0[1]) - (alphas_edge1[0]*(b1[1]-b0[1]) +alphas_edge1[1]*(b2[1]-b0[1])));
        //printf("test_bebug y %10.16e ( alpha 1 %10.16e alpha2 %10.16e  \n",test_debug,alphas_edge1[0],alphas_edge1[1]);
        //
        assert(fabs(test_debug) < 1e-5);
    }
    if(alpha1_elements[1] ==3)
    {
        test_debug = ((m1[0]-b0[0]) - (alphas_edge1[0]*(b1[0]-b0[0]) +alphas_edge1[1]*(b3[0]-b0[0])));
        //printf("test_bebug x  %10.16e alpha1 %10.16e alpha2 %10.16e \n",test_debug,alphas_edge1[0],alphas_edge1[1]);
        //assert(((m1[0]-b0[0]) - (alphas_edge1[0]*(b1[0]-b0[0]) +alphas_edge1[1]*(b3[0]-b0[0]))) < EPSILON);

        assert(fabs(test_debug) < 1e-10);
        
      
        test_debug = ((m1[1]-b0[1]) - (alphas_edge1[0]*(b1[1]-b0[1]) +alphas_edge1[1]*(b3[1]-b0[1])));
        //printf("test_bebug y  %10.16e alpha1 %10.16e alpha2 %10.16e \n",test_debug,alphas_edge1[0],alphas_edge1[1]);   
        assert(fabs(test_debug) < 1e-10);
    }    //compute tilde_uh_m1_K0 and delta_u_bar_m1_K 



    for(s_dim=0;s_dim<SYSTEM_DIM;s_dim++)
    {
        tilde_uh_m1_K0[s_dim] = uhm1_K0[s_dim] - uh_K0_bar[s_dim];
       

        if(alpha1_elements[1] ==2)
        {
            delta_u_bar_m1_K0[s_dim] = alphas_edge1[0]*(uh_K1_bar[s_dim] - uh_K0_bar[s_dim]) + alphas_edge1[1]*(uh_K2_bar[s_dim] - uh_K0_bar[s_dim]);
        }
        else if(alpha1_elements[1] ==3)
        {
            delta_u_bar_m1_K0[s_dim] = alphas_edge1[0]*(uh_K1_bar[s_dim] - uh_K0_bar[s_dim]) + alphas_edge1[1]*(uh_K3_bar[s_dim] - uh_K0_bar[s_dim]);
        }
        else
        {
            fprintf(stderr,"Error in edge1 delta computation \n");
            exit(1);
        }
        
        //printf("s_dim %d tile_uh_m1_K0 %lf \n",s_dim,tilde_uh_m1_K0[s_dim]);
        //printf("delta_ubar %lf \n",delta_u_bar_m1_K0[s_dim]);
    }
    //printf("==========================================================>\n");
    //getchar();

    //compute delta_1
    for(s_dim=0;s_dim<SYSTEM_DIM;s_dim++)
    {
        delta_1[s_dim]= min_mod_bar(tilde_uh_m1_K0[s_dim],(*user_param).nu*delta_u_bar_m1_K0[s_dim],user_param);
        //printf("tilde uh %lf mu_delta %lf \n",tilde_uh_m1_K0[s_dim],(*user_param).nu*delta_u_bar_m1_K0[s_dim]);
        //printf("delta_1 %lf \n",delta_1[s_dim]);//getchar();
    } 
    compute_edge_alphas(m2,b0,b1,b2,b3,2,alphas_edge2,alpha2_elements);
    if(alpha2_elements[1] ==1)
    {
        assert(fabs((m2[0]-b0[0]) - (alphas_edge2[0]*(b2[0]-b0[0]) +alphas_edge2[1]*(b1[0]-b0[0]))) < EPSILON);
        assert(fabs((m2[1]-b0[1]) - (alphas_edge2[0]*(b2[1]-b0[1]) +alphas_edge2[1]*(b1[1]-b0[1]))) < EPSILON);
    }
    if(alpha2_elements[1] ==3)
    {
        assert(fabs((m2[0]-b0[0]) - (alphas_edge2[0]*(b2[0]-b0[0]) +alphas_edge2[1]*(b3[0]-b0[0]))) < EPSILON);
        assert(fabs((m2[1]-b0[1]) - (alphas_edge2[0]*(b2[1]-b0[1]) +alphas_edge2[1]*(b3[1]-b0[1]))) < EPSILON);
    }
    //edge 2
    for(s_dim=0;s_dim<SYSTEM_DIM;s_dim++)
    {
        tilde_uh_m2_K0[s_dim] = uhm2_K0[s_dim] - uh_K0_bar[s_dim];
        

        if(alpha2_elements[1] ==1)
        {
            delta_u_bar_m2_K0[s_dim] = alphas_edge2[0]*(uh_K2_bar[s_dim] - uh_K0_bar[s_dim]) + alphas_edge2[1]*(uh_K1_bar[s_dim] - uh_K0_bar[s_dim]);
            //printf("uh_K2_bar[s_dim] %lf \n",uh_K2_bar[s_dim]);
        }
        else if(alpha2_elements[1] ==3)
        {
            delta_u_bar_m2_K0[s_dim] = alphas_edge2[0]*(uh_K2_bar[s_dim] - uh_K0_bar[s_dim]) + alphas_edge2[1]*(uh_K3_bar[s_dim] - uh_K0_bar[s_dim]);
            //printf("uh_K3_bar[s_dim] %lf \n",uh_K3_bar[s_dim]);
        }
        else
        {
            fprintf(stderr,"Error in edge2 delta computation \n");
            exit(1);
        }
    }
    //compute delta_2
    for(s_dim=0;s_dim<SYSTEM_DIM;s_dim++)
    {
        delta_2[s_dim]= min_mod_bar(tilde_uh_m2_K0[s_dim],(*user_param).nu*delta_u_bar_m2_K0[s_dim],user_param);
        //printf("tilde uh %lf mu_delta %lf \n",tilde_uh_m2_K0[s_dim],(*user_param).nu*delta_u_bar_m2_K0[s_dim]);
        //printf("delta_2 %lf \n",delta_2[s_dim]);//getchar();
    }    
    
    //printf("alphas_edge2 %lf %lf alpha_elements %d %d \n",alphas_edge2[0],alphas_edge2[1],alpha2_elements[0],alpha2_elements[1]);getchar();

    //printf("============================>m3 before compute alphas %lf %lf \n",m3[0],m3[1]);
    compute_edge_alphas(m3,b0,b1,b2,b3,3,alphas_edge3,alpha3_elements);

    if(alpha3_elements[1] ==1)
    {
        test_debug = ((m3[0]-b0[0]) - (alphas_edge3[0]*(b3[0]-b0[0]) +alphas_edge3[1]*(b1[0]-b0[0])));
        //printf("test_bebug x  %10.16e alpha1 %10.16e alpha2 %10.16e \n",test_debug,alphas_edge3[0],alphas_edge3[1]);
        assert(fabs((m3[0]-b0[0]) - (alphas_edge3[0]*(b3[0]-b0[0]) +alphas_edge3[1]*(b1[0]-b0[0]))) < 1e-10);
       
        test_debug = ((m3[1]-b0[1]) - (alphas_edge3[0]*(b3[1]-b0[1]) +alphas_edge3[1]*(b1[1]-b0[1])));
        //printf("test_bebug y  %10.16e alpha1 %10.16e alpha2 %10.16e \n",test_debug,alphas_edge3[0],alphas_edge3[1]);
        assert(fabs((m3[1]-b0[1]) - (alphas_edge3[0]*(b3[1]-b0[1]) +alphas_edge3[1]*(b1[1]-b0[1]))) < 1e-10);
    }
    if(alpha3_elements[1] ==2)
    {
         test_debug = ((m3[0]-b0[0]) - (alphas_edge3[0]*(b3[0]-b0[0]) +alphas_edge3[1]*(b2[0]-b0[0])));
         //printf("test_bebug x %10.16e \n",test_debug);        
         assert(fabs(test_debug) < 1e-10);
         test_debug = ((m3[1]-b0[1]) - (alphas_edge3[0]*(b3[1]-b0[1]) +alphas_edge3[1]*(b2[1]-b0[1])));
         //printf("test_bebug y  %10.16e alpha1 %10.16e alpha2 %10.16e \n",test_debug,alphas_edge3[0],alphas_edge3[1]);       
         assert(fabs(test_debug) < 1e-10);

    } 
    for(s_dim=0;s_dim<SYSTEM_DIM;s_dim++)
    {
        tilde_uh_m3_K0[s_dim] = uhm3_K0[s_dim] - uh_K0_bar[s_dim];

        if(alpha3_elements[1] ==1)
        {
            delta_u_bar_m3_K0[s_dim] = alphas_edge3[0]*(uh_K3_bar[s_dim] - uh_K0_bar[s_dim]) + alphas_edge3[1]*(uh_K1_bar[s_dim] - uh_K0_bar[s_dim]);
            //printf("uh_K3_bar[s_dim]-uh_K0_bar[s_dim]  %10.16e \n",(uh_K3_bar[s_dim] - uh_K0_bar[s_dim]));getchar();
        }
        else if(alpha3_elements[1] ==2)
        {
            delta_u_bar_m3_K0[s_dim] = alphas_edge3[0]*(uh_K3_bar[s_dim] - uh_K0_bar[s_dim]) + alphas_edge3[1]*(uh_K2_bar[s_dim] - uh_K0_bar[s_dim]);
        }
        else
        {
            fprintf(stderr,"Error in edge2 delta computation \n");
            exit(1);
        }

    }
    //compute delta_3
    for(s_dim=0;s_dim<SYSTEM_DIM;s_dim++)
    {
        delta_3[s_dim]= min_mod_bar(tilde_uh_m3_K0[s_dim],(*user_param).nu*delta_u_bar_m3_K0[s_dim],user_param);

        /*if(tilde_uh_m2_K0[s_dim] > 1e-6)
        {
             printf("tilde uh %10.16e delta %10.16e \n",tilde_uh_m3_K0[s_dim],delta_u_bar_m3_K0[s_dim]);
             printf("delta_3 %10.16e \n",delta_3[s_dim]);getchar();
        }*/
    }
    //printf("alphas_edge3 %lf %lf alpha_elements %d %d \n",alphas_edge3[0],alphas_edge3[1],alpha3_elements[0],alpha3_elements[1]);getchar();


    
  
   double t1,t2,t3;

  if((*user_param).flag_trans)
  {

      vec_m1minusb0[0] = m1[0]-b0[0];
      vec_m1minusb0[1] = m1[1]-b0[1];
      normalize_vector(2,vec_m1minusb0);


      vec_m2minusb0[0] = m2[0]-b0[0];
      vec_m2minusb0[1] = m2[1]-b0[1]; 
      normalize_vector(2,vec_m2minusb0);


      vec_m3minusb0[0] = m3[0]-b0[0];
      vec_m3minusb0[1] = m3[1]-b0[1];
      normalize_vector(2,vec_m3minusb0);

      transf_to_char_var(uh_K0_bar,tilde_uh_m1_K0,vec_m1minusb0,user_param,transf_tilde_uh_m1_K0);
      transf_to_char_var(uh_K0_bar,tilde_uh_m2_K0,vec_m2minusb0,user_param,transf_tilde_uh_m2_K0);
      transf_to_char_var(uh_K0_bar,tilde_uh_m3_K0,vec_m3minusb0,user_param,transf_tilde_uh_m3_K0);

      
      
      if(tilde_uh_m1_K0[0] > 1e-6)
      {
          //printf("transf_tilde_uh_m1_K0 %10.6e %10.6e %10.6e \n",transf_tilde_uh_m1_K0[0],transf_tilde_uh_m1_K0[1],transf_tilde_uh_m1_K0[2]);
          //printf("====>tilde_uh_m1_K0 %10.6e %10.6e %10.6e \n",tilde_uh_m1_K0[0],tilde_uh_m1_K0[1],tilde_uh_m1_K0[2]);
      }

      
      transf_to_char_var(uh_K0_bar,delta_u_bar_m1_K0,vec_m1minusb0,user_param,transf_delta_u_bar_m1_K0);
      transf_to_char_var(uh_K0_bar,delta_u_bar_m2_K0,vec_m2minusb0,user_param,transf_delta_u_bar_m2_K0);
      transf_to_char_var(uh_K0_bar,delta_u_bar_m3_K0,vec_m3minusb0,user_param,transf_delta_u_bar_m3_K0);

      for(s_dim=0;s_dim <SYSTEM_DIM;s_dim++)
      {

        
          transf_delta_1[s_dim]= min_mod_bar(transf_tilde_uh_m1_K0[s_dim],(*user_param).nu*transf_delta_u_bar_m1_K0[s_dim],user_param); 
          transf_delta_2[s_dim]= min_mod_bar(transf_tilde_uh_m2_K0[s_dim],(*user_param).nu*transf_delta_u_bar_m2_K0[s_dim],user_param);
          transf_delta_3[s_dim]= min_mod_bar(transf_tilde_uh_m3_K0[s_dim],(*user_param).nu*transf_delta_u_bar_m3_K0[s_dim],user_param);
      }
      transf_from_char_var(uh_K0_bar,transf_delta_1,vec_m1minusb0,user_param,delta_1);
      transf_from_char_var(uh_K0_bar,transf_delta_2,vec_m2minusb0,user_param,delta_2);
      transf_from_char_var(uh_K0_bar,transf_delta_3,vec_m3minusb0,user_param,delta_3);

      for(s_dim=0;s_dim<SYSTEM_DIM;s_dim++)
      {
          if(fabs(transf_delta_1[s_dim]+transf_delta_2[s_dim]+transf_delta_3[s_dim]) < EPSILON)
          {
              if((fabs((transf_tilde_uh_m1_K0[s_dim]-transf_delta_1[s_dim])) > 1e-10)||
                 (fabs((transf_tilde_uh_m2_K0[s_dim]-transf_delta_2[s_dim])) > 1e-10)||
                 (fabs((transf_tilde_uh_m3_K0[s_dim]-transf_delta_3[s_dim])) > 1e-10))
              {
                 
                  project_linear_to_dg_quad(K0,mesh_element,mesh_node,user_param,uh_K0_bar[s_dim],delta_1[s_dim],delta_2[s_dim],delta_3[s_dim],limited_dg_coefs);
                  grad_x = 2.0*(uhm2_K0[s_dim]-uhm3_K0[s_dim]);
                  grad_y = 2.0*(uhm2_K0[s_dim]-uhm1_K0[s_dim]);

                  limit_grad_x = 2.0*(limited_dg_coefs[1]-limited_dg_coefs[2]);
                  limit_grad_y = 2.0*(limited_dg_coefs[1]-limited_dg_coefs[0]);

                  limit_grad = sqrt(pow(limit_grad_x,2.0) + pow(limit_grad_y,2.0));
                  grad = sqrt(pow(grad_x,2.0) + pow(grad_y,2.0));
                  diff = grad-limit_grad;
                  //printf("limit_grad %lf grad %lf diff %10.16e \n",limit_grad,grad,diff);

                  //if(fabs(diff) > EPSILON)
                     //assert(grad-limit_grad > EPSILON);

                  for(idofs=0;idofs<Nloc;idofs++)
                  {
                    //printf("uh_limit before @idofs %d %lf \n",idofs,uh_limit[SYSTEM_DIM*(K0-1)*Nloc + (Nloc*s_dim)+idofs]);
                    uh_limit[SYSTEM_DIM*(K0-1)*Nloc + (Nloc*s_dim)+idofs] = limited_dg_coefs[idofs];
                    //printf(" after %lf \n",limited_dg_coefs[idofs]);
                  }
              }//if limiting makes diff
          }
          else
          {
                //printf("In pos-neg! sum_dim!=0!\n");//getchar();
                //printf("delta1 %10.16e  delta2 %10.16e delta3 %10.16e \n",transf_delta_1[s_dim],transf_delta_2[s_dim],transf_delta_3[s_dim]); //getchar();
                pos = max_fun(0.0,transf_delta_1[s_dim]) + max_fun(0.0,transf_delta_2[s_dim]) + max(0.0,transf_delta_3[s_dim]);
                neg = max_fun(0.0,-1.0*transf_delta_1[s_dim]) + max_fun(0.0,-1.0*transf_delta_2[s_dim]) + max_fun(0.0,-1.0*transf_delta_3[s_dim]);
                //printf("pos %10.16e neg %10.16e \n",pos,neg); //getchar();
                theta_plus = min_fun(1.0, neg/pos);
                theta_minus= min_fun(1.0, pos/neg);
                //printf("theta_plus %10.16e theta_minus %10.16e \n",theta_plus,theta_minus);

                tdelta_1_hat[s_dim] = theta_plus*max_fun(0.0,transf_delta_1[s_dim]) - theta_minus*max_fun(0.0,-1.0*transf_delta_1[s_dim]);
                tdelta_2_hat[s_dim] = theta_plus*max_fun(0.0,transf_delta_2[s_dim]) - theta_minus*max_fun(0.0,-1.0*transf_delta_2[s_dim]);
                tdelta_3_hat[s_dim] = theta_plus*max_fun(0.0,transf_delta_3[s_dim]) - theta_minus*max_fun(0.0,-1.0*transf_delta_3[s_dim]);
          }
      }//s_dim



      transf_from_char_var(uh_K0_bar,tdelta_1_hat,vec_m1minusb0,user_param,delta_1);
      transf_from_char_var(uh_K0_bar,tdelta_2_hat,vec_m2minusb0,user_param,delta_2);
      transf_from_char_var(uh_K0_bar,tdelta_3_hat,vec_m3minusb0,user_param,delta_3);
      
      for(s_dim=0;s_dim< SYSTEM_DIM;s_dim++)
      {

        
          if(fabs(transf_delta_1[s_dim]+transf_delta_2[s_dim]+transf_delta_3[s_dim]) > EPSILON)
          {
              if((fabs((transf_tilde_uh_m1_K0[s_dim]-tdelta_1_hat[s_dim])) > 1e-10)||
                (fabs((transf_tilde_uh_m2_K0[s_dim]-tdelta_2_hat[s_dim])) > 1e-10)||
                (fabs((transf_tilde_uh_m3_K0[s_dim]-tdelta_3_hat[s_dim])) > 1e-10))
            
              {

                      project_linear_to_dg_quad(K0,mesh_element,mesh_node,user_param,uh_K0_bar[s_dim],delta_1[s_dim],delta_2[s_dim],delta_3[s_dim],limited_dg_coefs);

                      grad_x = 2.0*(uhm2_K0[s_dim]-uhm3_K0[s_dim]);
                      grad_y = 2.0*(uhm2_K0[s_dim]-uhm1_K0[s_dim]);

                      limit_grad_x = 2.0*(limited_dg_coefs[1]-limited_dg_coefs[2]);
                      limit_grad_y = 2.0*(limited_dg_coefs[1]-limited_dg_coefs[0]);

                      limit_grad = sqrt(pow(limit_grad_x,2.0) + pow(limit_grad_y,2.0));
                      grad = sqrt(pow(grad_x,2.0) + pow(grad_y,2.0));
                      diff = grad-limit_grad;
                      //printf("limit_grad %lf grad %lf diff %10.16e \n",limit_grad,grad,diff);

                      //if(fabs(diff) > EPSILON)
                         //assert(grad-limit_grad > EPSILON);

                      for(idofs=0;idofs<Nloc;idofs++)
                      {
                        //printf("uh_limit before @idofs %d %lf \n",idofs,uh_limit[SYSTEM_DIM*(K0-1)*Nloc + (Nloc*s_dim)+idofs]);
                        uh_limit[SYSTEM_DIM*(K0-1)*Nloc + (Nloc*s_dim)+idofs] = limited_dg_coefs[idofs];
                        //printf(" after %lf \n",limited_dg_coefs[idofs]);
                      }
              }
          }
      }//s_dim

  }
  else
  {
      for(s_dim=0;s_dim<SYSTEM_DIM;s_dim++)
        {
            

            if(fabs(delta_1[s_dim]+delta_2[s_dim]+delta_3[s_dim]) < 1e-10)
            {
                //printf("sum delta_i =0 !\n");

                t1 =uh_K0_bar[s_dim] + delta_1[s_dim];
                t2 =uh_K0_bar[s_dim] + delta_2[s_dim];
                t3 =uh_K0_bar[s_dim] + delta_3[s_dim];
                //printf("s_dim %d t1 %lf t2 %lf t3 %lf \n",s_dim,t1,t2,t3);
                //printf("s_dim %d uhm1 %lf uhm2 %lf uhm3 %lf \n",s_dim,uhm1_K0[s_dim],uhm2_K0[s_dim],uhm3_K0[s_dim]);//getchar();

             

                project_linear_to_dg_quad(K0,mesh_element,mesh_node,user_param,uh_K0_bar[s_dim],delta_1[s_dim],delta_2[s_dim],delta_3[s_dim],limited_dg_coefs);
                /*if((fabs(limited_dg_coefs[0]- uhm1_K0[s_dim]) > 1e-10) ||
                   (fabs(limited_dg_coefs[1]- uhm2_K0[s_dim]) > 1e-10) ||
                   (fabs(limited_dg_coefs[2]- uhm3_K0[s_dim]) > 1e-10))*/

                if((fabs((tilde_uh_m1_K0[s_dim]-delta_1[s_dim])) > 1e-10)||
                 (fabs((tilde_uh_m2_K0[s_dim]-delta_2[s_dim])) > 1e-10)||
                 (fabs((tilde_uh_m3_K0[s_dim]-delta_3[s_dim])) > 1e-10))
                {
                      //printf("K0 %d limited\n",K0);
                      grad_x = 2.0*(uhm2_K0[s_dim]-uhm3_K0[s_dim]);
                      grad_y = 2.0*(uhm2_K0[s_dim]-uhm1_K0[s_dim]);

                      limit_grad_x = 2.0*(limited_dg_coefs[1]-limited_dg_coefs[2]);
                      limit_grad_y = 2.0*(limited_dg_coefs[1]-limited_dg_coefs[0]);

                      limit_grad = sqrt(pow(limit_grad_x,2.0) + pow(limit_grad_y,2.0));
                      grad = sqrt(pow(grad_x,2.0) + pow(grad_y,2.0));
                      diff = grad-limit_grad;
                     //printf("limit_grad %lf grad %lf diff %10.16e \n",limit_grad,grad,diff);

                    // if(fabs(diff) > EPSILON)
                     //    assert(grad-limit_grad > EPSILON);



                    mesh_element[K0].limited =1;
                    project_linear_to_dg_quad(K0,mesh_element,mesh_node,user_param,uh_K0_bar[s_dim],delta_1[s_dim],delta_2[s_dim],delta_3[s_dim],limited_dg_coefs);

                    //project_linear_to_dg_mon(K0,mesh_element,mesh_node,user_param,uh_K0_bar[s_dim],delta_1[s_dim],delta_2[s_dim],delta_3[s_dim],limited_dg_coefs1);
                    /*for(idofs=0;idofs<Nloc;idofs++)
                    {
                        printf("quad %lf linear %lf \n",limited_dg_coefs[idofs],limited_dg_coefs1[idofs]);
                        assert(fabs(limited_dg_coefs[idofs]-limited_dg_coefs1[idofs]) < EPSILON);
                    }*/

                    for(idofs=0;idofs<Nloc;idofs++)
                    {

                        //printf("uh_limit before @idofs %d %lf \n",idofs,uh_limit[SYSTEM_DIM*(K0-1)*Nloc + (Nloc*s_dim)+idofs]);
                        //if((s_dim==1))
                            uh_limit[SYSTEM_DIM*(K0-1)*Nloc + (Nloc*s_dim)+idofs] = limited_dg_coefs[idofs];
                        //printf(" after %lf \n",limited_dg_coefs[idofs]);

                    }
                    //getchar();
                }
                else
                {
                    //printf("element %d K0 NOT limited s_dim %d  \n",K0,s_dim);//getchar();
                }
     
            }
            else
            {

                //printf("In pos-neg! sum_dim!=0!\n");//getchar();
                //printf("delta1 %10.16e  delta2 %10.16e delta3 %10.16e \n",delta_1[s_dim],delta_2[s_dim],delta_3[s_dim]); //getchar();
                pos = max_fun(0.0,delta_1[s_dim]) + max_fun(0.0,delta_2[s_dim]) + max(0.0,delta_3[s_dim]);
                neg = max_fun(0.0,-1.0*delta_1[s_dim]) + max_fun(0.0,-1.0*delta_2[s_dim]) + max_fun(0.0,-1.0*delta_3[s_dim]);
                //printf("pos %10.16e neg %10.16e \n",pos,neg); //getchar();
                theta_plus = min_fun(1.0, neg/pos);
                theta_minus= min_fun(1.0, pos/neg);
                //printf("theta_plus %10.16e theta_minus %10.16e \n",theta_plus,theta_minus);

                delta_1_hat = theta_plus*max_fun(0.0,delta_1[s_dim]) - theta_minus*max_fun(0.0,-1.0*delta_1[s_dim]);
                delta_2_hat = theta_plus*max_fun(0.0,delta_2[s_dim]) - theta_minus*max_fun(0.0,-1.0*delta_2[s_dim]);
                delta_3_hat = theta_plus*max_fun(0.0,delta_3[s_dim]) - theta_minus*max_fun(0.0,-1.0*delta_3[s_dim]);
                
                
                //printf("delta_1_hat %10.16e delta_2_hat %10.16e delta_3_hat %10.16e \n",delta_1_hat,delta_2_hat,delta_3_hat);///getchar();
                //printf("uh_K0_bar %lf s_dim %d \n",uh_K0_bar[s_dim],s_dim);//getchar();


                project_linear_to_dg_quad(K0,mesh_element,mesh_node,user_param,uh_K0_bar[s_dim],delta_1_hat,delta_2_hat,delta_3_hat,limited_dg_coefs);
               
                /*if((fabs(limited_dg_coefs[0]- uhm1_K0[s_dim]) > 1e-9) ||
                   (fabs(limited_dg_coefs[1]- uhm2_K0[s_dim]) > 1e-9) ||
                   (fabs(limited_dg_coefs[2]- uhm3_K0[s_dim]) > 1e-9))*/
                if((fabs((tilde_uh_m1_K0[s_dim]-delta_1_hat)) > 1e-10)||
                 (fabs((tilde_uh_m2_K0[s_dim]-delta_2_hat)) > 1e-10)||
                 (fabs((tilde_uh_m3_K0[s_dim]-delta_3_hat)) > 1e-10))
                {

                      //printf("K0 %d limited\n",K0);

                      grad_x = 2.0*(uhm2_K0[s_dim]-uhm3_K0[s_dim]);
                      grad_y = 2.0*(uhm2_K0[s_dim]-uhm1_K0[s_dim]);

                      limit_grad_x = 2.0*(limited_dg_coefs[1]-limited_dg_coefs[2]);
                      limit_grad_y = 2.0*(limited_dg_coefs[1]-limited_dg_coefs[0]);

                      limit_grad = sqrt(pow(limit_grad_x,2.0) + pow(limit_grad_y,2.0));
                      grad = sqrt(pow(grad_x,2.0) + pow(grad_y,2.0));



                      diff = grad-limit_grad;
                      //printf("limit_grad %lf grad %lf diff %10.16e \n",limit_grad,grad,diff);

                      //if(fabs(diff)> EPSILON)
                       //   assert(grad-limit_grad > EPSILON);




                    mesh_element[K0].limited=1;
                    project_linear_to_dg_quad(K0,mesh_element,mesh_node,user_param,uh_K0_bar[s_dim],delta_1_hat,delta_2_hat,delta_3_hat,limited_dg_coefs);

                    //project_linear_to_dg_mon(K0,mesh_element,mesh_node,user_param,uh_K0_bar[s_dim],delta_1_hat,delta_2_hat,delta_3_hat,limited_dg_coefs1);

                    /*for(idofs=0;idofs<Nloc;idofs++)
                    {

                        printf("quad %lf linear %lf \n",limited_dg_coefs[idofs],limited_dg_coefs1[idofs]);
                        assert(fabs(limited_dg_coefs[idofs]-limited_dg_coefs1[idofs]) < EPSILON);
                    }*/
                    



                    for(idofs=0;idofs<Nloc;idofs++)
                    {
                        //printf("uh_limit before %lf \n",uh_limit[SYSTEM_DIM*(K0-1)*Nloc + (Nloc*s_dim)+idofs]);
                        
                        //if(s_dim==1)
                            uh_limit[SYSTEM_DIM*(K0-1)*Nloc + (Nloc*s_dim)+idofs] = limited_dg_coefs[idofs];

                        //printf("%lf \n",limited_dg_coefs[idofs]);

                        //printf("uh_limit after %lf \n",uh_limit[SYSTEM_DIM*(K0-1)*Nloc + (Nloc*s_dim)+idofs]);
                    }




                   // getchar();
                }
                else
                {
                    //printf("element %d K0 NOT limited s_dim %d  \n",K0,s_dim);//getchar();
                }


            }
            
        }//s_dim
  }

    double dg_sol_av[SYSTEM_DIM];
    dg_quadratic_basis_cell_average(K0,mesh_element,mesh_node,uh_limit,dg_sol_av);
    //assert(fabs(dg_sol_av[2]-dg_sol_av1[2]) < EPSILON);
    //printf("------------------done with K0 (%d)  cell average %lf %lf %lf  \n",K0,dg_sol_av[0],dg_sol_av[1],dg_sol_av[2]);///getchar();
   free(limited_dg_coefs1);
   free(limited_dg_coefs);


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
        quad_sol_linear_phi(E,mesh_element,mesh_node,user_param,uh_tn_sol,s_dim,quad_sol_linear_phi_vec);
        //solve linear system
        qrsolv(A_loc_QR,linear_NLOC,Q_coeff,R_coeff,quad_sol_linear_phi_vec);
        for(idofs=0;idofs<linear_NLOC;idofs++)
        {
            uh_linear[SYSTEM_DIM*(E-1)*linear_NLOC + (s_dim*linear_NLOC) +idofs] = quad_sol_linear_phi_vec[idofs];
            //assert(fabs(quad_sol_linear_phi_vec[idofs]-uh_tn_sol[SYSTEM_DIM*(E-1)*linear_NLOC + (s_dim*linear_NLOC) +idofs]) < EPSILON);
        }

    }

    free(Q_coeff);
    free(R_coeff);
    free(linear_Aloc_matrix);
    for(i=0;i<linear_NLOC+1;i++)
        free(A_loc_QR[i]);

    

}

void limiting_linear_basis_functions(double x,
                                    double y,
                                    double* phi_linear_vec)
{
   
    phi_linear_vec[0] = 1.0 - 2.0*y;
    phi_linear_vec[1] = 2.0*x+2.0*y-1.0;
    phi_linear_vec[2] = 1.0-2.0*x;
    
}

void quad_sol_linear_phi(int E,
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

void linear_sol_phi(int E,
                         element* mesh_element,
                         node* mesh_node,
                         parameters* user_param,
                         double uh_K0_bar,
                         double delta1,
                         double delta2,
                         double delta3,
                         double* linear_sol_phi_vec)
{
    int idofs=0;
    int k=0;
    int Nloc=0;
    double linear_sol_val=0.0;

    Nloc = (mesh_element[E].degree+1)*(mesh_element[E].degree+2)/2;
    init_zero_d(linear_sol_phi_vec,Nloc);

    for(k=1;k<ngpts+1;k++) 
    {
        //get linear solution at gauss points
        linear_sol_val =limiting_linear_function_val(uh_K0_bar,delta1,delta2,delta3,mesh_element[E].el_gpts_ref_x[k],mesh_element[E].el_gpts_ref_y[k]);
        //printf("uk0_bar %lf \n",uh_K0_bar); //getchar();
        //printf("linear_sol_val %lf \n",linear_sol_val);
        for(idofs =0; idofs < Nloc;idofs++) 
        {
            linear_sol_phi_vec[idofs]+=mesh_element[E].det*(linear_sol_val*mesh_element[E].el_gpts_basis[(k-1)*Nloc+idofs])*mesh_element[E].el_gpts_w[k];
        }
    }//loop over quad points

   
}

double limiting_linear_function_val(double uh_K0_bar,
                                    double delta1,
                                    double delta2,
                                    double delta3,
                                    double xhat,
                                    double yhat)
{
    double limited_sol_val =0.0;

    limited_sol_val = uh_K0_bar + delta1*(1.0-2*yhat) + delta2*(2.0*xhat + 2.0*yhat -1.0) + delta3*(1.0-2.0*xhat);

    return limited_sol_val;
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


void limiting_linear_cell_average(int E,
                                  element* mesh_element,
                                  node* mesh_node,
                                  double* uh_tn_linear,
                                  double* uh_linear_average)
{
    int idofs=0;
    int k=0;
    int s_dim=0;
    double linear_approx_sol[SYSTEM_DIM];
    init_zero_d(uh_linear_average,SYSTEM_DIM);
    for(k=1;k<ngpts+1;k++)
    {
        init_zero_d(linear_approx_sol,SYSTEM_DIM);
        get_approx_sol_linear_basis(E,mesh_element,mesh_node,uh_tn_linear,mesh_element[E].el_gpts_ref_x[k],mesh_element[E].el_gpts_ref_y[k],linear_approx_sol);
        
        for(s_dim=0;s_dim<SYSTEM_DIM;s_dim++)
        {
            uh_linear_average[s_dim] += (1.0/(0.5*mesh_element[E].det))*(mesh_element[E].det*linear_approx_sol[s_dim]*mesh_element[E].el_gpts_w[k]);
        }
    }
}

void get_approx_sol_linear_basis(int E,
                                 element* mesh_element,
                                 node* mesh_node,
                                 double* uh_tn_linear,
                                 double xhat,
                                 double yhat,
                                 double* uh_linear_sol)
{

    int idofs =0;
    int k=0; 
    int s_dim=0;
    double temp =0.0;
    double linear_phi_vec[linear_NLOC];


    init_zero_d(linear_phi_vec,linear_NLOC);
    limiting_linear_basis_functions(xhat,yhat,linear_phi_vec);
    for(s_dim=0;s_dim <SYSTEM_DIM;s_dim++)
    {
        temp=0.0;
        for(idofs =0;idofs<linear_NLOC;idofs++)
        {
            temp+= uh_tn_linear[SYSTEM_DIM*linear_NLOC*(E-1) + (s_dim*linear_NLOC) +idofs]*linear_phi_vec[idofs];
        }
        uh_linear_sol[s_dim] = temp;
    }
}

void dg_quadratic_basis_cell_average(int E,
                                     element* mesh_element,
                                     node* mesh_node,
                                     double* uh_tn_sol,
                                     double* uh_quad_average)
{
    int idofs=0;
    int k=0;
    int s_dim=0;
    double approx_sol[SYSTEM_DIM];

    init_zero_d(uh_quad_average,SYSTEM_DIM);
    
    for(k=1;k<ngpts+1;k++)
    {
        init_zero_d(approx_sol,SYSTEM_DIM);
        get_approx_solution(uh_tn_sol,E,mesh_element,mesh_node,mesh_element[E].el_gpts_ref_x[k],mesh_element[E].el_gpts_ref_y[k],approx_sol);
        
        for(s_dim=0;s_dim<SYSTEM_DIM;s_dim++)
        {
            uh_quad_average[s_dim] += (1.0/(0.5*mesh_element[E].det))*(mesh_element[E].det*approx_sol[s_dim]*mesh_element[E].el_gpts_w[k]);
        }
    }
}

void get_barycenter(int E,
                    element* mesh_element,
		            node* mesh_node,
		            double* bary_coords
                    )
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
}


/*returns a positive combination of alpha
 * and corresponding elements*/
void compute_edge_alphas(double* mi,
                         double* b0,
                         double* b1,
                         double* b2,
                         double* b3,
                         int iedge_index,
                         double* alpha_edgei,
                         int* alpha_neighbours)
{
    double alpha1=0.0;
    double alpha2=0.0;

    /*printf("iedge_index %d mi (%lf %lf ) \n",iedge_index,mi[0],mi[1]);
    printf("b0 (%lf %lf ) \n",b0[0],b0[1]);
    printf("b1 (%lf %lf)  b2 (%lf %lf )  b3 (%lf %lf ) \n",b1[0],b1[1],b2[0],b2[1],b3[0],b3[1]);*/

    switch(iedge_index)
    {
        case(1):
            {
                //try K1 and K2
                alpha1 = ((mi[1]-b0[1])*(b2[0]-b0[0]) - (mi[0] - b0[0])*(b2[1]-b0[1]))/((b1[1]-b0[1])*(b2[0]-b0[0]) - (b1[0]-b0[0])*(b2[1]-b0[1]));

                if((!(isnan(alpha1))) & (!isinf(alpha1)) & (fabs(alpha1) < EPSILON))
                {
                    //printf("alpha1 zero correct %10.16e \n",alpha1);
                    alpha1 = fabs(alpha1);
                }
                //printf("alpha1 with K1 and K2 %lf \n",alpha1);
                if((alpha1 <0) || (isnan(alpha1)) || isinf(alpha1))
                {
                    //try K1 and K3
                    alpha1 = ((mi[1]-b0[1])*(b3[0]-b0[0]) - (mi[0] - b0[0])*(b3[1]-b0[1]))/((b1[1]-b0[1])*(b3[0]-b0[0]) - (b1[0]-b0[0])*(b3[1]-b0[1]));
                    if((alpha1 <0) || (isnan(alpha1))|| (isinf(alpha1)))
                    {
                        fprintf(stderr,"Non-negative alpha1 not found! \n");
                        exit(1);
                    }
                    else
                    {
                        alpha2 =  ((mi[0]-b0[0]) - alpha1*(b1[0]-b0[0]))/(b3[0]-b0[0]);
                        if((!(isnan(alpha2))) & (!(isinf(alpha2))) & (fabs(alpha2) < EPSILON))
                        {
                            //printf("alpha2 zero correct %10.16e \n",alpha2);
                            alpha2 =fabs(alpha2);
                        }
                        if((alpha2 <0) || (isnan(alpha2)) || (isinf(alpha2)))
                        {
                            fprintf(stderr,"Non-negative alpha2 not found edge1 ! \n");
                            exit(1);
                        }
                        else
                        {
                            //K1 and K3 give non-negative alpha
                            alpha_neighbours[0] = 1; //edge 1
                            alpha_neighbours[1] = 3; //edge 3
                        }
                    }
                }
                else
                {                    
                    alpha2 =  ((mi[0]-b0[0]) - alpha1*(b1[0]-b0[0]))/(b2[0]-b0[0]);
                    if((!(isnan(alpha2))) & (!(isinf(alpha2))) & (fabs(alpha2) < EPSILON))
                    {
                     
                        //printf("alpha2 zero correct %10.16e \n",alpha2);
                        alpha2 =fabs(alpha2);
                    }
                    //printf("alpha 2  with K1 and K2 %10.16e \n",alpha2);
                    if((alpha2 <0) || (isnan(alpha2)) || (isinf(alpha2)))
                    {
                        //try K1 and K3
                        alpha1 = ((mi[1]-b0[1])*(b3[0]-b0[0]) - (mi[0] - b0[0])*(b3[1]-b0[1]))/((b1[1]-b0[1])*(b3[0]-b0[0]) - (b1[0]-b0[0])*(b3[1]-b0[1]));
                        //printf("m1-b0 %lf \n",mi[0]-b0[0]);
                        //printf("b1-b0 %lf \n",b1[0]-b0[0]);
                        //printf("alpha1 with K1 and K3 %10.16e \n",alpha1);
                        if((!(isnan(alpha1))) & (!isinf(alpha1)) & (fabs(alpha1) < EPSILON))
                        {
                            //printf("alpha1 zero correct %10.16e \n",alpha1);
                            alpha1=fabs(alpha1);
                        }                        
                        if((alpha1 <0) || (isnan(alpha1))|| (isinf(alpha1)))
                        {
                            fprintf(stderr,"Non-negative alpha1 not found with K1 and K3 edge 1! \n");
                            exit(1);
                        }
                        else
                        {
                            alpha2 =  ((mi[0]-b0[0]) - alpha1*(b1[0]-b0[0]))/(b3[0]-b0[0]);
                            //printf("b3[0] %lf b0[0] %lf \n",b3[0],b0[0]); 
                            if((!(isnan(alpha2))) & (!(isinf(alpha2))) & (fabs(alpha2) < EPSILON))
                            {
                                alpha2 = fabs(alpha2);
                            }
                            //printf("alpha 2 with K1 and K3  %10.16e \n",alpha2); //getchar();
                            if((alpha2 <0) || (isnan(alpha2)) || isinf(alpha2))
                            {
                                fprintf(stderr,"Non-negative alpha2 not found edge 1! \n");
                                exit(1);
                            }
                            else
                            {
                                //K1 and K3 give non-negative alpha
                                alpha_neighbours[0] = 1; //edge 1
                                alpha_neighbours[1] = 3; //edge 3
                                //printf("alpha_neighbours[0] %d alpha_neighbours[1] %d\n",alpha_neighbours[0],alpha_neighbours[1]);
                            }
                        }
                    }
                    else
                    {

                        //K1 and K2 give non-negative alphas
                        alpha_neighbours[0] = 1; //edge 1
                        alpha_neighbours[1] = 2; //edge 2
                    }
                }

                break;
            }//case 1: edge 1
        case(2):
            {
                //try K2 and K1
                alpha1 = ((mi[1]-b0[1])*(b1[0] -b0[0])-(mi[0]-b0[0])*(b1[1]-b0[1]))/((b2[1]-b0[1])*(b1[0]-b0[0])-(b2[0]-b0[0])*(b1[1]-b0[1]));
                if((!(isnan(alpha1))) & (!isinf(alpha1)) & (fabs(alpha1) < EPSILON)) 
                {
                    alpha1 =fabs(alpha1);
                }
                //printf("alpha 1 edge 2 with K2 and K1  %10.16e \n",alpha1); 
                if((alpha1 <0) || (isnan(alpha1)) || isinf(alpha1)) //K2 and K1 don't work
                {
                    //try K2 and K3
                    alpha1 = ((mi[1]-b0[1])*(b3[0] -b0[0])-(mi[0]-b0[0])*(b3[1]-b0[1]))/((b2[1]-b0[1])*(b3[0]-b0[0])-(b2[0]-b0[0])*(b3[1]-b0[1]));
                    if((alpha1 <0) || (isnan(alpha1)) || (isinf(alpha1)))
                    {
                        fprintf(stderr,"Non-negative alpha1 not found! \n");
                        exit(1);
                    }
                    else //alpha1 found with K2 and K3 
                    {
                        alpha2 = ((mi[0]-b0[0]) - alpha1*(b2[0]-b0[0]))/(b3[0]-b0[0]);
                         if((alpha2 <0) || (isnan(alpha2)) || (isinf(alpha2)))
                         {
                             fprintf(stderr,"Non-negative alpha2 not found! \n");
                             exit(1);
                         }
                         else
                         {
                             alpha_neighbours[0] = 2;
                             alpha_neighbours[1] = 3;
                         }
                    }
                }
                else
                {
                    //using K2 and K1 (alpha 2)
                    alpha2 = ((mi[0]-b0[0]) - alpha1*(b2[0]-b0[0]))/(b1[0]-b0[0]); 
                    if((!(isnan(alpha2))) & (!isinf(alpha2)) & (fabs(alpha2) < EPSILON)) 
                    {
                        alpha2 = fabs(alpha2);
                    }
                    //printf("alpha2  in edge 2 %10.16e \n",alpha2); 
                    if((alpha2 <0) || (isnan(alpha2)) || (isinf(alpha2)))
                    {
                        //try K2 and K3
                        alpha1 = ((mi[1]-b0[1])*(b3[0] -b0[0])-(mi[0]-b0[0])*(b3[1]-b0[1]))/((b2[1]-b0[1])*(b3[0]-b0[0])-(b2[0]-b0[0])*(b3[1]-b0[1]));
                        if((!(isnan(alpha1))) & (!isinf(alpha1)) & (fabs(alpha1) < EPSILON)) 
                        {
                            alpha1 = fabs(alpha1);
                        }
                        //printf("alpha1 in edge 2 using K2 and K3 %lf \n",alpha1);
                        if((alpha1 <0) || (isnan(alpha1)) || (isinf(alpha1)))
                        {
                            fprintf(stderr,"Non-negative alpha1 not found! \n");
                            exit(1);
                        }
                        else
                        {
                            alpha2 = ((mi[0]-b0[0]) - alpha1*(b2[0]-b0[0]))/(b3[0]-b0[0]);

                            if((!(isnan(alpha2))) & (!isinf(alpha2)) & (fabs(alpha2) < EPSILON)) 
                            {
                                alpha2 = fabs(alpha2);
                            }
                            //printf("m2 ( %lf  %lf ) b0 (%lf %lf ) b2 (%lf %lf ) b3 (%lf %lf ) \n",mi[0],mi[1],b0[0],b0[1],b2[0],b2[1],b3[0],b3[1]); getchar();
                            //printf("alpha2  in edge 2 using K2 and K3 %10.16e \n",alpha2); //getchar();
                            if((alpha2 <0) || (isnan(alpha2)) || (isinf(alpha2)))
                            {
                                 fprintf(stderr,"Non-negative alpha2 not found! \n");
                                 exit(1);
                            }
                            else
                            {

                                alpha_neighbours[0] = 2;
                                alpha_neighbours[1] = 3;
                            }
                        }
                    }
                    else
                    {
                        alpha_neighbours[0]  = 2;
                        alpha_neighbours[1]  = 1;
                    }//alphas with K2 and K1
                }
                break;
            }
        case(3):
            {
                //try K3 and K1
                //printf("iedge_index %d m3 (%lf %lf ) b0 (%lf %lf ) b3 (%lf %lf ) b1 (%lf %lf ) \n",iedge_index,mi[0],mi[1],b0[0],b0[1],b3[0],b3[1],b1[0],b1[1]);getchar();
                alpha1 = ((mi[1]-b0[1])*(b1[0]-b0[0])-(mi[0]-b0[0])*(b1[1]-b0[1]))/((b3[1]-b0[1])*(b1[0]-b0[0])-(b3[0]-b0[0])*(b1[1]-b0[1]));
                if((!(isnan(alpha1))) & (!isinf(alpha1)) & (fabs(alpha1) < EPSILON)) 
                {
                    alpha1 = fabs(alpha1);
                }
                //printf("alpha1 in edge 3  K3 and K1 %10.16e \n",alpha1);

                if((alpha1 <0) || (isnan(alpha1)) || (isinf(alpha1)))
                {
                    //try K3 and K2
                    alpha1 = ((mi[1]-b0[1])*(b2[0]-b0[0])-(mi[0]-b0[0])*(b2[1]-b0[1]))/((b3[1]-b0[1])*(b2[0]-b0[0])- (b3[0]-b0[0])*(b2[1]-b0[1]));
                    if((!(isnan(alpha1))) & (!isinf(alpha1)) & (fabs(alpha1) < EPSILON)) 
                    {
                        alpha1 = fabs(alpha1);
                    }
                    //printf("alpha 1 using K3 and K2 %10.16e \n",alpha1); getchar();
                    if((alpha1 <0) || (isnan(alpha1)) || (isinf(alpha1)))
                    {
                        fprintf(stderr,"Non-negative alpha1 not found ! \n");
                        exit(1);
                    }
                    else
                    {
                        //using K3 and K2 with K0
                        alpha2 = ((mi[0]-b0[0]) - alpha1*(b3[0]-b0[0]))/(b2[0]-b0[0]);
                        if((!(isnan(alpha2))) & (!isinf(alpha2)) & (fabs(alpha2) < EPSILON)) 
                        {
                            alpha2 = fabs(alpha2);
                        }
                        //printf("alpha 2 in edge 3 %10.16e \n",alpha2);
                        if((alpha2<0)||(isnan(alpha2)) ||(isinf(alpha2)))
                        {
                            fprintf(stderr,"Non-negative alpha2 not found !\n");
                            exit(1);
                        }
                        else
                        {
                            alpha_neighbours[0] = 3;
                            alpha_neighbours[1] = 2;
                        }
                    }
                }
                else
                {
                    alpha2 = ((mi[0]-b0[0]) - alpha1*(b3[0]-b0[0]))/(b1[0]-b0[0]);
                    if((!(isnan(alpha2))) & (!isinf(alpha2)) & (fabs(alpha2) < EPSILON)) 
                    {
                        alpha2 = fabs(alpha2);
                    }
                    //printf("alpha2 using K3 and K1 %10.16e \n",alpha2); getchar();
                    if((alpha2 < 0) || (isnan(alpha2))|| (isinf(alpha2)))
                    {
                        //try K3 and K2
                        alpha1 = ((mi[1]-b0[1])*(b2[0]-b0[0])-(mi[0]-b0[0])*(b2[1]-b0[1]))/((b3[1]-b0[1])*(b2[0]-b0[0])- (b3[0]-b0[0])*(b2[1]-b0[1]));
                        if((!(isnan(alpha1))) & (!isinf(alpha1)) & (fabs(alpha1) < EPSILON)) 
                        {
                            alpha1 = fabs(alpha1);
                        }
                        if((alpha1 <0) || (isnan(alpha1)) || (isinf(alpha1)))
                        {
                            fprintf(stderr,"Non-negative alpha1 not found ! \n");
                            exit(1);
                        }
                        else
                        {
                            //using K3 and K2 with K0
                            alpha2 = ((mi[0]-b0[0]) - alpha1*(b3[0]-b0[0]))/(b2[0]-b0[0]);
                            if((!(isnan(alpha2))) & (!isinf(alpha2)) & (fabs(alpha2) < EPSILON))
                            {
                                alpha2 = fabs(alpha2);
                            }
                            if((alpha2<0)||(isnan(alpha2))|| (isinf(alpha2)))
                            {
                                fprintf(stderr,"Non-negative alpha2 not found !\n");
                                exit(1);
                            }
                            else
                            {
                                alpha_neighbours[0] = 3;
                                alpha_neighbours[1] = 2;
                            }
                        }
                    }
                    else
                    {
                        alpha_neighbours[0] = 3;
                        alpha_neighbours[1] = 1;
                    }
                }
                break;
            }//case 3
    }//switch edge index

    alpha_edgei[0] = alpha1;
    alpha_edgei[1] = alpha2;

}

void project_linear_to_dg_quad(int E,
                            element* mesh_element,
                            node* mesh_node,
                            parameters* user_param,
                            double uh_K0_bar,
                            double delta1,
                            double delta2,
                            double delta3,
                            double* dg_loc_coefs 
                            )
{
    double* Aloc_matrix = NULL;
    double* linear_sol_phi_vec =NULL;
    // QR-decomposition and solving
    double** A_loc_QR =NULL; 
    double* Q_coeff = NULL; 
    double* R_coeff=NULL; 
    int sing_decomp;
    int i,j;
    int idofs=0;
    int Nloc =0;

    Nloc = (mesh_element[E].degree+1)*(mesh_element[E].degree+2)/2;
   

    Q_coeff = (double*)malloc((Nloc+1)*sizeof(double));
    R_coeff = (double*)malloc((Nloc+1)*sizeof(double));
    A_loc_QR = (double**)malloc((Nloc+1)*sizeof(double*));
    Aloc_matrix = (double*)malloc((Nloc*Nloc)*sizeof(double));
    linear_sol_phi_vec = (double*)malloc(Nloc*sizeof(double));

    for(i=0;i<Nloc+1;i++)
    A_loc_QR[i] = (double*)malloc((Nloc+1)*sizeof(double));


    init_zero_d(Aloc_matrix,Nloc*Nloc);
    init_zero_d(linear_sol_phi_vec,Nloc);
    init_zero_m(A_loc_QR,Nloc+1,Nloc+1);

    Alocal_mat(E,mesh_element,mesh_node,user_param,Aloc_matrix);

    assemble_qr(Aloc_matrix,Nloc,Q_coeff,R_coeff,A_loc_QR,&sing_decomp);
    init_zero_d(linear_sol_phi_vec,Nloc);
    linear_sol_phi(E,mesh_element,mesh_node,user_param,uh_K0_bar,delta1,delta2,delta3,linear_sol_phi_vec);
    
    //solve linear system
    qrsolv(A_loc_QR,Nloc,Q_coeff,R_coeff,linear_sol_phi_vec);
    for(idofs=0;idofs<Nloc;idofs++)
    {
        dg_loc_coefs[idofs] = linear_sol_phi_vec[idofs];
    }
    free(Q_coeff);
    free(R_coeff);
    free(Aloc_matrix);
    free(linear_sol_phi_vec); 
    for(i=0;i<Nloc+1;i++)
        free(A_loc_QR[i]);

    free(A_loc_QR);



}
void  get_edge_mid_point(edge* mesh_edge,
                         int iedge,
                         node* mesh_node,
                         double* mid_point)
{
    int node_a,node_b;
    double x1,y1,x2,y2;
    node_a = mesh_edge[iedge].vertex[1];
    node_b = mesh_edge[iedge].vertex[2];

    x1 = mesh_node[node_a].coord[0];
    y1 = mesh_node[node_a].coord[1];

    x2 = mesh_node[node_b].coord[0];
    y2 = mesh_node[node_b].coord[1];


    mid_point[0] = 0.5*(x1+x2);
    mid_point[1] = 0.5*(y1+y2);
}


double min_mod_bar(double q1,
                   double q2,
                   parameters* user_param)
{

    double return_val =0.0;

    if(fabs(q1) <= (*user_param).mdx2)
    {
        return_val = q1;
    }
    else
    {
        return_val = min_mod(q1,q2);
    }

    return return_val;
}
double min_mod(double q1,
               double q2)
{
    double return_val=0.0;

    if((q1 <=0) && (q2 <=0))
    {
        return_val = -1.0*min_fun(fabs(q1),fabs(q2));
    }
    else if((q1 >=0) && (q2 >=0))
    {
        return_val = 1.0*min_fun(fabs(q1),fabs(q2));
    }
    else
    {
        return_val =0.0;
    }

    return return_val;
}

double min_fun(double a,
           double b)
{
    double return_val =0.0;

    if(a < b)
    {
        return_val =a;
    }
    else if(b < a)
    {
        return_val = b;
    }
    else
    {
        return_val = a;
    }
    return return_val;
}
double max_fun(double a,
              double b)
{
    double return_val=0.0;
    if(a > b)
    {
        return_val = a;
    }
    else if(b >a)
    {
        return_val =b;
    }
    else
    {
        return_val = a;
    }
    return return_val;
}

void project_linear_to_dg_mon(int E,
                            element* mesh_element,
                            node* mesh_node,
                            parameters* user_param,
                            double uh_K0_bar,
                            double delta1,
                            double delta2,
                            double delta3,
                            double* dg_loc_coefs 
                            )
{
    double* Aloc_matrix = NULL;
    double* linear_sol_phi_vec =NULL;
    // QR-decomposition and solving
    double** A_loc_QR =NULL; 
    double* Q_coeff = NULL; 
    double* R_coeff=NULL; 
    int sing_decomp;
    int i,j;
    int idofs=0;

   

    Q_coeff = (double*)malloc((linear_NLOC+1)*sizeof(double));
    R_coeff = (double*)malloc((linear_NLOC+1)*sizeof(double));
    A_loc_QR = (double**)malloc((linear_NLOC+1)*sizeof(double*));
    Aloc_matrix = (double*)malloc((linear_NLOC*linear_NLOC)*sizeof(double));
    linear_sol_phi_vec = (double*)malloc(linear_NLOC*sizeof(double));

    for(i=0;i<linear_NLOC+1;i++)
    A_loc_QR[i] = (double*)malloc((linear_NLOC+1)*sizeof(double));


    init_zero_d(Aloc_matrix,linear_NLOC*linear_NLOC);
    init_zero_d(linear_sol_phi_vec,linear_NLOC);
    init_zero_m(A_loc_QR,linear_NLOC+1,linear_NLOC+1);

    Aloc_mat_linear_mon(E,mesh_element,mesh_node,user_param,Aloc_matrix);

    assemble_qr(Aloc_matrix,linear_NLOC,Q_coeff,R_coeff,A_loc_QR,&sing_decomp);
    init_zero_d(linear_sol_phi_vec,linear_NLOC);
    linear_sol_phi(E,mesh_element,mesh_node,user_param,uh_K0_bar,delta1,delta2,delta3,linear_sol_phi_vec);
    
    //solve linear system
    qrsolv(A_loc_QR,linear_NLOC,Q_coeff,R_coeff,linear_sol_phi_vec);
    init_zero_d(dg_loc_coefs,6);
    for(idofs=0;idofs<linear_NLOC;idofs++)
    {
        dg_loc_coefs[idofs] = linear_sol_phi_vec[idofs];
    }
     
    free(Q_coeff);
    free(R_coeff);
    free(Aloc_matrix);
    free(linear_sol_phi_vec);
    for(i=0;i<linear_NLOC+1;i++)
        free(A_loc_QR[i]);

    free(A_loc_QR);

}

void Aloc_mat_linear_mon(int E,
                     element* mesh_element,
                     node* mesh_node,
                     parameters* user_param,
                     double* linear_Aloc_matrix)
{
    int idofs=0;
    int jdofs=0;
    int k=0;
    double* linear_phi_vec = NULL;

    
    linear_phi_vec = (double*)malloc(linear_NLOC*sizeof(double));

    init_zero_d(linear_Aloc_matrix,linear_NLOC*linear_NLOC);
    for(k=1;k<ngpts+1;k++) 
    {
        //evaluate linear basis functions at gauss points
        init_zero_d(linear_phi_vec,linear_NLOC);
        //init_monomial_basis(E,mesh_element,mesh_node,1,mesh_element[E].el_gpts_ref_x[k],mesh_element[E].el_gpts_ref_y[k],linear_phi_vec);

        init_nodal_basis(E,mesh_element,mesh_node,1,mesh_element[E].el_gpts_ref_x[k],mesh_element[E].el_gpts_ref_y[k],linear_phi_vec);
        for(idofs =1; idofs < linear_NLOC+1;idofs++) 
        {
            for(jdofs =1;jdofs < linear_NLOC+1;jdofs++) 
            {
                linear_Aloc_matrix(idofs,jdofs,linear_NLOC) +=mesh_element[E].det*(linear_phi_vec[idofs-1]*linear_phi_vec[jdofs-1])*mesh_element[E].el_gpts_w[k];
            }//j
        }//
    }//loop over quad points

    free(linear_phi_vec);
}


double limiting_linear_function_val(double uh_K0_bar,
                                    double delta1,
                                    double delta2,
                                    double delta3,
                                    double xhat,
                                    double yhat)
{
    double limited_sol_val =0.0;

    limited_sol_val = uh_K0_bar + delta1*(1.0-2*yhat) + delta2*(2.0*xhat + 2.0*yhat -1.0) + delta3*(1.0-2.0*xhat);

    return limited_sol_val;
}

