#include"mesh_utils.h"

/*--------------------------------------------------------------------------------------
 * init_mesh_data_structure
 * initializes the mesh data structure 
 *-------------------------------------------------------------------------------------*/

void init_mesh_data_structure(
		element* mesh_element,
		node* mesh_node,
		edge* mesh_edge,
		int num_elements,
		int current_node_count,
		int *num_edges,
		int *active_element_count)
{
    int i,j,k,id;
    int nodes[4];
    int face_node_a, face_node_b;
    int edge_count =0;
    int element_faces[7];
    int edge_index;
    int n_nodes_count=0;
    new_edge_check* edge_matrix= NULL;
    int existing_edges=0;
    existing_edges = *num_edges;
    edge_matrix = (new_edge_check*)malloc((current_node_count+1)*sizeof(new_edge_check));
    for(i=1;i < current_node_count+1;i++) 
    {
        edge_matrix[i].adjacent_count =0;
        for(j=1;j <10;j++) 
        {
            edge_matrix[i].adjacent_nodes[j] = 0;
            edge_matrix[i].connecting_edge[j]=0;
        }
    }
    //loop over the elements in darcy region
    for(id=1;id < num_elements+1;id++) 
    {
        //if element has not been refined
        if(mesh_element[id].refined == 0) 
        {
            *active_element_count +=1;
            mesh_element[id].elmt_index = *active_element_count;

            //loop over nodes of current element
            for(j=1;j<4;j++) 
            {
                nodes[j]= mesh_element[id].vertex[j];
            }
            //label each edge of the current element
            
            //store the coordinates of the edges in order
            //edge 1
            element_faces[1]=  nodes[1];
            element_faces[2] = nodes[2];
            //edge 2
            element_faces[3] = nodes[2];
            element_faces[4] = nodes[3];
            //edge 3
            element_faces[5] = nodes[3];
            element_faces[6] = nodes[1];

            //initialize faces of current element
            for(k=1;k<4;k++) 
            {
                face_node_a = element_faces[2*k-1];
                face_node_b = element_faces[2*k];
                
                //new edge
                new_edge(edge_matrix,face_node_a,face_node_b,&edge_index);

                if(edge_index == 0) 
                {
                    edge_count = edge_count+1;

                    mesh_edge[edge_count+existing_edges].vertex[1] = face_node_a;
                    mesh_edge[edge_count+existing_edges].vertex[2] = face_node_b;
                    mesh_edge[edge_count+existing_edges].neighbour[1] = id;
		    /*mesh_edge[edge_count+existing_edges].neighbour[2] = 0;*/
                    
                    n_nodes_count = edge_matrix[face_node_a].adjacent_count;
                    edge_matrix[face_node_a].adjacent_nodes[n_nodes_count+1] = face_node_b;
                    edge_matrix[face_node_a].connecting_edge[n_nodes_count+1] = edge_count;
                    edge_matrix[face_node_a].adjacent_count = n_nodes_count+1;

                    n_nodes_count = edge_matrix[face_node_b].adjacent_count;
                    edge_matrix[face_node_b].adjacent_nodes[n_nodes_count+1] = face_node_a;
                    edge_matrix[face_node_b].connecting_edge[n_nodes_count+1] = edge_count;
                    edge_matrix[face_node_b].adjacent_count = n_nodes_count+1;
                    
                    mesh_element[id].edge[k] = existing_edges+edge_count;
                }
                //else the edge has already been labelled assign second neigbour
                else 
                {
                    mesh_edge[existing_edges+ edge_index].neighbour[2] = id;
                    mesh_element[id].edge[k] = existing_edges+ edge_index;
                }
            }//loop over faces of element
        }//unrefined
    }//loop over elements
    *num_edges = edge_count;
    for(i=1;i <edge_count+1;i++) 
    {
	if((mesh_edge[i].neighbour[1] >0)&& (mesh_edge[i].neighbour[2] >0))
        {
	    mesh_edge[i].edge_type = INTERIOR;
        }
        else 
        {
            mesh_edge[i].edge_type = EXTERIOR;
        }

    }
    free(edge_matrix);
}

void check_symmetry_edges(element* mesh_element,
		                  node* mesh_node,
		                  edge* mesh_edge,
		                  int Nelts)
{
    int E;
    int Es;
    int i=0;
    int E_iedge=0;
    int Es_edge_index =0;
    int temp=0;
    int k=0;
    int verbose =0;

    
    for(E=1;E<Nelts+1;E++)
    {
        Es = find_symmetric_element(mesh_element,Nelts,mesh_node,E);

        
        if(E < Es)
        {
            if(verbose)
            {
                printf("old edge ordering for Es = %d \n",Es);
                for(i=1;i<4;i++)
                {
                    printf("i %d edge %d \n",i,mesh_element[Es].edge[i]);
                }
            }


            for(i=1;i<4;i++)
            {
                E_iedge = mesh_element[E].edge[i]; 
                Es_edge_index = find_symmetric_edge(mesh_edge,mesh_element,mesh_node,E,Es,E_iedge);
                assert(Es_edge_index);
                if(!(Es_edge_index == i))
                {
                    if(verbose)
                    {
                        printf("E %d iedge %d Es %d Es_edge_index %d \n",E,i,Es,Es_edge_index);
                        //getchar();
                    }
                    temp = mesh_element[Es].edge[i];
                    mesh_element[Es].edge[i] =  mesh_element[Es].edge[Es_edge_index];
                    mesh_element[Es].edge[Es_edge_index]= temp;
                }
                if(verbose)
                {
                    printf("new edge order in Es = %d \n",Es);
                    for(k=1;k<4;k++)
                    {
                        printf("i %d edge %d \n",k,mesh_element[Es].edge[k]);
                    }
                    //getchar();
                }
            }

        }
    }
}

/*-------------------------------------
 * initialize a double array to zero
 * ------------------------------------*/
void init_zero_d(double* vec, //inout:vector
		 int length   //in: length
		 )
{
    int i;
    for(i=0;i < length;i++) {
	vec[i] = 0.0;
    }
}

/*-------------------------------------
 * initialize a double matrix to zero
 * ------------------------------------*/
void init_zero_m(
		double** vec, //inout:vector
		int l1,   //in: l1
		int l2   //in: l2
		)
{
    int i,j;
    for(i=0;i <l1;i++) 
	for(j=0;j <l2;j++) 
	    vec[i][j] = 0.0;
}

/*-----------------------------------
 * initialize an integer array to zero
 * ----------------------------------*/
void init_zero_int(
		int* vec, 
		int length
		)
{
    int i;
    for(i=0;i <length;i++) {
	vec[i] = 0;
    }
}

void new_edge(new_edge_check* edge_matrix,
	      int node_a,
	      int node_b,
	      int* edge)
{
    int i=1;
    int test_node;
    *edge = 0;
    while(edge_matrix[node_a].adjacent_nodes[i] > 0) {
	test_node = edge_matrix[node_a].adjacent_nodes[i];
	if(test_node == node_b) {
	    *edge = edge_matrix[node_a].connecting_edge[i];
	    break;
	}
	i++;
    }
}

/*---------------------------------------
 * get_Element_data
 * -returns det, and inv_BE_T
 *  ------------------------------------*/
double get_Element_data(
        double* inv_BE_T,      //out: inverse of mapping from ref to phys
	    element* mesh_element, //in: element data structure
		node* mesh_node,       //in: node data structure
		int E                 //in: current element
		)
{
    //3 nodes per element
    int n1,n2,n3;
    //coordinates of each node
    double x1,y1,x2,y2,x3,y3;
    double det;

    n1 = mesh_element[E].vertex[1];
    n2 = mesh_element[E].vertex[2];
    n3 = mesh_element[E].vertex[3];
    x1 = mesh_node[n1].coord[0];
    y1 = mesh_node[n1].coord[1];
    x2 = mesh_node[n2].coord[0];
    y2 = mesh_node[n2].coord[1];
    x3 = mesh_node[n3].coord[0];
    y3 = mesh_node[n3].coord[1];

    //fill in matrix BE transpose inverse
    inv_BE_T(1,1,2) = (1.0/mesh_element[E].det)*(y3-y1);
    inv_BE_T(1,2,2) = (1.0/mesh_element[E].det)*(y1-y2);
    inv_BE_T(2,1,2) = (1.0/mesh_element[E].det)*(x1-x3);
    inv_BE_T(2,2,2) = (1.0/mesh_element[E].det)*(x2-x1);

    return mesh_element[E].det;
}//get_Element_data

/************************************************************************
 * computes the transformation from reference element to physical element
 * input: elements, nodes, E, x and y coordinates on reference element
 * output: array phys_coords 
 * **********************************************************************/
void map_to_physical_element(
		element* mesh_element, 
		node* mesh_node, 
		int E, 
		double* phys_coords, 
		double x_hat, 
		double y_hat)
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
     phys_coords[0] = x1 + x_hat*(x2-x1) + y_hat*(x3-x1);
     phys_coords[1] = y1 + x_hat*(y2-y1) + y_hat*(y3-y1);
}

/*******************************************************************************************
 * viod map_to_reference
 * in: mesh_element- array of elements
 * in: mesh_node - array of mesh_nodes
 * in: iface - the index of the current element
 * in: x_phys, y_phys - coordinates on the physical element
 * out: ref_coords - coordinates on reference element corresponding to physical coordinates
 * ******************************************************************************************/
void map_to_reference_element(element* mesh_element, 
                              node* mesh_node, 
                              int E, 
                              double* ref_coords, 
                              double x_phys, 
                              double y_phys)
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

     //if((E == 2) || (E ==24))
      //   printf("E = %d\t n1 = %d\t n2 = %d\t n3 = %d\n",E,n1,n2,n3);
     ref_coords[0] = 1.0/mesh_element[E].det*((y3-y1)*(x_phys-x1) + (x1-x3)*(y_phys-y1));
     ref_coords[1] = 1.0/mesh_element[E].det*((y1-y2)*(x_phys-x1) + (x2-x1)*(y_phys-y1));
}

void load_mesh_data_structure(FILE* meshfile,
			      node* mesh_node,
			      element* mesh_element,
			      int p_degree,
			      int* nnodes,
			      int* nelts,
                           parameters* user_param)
{
    int i=0;
    int node1=0;
    int node2=0;
    int node3=0;
    int num_nodes=0;
    int num_elts=0;
    double x_coord=0.0;
    double y_coord=0.0;
    int result=0;
    int E=0;
    char filename[100];
    string s_buffer;
    double d_buffer;

    //Init infilestream
    filename[0] = '\0';
    strcat(filename,(*user_param).meshfilename);

    ifstream infile(filename);
    double x1,y1,x2,y2,x3,y3,det;

    int verbose = 0;

    infile>>s_buffer;

    //check format of meshfile
    if(s_buffer=="$MeshFormat")
    {
        verbose=0;
        while(s_buffer!="$Nodes")
        {
            infile>>s_buffer;
        }
        infile>>num_nodes;
        for(i=1; i<(num_nodes+1);i++)
        { 
	     //read coordinates for each node  
            infile>>d_buffer; //skip nodeID
            infile>>x_coord;
            infile>>y_coord;
	     mesh_node[i].coord[0]= x_coord;
	     mesh_node[i].coord[1]= y_coord;
            infile>>d_buffer; //skip z_coord
        }
        while(s_buffer!="$Elements")
        {
            infile>>s_buffer;
        }
        infile>>num_elts;
        //fill in the element structure with the connectivity and domain 
        int element_type=0; //gmsh element-type
        int num_coords=0;   //number of nodes of element
        int num_tags=0;
        int E_count=0; // counter for triangles in meshfile	   
        for(E=1; E<(num_elts+1);E++) 
        {
            infile>>d_buffer; //skip elementID
            infile>>element_type;
            switch(element_type)
	     {
                case(1):
		  {
		      num_coords=2;
		      break;
		  }
                case(2):
		  {
		      break;
		  }
                case(15):
		  {
		      num_coords=1;
		      break;
		  }
                default:
		  {
		      cout<<"Unknown element_type ("<<element_type<<") in gmsh-meshfile"<<endl;
		      exit(1);
		  }
            }
            if(element_type!=2) // read only triangles
            {
                //skip tags and coords
                infile>>num_tags;
                for(int i=0;i<num_tags+num_coords;i++)
                {
                    infile>>d_buffer;
                } 
            }
            else 
            {
                E_count++;
                //skip tags
                infile>>num_tags;
                for(int i=0;i<num_tags;i++)
                {
                    infile>>d_buffer;
                } 
                infile>>node1;  
                infile>>node2;
                infile>>node3;
	         mesh_element[E_count].vertex[1]= node1;
	         mesh_element[E_count].vertex[2]= node2;
	         mesh_element[E_count].vertex[3]= node3;
	         mesh_element[E_count].domain =1;
	         mesh_element[E_count].refined = 0;
	         mesh_element[E_count].degree = p_degree;

                // Calculate det 
	         x1 = mesh_node[node1].coord[0];
	         y1 = mesh_node[node1].coord[1];
	         x2 = mesh_node[node2].coord[0];
	         y2 = mesh_node[node2].coord[1];
	         x3 = mesh_node[node3].coord[0];
	         y3 = mesh_node[node3].coord[1];

                det = (x2-x1)*(y3-y1)-(y2-y1)*(x3-x1);
	         if (det<=0) 
                {
	             mesh_element[E_count].vertex[1] = node3;
	             mesh_element[E_count].vertex[2] = node2;
	             mesh_element[E_count].vertex[3] = node1;

	             node1 = mesh_element[E_count].vertex[1];
	             node2 = mesh_element[E_count].vertex[2];
	             node3 = mesh_element[E_count].vertex[3];
	             x1 = mesh_node[node1].coord[0];
	             y1 = mesh_node[node1].coord[1];
	             x2 = mesh_node[node2].coord[0];
	             y2 = mesh_node[node2].coord[1];
	             x3 = mesh_node[node3].coord[0];
	             y3 = mesh_node[node3].coord[1];

	             if (verbose == 1)
                    { 
		          printf("WARNING::Element %d may have wrong ordering of nodes\n", E_count);
		          printf("node1 %d node2 %d node3 %d \n",node1,node2,node3);
		          printf("x1 %lf y1 %lf \n",x1,y1);
		          printf("x2 %lf y2 %lf \n",x2,y2);
		          printf("x3 %lf y3 %lf \n",x3,y3);

	             }
	             det = (x2-x1)*(y3-y1)-(y2-y1)*(x3-x1);

	             if (det<=0)
	             {
                     exit(1);
		          mesh_element[E_count].vertex[1] = node2;
		          mesh_element[E_count].vertex[2] = node3;
		          mesh_element[E_count].vertex[3] = node1;

		          node1 = mesh_element[E_count].vertex[1];
		          node2 = mesh_element[E_count].vertex[2];
		          node3 = mesh_element[E_count].vertex[3];
		          x1 = mesh_node[node1].coord[0];
		          y1 = mesh_node[node1].coord[1];
		          x2 = mesh_node[node2].coord[0];
		          y2 = mesh_node[node2].coord[1];
		          x3 = mesh_node[node3].coord[0];
		          y3 = mesh_node[node3].coord[1];

		          if (verbose == 1)
                        { 
		              printf("WARNING::Element %d may have wrong ordering of nodes\n", E_count);
		              printf("node1 %d node2 %d node3 %d \n",node1,node2,node3);
		              printf("x1 %lf y1 %lf \n",x1,y1);
		              printf("x2 %lf y2 %lf \n",x2,y2);
		              printf("x3 %lf y3 %lf \n",x3,y3);
		              getchar();
		              det = (x2-x1)*(y3-y1)-(y2-y1)*(x3-x1);	
		          }
	             }
                }

	         mesh_element[E_count].det = fabs(det);
            }          
        }//loop over all E
        *nnodes = num_nodes;
        *nelts = E_count;    
    }
    else
    {
        result =fscanf(meshfile, "%d",&num_nodes);
        if(result == 0)
	 printf("Read no number of number of nodes \n");

        for(i=1; i<(num_nodes+1);i++)
        {
	    //read in the different values for each of the node parameters
	    result =fscanf(meshfile,"%lf %lf", &x_coord, &y_coord);
	    if(result==0)
	    {
	        printf("Read no node information list at node %d \n",i);
	    }    
	    mesh_node[i].coord[0]= x_coord;
	    mesh_node[i].coord[1]= y_coord;
        }
        result =fscanf(meshfile, "%d",&num_elts);
        if(result == 0)
        {
	    printf("Read no number of elements \n");
        }

        //fill in the element structure with the connectivity and domain 	
        for(E=1; E<(num_elts+1);E++) 
        {
	    result =fscanf(meshfile, "%d %d %d", &node1, &node2, &node3);
	    if(result ==0)
	    {
	        printf("Read no element list at %d \n",i);
	    }
	    mesh_element[E].vertex[1]= node1;
	    mesh_element[E].vertex[2]= node2;
	    mesh_element[E].vertex[3]= node3;
	    mesh_element[E].domain =1;
	    mesh_element[E].refined = 0;
	    mesh_element[E].degree = p_degree;

	    // Calculate det 
	    x1 = mesh_node[node1].coord[0];
	    y1 = mesh_node[node1].coord[1];
	    x2 = mesh_node[node2].coord[0];
	    y2 = mesh_node[node2].coord[1];
	    x3 = mesh_node[node3].coord[0];
	    y3 = mesh_node[node3].coord[1];


	    det = (x2-x1)*(y3-y1)-(y2-y1)*(x3-x1);
	    if (det<=0) {
	        mesh_element[E].vertex[1] = node3;
	        mesh_element[E].vertex[2] = node2;
	        mesh_element[E].vertex[3] = node1;

	        node1 = mesh_element[E].vertex[1];
	        node2 = mesh_element[E].vertex[2];
	        node3 = mesh_element[E].vertex[3];
	        x1 = mesh_node[node1].coord[0];
	        y1 = mesh_node[node1].coord[1];
	        x2 = mesh_node[node2].coord[0];
	        y2 = mesh_node[node2].coord[1];
	        x3 = mesh_node[node3].coord[0];
	        y3 = mesh_node[node3].coord[1];

	        if (verbose == 1){ 
		    printf("WARNING::Element %d may have wrong ordering of nodes\n", E);
		    printf("node1 %d node2 %d node3 %d \n",node1,node2,node3);
		    printf("x1 %lf y1 %lf \n",x1,y1);
		    printf("x2 %lf y2 %lf \n",x2,y2);
		    printf("x3 %lf y3 %lf \n",x3,y3);
		    //getchar();
	        }
	        det = (x2-x1)*(y3-y1)-(y2-y1)*(x3-x1);

	        if (det<=0)
	        {
		    mesh_element[E].vertex[1] = node2;
		    mesh_element[E].vertex[2] = node3;
		    mesh_element[E].vertex[3] = node1;

		    node1 = mesh_element[E].vertex[1];
		    node2 = mesh_element[E].vertex[2];
		    node3 = mesh_element[E].vertex[3];
		    x1 = mesh_node[node1].coord[0];
		    y1 = mesh_node[node1].coord[1];
		    x2 = mesh_node[node2].coord[0];
		    y2 = mesh_node[node2].coord[1];
		    x3 = mesh_node[node3].coord[0];
		    y3 = mesh_node[node3].coord[1];

		    if (verbose == 1){ 
		        printf("WARNING::Element %d may have wrong ordering of nodes\n", E);
		        printf("node1 %d node2 %d node3 %d \n",node1,node2,node3);
		        printf("x1 %lf y1 %lf \n",x1,y1);
		        printf("x2 %lf y2 %lf \n",x2,y2);
		        printf("x3 %lf y3 %lf \n",x3,y3);
		        getchar();
		        det = (x2-x1)*(y3-y1)-(y2-y1)*(x3-x1);	
		    }
	        }
	    }

	    mesh_element[E].det = fabs(det);
        }

        *nnodes = num_nodes;
        *nelts = num_elts;
    }
}

/*void make_boundary_periodic
 * This assumes that the computational domain
 * is rectangular/square
 */
void make_boundary_periodic(edge* mesh_edge,
                            node* mesh_node,
                            element* mesh_element,
                            int num_edges,
                            parameters* user_param
                            )
{
    int iedge=0;
    int node_a=0;
    int node_b=0;
    int found =0;
    double xcoord_na =0.0;
    double ycoord_na =0.0;
    double xcoord_nb =0.0;
    double ycoord_nb =0.0;

    //loop over the edges and classify them TOP, BOTTOM, LEFT, RIGHT BOUNDARY
    for(iedge=1;iedge<num_edges+1;iedge++)
    {
        if(mesh_edge[iedge].edge_type == EXTERIOR)
        {
            node_a = mesh_edge[iedge].vertex[1];
            xcoord_na = mesh_node[node_a].coord[0];
            ycoord_na = mesh_node[node_a].coord[1];

            node_b = mesh_edge[iedge].vertex[2];
            xcoord_nb = mesh_node[node_b].coord[0];
            ycoord_nb = mesh_node[node_b].coord[1];

            if(fabs(ycoord_nb - ycoord_na)< EPSILON)
            {
                mesh_edge[iedge].slope = HORIZONTAL;

                if(fabs(ycoord_nb) < EPSILON)
                {
                    mesh_edge[iedge].boundary_side = BOTTOM;
                }
                else if(fabs(ycoord_nb - (*user_param).periodic_domain_size_y) < EPSILON)
                {
                    mesh_edge[iedge].boundary_side = TOP;
                }
                else
                {
                    fprintf(stderr,"Error in boundary edge classification on horizontal boundary\n");
                    exit(1);
                }
            }
            else if(fabs(xcoord_nb - xcoord_na) < EPSILON)
            {
                mesh_edge[iedge].slope = VERTICAL;
                if(fabs(xcoord_nb) < EPSILON)
                {
                    mesh_edge[iedge].boundary_side = LEFT;
                }
                else if(fabs(xcoord_nb- (*user_param).periodic_domain_size_x) < EPSILON)
                {
                    mesh_edge[iedge].boundary_side = RIGHT;
                }
                else
                {
                     fprintf(stderr,"Error in boundary edge classification on vertical boundary\n");
                     exit(1);
                }
            }
            else
            {
                fprintf(stderr,"Error in edge slopes in mesh utils make_boundary_periodic. Domain not rectangular?\n");
                exit(1);
            }

        }
    }//loop over edges

    //loop over edges and search for corresponding edge on periodic side
    for(iedge =1;iedge<num_edges+1;iedge++)
    {
        if(mesh_edge[iedge].edge_type == EXTERIOR)
        {
            found = search_periodic_edge(iedge,mesh_edge,mesh_node,mesh_element,num_edges);
            assert(found >0);
        }
    }

}


int search_periodic_edge(int current_edge,
                         edge* mesh_edge,
                         node* mesh_node,
                         element* mesh_element,
                         int num_edges
                         )
{
   int iedge=0;
   int found =0;
   int node_a_current_edge=0;
   int node_b_current_edge=0;
   int node_a_iedge=0;
   int node_b_iedge=0;
   int E_current_edge=0;
   int E_iedge=0;
   double diff=0.0;

   node_a_current_edge = mesh_edge[current_edge].vertex[1];
   node_b_current_edge = mesh_edge[current_edge].vertex[2];
   E_current_edge = mesh_edge[iedge].neighbour[1];

   iedge=1;

   while((found ==0) &&(iedge < (num_edges+1)))
   {
       if((mesh_edge[iedge].edge_type == EXTERIOR) && (!(iedge == current_edge)))
       {
           node_a_iedge = mesh_edge[iedge].vertex[1];
           node_b_iedge = mesh_edge[iedge].vertex[2];
           if(mesh_edge[iedge].slope == mesh_edge[current_edge].slope)
           {
               if((mesh_edge[iedge].slope == HORIZONTAL))
               {
                   if(((fabs(mesh_node[node_a_current_edge].coord[0]-mesh_node[node_a_iedge].coord[0])< EPSILON) &&
                      (fabs(mesh_node[node_b_current_edge].coord[0]-mesh_node[node_b_iedge].coord[0])< EPSILON)) ||
                      ((fabs(mesh_node[node_a_current_edge].coord[0]-mesh_node[node_b_iedge].coord[0])< EPSILON) &&
                      (fabs(mesh_node[node_b_current_edge].coord[0]-mesh_node[node_a_iedge].coord[0])< EPSILON)))
                   {
                       //printf("found! \n");
                       //getchar();
                       found =iedge;
                       mesh_edge[iedge].periodic_found = found;

                       mesh_edge[current_edge].periodic_nedge = found;
                       mesh_edge[current_edge].neighbour[2] = mesh_edge[iedge].neighbour[1];
                       mesh_edge[iedge].neighbour[2] = mesh_edge[current_edge].neighbour[1];
                   }
                   else
                   {
                       found =0;
                   }
               }
               else if((mesh_edge[iedge].slope == VERTICAL))
               {
                   node_a_iedge = mesh_edge[iedge].vertex[1];
                   node_b_iedge = mesh_edge[iedge].vertex[2];

                   if(((fabs(mesh_node[node_a_current_edge].coord[1]-mesh_node[node_a_iedge].coord[1])< EPSILON) &&
                      (fabs(mesh_node[node_b_current_edge].coord[1]-mesh_node[node_b_iedge].coord[1])< EPSILON)) ||
                      ((fabs(mesh_node[node_a_current_edge].coord[1]-mesh_node[node_b_iedge].coord[1])< EPSILON) &&
                      (fabs(mesh_node[node_b_current_edge].coord[1]-mesh_node[node_a_iedge].coord[1])< EPSILON)))
                   {
                       found =iedge;
                       mesh_edge[iedge].periodic_found = found;
                       mesh_edge[current_edge].periodic_nedge = found;
                       mesh_edge[current_edge].neighbour[2] = mesh_edge[iedge].neighbour[1];
                       mesh_edge[iedge].neighbour[2] = mesh_edge[current_edge].neighbour[1];

                   }
                   else
                   {
                       found =0;
                   }
               }
           }
           else 
           {
               found =0;
           }
       }
       else
       {
           found =0;
       }

       iedge++;
   }

   return found;
}

void label_reflective_boundary(int num_edges,
                               edge* mesh_edge,
                               node* mesh_node,
                               parameters* user_param)
{
    int iedge=0;
    int node_a =0;
    int node_b =0;
    double xcoord_na =0.0;
    double ycoord_na =0.0;
    double xcoord_nb =0.0;
    double ycoord_nb =0.0;

    for(iedge =1;iedge<num_edges+1;iedge++)
    {
        if(mesh_edge[iedge].edge_type == EXTERIOR)
        {
            node_a = mesh_edge[iedge].vertex[1];
            node_b = mesh_edge[iedge].vertex[2];

            xcoord_na = mesh_node[node_a].coord[0];
            ycoord_na = mesh_node[node_a].coord[1];

            xcoord_nb = mesh_node[node_b].coord[0];
            ycoord_nb = mesh_node[node_b].coord[1];

            switch ((*user_param).test_case) 
            {
                case(10): 
                    {
                        if((xcoord_na == 1.0) && (xcoord_nb == 1.0)) 
                       {
                           mesh_edge[iedge].reflective=1;
                       }
                        else 
                        {
                            mesh_edge[iedge].reflective=0;
                        }
                        break;
                    }
                case(11): 
                    {
                        if(((ycoord_na == 1.0) && (ycoord_nb == 1.0)) || ((ycoord_na == 0.0) && (ycoord_nb == 0.0))) 
                        {
                            mesh_edge[iedge].reflective=1;
                        }
                        else 
                        {
                            mesh_edge[iedge].reflective=0;
                        }
                        break;
                    }
                case(17): 
                    {
                        if(((xcoord_na ==0.0) && (xcoord_nb==0.0)) || ((xcoord_na ==11.0) &&(xcoord_nb==11.0)))
                        {
                            mesh_edge[iedge].reflective=0;
                        }
                        else
                        {
                            mesh_edge[iedge].reflective=1;
                        }
                        break;
                    }
                case(26):
                    {
                        if (((fabs(ycoord_na) < EPSILON) && (fabs(ycoord_nb) < EPSILON)) || ((fabs(ycoord_na -6.0) < EPSILON) && (fabs(ycoord_nb-6.0) < EPSILON)))
                        {
                            mesh_edge[iedge].reflective=1;
                        }
                        else
                        {
                            mesh_edge[iedge].reflective=0;
                        }
                        break;
                    }
                default: 
                    {
                        mesh_edge[iedge].reflective = 0;
                        break;
                    }
            }//switch
        }//ext
    }//loop over edges
}//label edges


void check_symmetry_of_gpts(element* mesh_element,
                            node* mesh_node,
                            edge* mesh_edge,
                            int global_element_count,
                            int global_edge_count)
{
    int E=0;
    int Es=0;
    int k=0;
    int j=0;
    int found =0;
    double temp_gpts_data[5];
    double temp_basis_val[NLOC];
    int i=0;
    int E_iedge = 0;
    int Es_iedge =0;
    int Es_edge_index=0;
    int verbose=1;
    //check symmetry of volume gpts

    for(E=1;E<global_element_count+1;E++)
    {
        Es = find_symmetric_element(mesh_element,global_element_count,mesh_node,E);
        printf("E %d Es %d \n",E,Es);
        //getchar();
        if(E < Es)
        {
            j=1;

            for(k=1;k<ngpts+1;k++)
            {
                found=0;
                j=1;
                while(j < (ngpts+1) && (found ==0))
                {
                    if((fabs(mesh_element[E].el_gpts_y[k] - (1.0 - mesh_element[Es].el_gpts_x[j])) < EPSILON) &&
                       (fabs(mesh_element[Es].el_gpts_y[j] - (1.0 - mesh_element[E].el_gpts_x[k])) < EPSILON))
                    {
                        found =j;
                    }
                    j++;
                }
                //printf("k = %d and found = %d \n",k,found);
                if(!(found == k))
                {
                    init_zero_d(temp_gpts_data,5);
                    //printf("found = %d in Es but k = %d in E \n",found,k);

                    temp_gpts_data[0] = mesh_element[Es].el_gpts_x[k];  // phys_coords[0];
                    temp_gpts_data[1] = mesh_element[Es].el_gpts_y[k];   // phys_coords[1];
                    temp_gpts_data[2] = mesh_element[Es].el_gpts_ref_x[k];//gpx[k];
                    temp_gpts_data[3] = mesh_element[Es].el_gpts_ref_y[k];//gpy[k];
                    temp_gpts_data[4] = mesh_element[Es].el_gpts_w[k];    //w_el[k];
                    
                    mesh_element[Es].el_gpts_x[k] = mesh_element[Es].el_gpts_x[found]; 
                    mesh_element[Es].el_gpts_y[k] = mesh_element[Es].el_gpts_y[found];
                    mesh_element[Es].el_gpts_ref_x[k] = mesh_element[Es].el_gpts_ref_x[found];
                    mesh_element[Es].el_gpts_ref_y[k] = mesh_element[Es].el_gpts_ref_y[found];
                    mesh_element[Es].el_gpts_w[k] = mesh_element[Es].el_gpts_w[found];

                    mesh_element[Es].el_gpts_x[found] = temp_gpts_data[0];
                    mesh_element[Es].el_gpts_y[found] = temp_gpts_data[1];
                    mesh_element[Es].el_gpts_ref_x[found] = temp_gpts_data[2];
                    mesh_element[Es].el_gpts_ref_y[found] = temp_gpts_data[3];
                    mesh_element[Es].el_gpts_w[found] = temp_gpts_data[4];

                    init_zero_d(temp_basis_val,NLOC);
                    for(i=0;i<NLOC;i++)
                    {
                        temp_basis_val[i] = mesh_element[Es].el_gpts_basis[(k-1)*NLOC+i];
                    }
                    for(i=0;i<NLOC;i++)
                    {
                        mesh_element[Es].el_gpts_basis[(k-1)*NLOC+i] = mesh_element[Es].el_gpts_basis[(found-1)*NLOC+i];
                    }
                    for(i=0;i<NLOC;i++)
                    {
                        mesh_element[Es].el_gpts_basis[(found-1)*NLOC+i] = temp_basis_val[i];
                    }

                    //getchar();
                }
            }//loop over E gpts
        }

    }//loop over elements

    //check symmetry of edges
    
    for(E=1; E < global_element_count +1;E++)
    {
        Es = find_symmetric_element(mesh_element,global_element_count,mesh_node,E);
        if(E < Es)
        {
            for(i=1;i<4;i++)
            {
                E_iedge = mesh_element[E].edge[i];
                Es_edge_index = find_symmetric_edge(mesh_edge,mesh_element,mesh_node,E,Es,E_iedge);
                assert(Es_edge_index == i);
            }
        }
    }

    /*if(verbose)
    {
        for(E=1;E<global_element_count+1;E++)
        {
            Es = find_symmetric_element(mesh_element,global_element_count,mesh_node,E);
            if(E < Es)
            {
               printf("gpx  %10.12lf %10.12lf \n",mesh_element[Es].el_gpts_x[k] = mesh_element[Es].el_gpts_x[found]; 
               mesh_element[Es].el_gpts_y[k] = mesh_element[Es].el_gpts_y[found];
               mesh_element[Es].el_gpts_ref_x[k] = mesh_element[Es].el_gpts_ref_x[found];
               mesh_element[Es].el_gpts_ref_y[k] = mesh_element[Es].el_gpts_ref_y[found];
               mesh_element[Es].el_gpts_w[k] = mesh_element[Es].el_gpts_w[found];*/


            

}//check_symmetry_of_gpts
