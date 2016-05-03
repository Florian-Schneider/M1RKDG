#include"refine.h"
/*
This modules refines 2D triangles
Each cell will be associated with a radiative transfer property
which will replace the permeability variable
*/
void refine_mesh(
	element* mesh_element,          //mesh element data structure
	node* mesh_node,		//mesh node data structure
	edge* mesh_edge,		//mesh edge data structure
	int *global_node_count,	        //total number of nodes
	int *global_element_count)      //total number of elements
{
    int E=0;
    int i=0;	
    int edge1=0;
    int edge2=0;
    int edge3=0;
    int refined_vec[7]; 
    double mid_point[3];
    int node_a =0;
    int node_b =0;
    int refined_node_count=0;
    int rcount = 0;
    int num_elts=0;
    int num_nodes_coarse=0;
    int element_refine_count = 0;
    
    num_elts = *global_element_count;
    num_nodes_coarse = *global_node_count;

    init_zero_int(refined_vec,7);
    init_zero_d(mid_point,3);

    for(E=1; E < num_elts +1; E++) 
    {
	    if(mesh_element[E].refined ==0) 
	    {
		    init_zero_int(refined_vec,7);
		    //edge1
		    edge1 = mesh_element[E].edge[1];
		    node_a = mesh_edge[edge1].vertex[1];
		    node_b = mesh_edge[edge1].vertex[2];
		    //edge1 has not been bisected 
		    if(mesh_edge[edge1].refined==0) 
		    {
			    refined_node_count++;
			    get_mid_point(mesh_edge,edge1,mesh_node,mid_point);
			    mesh_node[num_nodes_coarse+refined_node_count].coord[0] = mid_point[1];
		    	    mesh_node[num_nodes_coarse+refined_node_count].coord[1] = mid_point[2];

			    if(mesh_edge[edge1].edge_type == INTERIOR) 
			    {
				     mesh_node[num_nodes_coarse+refined_node_count].node_type = 1;
			    }//if interior
			   else if(mesh_edge[edge1].edge_type == EXTERIOR) 
			   {
				    mesh_node[num_nodes_coarse+refined_node_count].node_type = -1;
			   }//else exterior

			   refined_vec[1] = node_a;
			   refined_vec[2] = num_nodes_coarse+refined_node_count;
			   refined_vec[3] = node_b;
			   mesh_edge[edge1].midpoint = num_nodes_coarse+refined_node_count;
			   mesh_edge[edge1].refined = 1;
		    }//if refined==0
	    	   else if(mesh_edge[edge1].refined ==1) 
		   {
			   refined_vec[1] = node_b;
			   refined_vec[2] = mesh_edge[edge1].midpoint;
			   refined_vec[3] = node_a;
		   }
		    //edge2
		   edge2 = mesh_element[E].edge[2];
		   node_a = mesh_edge[edge2].vertex[1];
		   node_b = mesh_edge[edge2].vertex[2];
	    	   if(mesh_edge[edge2].refined ==0) 
		   {
			   refined_node_count++;
			   get_mid_point(mesh_edge,edge2,mesh_node,mid_point);
			   mesh_node[num_nodes_coarse+refined_node_count].coord[0] = mid_point[1];
			   mesh_node[num_nodes_coarse+refined_node_count].coord[1] = mid_point[2];
			   refined_vec[4] = num_nodes_coarse+refined_node_count;
			   refined_vec[5] = node_b;
			   mesh_edge[edge2].midpoint= num_nodes_coarse+refined_node_count;
			   mesh_edge[edge2].refined = 1;

			  if(mesh_edge[edge2].edge_type == INTERIOR) 
			  {
				   mesh_node[num_nodes_coarse+refined_node_count].node_type = 1;
			  }
			 else if(mesh_edge[edge2].edge_type == EXTERIOR) 
			 {
				 mesh_node[num_nodes_coarse+refined_node_count].node_type = -1;
			 }
		   }
		   else if(mesh_edge[edge2].refined ==1) 
		   {
			   refined_vec[4] = mesh_edge[edge2].midpoint;
			   refined_vec[5] = node_a;
		   }

		    //edge3
		   edge3 = mesh_element[E].edge[3];
		   node_a = mesh_edge[edge3].vertex[1];
		   node_b = mesh_edge[edge3].vertex[2];
	    	   if(mesh_edge[edge3].refined ==0) 
		   {
			   refined_node_count++;
			   get_mid_point(mesh_edge,edge3,mesh_node,mid_point);
			   mesh_node[num_nodes_coarse+refined_node_count].coord[0] = mid_point[1];
			   mesh_node[num_nodes_coarse+refined_node_count].coord[1] = mid_point[2];

			   refined_vec[6] = num_nodes_coarse+refined_node_count;
		           mesh_edge[edge3].midpoint = num_nodes_coarse+refined_node_count;

			  if(mesh_edge[edge3].edge_type == INTERIOR) 
			  {
				  mesh_node[num_nodes_coarse+refined_node_count].node_type = 1;
			  }
			  else if(mesh_edge[edge3].edge_type == EXTERIOR) 
			  {
				   mesh_node[num_nodes_coarse+refined_node_count].node_type = -1;
			  }
			  mesh_edge[edge3].refined = 1;
		   }
	   	   else if(mesh_edge[edge3].refined ==1) 
		   {
			   refined_vec[6] = mesh_edge[edge3].midpoint;
		   }

	    	   //children of element E
	    	   for(i=1;i < 5;i++) 
		   {
			   mesh_element[E].children[i] = num_elts+element_refine_count+i;
		   }

		    //child 1
		    mesh_element[num_elts+element_refine_count+1].vertex[1] = refined_vec[1];
		    mesh_element[num_elts+element_refine_count+1].vertex[2] = refined_vec[2];
		    mesh_element[num_elts+element_refine_count+1].vertex[3] = refined_vec[6];
		    mesh_element[num_elts+element_refine_count+1].domain = mesh_element[E].domain;
		    mesh_element[num_elts+element_refine_count+1].parent = E;
		    mesh_element[num_elts+element_refine_count+1].refined = 0;
		    mesh_element[num_elts+element_refine_count+1].degree = mesh_element[E].degree;
		    //child 2
		    mesh_element[num_elts+element_refine_count+2].vertex[1] = refined_vec[2];
		    mesh_element[num_elts+element_refine_count+2].vertex[2] = refined_vec[3];
		    mesh_element[num_elts+element_refine_count+2].vertex[3] = refined_vec[4];
		    mesh_element[num_elts+element_refine_count+2].domain = mesh_element[E].domain;
		    mesh_element[num_elts+element_refine_count+2].parent = E;
		    mesh_element[num_elts+element_refine_count+2].refined = 0;
		    mesh_element[num_elts+element_refine_count+2].degree = mesh_element[E].degree;

		    //child 3
		    mesh_element[num_elts+element_refine_count+3].vertex[1] = refined_vec[2];
		    mesh_element[num_elts+element_refine_count+3].vertex[2] = refined_vec[4];
		    mesh_element[num_elts+element_refine_count+3].vertex[3] = refined_vec[6];
		    mesh_element[num_elts+element_refine_count+3].domain = mesh_element[E].domain;
		    mesh_element[num_elts+element_refine_count+3].parent = E;
		    mesh_element[num_elts+element_refine_count+3].refined = 0;
		    mesh_element[num_elts+element_refine_count+3].degree = mesh_element[E].degree;


		    //child 4
		    mesh_element[num_elts+element_refine_count+4].vertex[1] = refined_vec[4];
		    mesh_element[num_elts+element_refine_count+4].vertex[2] = refined_vec[5];
		    mesh_element[num_elts+element_refine_count+4].vertex[3] = refined_vec[6];
		    mesh_element[num_elts+element_refine_count+4].domain = mesh_element[E].domain;
		    mesh_element[num_elts+element_refine_count+4].parent = E;
		    mesh_element[num_elts+element_refine_count+4].refined = 0;
		    mesh_element[num_elts+element_refine_count+4].degree = mesh_element[E].degree;

		    mesh_element[E].refined=1;
		    element_refine_count  = element_refine_count+4;
		    rcount++;
	    }//if unrefined
    }
    *global_node_count  = *global_node_count+refined_node_count;
    *global_element_count = *global_element_count+ element_refine_count;
}

void  get_mid_point(
		edge* mesh_edge,
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

    mid_point[1] = 0.5*(x1+x2);
    mid_point[2] = 0.5*(y1+y2);
}

void set_neighbour(edge* mesh_edge,
		   int iedge,
		   int E)
{
    if(mesh_edge[iedge].neighbour[1] == 0) 
    {
	mesh_edge[iedge].neighbour[1] = E;
    }
    else 
    {
	mesh_edge[iedge].neighbour[2] = E;
    }
}


/****************************************
 * reindex the free nodes for solve routine
 * ****************************************/
void reindex_mesh(
		element* mesh_element,
		node* mesh_node,
		int num_elts,
		int* active_elmt_count)
{
    int E=0;
    int acounter=0;

    for(E=1;E < num_elts+1;E++) 
    {
	if(mesh_element[E].refined == 0) 
	{
	    acounter++;
	    mesh_element[E].elmt_index = acounter;
	}
    }
    *active_elmt_count =acounter;
}
