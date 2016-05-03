/*
Monomial basis
Prince Chidyagwai
Philipp Monreal
*/
#include"dg_monomial_basis.h"

/* void init_monomial_basis
 * calculates value and derivatives of basis functions
 *in: deg : polynomial degree
 *in: xhat yhat : x and y coordinates of point
 *in: E : elment number
 *in: mesh_element
 *in: mesh_node
 *out: phi_vec, grad_phi_mat ==> value and derivates of basis functions
 */
void init_monomial_basis(int E,
			 element* mesh_element, 
			 node* mesh_node,
			 int deg, 
			 double xhat, 
			 double yhat,
			 double* phi_vec)
{
    //local variables
    int i,j,ii, Nloc;
    double det_BE_T;
    double* ssn = NULL;
    double* ttn = NULL;

    ssn  = (double*)malloc((deg+1)*sizeof(double));
    ttn  = (double*)malloc((deg+1)*sizeof(double));

    Nloc = ((deg+1)*(deg+2))/2;

    init_zero_d(ssn,deg+1);
    init_zero_d(ttn,deg+1);

    det_BE_T = mesh_element[E].det;

    //double* inv_BE_T = NULL;
    //inv_BE_T = (double*)malloc(4*sizeof(double));
    //init_zero_d(inv_BE_T,4);
    //det_BE_T = get_Element_data(inv_BE_T,mesh_element,mesh_node,E);
    //if ((det_BE_T - mesh_element[E].det) > 1e-16) 
   // {
	//printf("det(store) = %lf ",mesh_element[E].det);
	//getchar();
    //}

    //1.0
    ssn[0] = 1.0;
    ttn[0] = 1.0;

    if (deg > 0){
	//x
	ssn[1] = xhat;
	ttn[1] = yhat;

	for(i=2;i < deg+1;i++) {
	    ssn[i] = xhat*ssn[i-1];
	    ttn[i] = yhat*ttn[i-1];
	}
    }
    //values of basis functions
    ii=0;
    for(i=0;i<deg+1;i++) {
	for(j=0;j <i+1;j++) {
	    phi_vec[ii] = ssn[i-j]*ttn[j];
	    ii = ii+1;
	}//j
    }//i

    //free(inv_BE_T);
    free(ssn);
    free(ttn);
}

/* void init_monomial_basis_deriv
 * calculates value and derivatives of basis functions
 *in: deg : polynomial degree
 *in: xhat yhat : x and y coordinates of point
 *in: E : elment number
 *in: mesh_element
 *in: mesh_node
 *out: phi_vec, grad_phi_mat ==> value and derivates of basis functions
 */
void init_monomial_basis_deriv(int E,
			 element* mesh_element, 
			 node* mesh_node,
			 int deg, 
			 double xhat, 
			 double yhat,
			 double* phi_vec, 
			 double* grad_phi_mat)
{
    //local variables
    int i,j,ii, Nloc;
    double det_BE_T,temp;
    double* inv_BE_T = NULL;
    double* ssn = NULL;
    double* ttn = NULL;
    double* dssn = NULL;
    double* dttn = NULL;

    ssn  = (double*)malloc((deg+1)*sizeof(double));
    ttn  = (double*)malloc((deg+1)*sizeof(double));
    dssn = (double*)malloc((deg+1)*sizeof(double));
    dttn = (double*)malloc((deg+1)*sizeof(double));
    
    inv_BE_T = (double*)malloc(2*2*sizeof(double));
    Nloc = ((deg+1)*(deg+2))/2;

    init_zero_d(ssn,deg+1);
    init_zero_d(ttn,deg+1);
    init_zero_d(dssn,deg+1);
    init_zero_d(dttn,deg+1);
    init_zero_d(inv_BE_T,4);

    det_BE_T = get_Element_data(inv_BE_T,mesh_element,mesh_node,E);

    //1.0
    ssn[0] = 1.0;
    ttn[0] = 1.0;
    dssn[0] = 0.0;
    dttn[0] = 0.0;

    if (deg > 0){
	//x
	ssn[1] = xhat;
	ttn[1] = yhat;
	dssn[1] = 1.0;
	dttn[1] = 1.0;

	for(i=2;i < deg+1;i++) {
	    ssn[i] = xhat*ssn[i-1];
	    ttn[i] = yhat*ttn[i-1];
	    dssn[i] = xhat*dssn[i-1];
	    dttn[i] = yhat*dttn[i-1];
	}
	//adjust constants for derivatives
	for(i=2; i < deg+1;i++) {
	    dssn[i] *= i;
	    dttn[i] *= i;
	}
    }
    //values of basis functions
    ii=0;
    for(i=0;i<deg+1;i++) {
	for(j=0;j <i+1;j++) {
	    phi_vec[ii] = ssn[i-j]*ttn[j];
	    ii = ii+1;
	}//j
    }//i

    //derivatives of basis functions
    ii =1;
    for(i=0;i < deg+1 ;i++) {
	for(j=0;j <i+1;j++) {
	    //local derivatives
	    grad_phi_mat(1,ii,2) = dssn[i-j]*ttn[j];
	    grad_phi_mat(2,ii,2) = ssn[i-j]*dttn[j];

	    //global derivatives
	    temp = grad_phi_mat(1,ii,2)*inv_BE_T(1,1,2) + grad_phi_mat(2,ii,2)*inv_BE_T(1,2,2);
	    grad_phi_mat(2,ii,2) = grad_phi_mat(1,ii,2)*inv_BE_T(2,1,2) + grad_phi_mat(2,ii,2)*inv_BE_T(2,2,2);
	    grad_phi_mat(1,ii,2) = temp;
	    ii = ii+1;
	}//j
    }//i
    free(ssn);
    free(ttn);
    free(dssn);
    free(dttn);
    free(inv_BE_T);
}

void init_monomial_basis_limit(int E,
                              element* mesh_element,
                              node* mesh_node,
                              int deg,
                              double xhat,
                              double yhat,
                              double* phi_vec)
{
    if(deg ==1)
    {
        phi_vec[0] = 1.0-2.0*yhat;
        phi_vec[1] = 2.0*xhat+2.0*yhat-1.0;
        phi_vec[2] = 1.0-2.0*xhat;
    }

}

