#include"limiting_gauss_pts.h"
void limiting_gpts(double* inner_pts_x,
                   double* inner_pts_y,
                   double* w_inner_pts,
                   double* edge1_pts_x,
                   double* edge1_pts_y,
                   double* edge2_pts_x,
                   double* edge2_pts_y,
                   double* edge3_pts_x,
                   double* edge3_pts_y,
                   double* w_edge_pts)
{
    
    int k=0;

    //inner points
    inner_pts_x[0] =  0.5- (sqrt(15.0)/10.0);
    inner_pts_y[0] =  (sqrt(15.0)/20.0) + (1.0/4.0);
    w_inner_pts[0] =  (sqrt(15.0)/81.0) +(5.0/81.0);

    inner_pts_x[1] =  (sqrt(15.0)/20.0) + (1.0/4.0);
    inner_pts_y[1] =  (1.0/2.0) - (sqrt(15.0)/10.0);
    w_inner_pts[1] =  (sqrt(15.0)/81.0) + (5.0/81.0);

    inner_pts_x[2] =  (sqrt(15.0)/20.0) + (1.0/4.0);
    inner_pts_y[2] =  (sqrt(15.0)/20.0) + (1.0/4.0);
    w_inner_pts[2] =  (sqrt(15.0)/81.0) + (5.0/81.0);

    inner_pts_x[3] = 0.5;
    inner_pts_y[3] = 0.25;
    w_inner_pts[3] = 8.0/81.0;

    inner_pts_x[4] = 0.25;
    inner_pts_y[4] = 0.50;
    w_inner_pts[4] = 8.0/81.0;

    inner_pts_x[5] = 0.25;
    inner_pts_y[5] = 0.25;
    w_inner_pts[5] = 8.0/81.0;

    inner_pts_x[6] = (sqrt(15.0)/10.0) + 0.5;
    inner_pts_y[6] = 0.25 - (sqrt(15.0)/20.0);
    w_inner_pts[6] = (5.0/81.0) - (sqrt(15.0)/81.0);

    inner_pts_x[7] = 0.25- (sqrt(15.0)/20.0);
    inner_pts_y[7] = (sqrt(15.0))/10.0 + 0.50;
    w_inner_pts[7] = (5.0/81.0) - (sqrt(15.0)/81.0);

    inner_pts_x[8] = 0.25 - (sqrt(15.0)/20.0);
    inner_pts_y[8] = 0.25 - (sqrt(15.0)/20.0);
    w_inner_pts[8] = (5.0/81.0) - (sqrt(15.0)/81.0);


   /* for(k=0;k<9;k++)
    {
        w_inner_pts[k] = 0.5*w_inner_pts[k];
    }*/

    
  
    //edge1 
    edge1_pts_x[0] = (sqrt(15.0)/10.0) + 0.5;
    edge1_pts_y[0] = 0.0;
    
    edge1_pts_x[1] = 0.5;
    edge1_pts_y[1] = 0.0;

    edge1_pts_x[2] = 0.5 - (sqrt(15.0)/10.0);
    edge1_pts_y[2] = 0.0;



    //edge2
    
    edge2_pts_x[0] = 0.5 - (sqrt(15.0)/10.0);
    edge2_pts_y[0] = (sqrt(15.0)/10.0)+ 0.5;

    edge2_pts_x[1] = 0.5;
    edge2_pts_y[1] = 0.5;

    edge2_pts_x[2] = (sqrt(15.0)/10.0) + 0.5;
    edge2_pts_y[2] = 0.5 - (sqrt(15.0)/10.0);


 
    
    //edge3
    edge3_pts_x[0] = 0.0;
    edge3_pts_y[0] = 0.5 - (sqrt(15.0)/10.0);

    edge3_pts_x[1] = 0.0;
    edge3_pts_y[1] = 0.5;

    edge3_pts_x[2] = 0.0;
    edge3_pts_y[2] = (sqrt(15.0)/10.0)+0.5;

    //edge weights

  
    w_edge_pts[0] = (5.0/162.0);
    w_edge_pts[1] = (4.0/81.0);
    w_edge_pts[2] = (5.0/162.0);

}
