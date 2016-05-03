#include <complex>
using namespace std;

typedef complex<double> dcmplx;

double e = E[j];
            if (e<=0)
                e=1.0e-10;
            
            double fx = Fx[j]/e;
            double fy = Fy[j]/e;
            if (fx*fx+fy*fy>1.0-1.0e-5)
                {
                    double tmp = sqrt(fx*fx+fy*fy);
                    fx = fx/tmp*(1.0-1.0e-5);
                    fy = fy/tmp*(1.0-1.0e-5);
                
                }
            
            
            if (abs(fx)<1.0e-12)
                {
                    fx = 1.0e-12*sqrt(0.5); //disturb f   
                }
            if (abs(fy)<1.0e-12)
                {
                    fy = -1.0e-12*sqrt(0.5); //disturb f   
                }
            double n = fx*fx+fy*fy; //norm(f)^2
            double f = sqrt(n);

            double xi = sqrt(4.0-3.0*n);
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
            
            
            dcmplx cmplxi(0.0,1.0);
            
            dcmplx t2 = dPxydFx*dPxydFx;
            dcmplx t3 = dPxydE*dPxydE;
            dcmplx t4 = dPxydE*dPxxdFx*dPxydFx;
            dcmplx t27 = dPxxdE*t2;
            dcmplx t5 = t3+t4-t27;
            dcmplx t6 = 1.0/t5;
            dcmplx t7 = dPxxdFx+dPxydFy;
            dcmplx t8 = dPxydE*dPxxdFy*(1.0/2.0);
            dcmplx t9 = t7*t7;
            dcmplx t16 = t7*t9*(1.0/2.7E1);
            dcmplx t17 = dPxxdFx*dPxydFy;
            dcmplx t18 = dPxxdFy*dPxydFx;
            dcmplx t19 = dPxxdE-t17+t18;
            dcmplx t20 = t7*t19*(1.0/6.0);
            dcmplx t21 = dPxxdE*dPxydFy*(1.0/2.0);
            dcmplx t10 = t8+t16+t20-t21;
            dcmplx t11 = dPxxdE*(1.0/3.0);
            dcmplx t12 = dPxxdFy*dPxydFx*(1.0/3.0);
            dcmplx t13 = t9*(1.0/9.0);
            dcmplx t23 = dPxxdFx*dPxydFy*(1.0/3.0);
            dcmplx t14 = t11+t12+t13-t23;
            dcmplx t15 = t14*t14;
            dcmplx t22 = t10*t10;
            dcmplx t30 = t14*t15;
            dcmplx t24 = t22-t30;
            dcmplx t25 = sqrt((t24));
            dcmplx t28 = dPxxdFx*(1.0/3.0);
            dcmplx t29 = dPxydFy*(1.0/3.0);
            dcmplx t31 = t8+t16+t20-t21+t25;
            dcmplx t32 = pow(t31,1.0/3.0);
            dcmplx t26 = t28+t29+t32+t14*1.0/pow(t8+t16+t20+t25-dPxxdE*dPxydFy*(1.0/2.0),1.0/3.0);
            dcmplx t33 = dPxydE*dPxydFy;
            dcmplx t34 = dPxxdFx*dPxydFx*dPxydFy;
            dcmplx t43 = dPxxdFy*t2;
            dcmplx t35 = t33+t34-t43;
            dcmplx t36 = dPxxdFx*dPxydFx;
            dcmplx t37 = dPxydFx*dPxydFy;
            dcmplx t38 = dPxydE+t36+t37;
            dcmplx t39 = 1.0/pow(t31,1.0/3.0);
            dcmplx t40 = t14*t39;
            dcmplx t41 = sqrt(3.0);
            dcmplx t44 = t41*(t32-t40)*-5.0E-1*cmplxi;
            dcmplx t45 = t14*t39*(1.0/2.0);
            dcmplx t46 = t32*(1.0/2.0);
            dcmplx t42 = t28+t29-t44-t45-t46;
            dcmplx t47 = t28+t29+t44-t45-t46;
            dcmplx t48 = t28+t29+t32+t40;
            dcmplx t49 = dPxxdE*dPxydFy;
            dcmplx t55 = dPxydE*dPxxdFy;
            dcmplx t50 = t49-t55;
            dcmplx t51 = dPxydFx*t6*t50;
            dcmplx t52 = dPxxdE*dPxydFx;
            dcmplx t53 = t33+t52;
            dcmplx t54 = -t28-t29+t44+t45+t46;
            dcmplx t56 = t47*t47;
            double Vx11 = real(-t6*t35-dPxydFx*t6*(t26*t26)+t6*t38*t48);
            double Vx12 = real(-t6*t35-dPxydFx*t6*(t42*t42)+t6*t38*(t28+t29-t32*(1.0/2.0)-t14*t39*(1.0/2.0)+t41*(t32-t40)*5.0E-1*cmplxi));
            double Vx13 = real(-t6*t35-dPxydFx*t6*t56+t6*t38*t47);
            double Vx21 = real(t51+dPxydE*t6*(t48*t48)-t6*t48*t53);
            double Vx22 = real(t51+dPxydE*t6*(t54*t54)+t6*t53*t54);
            double Vx23 = real(t51+dPxydE*t6*t56-t6*t47*t53);
            double Vx31 = 1.0;
            double Vx32 = 1.0;
            double Vx33 = 1.0;
           
            

            double n1 = sqrt(Vx11*Vx11+Vx21*Vx21+Vx31*Vx31);
            double n2 = sqrt(Vx12*Vx12+Vx22*Vx22+Vx32*Vx32);
            double n3 = sqrt(Vx13*Vx13+Vx23*Vx23+Vx33*Vx33);
            

            
            if (n>1.0e-14)
            {
                Vx11 = Vx11/n1;
                Vx21 = Vx21/n1;
                Vx31 = Vx31/n1;
                
                Vx12 = Vx12/n2;
                Vx22 = Vx22/n2;
                Vx32 = Vx32/n2;
                
                Vx13 = Vx13/n3;
                Vx23 = Vx23/n3;
                Vx33 = Vx33/n3;
                
            }
            else
            {
                Vx11 = 0.0;
                Vx12 = -sqrt(3.0);
                Vx13 = sqrt(3.0);
                Vx21 = 0.0;
                Vx22 = 1.0;
                Vx23 = 1.0;
                Vx31 = 1.0;
                Vx32 = 0.0;
                Vx33 = 0.0;
                
            }
            

            
            
            t2 = dPxydFx+dPyydFy;
            t3 = dPxydE*dPyydFx*(1.0/2.0);
            t4 = t2*t2;
            t11 = t2*t4*(1.0/2.7E1);
            t12 = dPxydFx*dPyydFy;
            t13 = dPxydFy*dPyydFx;
            t14 = dPyydE-t12+t13;
            t15 = t2*t14*(1.0/6.0);
            t16 = dPxydFx*dPyydE*(1.0/2.0);
            t5 = t3+t11+t15-t16;
            t6 = dPyydE*(1.0/3.0);
            t7 = dPxydFy*dPyydFx*(1.0/3.0);
            t8 = t4*(1.0/9.0);
            t18 = dPxydFx*dPyydFy*(1.0/3.0);
            t9 = t6+t7+t8-t18;
            t10 = t9*t9;
            t17 = t5*t5;
            t28 = t9*t10;
            t19 = t17-t28;
            t20 = sqrt(t19);
            t26 = dPxydFx*(1.0/3.0);
            t27 = dPyydFy*(1.0/3.0);
            t29 = t3+t11+t15-t16+t20;
            t30 = pow(t29,1.0/3.0);
            t21 = t26+t27+t30+t9*1.0/pow(t3+t11+t15+t20-dPxydFx*dPyydE*(1.0/2.0),1.0/3.0);
            t22 = dPxydE*dPyydFx;
            t25 = dPxydFx*dPyydE;
            t23 = t22-t25;
            t24 = 1.0/t23;
            t31 = 1.0/pow(t29,1.0/3.0);
            t32 = t9*t31;
            t34 = sqrt(3.0);
            t35 = t34*(t30-t32)*-5.0E-1*cmplxi;
            t36 = t9*t31*(1.0/2.0);
            t37 = t30*(1.0/2.0);
            t33 = t26+t27-t35-t36-t37;
            t38 = t26+t27+t35-t36-t37;
            t39 = 1.0/dPyydFx;
            t40 = t26+t27+t30+t32;
            t41 = dPyydE*dPyydE;
            t42 = dPxydFy*dPyydE*dPyydFx;
            t48 = dPxydE*dPyydFx*dPyydFy;
            t43 = t41+t42-t48;
            t44 = t24*t39*t43;
            t45 = dPyydE*dPyydFy;
            t46 = t22+t45;
            t47 = -t26-t27+t35+t36+t37;
            t49 = t38*t38;
            double Vy11 = real(-t14*t24+(t21*t21)*t24-t2*t24*t40);
            double Vy12 = real(-t14*t24+t24*(t33*t33)+t2*t24*t47);
            double Vy13 = real(-t14*t24+t24*t49-t2*t24*t38);
            double Vy21 = real(t44+t24*t39*t40*t46-dPyydE*t24*t39*(t40*t40));
            double Vy22 = real(t44-t24*t39*t46*t47-dPyydE*t24*t39*(t47*t47));
            double Vy23 = real(t44-dPyydE*t24*t39*t49+t24*t38*t39*t46);
            double Vy31 = 1.0;
            double Vy32 = 1.0;
            double Vy33 = 1.0;

            

            
            
            n1 = sqrt(Vy11*Vy11+Vy21*Vy21+Vy31*Vy31);
            n2 = sqrt(Vy12*Vy12+Vy22*Vy22+Vy32*Vy32);
            n3 = sqrt(Vy13*Vy13+Vy23*Vy23+Vy33*Vy33);
            

            
            if (n>1.0e-14)
            {
                Vy11 = Vy11/n1;
                Vy21 = Vy21/n1;
                Vy31 = Vy31/n1;
                
                Vy12 = Vy12/n2;
                Vy22 = Vy22/n2;
                Vy32 = Vy32/n2;
                
                Vy13 = Vy13/n3;
                Vy23 = Vy23/n3;
                Vy33 = Vy33/n3;
                
            }
            else
            {
                Vy11 = 0.0;
                Vy12 = -sqrt(3.0);
                Vy13 = sqrt(3.0);
                Vy21 = 1.0;
                Vy22 = 0.0;
                Vy23 = 0.0;
                Vy31 = 0.0;
                Vy32 = 1.0;
                Vy33 = 1.0;
                
            }
            
                    
            double Vxinv11 = (Vx22*Vx33-Vx23*Vx32)/(Vx11*(Vx22*Vx33-Vx23*Vx32)-Vx12*Vx21*Vx33+Vx12*Vx23*Vx31+Vx13*Vx21*Vx32-Vx13*Vx22*Vx31);
            double Vxinv12 = -(Vx12*Vx33-Vx13*Vx32)/(Vx11*(Vx22*Vx33-Vx23*Vx32)-Vx12*Vx21*Vx33+Vx12*Vx23*Vx31+Vx13*Vx21*Vx32-Vx13*Vx22*Vx31);
            double Vxinv13 = (Vx12*Vx23-Vx13*Vx22)/(Vx11*(Vx22*Vx33-Vx23*Vx32)-Vx12*Vx21*Vx33+Vx12*Vx23*Vx31+Vx13*Vx21*Vx32-Vx13*Vx22*Vx31);
            double Vxinv21 = -(Vx21*Vx33-Vx23*Vx31)/(Vx11*(Vx22*Vx33-Vx23*Vx32)-Vx12*Vx21*Vx33+Vx12*Vx23*Vx31+Vx13*Vx21*Vx32-Vx13*Vx22*Vx31);
            double Vxinv22 = (Vx11*Vx33-Vx13*Vx31)/(Vx11*(Vx22*Vx33-Vx23*Vx32)-Vx12*Vx21*Vx33+Vx12*Vx23*Vx31+Vx13*Vx21*Vx32-Vx13*Vx22*Vx31);
            double Vxinv23 = -(Vx11*Vx23-Vx13*Vx21)/(Vx11*(Vx22*Vx33-Vx23*Vx32)-Vx12*Vx21*Vx33+Vx12*Vx23*Vx31+Vx13*Vx21*Vx32-Vx13*Vx22*Vx31);
            double Vxinv31 = (Vx21*Vx32-Vx22*Vx31)/(Vx11*(Vx22*Vx33-Vx23*Vx32)-Vx12*Vx21*Vx33+Vx12*Vx23*Vx31+Vx13*Vx21*Vx32-Vx13*Vx22*Vx31);
            double Vxinv32 = -(Vx11*Vx32-Vx12*Vx31)/(Vx11*(Vx22*Vx33-Vx23*Vx32)-Vx12*Vx21*Vx33+Vx12*Vx23*Vx31+Vx13*Vx21*Vx32-Vx13*Vx22*Vx31);
            double Vxinv33 = (Vx11*Vx22-Vx12*Vx21)/(Vx11*(Vx22*Vx33-Vx23*Vx32)-Vx12*Vx21*Vx33+Vx12*Vx23*Vx31+Vx13*Vx21*Vx32-Vx13*Vx22*Vx31);

            
            double Vyinv11 = (Vy22*Vy33-Vy23*Vy32)/(Vy11*(Vy22*Vy33-Vy23*Vy32)-Vy12*Vy21*Vy33+Vy12*Vy23*Vy31+Vy13*Vy21*Vy32-Vy13*Vy22*Vy31);
            double Vyinv12 = -(Vy12*Vy33-Vy13*Vy32)/(Vy11*(Vy22*Vy33-Vy23*Vy32)-Vy12*Vy21*Vy33+Vy12*Vy23*Vy31+Vy13*Vy21*Vy32-Vy13*Vy22*Vy31);
            double Vyinv13 = (Vy12*Vy23-Vy13*Vy22)/(Vy11*(Vy22*Vy33-Vy23*Vy32)-Vy12*Vy21*Vy33+Vy12*Vy23*Vy31+Vy13*Vy21*Vy32-Vy13*Vy22*Vy31);
            double Vyinv21 = -(Vy21*Vy33-Vy23*Vy31)/(Vy11*(Vy22*Vy33-Vy23*Vy32)-Vy12*Vy21*Vy33+Vy12*Vy23*Vy31+Vy13*Vy21*Vy32-Vy13*Vy22*Vy31);
            double Vyinv22 = (Vy11*Vy33-Vy13*Vy31)/(Vy11*(Vy22*Vy33-Vy23*Vy32)-Vy12*Vy21*Vy33+Vy12*Vy23*Vy31+Vy13*Vy21*Vy32-Vy13*Vy22*Vy31);
            double Vyinv23 = -(Vy11*Vy23-Vy13*Vy21)/(Vy11*(Vy22*Vy33-Vy23*Vy32)-Vy12*Vy21*Vy33+Vy12*Vy23*Vy31+Vy13*Vy21*Vy32-Vy13*Vy22*Vy31);
            double Vyinv31 = (Vy21*Vy32-Vy22*Vy31)/(Vy11*(Vy22*Vy33-Vy23*Vy32)-Vy12*Vy21*Vy33+Vy12*Vy23*Vy31+Vy13*Vy21*Vy32-Vy13*Vy22*Vy31);
            double Vyinv32 = -(Vy11*Vy32-Vy12*Vy31)/(Vy11*(Vy22*Vy33-Vy23*Vy32)-Vy12*Vy21*Vy33+Vy12*Vy23*Vy31+Vy13*Vy21*Vy32-Vy13*Vy22*Vy31);
            double Vyinv33 = (Vy11*Vy22-Vy12*Vy21)/(Vy11*(Vy22*Vy33-Vy23*Vy32)-Vy12*Vy21*Vy33+Vy12*Vy23*Vy31+Vy13*Vy21*Vy32-Vy13*Vy22*Vy31);
