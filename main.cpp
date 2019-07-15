/*this program minimizes a function (which can be defined in the function 
 * myfunc) within a cube. If you change the function you
 * may need to change the vales of the initial "t", "alpha", and "beta". Without
 * proper adjustment this code may not give the proper minima.
 * This program, as there are no methods included to account for local minima 
 * can get stuck on a local minima without reaching the global minima.
 * This program also only goes to the edges of the square, not along them, 
 * although I have ideas of how to implement that in the future. We use Newtons
 * method to minimize the function. Also due to the constraints of Newton's 
 * method this code can get stuck at inflection points
 */


#include <cstdlib>
#include <iostream>
#include <cmath>


using namespace std;

//we create our function, gradient, and constraints (defined below)
double myfunc(double x, double y, double z);
double gradx(double x, double y, double z, double h);
double grady(double x, double y, double z, double h);
double gradz(double x, double y, double z, double h);
double H_xx(double x, double y, double z, double h);
double H_xy(double x, double y, double z, double h);
double H_xz(double x, double y, double z, double h);
double H_yy(double x, double y, double z, double h);
double H_yz(double x, double y, double z, double h);
double H_zz(double x, double y, double z, double h);
double H_inv_11(double x, double y, double z, double h);
double H_inv_12(double x, double y, double z, double h);
double H_inv_13(double x, double y, double z, double h);
double H_inv_22(double x, double y, double z, double h);
double H_inv_23(double x, double y, double z, double h);
double H_inv_33(double x, double y, double z, double h);
double H_det(double x, double y, double z, double h);
double step_x(double x, double y, double z, double h);
double step_y(double x, double y, double z, double h);
double step_z(double x, double y, double z, double h);

double constraints(double x, double y, double z);


int main(int argc, char** argv) {
    
    //first we define our starting point (depending on the function be careful 
    //to make sure the gradient isnt zero at the starting point
    double loc_i[3];
    double loc_f[3] = {0.25,0.25,0.25};
    
    //We initialize our gradient
    double grad[3];
    double del[3];
    double h = 0.0000001;
    
    //create variables for our backtracking line search
    double alpha = 0.25;
    double beta = 0.5;
    double t=1;
    double f_step;
    double f_orig;
    double f_bls_diff;
    double grad_del;
    
    //and we create variables for our stopping conditions
    double constrain_break = 0;
    double eps = 0.0000001;
    double f_diff = 1;
    
    constrain_break = constraints(loc_f[0],loc_f[1],loc_f[2]); 
    if(constrain_break == 1){
        cout << "error: initial point is outside the cube" << endl;
    }
    
    //Then we create the while loop to implement the gradient descent technique
    while(f_diff >= eps && constrain_break == 0){
        //first we make our initial location the final location from the last step
        for(int i=0; i<3; i++){
            loc_i[i] = loc_f[i];
        }
        
        //Then we find our step directions by finding the gradient
        grad[0] = gradx(loc_i[0], loc_i[1], loc_i[2], h);
        grad[1] = grady(loc_i[0], loc_i[1], loc_i[2], h);
        grad[2] = gradz(loc_i[0], loc_i[1], loc_i[2], h);
        
        del[0] = step_x(loc_i[0], loc_i[1], loc_i[2], h);
        del[1] = step_y(loc_i[0], loc_i[1], loc_i[2], h);
        del[2] = step_z(loc_i[0], loc_i[1], loc_i[2], h);
        
        /*
        del[0] = del[0]/(pow(del[0],2)+ pow(del[1],2) + pow(del[2],2));
        del[1] = del[1]/(pow(del[0],2)+ pow(del[1],2) + pow(del[2],2));
        del[2] = del[2]/(pow(del[0],2)+ pow(del[1],2) + pow(del[2],2));
        */
        
        
        //then we implement the backtracking line search 
        t=0.01;
        f_step = myfunc(loc_i[0] + t*del[0], loc_i[1] + t*del[1], loc_i[2] + t*del[2]);
        f_orig = myfunc(loc_i[0], loc_i[1], loc_i[2]);
        grad_del = 0;
        for(int i=0; i<3; i++){
            grad_del += grad[i]*del[i];
        }
        f_bls_diff = f_step - f_orig - alpha*t*grad_del;
        
        
        while(f_bls_diff > 0){
            t *= beta;
            f_step = myfunc(loc_i[0] + t*del[0], loc_i[1] + t*del[1], loc_i[2] + t*del[2]);
            f_orig = myfunc(loc_i[0], loc_i[1], loc_i[2]);
            grad_del = 0;
            
            for(int i=0; i<3; i++){
                grad_del += grad[i]*del[i];
            }
            f_bls_diff = f_step - f_orig - alpha*t*grad_del;
            
        }
        //and we advance the step
        for(int i=0; i<3; i++){
            loc_f[i] = loc_i[i]+ t*del[i];
        }  
        f_diff = abs(myfunc(loc_i[0], loc_i[1], loc_i[2]) - myfunc(loc_f[0], loc_f[1], loc_f[2]));
        constrain_break = constraints(loc_f[0],loc_f[1],loc_f[2]);          
    }
    
    cout << "The value of the minimized function is f=" << f_orig << endl;
    cout << "the location is x=" << loc_i[0] << "," << "\t" << "y=" << loc_i[1] << "," << "\t" << "z=" << loc_i[2] << endl;
    
    return 0;
}


//we define our functions
double myfunc(double x, double y, double z){
    double f = pow(x,2)+ pow(y,2) + pow(z,2);    
    return f;
}
double gradx(double x, double y, double z, double h){
    double x_deriv = (myfunc(x+h,y,z) - myfunc(x,y,z))/h;
    return x_deriv;
}
double grady(double x, double y, double z, double h){
    double y_deriv = (myfunc(x,y+h,z) - myfunc(x,y,z))/h;
    return y_deriv;
}
double gradz(double x, double y, double z, double h){
    double z_deriv = (myfunc(x,y,z+h) - myfunc(x,y,z))/h;
    return z_deriv;
}
double H_xx(double x, double y, double z, double h){
    double xx_deriv = (gradx(x+h,y,z,h) - gradx(x,y,z,h))/h;
    return xx_deriv;
}
double H_xy(double x, double y, double z, double h){
    double xy_deriv = (gradx(x,y+h,z,h) - gradx(x,y,z,h))/h;
    return xy_deriv;
}
double H_xz(double x, double y, double z, double h){
    double xz_deriv = (gradx(x,y,z+h,h) - gradx(x,y,z,h))/h;
    return xz_deriv;
}
double H_yy(double x, double y, double z, double h){
    double yy_deriv = (grady(x,y+h,z,h) - grady(x,y,z,h))/h;
    return yy_deriv;
}
double H_yz(double x, double y, double z, double h){
    double yz_deriv = (grady(x,y,z+h,h) - grady(x,y,z,h))/h;
    return yz_deriv;
}
double H_zz(double x, double y, double z, double h){
    double zz_deriv = (gradz(x,y,z+h,h) - gradz(x,y,z,h))/h;
    return zz_deriv;
}
double H_det(double x, double y, double z, double h){
    double h_det = -1.*H_xx(x,y,z,h)*H_yy(x,y,z,h)*H_zz(x,y,z,h) + H_xx(x,y,z,h)*H_yz(x,y,z,h)*H_yz(x,y,z,h) + H_xy(x,y,z,h)*H_xy(x,y,z,h)*H_zz(x,y,z,h) - 2.*H_xy(x,y,z,h)*H_xz(x,y,z,h)*H_yz(x,y,z,h) + H_xz(x,y,z,h)*H_xz(x,y,z,h)*H_yy(x,y,z,h); 
    return h_det;
}
double H_inv_11(double x, double y, double z, double h){
    double h_inv_11 = (1./H_det(x,y,z,h))*(H_yz(x,y,z,h)*H_yz(x,y,z,h)-H_yy(x,y,z,h)*H_zz(x,y,z,h));
    return h_inv_11;
}
double H_inv_12(double x, double y, double z, double h){
    double h_inv_12 = (1./H_det(x,y,z,h))*(H_xy(x,y,z,h)*H_zz(x,y,z,h)-H_xz(x,y,z,h)*H_yz(x,y,z,h));
    return h_inv_12;
}
double H_inv_13(double x, double y, double z, double h){
    double h_inv_13 = (1./H_det(x,y,z,h))*(H_xz(x,y,z,h)*H_yy(x,y,z,h)-H_xy(x,y,z,h)*H_yz(x,y,z,h));
    return h_inv_13;
}
double H_inv_22(double x, double y, double z, double h){
    double h_inv_22 = (1./H_det(x,y,z,h))*(H_xz(x,y,z,h)*H_xz(x,y,z,h)-H_xx(x,y,z,h)*H_zz(x,y,z,h));
    return h_inv_22;
}
double H_inv_23(double x, double y, double z, double h){
    double h_inv_23 = (1./H_det(x,y,z,h))*(H_xx(x,y,z,h)*H_yz(x,y,z,h)-H_xy(x,y,z,h)*H_xz(x,y,z,h));
    return h_inv_23;
}
double H_inv_33(double x, double y, double z, double h){
    double h_inv_33 = (1./H_det(x,y,z,h))*(H_xy(x,y,z,h)*H_xy(x,y,z,h)-H_xx(x,y,z,h)*H_yy(x,y,z,h));
    return h_inv_33;
}
double step_x(double x, double y, double z, double h){
    double x_step = -1.*(H_inv_11(x,y,z,h)*gradx(x,y,z,h) + H_inv_12(x,y,z,h)*grady(x,y,z,h) + H_inv_13(x,y,z,h)*gradz(x,y,z,h));
    return x_step;
}
double step_y(double x, double y, double z, double h){
    double y_step = -1.*(H_inv_12(x,y,z,h)*gradx(x,y,z,h) + H_inv_22(x,y,z,h)*grady(x,y,z,h) + H_inv_23(x,y,z,h)*gradz(x,y,z,h));
    return y_step;
}
double step_z(double x, double y, double z, double h){
    double z_step = -1.*(H_inv_13(x,y,z,h)*gradx(x,y,z,h) + H_inv_23(x,y,z,h)*grady(x,y,z,h) + H_inv_33(x,y,z,h)*gradz(x,y,z,h));
    return z_step;
}
double constraints(double x, double y, double z){
    double satisfied_count = 0;
    double sat_or_not;
    if(x < 1){
        satisfied_count++;
    }
    if(x > -1){
        satisfied_count++;
    }
    if(y < 1){
        satisfied_count++;
    }
    if(y > -1){
        satisfied_count++;
    }
    if(z < 1){
        satisfied_count++;
    }
    if(z > -1){
        satisfied_count++;
    }
    if(satisfied_count == 6){
        sat_or_not = 0;
    }
    else {
        sat_or_not = 1;
    }
    return sat_or_not;
}