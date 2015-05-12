// to run:
//    g++  -L/usr/include -I/usr/include -L/usr/lib test.cpp   -lgsl -lgslcblas -lm  -o  test
// exact parameters as ExactIso 

#include <iostream>
#include <complex>
#include <stdio.h>
#include<vector> 
#include<math.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_matrix.h>
#include "HW_aniso.h"
using namespace std;


int main(int argc,char **args)
{
    cout<<"This is the anisotropic version of the shell model for HW";
    
    /****************/
    /* Declare      */
    /****************/
    int i,N,N_V;
    double ink, ing;
    double qIn;
    double alphbarIn,nuZFIn;
    int ZF;//0 if there is no zonal flow any other integer if there is.
    //complex *phi;
    int status;
    double Par, alph0, alphX, alphY,alph0_n,  alphX_n, alphY_n;
    complex<double> I (0.0,1.0);
    
    
    
    V V1;
    
    
    FILE *rs,*rsZF,*EFGH;
    clock_t bclc, eclc;
    bclc = clock();
    
    EFGH = fopen("efgh.txt","w");
    rsZF=fopen("rsltsZF.txt","w");
    rs=fopen("rslts.txt","w");
    /****************/
    /* initialize   */
    /****************/
    
    FILE *fInp;
    /*** values from input file ****/
    
    
    double kmax;
    fInp = fopen("INPUT","r");   
    N =  (int) get_data_NAME(fInp,"Number of Shells");
    ink = get_data_NAME(fInp,"k min"); 
    kmax = get_data_NAME(fInp,"k max"); 
    ing = exp(log(kmax/ink)/((double) N-1));
    fclose(fInp);
    
    
    /** System **/
    qIn = 0.01;
    alphbarIn=1.; nuZFIn = 2.e2;
    ZF = 0;
    
    Par = ing*ing;
    alph0= Par/2;   alphX=Par; alphY=Par; 
    alph0_n = Par/2;  alphX_n=Par; alphY_n=Par;
    
    N_V = 2+6*N;
    
    
    V1.setV(N, ing, ink, qIn, alphbarIn,nuZFIn,ZF);
    
    V1.setCoef_HW();
    
    V1.set_alph_HW(alph0,  alphX, alphY, alph0_n,alphX_n, alphY_n);
    //////////////////////////////////
    ////////////////////////////////////
    
    /****************/
    /*    GSL       */ 
    /****************/
    
    
    /** Integrator */
    
    
    int dimension = sizeof(complex<double>)*N_V;
    
    double eps_abs = 1.0e-21;	// absolute error requested 
    double eps_rel = 1.e-7;	// relative error requested 
    
    // define the type of routine for making steps: 
    const gsl_odeiv2_step_type *type_ptr = gsl_odeiv2_step_rk8pd;
    
    
    gsl_odeiv2_step *step_ptr = gsl_odeiv2_step_alloc (type_ptr, dimension);
    gsl_odeiv2_control *control_ptr  = gsl_odeiv2_control_standard_new (eps_abs, eps_rel, 1.0,0.0);
    gsl_odeiv2_evolve *evolve_ptr = gsl_odeiv2_evolve_alloc (dimension);
    
    gsl_odeiv2_system my_system;
    
    double t, t_next,delta_t;		// current and next independent variable 
    double tmin, tmax; // range of t and step size for output 
    
    
    double h = 1.e-12;		// starting step size for ode solver 
    
    complex<double> *y0;
    int szy;
    
    
    szy = sizeof(complex<double>)*N_V;
    double y[szy]; 
    
    y0 = (complex<double>*)&y[0];
    for(i=0;i<N_V;i++){
        y0[i] = 1.e-10*exp(rand()*2.*M_PI*I);
    }
    
    
    // load values into the my_system structure 
    my_system.function = func_hmdsi;	// the right-hand-side functions dy[i]/dt 
    my_system.jacobian = NULL;	// the Jacobian df[i]/dy[j] 
    my_system.dimension = dimension;	// number of diffeq's 
    my_system.params = &V1;	// parameters to pass to rhs and jacobian 
    
    tmin = 0;
    
    fInp = fopen("INPUT","r");   
    tmax = get_data_NAME(fInp,"t max");
    fclose(fInp);
    
    
    
    
    
    fInp = fopen("INPUT","r");   
    delta_t = get_data_NAME(fInp,"dt");
    fclose(fInp);
    
    
    
    
    t=tmin;
    complex<double> *phiZF = (complex<double> *)&y[0];
    
    complex<double> *phiY = (complex<double> *)&y[2*V1.ipYB];
    complex<double> *phi0 = (complex<double> *)&y[2*V1.ip0B]; 
    complex<double> *phiX = (complex<double> *)&y[2*V1.ipXB];
    
    complex<double> *nY = (complex<double> *)&y[2*V1.inYB];
    complex<double> *n0 = (complex<double> *)&y[2*V1.in0B];
    complex<double> *nX = (complex<double> *)&y[2*V1.inXB];
    
    complex<double> *nZF = (complex<double> *)&y[2*V1.inXE+2];
    
    
    for (t_next = tmin + delta_t; t_next <= tmax; t_next += delta_t)
    {
        while (t < t_next)	// evolve from t to t_next 
        {
            status = gsl_odeiv2_evolve_apply (evolve_ptr, control_ptr, step_ptr, &my_system, &t, t_next, &h, y);
            
            if (status != GSL_SUCCESS)
                break;
        }
        printf("%0.2f\n",t);
        std::scientific;
        std::cout <<  abs(phiZF[0])*abs(phiZF[0])<<"\t"<< abs(nZF[0])*abs(nZF[0]) <<"\n";
        fprintf(rsZF,"%e \t %e \n",abs(phiZF[0])*abs(phiZF[0]), abs(nZF[0])*abs(nZF[0]));
        
        for(i=0;i<N;i++){
            std::scientific;
            std::cout << V1.kn[i]<<"\t"<< abs(phiY[i])*abs(phiY[i]) <<"\t"<< abs(phi0[i])*abs(phi0[i])<<"\t"<< abs(phiX[i])*abs(phiX[i])<<"\t"<< abs(nY[i])*abs(nY[i])<<"\t"<< abs(n0[i])*abs(n0[i])<<"\t"<< abs(nX[i])*abs(nX[i])<< '\n';
            fprintf(rs,"%e \t %e\t %e \t %e \t %e \t %e \t %e \n", V1.kn[i], abs(phiY[i])*abs(phiY[i]),abs(phi0[i])*abs(phi0[i]), abs(phiX[i])*abs(phiX[i]), abs(nY[i])*abs(nY[i]), abs(n0[i])*abs(n0[i]), abs(nX[i])*abs(nX[i]));
            fprintf(EFGH,"%e \t %e\t %e \t %e \t %e \t %e \t %e \t %e\t %e \t %e \t %e \t %e \t %e \n", V1.kn[i], abs(phiY[i])*abs(phiY[i])*V1.kn[i],abs(phi0[i])*abs(phi0[i])*V1.kn[i], abs(phiX[i])*abs(phiX[i])*V1.kn[i], abs(nY[i])*abs(nY[i])/V1.kn[i], abs(n0[i])*abs(n0[i])/V1.kn[i], abs(nX[i])*abs(nX[i])/V1.kn[i],  imag(conj(phiY[i])*nY[i]),imag(conj(phi0[i])*n0[i]), imag(conj(phiX[i])*nX[i]),real(conj(phiY[i])*nY[i]),real(conj(phi0[i])*n0[i]), real(conj(phiX[i])*nX[i]));
            
        }
        
    } 
    
    
    fclose(rs);
    fclose(rsZF);
    fclose(EFGH);
    
    // all done; free up the gsl_odeiv stuff
    gsl_odeiv2_evolve_free (evolve_ptr);
    gsl_odeiv2_control_free (control_ptr);
    gsl_odeiv2_step_free (step_ptr);
    
    
    eclc = clock();
    return 0; 
    
    
}


