//
//  HW_iso.h
//  
//
//  Created by katy ghantous on 5/6/15.
//  Copyright 2015 __MyCompanyName__. All rights reserved.
//

#ifndef _HW_aniso_h
#define _HW_aniso_h



#include <iostream>
#include <complex>
#include <stdio.h>
#include<vector> 
#include<math.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_matrix.h>
using namespace std;




class V{
    
public:
    vector<double> kn;
    vector<double> phi_an,phi_bn,phi_cn, n_an, n_bn,n_cn;
    double g, k_0;
    double apy,bpy, cpy, ap0,bp0x,cp0x,bp0y,cp0y,apx,bpx,cpx;
    double any,bny, cny, an0,bn0x,cn0x,bn0y,cn0y,anx,bnx,cnx;
    vector<int> indx;
    
    double q;
    double alphbar,nuZF;
    int Nshls;
    int ZF;
    int ipYB,ipYE,ip0B,ip0E,ipXB,ipXE; 
    int inYB,inYE,in0B,in0E,inXB,inXE;
    
    void setV(int N, double gIn, double k_0i, double qIn,double alphbarIn,double nuZFIn,int ZFi);
    void set_alph_HW(double alph0, double alphX, double alphY,double alph0_n, double alphX_n, double alphY_n);
    void setCoef_HW();
    
    V();
    ~V();
    
};



int func_hmdsi(double t, const double y[], double dydt[], void *params);

double get_data_NAME(FILE *someFile, const char * nameData);


#endif
