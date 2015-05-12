//
//  HW_iso.cpp
//  
//
//  Created by katy ghantous on 5/6/15.
//  Copyright 2015 __MyCompanyName__. All rights reserved.
//

#include <iostream>
#include <complex>
#include <stdio.h>
#include<vector> 
#include<math.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_matrix.h>
#include "HW_aniso.h"
using namespace std;



V::V(){};



void V::setV(int N, double gIn, double k_0i,double qIn,double alphbarIn,double nuZFIn,int ZFi){
    
    int i;
    double n; 
    
    
    q= qIn;
    alphbar= alphbarIn ;nuZF = nuZFIn;
    ZF = ZFi;
    
    
    ipYB = 1; ipYE = N; 
    ip0B = N+1; ip0E = 2*N;
    ipXB = 2*N+1; ipXE = 3*N; 
    
    inYB = 3*N+1; inYE = 4*N;
    in0B = 4*N+1; in0E = 5*N;
    inXB = 5*N+1; inXE = 6*N;
    
    Nshls = N;
    g = gIn;
    k_0 = k_0i;
    
    
    
    for(i=0;i<N;i++){
        n = (double) i ;
        kn.push_back(pow(g,n)*k_0);
    }    
    
    indx.push_back(-1);//phi zonal flow
    for(i=0;i<6*N;i++){
        indx.push_back((int) i%N); 
    } 
    indx.push_back(-1); //n zonal flow
    
    
};


void V::set_alph_HW(double alph0, double alphX, double alphY,double alph0_n, double alphX_n, double alphY_n){
    
    
    apy = alphY; bpy = -2*alph0/g; cpy = 2*alph0/g/g;
    ap0 = alph0; bp0x = -alphX/2/g; cp0x = alphX/2/g/g; bp0y = -alphY/2/g; cp0y = alphY/2/g/g;
    apx = alphX; bpx = -2*alph0/g; cpx = 2*alph0/g/g;
    
    any = alphY_n; bny = -2*alph0_n/g; cny = 2*alph0_n/g/g;
    an0 = alph0_n; bn0x = -alphX_n/2/g; cn0x = alphX_n/2/g/g; bn0y = -alphY_n/2/g; cn0y = alphY_n/2/g/g;
    anx = alphX_n; bnx = -2*alph0_n/g; cnx = 2*alph0_n/g/g;
    
    
};

void V::setCoef_HW(){
    int i;
    for(i=0;i<Nshls;i++){
        
        phi_an.push_back(pow(kn[i],2)*(g*g-1)/pow(g,7));
        phi_bn.push_back(pow(kn[i],2)*(g*g*g*g-1)/pow(g,2));
        phi_cn.push_back(pow(kn[i],2)*(g*g-1)*pow(g,5));
        
        n_an.push_back(pow(kn[i],2)/pow(g,3));
        n_bn.push_back(pow(kn[i],2));
        n_cn.push_back(pow(kn[i],2)*pow(g,3));
        
    }
    
    
};


V::~V(){};



int func_hmdsi(double t, const double y[], double dydt[], void *params){
    
    
    V *VIn = (V*)params;
    FILE * fInp;
    
    
    
    complex<double> *dphidt=(complex<double> *)&dydt[0];
    
    
    
    int tid,nthreads;
    int i,n,N_V,fi,N,N_trs;
    complex<double> C_phi_phi, C_phi_n,Sm,disp,Z_fl,frcng;
    double musm,mubg,Fn,Fp,C2,vn;
    double musmFac,mubgFac;   
    
    complex<double> I (0.0,1.0);
    
    complex<double> SpZF (0.0,0.0);
    complex<double> SnZF (0.0,0.0);
    
    //-------------------------------------------------------
    //                  Setting up the equations    
    //-------------------------------------------------------
    
    
    /////   values from input file 
    
    fInp = fopen("INPUT","r");   
    N_trs = get_data_NAME(fInp,"Num of threads");  
    fclose(fInp);
    
    
    
    fInp = fopen("INPUT","r");   
    C2 = get_data_NAME(fInp,"C");
    vn = get_data_NAME(fInp,"kappa"); 
    Fp = get_data_NAME(fInp,"forcing");       
    fclose(fInp);
    
    fInp = fopen("INPUT","r");   
    musmFac = get_data_NAME(fInp,"musmFac");
    mubgFac = get_data_NAME(fInp,"mubgFac"); 
    fclose(fInp);
    
    
    //////
    
    
    Fn= 0.0;
    musm = musmFac*pow(VIn->kn[0],6); // 1.e-18;//instead of -13.. 10*pow(VIn->kn[0],6);
    mubg = mubgFac*pow(VIn->kn[VIn->Nshls-1],-4); //1.e-24; // instead of -15;//100*pow(VIn->kn[VIn->Nshls-1],-4);

    
    N = VIn->Nshls;
    N_V = 2+6*N;
    fi =(int) (log(1/VIn->kn[0])/log(VIn->kn[1]/VIn->kn[0]));
    
    
    
    complex<double> *phiZF = (complex<double> *)&y[0];
    
    complex<double> *phiY = (complex<double> *)&y[2*VIn->ipYB];
    complex<double> *phi0 = (complex<double> *)&y[2*VIn->ip0B]; 
    complex<double> *phiX = (complex<double> *)&y[2*VIn->ipXB];
    
    complex<double> *nY = (complex<double> *)&y[2*VIn->inYB];
    complex<double> *n0 = (complex<double> *)&y[2*VIn->in0B];
    complex<double> *nX = (complex<double> *)&y[2*VIn->inXB];
    
    complex<double> *nZF = (complex<double> *)&y[2*VIn->inXE+2];
    
    for(n=0;n<N-1;n++){
        SpZF = SpZF + pow(VIn->kn[n],3)*((conj(phiY[n])*conj(phi0[n+1])) + (conj(phi0[n])*conj(phiX[n+1])));
        SnZF = SnZF + VIn->kn[n]*((conj(phiY[n])*conj(n0[n+1]) - conj(phi0[n+1])*conj(nY[n])) + (conj(phi0[n])*conj(nX[n+1]) - conj(phiX[n+1])*conj(n0[n])));
    }
    
    
    
    
    
    /////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////
#pragma omp parallel shared (phiY,phiX,phi0,nX,nY,n0,dphidt,VIn,Sm,C2,vn,fi,N_V,Fp,Fn,musm,mubg) private(tid,i,n,C_phi_phi,C_phi_n,disp,Z_fl,frcng)  num_threads(N_trs)
    {
#pragma omp for 
        for(i=1;i<N_V-1;i++){ 
            dphidt[i]=0;
            
            n = VIn->indx[i]; // 1+i%N
            
            
            //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//
            //                  Potential
            /////////////////////////////////  PHI Y/////////////////////
            
            if((i>=(VIn->ipYB)) && (i<=(VIn->ipYE))){
                
                C_phi_phi = (((n>1) ? VIn->apy*VIn->phi_an[n]*(conj(phi0[n-2])*conj(phi0[n-1]))  :0 )
                             +( (n>0 &&n<(N-1)) ?   VIn->bpy*VIn->phi_bn[n]*(conj(phiX[n-1])*conj(phi0[n+1]))  : 0 )
                             +((n<(N-2)) ?  VIn->cpy*VIn->phi_cn[n]*(conj(phiX[n+1])*conj(phi0[n+2])) : 0)
                             );
                
                
                
                // dissipations 
                disp  =  -(musm*pow(VIn->kn[n],-6)+mubg*pow(VIn->kn[n],4))*phiY[n];
                // forcing terms
                frcng = C2*(phiY[n]-nY[n])/(-pow(VIn->kn[n],2));
                if(n==fi||n==fi+1)
                    frcng=frcng+Fp; 
                
                //zonal part  :alpbr*q/k*(k^2/g^2-q^2/g)*[phibr*ph0_(n-1)]                                                                                                                                        
                Z_fl = (n<(N-1) && VIn->ZF!=0 ?  VIn->alphbar*VIn->q/VIn->kn[n]*(VIn->kn[n]*VIn->kn[n]*(VIn->g*VIn->g) -VIn->q*VIn->q)*(conj(phiZF[0])*conj(phi0[n+1])) : 0); 
                
                dphidt[i] = C_phi_phi + frcng  + disp + Z_fl;
                
            }
            
            /////////////////////////////////  PHI 0 //////////////////////
            if((i>=(VIn->ip0B)) && (i<=(VIn->ip0E))){
                
                C_phi_phi = (((n>1) ? VIn->ap0*VIn->phi_an[n]*(conj(phiX[n-2])*conj(phiY[n-1]) +conj(phiY[n-2])*conj(phiX[n-1]))  :0 )
                             +( (n>0 &&n<(N-1)) ?  VIn->phi_bn[n]*(VIn->bp0x*conj(phi0[n-1])*conj(phiX[n+1]) +VIn->bp0y*conj(phi0[n-1])*conj(phiY[n+1])) : 0 )
                             +((n<(N-2)) ?    VIn->phi_cn[n]*(VIn->cp0x*conj(phi0[n+1])*conj(phiX[n+2]) + VIn->cp0y*conj(phi0[n+1])*conj(phiY[n+2]) )  : 0)
                             );
                
                // dissipations 
                disp  =  -(musm*pow(VIn->kn[n],-6)+mubg*pow(VIn->kn[n],4))*phi0[n];
                // forcing terms
                frcng = C2*(phi0[n]-n0[n])/(-pow(VIn->kn[n],2));
                if(n==fi||n==fi+1)
                    frcng=frcng+Fp; 
                
                //zonal part  alphbr*qk(1+k^2/g^2-q^2)/(2(1+k^2))*[phbr*phY_(n-1)]
                
                Z_fl = (n>0 && n<(N-1) && VIn->ZF!=0 ?  -VIn->alphbar*VIn->q/2.*(VIn->kn[n]*VIn->kn[n]/(VIn->g*VIn->g) -VIn->q*VIn->q)/(VIn->kn[n]*VIn->g)*(conj(phiZF[0])*conj(phiY[n-1])) + VIn->alphbar*VIn->q/2.*(VIn->kn[n]*VIn->kn[n]*(VIn->g*VIn->g) -VIn->q*VIn->q)/(VIn->kn[n])*(conj(phiZF[0])*conj(phiX[n+1])) : 0); 
                
                
                
                dphidt[i] = C_phi_phi + frcng  + disp + Z_fl;
                
            }
            
            
            /////////////////////////////////  PHI X //////////////////////
            if((i>=(VIn->ipXB)) && (i<=(VIn->ipXE))){
                
                C_phi_phi = (((n>1) ?  VIn->apx*VIn->phi_an[n]*(conj(phi0[n-2])*conj(phi0[n-1]))  :0 )
                             +( (n>0 &&n<(N-1)) ?  VIn->bpx*VIn->phi_bn[n]*(conj(phiY[n-1])*conj(phi0[n+1]))  : 0 )
                             +((n<(N-2)) ?   VIn->cpx*VIn->phi_cn[n]*(conj(phiY[n+1])*conj(phi0[n+2]))  : 0)
                             );
                
                
                // dissipations 
                disp  =  -(musm*pow(VIn->kn[n],-6)+mubg*pow(VIn->kn[n],4))*phiX[n];
                // forcing terms
                frcng = C2*(phiX[n]-nX[n])/(-pow(VIn->kn[n],2));
                if(n==fi||n==fi+1)
                    frcng=frcng+Fp; 
                //zonal part
                
                Z_fl = (n>0 && VIn->ZF!=0 ?  -VIn->alphbar*VIn->q/VIn->kn[n]*(VIn->kn[n]*VIn->kn[n]/(VIn->g*VIn->g) -VIn->q*VIn->q)/(VIn->kn[n]*VIn->g)*(conj(phiZF[0])*conj(phi0[n-1])) : 0); 
                dphidt[i] = C_phi_phi + frcng  + disp + Z_fl;
                
            }
            
            
            
            //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%//
            //                  Denisty 
            ///////////////////////////////////  NY /////////////////////
            if((i>=(VIn->inYB)) && (i<=(VIn->inYE))){
                C_phi_n = (((n>1) ?   VIn->any*VIn->n_an[n]*(conj(n0[n-1])*conj(phi0[n-2])-conj(phi0[n-1])*conj(n0[n-2])) :0 )
                           +( (n>0 &&n<(N-1)) ?  VIn->bny*VIn->n_bn[n]*(conj(phiX[n-1])*conj(n0[n+1])-conj(phi0[n+1])*conj(nX[n-1])) : 0 )
                           +((n<(N-2)) ?   VIn->cny*VIn->n_cn[n]*(conj(phiX[n+1])*conj(n0[n+2])-conj(phi0[n+2])*conj(nX[n+1])) : 0)
                           );
                
                
                // dissipations 
                disp  =  -(musm*pow(VIn->kn[n],-6)+mubg*pow(VIn->kn[n],4))*nY[n];    
                // forcing terms
                frcng = C2*(phiY[n]-nY[n]) - 0.7653*I*vn*VIn->kn[n]*phiY[n];
                
                if(n==fi||n==fi+1)
                    frcng=frcng+Fn;
                //zonal part
                Z_fl =((n>0 && VIn->ZF!=0) ? VIn->alphbar*VIn->q*VIn->kn[n]*( (conj(phiZF[0])*conj(n0[n+1])-conj(nZF[0])*conj(phi0[n+1]))) : 0);
                
                dphidt[i] =   C_phi_n + frcng  + disp + Z_fl;
            }
            
            ///////////////////////////////////  N0 /////////////////////
            if((i>=(VIn->in0B)) && (i<=(VIn->in0E))){
                C_phi_n = (((n>1) ?  VIn->an0*VIn->n_an[n]*(conj(nY[n-1])*conj(phiX[n-2])-conj(phiY[n-1])*conj(nX[n-2])+(conj(nX[n-1])*conj(phiY[n-2])-conj(phiX[n-1])*conj(nY[n-2]))) :0 )
                           +( (n>0 &&n<(N-1)) ?  VIn->n_bn[n]*(VIn->bn0x*(conj(phi0[n-1])*conj(nX[n+1])-conj(phiX[n+1])*conj(n0[n-1]))+VIn->bn0y*((conj(phi0[n-1])*conj(nY[n+1])-conj(phiY[n+1])*conj(n0[n-1])))) : 0 )
                           +((n<(N-2)) ?  VIn->n_cn[n]*(VIn->cn0x*(conj(phi0[n+1])*conj(nX[n+2])-conj(phiX[n+2])*conj(n0[n+1]))+VIn->cn0y*((conj(phi0[n+1])*conj(nY[n+2])-conj(phiY[n+2])*conj(n0[n+1])))) : 0)
                           );
                
                
                // dissipations 
                disp  =  -(musm*pow(VIn->kn[n],-6)+mubg*pow(VIn->kn[n],4))*n0[n];    
                // forcing terms
                frcng = C2*(phi0[n]-n0[n]) - 0.5411*I*vn*VIn->kn[n]*phi0[n];
                
                
                if(n==fi||n==fi+1)
                    frcng=frcng+Fn;
                
                //zonal part
                Z_fl =((n>0 && n<(N-1) && VIn->ZF!=0) ? -VIn->alphbar*VIn->q*VIn->kn[n]*(conj(phiZF[0])*conj(nY[n-1])-conj(nZF[0])*conj(phiY[n-1]))/VIn->g/2. + VIn->alphbar*VIn->q*VIn->kn[n]*(conj(phiZF[0])*conj(nX[n+1])-conj(nZF[0])*conj(phiX[n+1]))/2.  : 0);
                
                dphidt[i] =  C_phi_n + frcng  + disp + Z_fl;
                
            }
            
            ///////////////////////////////////  NX /////////////////////
            if((i>=(VIn->inXB)) && (i<=(VIn->inXE))){
                
                C_phi_n = (((n>1) ?  VIn->anx*VIn->n_an[n]*(conj(n0[n-1])*conj(phi0[n-2])-conj(phi0[n-1])*conj(n0[n-2]))  :0 )
                           +( (n>0 &&n<(N-1)) ?  VIn->bnx*VIn->n_bn[n]*(conj(phiY[n-1])*conj(n0[n+1])-conj(phi0[n+1])*conj(nY[n-1])) : 0 )
                           +((n<(N-2)) ?  VIn->cnx*VIn->n_cn[n]*(conj(phiY[n+1])*conj(n0[n+2])-conj(phi0[n+2])*conj(nY[n+1])) : 0)
                           );
                
                
                
                // dissipations 
                disp  =  -(musm*pow(VIn->kn[n],-6)+mubg*pow(VIn->kn[n],4))*nX[n];    
                // forcing terms
                frcng = C2*(phiX[n]-nX[n]);
                
                if(n==fi||n==fi+1)
                    frcng=frcng+Fn;
                
                //zonal part
                Z_fl =((n>0 && VIn->ZF!=0) ? -VIn->alphbar*VIn->q*VIn->kn[n]*(conj(phiZF[0])*conj(n0[n-1])-conj(nZF[0])*conj(phi0[n-1]))/VIn->g  : 0);                                                                                                      
                
                dphidt[i] =  C_phi_n + frcng  + disp + Z_fl;
            }
            //  cout<<t<<"\t"<<i<<"\t"<<n<<"\t"<<abs(dphidt[i])*abs(dphidt[i])<<"\n";      
            
        }
        
    }
    if(VIn->ZF!=0){
        dphidt[0] =  -VIn->alphbar*(VIn->g*VIn->g-1)/VIn->q*SpZF -VIn->nuZF*phiZF[0];//VIn->alphbar*(VIn->g*VIn->q-1)/VIn->q*SpZF  -VIn->nuZF*phiZF[0];                                  
        dphidt[N_V-1] = -VIn->alphbar*VIn->q*SnZF;// no -VIn->nuZF*x[i] for n;                                                                                                                         
    }
    
    
    return GSL_SUCCESS;
};






///////////////////////////  INPUT DATA /////////////////////////

double get_data_NAME(FILE *someFile, const char * nameData)
{
    double vlu;
    int dataExists;
    char dat;
    int dg,i;
    char data[20]; // = malloc(20*sizeof(char));
    
    dataExists = 0;
    dat = fgetc(someFile);
    
    while(dat != EOF && dataExists == 0)
    {
        i=0;     
        dat = fgetc(someFile);
        
        if(dat == nameData[i])
        {
            while(dat == nameData[i])
            {	    
                dat = fgetc(someFile);
                i++;
            }
            
            if(nameData[i] == '\0')
            { 
                dataExists = 1;		  
            }
        }
    }
    
    dg =0;  
    dat = fgetc(someFile);
    
    while((dat == '\t') || (dat == ' '))
    {
        dat = fgetc(someFile);		 
    }
    while((dat != '\t') && (dat != ' ') && (dat != EOF) && (dat != '\0')  )
    {
        data[dg] = dat;
        dat = fgetc(someFile);		 
        dg++;
    }
    
    vlu =(double) atof(data);
    
    if(dataExists ==0){
        printf("\n You goofed! Make sure the correct spelling of %s exists in the INPUT file. \n \n",nameData);
        abort();
    }
    
    return vlu;
    
    
}







