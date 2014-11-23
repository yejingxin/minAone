#include <cstdlib>
#include <cmath>
#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>
#include <sstream>
#include <list>
#include <algorithm>
#include "RKA_double.h"

#define NSTAT 4
#define NPARA 19
#define NSTIM 1
#define NTIME 100000
#define NMEAS 1
#define DT 0.02

using namespace std;
void Func(double *x, double t, double *p, double *dx);

int main(int argc, char **argv){

    double *P = new double[NPARA+NSTIM];
    int i,j,k;
    double X[NSTAT];
    double Xout[NSTAT];
    const int WIDTH=15;
    const int PRECISION=7;

    int m_iter;
    double taus[1];
    taus[0] = DT;

    int mlistint[NMEAS] = {0};

    list<int> meas_list (mlistint,mlistint+NMEAS);
    P[0] = 120;
    P[1] = 50;
    P[2] = 20;
    P[3] = -77;
    P[4] = .3;
    P[5] = -54;
    P[6] = 0.8;
    P[7] = -40;
    P[8] = 0.06667;
    P[9] = .1;
    P[10] = .4;
    P[11] = -60;
    P[12] = -.06667;
    P[13] = 1;
    P[14] = 7;
    P[15] = -55;
    P[16] = .03333;
    P[17] = 1;
    P[18] = 5;
    X[0] = -64;
    X[1] = 0.05;
    X[2] = 0.6;
    X[3] = 0.4;
    ifstream stim[NSTIM];
    stim[0].open("current.dat");
    ofstream *meas = new ofstream[NMEAS];
    meas[0].open("measured.dat");
    //cout<<"here"<<endl;
    ofstream *outfile = new ofstream;
    outfile->open("allstates.dat");
    //(*outfile)<<setw(WIDTH)<<setprecision(PRECISION)<<0.0;
    for(i=0;i<NSTAT;i++){
        (*outfile)<<setw(WIDTH)<<setprecision(PRECISION)<<X[i];
        //cout<<setw(WIDTH)<<setprecision(PRECISION)<<X[i];
    }

    (*outfile)<<endl;
    //cout<<endl;
    for(i=0;i<NTIME;i++){
        if(X[0] != X[0]){
            cout<<"NAN error"<<endl;
            break;
        }
        for(j=0;j<NSTIM;j++){
            stim[j]>>P[NPARA+j];
        }
        rkA(NSTAT,X,i*DT,DT,Func,P,taus,Xout);

        for(j=0;j<NSTAT;j++){
            X[j]=Xout[j];
        }
        //(*outfile)<<setw(WIDTH)<<setprecision(PRECISION)<<(i+1)*DT;
        m_iter=0;
        for(j=0;j<NSTAT;j++){
            if(std::find( meas_list.begin(), meas_list.end(),j)!=meas_list.end()){
                meas[m_iter]<<X[j]<<endl;
                m_iter++;
            }
            (*outfile)<<setw(WIDTH)<<setprecision(PRECISION)<<X[j];

            //cout<<setw(WIDTH)<<setprecision(PRECISION)<<X[j];

            //cout<<(i+1)/RESO<<endl;
            }
        (*outfile)<<endl;
        //cout<<endl;
    }
    return 0;
}
void Func(double *x, double t, double *p, double *dx){
    dx[0] = p[0]*(x[1]*x[1]*x[1]*x[2])*(p[1]-x[0])+p[2]*x[3]*x[3]*x[3]*x[3]*(p[3]-x[0])+p[4]*(p[5]-x[0])+p[6]*p[NPARA+0];
    dx[1] = (0.5*(1+tanh((x[0]-p[7])*p[8])) - x[1])/(p[9]+p[10]*(1.0-tanh((x[0]-p[7])*p[8])*tanh((x[0]-p[7])*p[8])));
    dx[2] = (0.5*(1+tanh((x[0]-p[11])*p[12])) - x[2])/(p[13]+p[14]*(1.0-tanh((x[0]-p[11])*p[12])*tanh((x[0]-p[11])*p[12])));
    dx[3] = (0.5*(1+tanh((x[0]-p[15])*p[16])) - x[3])/(p[17]+p[18]*(1.0-tanh((x[0]-p[15])*p[16])*tanh((x[0]-p[15])*p[16])));
}

/*
    P[0] = 120;
    P[1] = 50;
    P[2] = 20;
    P[3] = -77;
    P[4] = .3;
    P[5] = -54;
    P[6] = 0.8;
    P[7] = -40;
    P[8] = 0.06667;
    P[9] = .1;
    P[10] = .4;
    P[11] = -60;
    P[12] = -.06667;
    P[13] = 1;
    P[14] = 7;
    P[15] = -55;
    P[16] = .03333;
    P[17] = 1;
    P[18] = 5;

*/

