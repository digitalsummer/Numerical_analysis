#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>

using namespace std;

//--------------------------------------------------------------

double f(double , double );
double real_f(double);
double *runge(double , double , int, double, ofstream &fout);
double *adams(double , double, int, double, ofstream &fout);
double *real(double, int, double, ofstream &fout);
double *lab(double , ofstream &fout);

//--------------------------------------------------------------

int main() {
    ofstream fout;
    fout.open("output.txt", ios::out );

    double h1=0.1, h2=0.05;
    int n1 = round(1/h1);
    int n2 = round(1/h2);

    double *E_runge_h1 = new double [n1];
    double *E_runge_h2 = new double [n2];
    double *e_runge = new double[n1];

    double h=0.05;

    for ( int i =0; i<5; i++){
        lab(h, fout);
        h=h+0.01;
    }

    fout<< " ВІДНОШЕННЯ ПОМИЛОК ПРИ Н=0.1 ТА Н=0.05"<<endl;
    E_runge_h1=lab(0.1, fout);
    E_runge_h2=lab(0.05, fout);
    for (int i=0; i<n1; i++){
        e_runge[i]=E_runge_h1[i]/E_runge_h2[2*i];
        fout<<e_runge[i]<<endl;
    }

    fout.close();
    return 0;
}

//----------------------------------------------------------------

double *lab(double h, ofstream &fout){

    int n=int(1/h);
    double x_0 = 0.0;
    double y_0 = 1.0;
    double *Y_real = new double [n];
    double *Y_runge = new double [n];
    double *Y_adams = new double [n];
    double *E_adams = new double [n];
    double *E_runge = new double [n];

    double *x_n1 = new double [n];

    fout<<endl<<"**************FOR H1 = "<< h <<"*****************"<< endl<<endl;
    fout << "-----------Runge method------------" << endl;
    Y_runge=runge(x_0, y_0, n, h, fout);
    fout << "-----------Adams-------------------" << endl;
    Y_adams=adams(x_0, y_0, n,h, fout);
    fout << "-----------Real--------------------" << endl;
    Y_real=real(x_0, n, h, fout);
    fout << "------------EPSILON----------------" << endl;
    x_n1[0]=x_0;
    fout << "------------epsilon Adams----------" << endl;
    for (int i=0; i<n; i++){
        E_adams[i]=fabs(Y_real[i]-Y_adams[i]);
        fout << E_adams[i] <<"   x=   "<<x_n1[i]<< endl;
        x_n1[i+1]=x_n1[i]+h;
    }
    fout << "-------------epsilon Runge---------" << endl;
    for (int i=0; i<n; i++){
        E_runge[i]=fabs(Y_real[i]-Y_runge[i]);
        fout << E_runge[i]<<"  x=  "<<x_n1[i] << endl;
        x_n1[i+1]=x_n1[i]+h;
    }
    return E_runge;
}
double real_f(double x) {
    return ((exp(x) + 8) / 9);
}
double f(double x, double y){
    return (y*exp(x)/(exp(x)+8));
}

double *runge(double x, double y, int n, double h, ofstream &fout){

    double *X=new double [n];
    double *Y=new double [n];
    int i;

    X[0]=x;
    Y[0]=y;

    for(i=0; i<=n; i++){
        double k1,k2,k3,k4;
        k1=f(X[i],Y[i]);
        k2=f(X[i]+h/2,Y[i]+h/2*k1);
        k3=f(X[i]+h/2,Y[i]+h/2*k2);
        k4=f(X[i]+h,Y[i]+h*k3);
        Y[i+1]=Y[i]+h*(k1+2*k2+2*k3+k4)/6;
        X[i+1]=X[i]+h;
        fout<< X[i]<< " = x "<<setprecision(9)<<Y[i]<<" = y "<<endl;
    }
    return Y;
}

double *adams(double x, double y, int n, double h, ofstream &fout){
    double *Y=new double [n], *X=new double [n];
    int i;
    for(i=0; i<=n; i++){
        Y[i]=0;
        X[i]=0;
    }
    Y[0]=y;
    X[0]=x;
    for(i=0; i<=3; i++){
        double k1,k2,k3,k4;
        k1=f(X[i],Y[i]);
        k2=f(X[i]+h/2,Y[i]+h/2*k1);
        k3=f(X[i]+h/2,Y[i]+h/2*k2);
        k4=f(X[i]+h,Y[i]+h*k3);
        Y[i+1]=Y[i]+h*(k1+2*k2+2*k3+k4)/6;
        X[i+1]=X[i]+h;
    }
    for(i=4; i<=n; i++){
        Y[i+1]=Y[i]+h*(55*f(X[i],Y[i])-59*f(X[i-1],Y[i-1])+37*f(X[i-2],Y[i-2])-9*f(X[i-3],Y[i-3]))/24;
        X[i+1]=X[i]+h;
    }
    for(i=0; i<=n; i++){
        fout<<"y "<<setprecision(9)<<Y[i]<<" x "<<X[i]<<endl;
    }
    return Y;
}

double *real(double x, int n, double h, ofstream &fout){
    double *X=new double [n];
    double *Y=new double [n];
    X[0]=x;
    int i;
    for ( i=0;i<=n;i++){
        Y[i]=real_f(X[i]);
        X[i+1]=X[i]+h;
        fout << "x = "<<setprecision(9)<< X[i]<<"  y =  "<<Y[i]<<endl;
    }
    return Y;
}
