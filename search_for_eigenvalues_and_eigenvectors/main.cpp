#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#define eps 0.0000000001

using namespace std;

//--------------------------------------

double **matr_L(double**, int );
double **mart_R_alternative(double**, int );
double **fout_matr(double**, int, ofstream & );
double **read_matr(int );
//double **out_matr(double**, int );
double **mult_matr(double**, double**, int, ofstream & );
double **zero_matr( int );
double **LR_m(double**, int, ofstream &);
double frebenius(double **a, int n);
double norm_vek ( double* , int );
double *mnojmatrvec(double** , double*, int );
double *mnojchvec(double , double* , int );
void scalar ( double**  , int , ofstream &);
double **create_matrix(int);
double **copy (int , double**, double**);
bool norma (double**, int );


//--------------------------------------

int main() {
    ofstream res("result.txt");
    int n = 4;
    double **A = create_matrix(n);

    A=zero_matr(4);
    A=read_matr(4);
    scalar(A, n, res);

    A=zero_matr(4);
    A=read_matr(4);
    LR_m(A,n,res);

    res.close();
    cout << endl;
}

//-------------------------------------------------------------------------


bool norma (double **A, int n){
    bool flag = true;
    double norm = 0.;
    for ( int i = 1; i < n ; i ++ ){
        for ( int j = 0; j < i; j ++){
            norm += A[i][j] * A[i][j];
        }
    }
    if (fabs(norm)<  eps){
        flag = false;
    }
    return flag;

}

double **LR_m(double **A, int n, ofstream& fout){
    double **L=create_matrix(n);
    double **E=create_matrix(n);
    double **B=create_matrix(n);
    double **R=create_matrix(n);
    E=zero_matr(n);
    for (int i = 0; i < n; i++) {
        E[i][i] = 1;
    }
    bool flag = true;
    B=read_matr(n);
    int r=1;
    do {

        fout_matr(B,n,fout);
        R=mart_R_alternative(B,n);
        fout<<endl<<"matr R"<<endl;
        fout_matr(R,n,fout);
        fout<<endl;
        L=matr_L(A,n);
        fout<<endl<<"matr L"<<endl;
        fout_matr(L,n,fout);
        cout<<endl;

        fout<< endl << endl << "mult matr" << endl << endl;
        A = mult_matr(R, L, n, fout);
        fout<< endl;
        B = copy(n,B,A);
        fout<< "iteration "<< r<< endl;
    r++;
            flag = norma (A,n);
        E = mult_matr(E, L, n, fout);
        fout << "frebenius norm  " << frebenius(R, n)<< endl<< endl;
    } while ((flag) && (r<100));


   // E = mult_matr(E, L, n, fout);

    double cc;
    double *v_v=new double [n];
    double *v = new double [n];
    double *vecnev = new double[n];
    for (int k = 0; k < n; k++) {
        v_v[k] = 1;
        for (int i = k + 1; i < n; i++)
            v_v[i] = 0;
        if (k > 0)
            for (int i = k - 1; i >= 0; i--) {
                cc = 0;
                for (int j = i; j < n; j++) {
                    cc += R[i][j + 1] * v_v[j + 1];
                }
                v_v[i] = cc / (R[k][k] - R[i][i]);
            }

        v = mnojmatrvec(E, v_v, n);
        fout << "Власний вектор_" << k + 1 << " = (";
        for (int i = 0; i < n; i++)
            fout << " " << v[i];
        fout << ")" << endl;

        fout << "Vector of residual_" << k + 1 << " = (";
        for (int i = 0; i < n; i++) {
            vecnev[i] = mnojmatrvec(A, v, n)[i] - mnojchvec(R[k][k], v, n)[i];
            fout << " " << vecnev[i];
        }
        fout << ")" << endl;
    }


    fout<< " r = "<< r;
}

double **zero_matr(int n){
    int i, j;
    double **A = create_matrix(n);

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            A[i][j] = 0;
        }
    }
    return A;
}
double **read_matr(int n) {
    ifstream L("matr.txt");
    double **A=create_matrix(n);
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            L >> A[i][j];
        }
    }
    L.close();
    return A;
}

double **copy (int n, double **A, double **B){
    for ( int i = 0; i < n; i++){
        for (int  j = 0; j < n; j++){
            A[i][j]=B[i][j];
        }
    }
    return A;
}

double **create_matrix(int n){
    double **B = new double *[n];
    for (int i = 0; i < n; i++){
        B[i] = new double [n];
    }

    return B;
}

double **mart_R_alternative(double **A, int n){
    double **B = create_matrix(n);
    int k;

    B=copy(n, B, A);

    for ( k = 1; k < n; k++) {
        for (int j = 0; j < n; j++){
            for (int i = k; i < n; i++  ) {
                B[i][j] = A[i][j] - A[k - 1][j] * A[i][k-1] / A[k - 1][k - 1];
            }
        }
        A=copy(n,A,B);
    }
    return B;
}

double **matr_L(double **A, int n) {
    int i,j;
    double **L = create_matrix(n);

    L = zero_matr(n);

    for (i = 0; i < n; i++){
        L[i][i]=1;
    }

    double c;
    for (int k = 0; k < n - 1; k++) {
        if (A[k][k] == 0)
            A[k][k] = eps;
        for (int i = k + 1; i < n; i++) {
            c = A[i][k] / A[k][k];
            for (int j = k; j < n; j++) {
                A[i][j] = A[i][j] - c * A[k][j];
                if (i > j)
                    L[i][j] = c;
            }
        }
    }

    return L;
}

//double **out_matr(double** A, int n){
//    for( int i =0; i<n; i++){
//        cout<< endl;
//        for ( int j=0; j<n; j++){
//            cout<< setw(15)<< A[i][j];
//        }
//    }
//    cout<<endl;
//}

double **fout_matr(double** A, int n, ofstream &fout){
    for( int i =0; i<n; i++){
        fout<< endl;
        for ( int j=0; j<n; j++){
            fout<< setw(15)<< A[i][j];
        }
    }
    fout<<endl;
}

double **mult_matr(double **A, double **B, int n, ofstream &fout) {
    double **C = create_matrix(n);

    for (int i = 0; i < n; i++)
    {
        for (int counter = 0; counter < n; counter++)
        {
            for (int j = 0; j < n; j++)
            {
                C[i][counter] += A[i][j] * B[j][counter];
            }
        }
    }
    fout_matr(C,n,fout);
    fout<<endl;
    return C;
}

double norm_vek ( double * v, int n){
    double norm = 0.;
    for ( int i =0 ; i < n ; i ++) norm += v[i]*v[i];
    return sqrt (norm);
}

double *mnojmatrvec(double **A, double *v, int n) {
    double *m = new double[n];
    for (int i = 0; i < n; i++) {
        m[i] = 0;
        for (int j = 0; j < n; j++)
            m[i] += A[i][j] * v[j];
    }
    return m;
}

double  scalar_v_v ( double * a, double * b, int n){
    double  v = 0.;
    for ( int i = 0; i < n; i++){
        v += a[i]*b[i];
    }
    return v;
}

double *mnojchvec(double ch, double *vec, int n) {
    double *x = new double[n];
    for (int i = 0; i < n; i++)
        x[i] = ch*vec[i];
    return x;
}

double frebenius(double **A, int n) {
    double norm = 0;
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            norm += pow(A[i][j], 2);
    norm = sqrt(norm);
    return norm;
}

double *get_resuidal ( double ** A, double * b, double  l, int n){
    double *res = new double [n];
    double *v1 = mnojmatrvec( A, b, n);
    double *v2 = mnojchvec( l, b, n);
    for (int i = 0; i < n; i++) res[i] = v1[i] - v2[i];
    return res;
}

void scalar ( double ** A , int n, ofstream &fout) {
    double *vl_vek = new double[n];
    for (int i = 0; i < n; i++) vl_vek[i] = 1;
    double *dop = new double[n];
    double lamda_n = 1124;
    double lamda_p;
    double c;
    double pohibka;
    int iteration = 0;
    fout << endl<< endl << "--------------------------------------" << endl;
    do {
        iteration++;
        fout << endl << "------------Ітерація------- " << iteration << endl;
        lamda_p = lamda_n;
        dop = mnojmatrvec(A, vl_vek, n);
        c = 1 / norm_vek(dop, n);
        vl_vek = mnojchvec(c, dop, n);
        lamda_n = scalar_v_v(vl_vek, mnojmatrvec(A, vl_vek, n), n);
        pohibka = fabs(fabs(lamda_p) - fabs(lamda_n));
        fout << " Найбыльше власне число  " << lamda_p << endl;
        fout << endl << "Власний вектор   ";
        for (int i = 0; i < n; i++){
            fout << setw(9) << vl_vek[i];
        }
    } while (pohibka > eps);
    double lamda_max = lamda_p;
    fout << " Найбыльше власне число  " << lamda_max;
    fout << endl << "Власний вектор  ";
    for (int i = 0; i < n; i++){
        fout << setw(9) << vl_vek[i];
    }

    double *res = get_resuidal(A, vl_vek, lamda_p, n);
    fout << endl << " Вектор невязкы для ( (" << lamda_p << "),";
    for (int j = 0; j < n; j++) {
        fout << setw(10) << vl_vek[j];
    }
    fout << ")      ";
    for (int j = 0; j < n; j++) {
        fout << setw(13) << res[j];
    }


    double **B = create_matrix(n);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            B[i][j] = A[i][j];
        }
    }
    int t = 0;
    for (int i = 0; i < n; i++) B[i][i] -= lamda_n;
    fout << endl << endl << "********** пошук мінімального власного числа *************" << endl << endl;
    do {t++;
        lamda_p = lamda_n;
        dop = mnojmatrvec(B, vl_vek, n);
        c = 1 / norm_vek(dop, n);
        vl_vek = mnojchvec(c, dop, n);
        lamda_n = scalar_v_v(vl_vek, mnojmatrvec(B, vl_vek, n), n);
        pohibka = fabs(fabs(lamda_p) - fabs(lamda_n));
        double lamda_min;
        lamda_min=(lamda_p + lamda_max);
        fout << endl << "Найменше власне число   " << lamda_min << endl;
        fout << endl << "Власний вектор   ";
        for (int i = 0; i < n; i++) {
            fout << setw(11) << vl_vek[i];
        }
    } while (pohibka > eps);
    fout<<endl<<"Ітерація "<< t <<endl;
    double pohibka_vlasne;
    pohibka_vlasne = (lamda_p + lamda_max);
    fout << endl << "Найменше власне число     " << pohibka_vlasne;
    fout << endl << "Власний вектор   ";
    for (int i = 0; i < n; i++) {
        fout << setw(11) << vl_vek[i];
    }
}