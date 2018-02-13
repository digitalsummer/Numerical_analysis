package com.company;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;

class  Hiperbolichna_zadacha {

    private double sigma = 0.0;
    private double x0 = 0.0;
    private double xn = 1.0;
    private int K = 100;
    private int N = 10;

    public void cleanfile (String file){

        File res = new File( file + ".txt");
        res.delete();
    }

    public double u0(double x){
        return (x*x*Math.cos(Math.PI*x));
    }

    public double u1(double t){
        return (u0(0)*Math.cos(Math.PI*t));
    }

    public double u2(double t){
        return (u0(xn)*Math.cos(Math.PI*t));
    }

    public double u_(double t, double x){
        return (u0(x)*Math.cos(Math.PI*t));
    }

    public double d2f_dx2(double x){
        return (2*Math.cos(Math.PI*x) - 2*x*Math.PI*Math.sin(Math.PI*x) - Math.PI*Math.PI*x*x*Math.cos(Math.PI*x) - 2*Math.PI*x*Math.sin(Math.PI*x));
    }

    public double F(double t, double x){
        return ( (Math.PI*Math.PI*(-u0(x))-d2f_dx2(x))*Math.cos(Math.PI*t) );
    }

    public double gauss(double f[], double a[][], int n, double y[]){
        double d =1;
        double aa[][];
        int ee[];
        ee = new int [n];
        aa = new double[n][n];
        for(int i = 0; i < n; i++){
            y[i] = f[i];
        }
        for (int i = 0; i < n; i++){
            for (int j = 0; j < n; j++){
                aa[i][j] = a[i][j];
            }
        }
        for (int i = 0; i < n; i++){
            ee[i] = i;
        }
        int qq =0;
        for (int i = 0; i < n; i++){
            if(qq == 0){
                double max = 0;
                int imax = 0;
                for (int j=i;j<n;j++)
                    if (max < Math.abs(aa[i][j])) {
                        max = Math.abs(aa[i][j]);
                        imax = j;
                    }
                if (max==0) {
                    qq = 1;
                    System.out.println("Oups");
                } else {
                    if (i!=imax) {
                        double temp;
                        for (int j=0;j<n;j++) {
                            temp = aa[j][i];
                            aa[j][i] = aa[j][imax];
                            aa[j][imax] = temp;
                        }
                        temp = ee[i];
                        ee[i] = ee[imax];
                        ee[imax] = (int)temp;
                        d *= -1;
                    }
                    d *= aa[i][i];
                    for (int j=0;j<n;j++)
                        if (j!=i){
                            double temp = aa[j][i];
                            for (int k = 0;k < n;k++)
                                aa[j][k] = aa[j][k] - (aa[i][k])*temp / (aa[i][i]);
                            y[j] = y[j] - y[i] * temp / aa[i][i];
                        }
                }
            }

        }
        if (qq == 1)
            d = 0;
        double[] temp;
        temp = new double [n];
        for (int i = 0;i < n;i++)
            temp[i] = y[i];
        for (int i = 0;i < n;i++)
            y[ee[i]] = temp[i] / aa[i][i];
        return d;
    }

    public void zadacha() {
        String filee = "file";
        String xxx  = "x";
        String norm = "norma";
        String error = "err";
        cleanfile(filee);
        cleanfile(error);
        cleanfile(norm);
        cleanfile(xxx);

        double tn = 2;
        double dx = (xn - x0) / N;
        double dt = tn / K;
        double a[][];
        double u[][];
        double f[];
        double y[];
        u = new double[K + 1][N + 1];
        a = new double[N + 1][N + 1];
        f = new double[N + 1];
        y = new double[N + 1];

        for (int i = 0; i <= N; i++)
            u[0][i] = u0(x0 + dx * i);
        for (int i = 1; i <= K; i++) {
            u[i][0] = u1(i * dt);
            u[i][N] = u2(i * dt);
        }
        for (int i = 1; i < N; i++)
            u[1][i] = u[0][i] + Math.pow(dt, 2) * ((u[0][i + 1] - 2 * u[0][i] + u[0][i - 1]) / Math.pow(dx, 2) + F(0, x0 + i * dx)) / 2;

        for (int k = 2; k <= K; k++) {
            for (int i = 0; i <= N; i++)
                for (int j = 0; j <= N; j++)
                    a[i][j] = 0;
            a[0][0] = 1;
            f[0] = u[k][0];
            a[N][N] = 1;
            f[N] = u[k][N];
            for (int i = 1; i < N; i++) {
                a[i][i] = 1 / Math.pow(dt, 2) + 2 * sigma / Math.pow(dx, 2);
                a[i][i - 1] = -sigma / Math.pow(dx, 2);
                a[i][i + 1] = -sigma / Math.pow(dx, 2);
                f[i] = (2 * u[k - 1][i] - u[k - 2][i]) / Math.pow(dt, 2) + (1 - 2 * sigma) * (u[k - 1][i + 1] - 2 * u[k - 1][i] + u[k - 1][i - 1]) / Math.pow(dx, 2) + sigma * (u[k - 2][i + 1] - 2 * u[k - 2][i] + u[k - 2][i - 1]) / Math.pow(dx, 2) + sigma * F(dt * k, x0 + i * dx) + (1 - 2 * sigma) * F(dt * (k - 1), x0 + i * dx) + sigma * F(dt * (k - 2), x0 + i * dx);
            }
            gauss(f, a, N + 1, y);
            for (int i = 1; i < N; i++)
                u[k][i] = y[i];
        }

        f = new double[K + 1];
        y = new double[K + 1];
        for (int i = 0; i <= K; i++) {
            f[i] = 0;
            y[i] = 0;
        }
        try {
            FileWriter writer = new FileWriter("file.txt", true);
            BufferedWriter file = new BufferedWriter(writer);

            FileWriter writer1 = new FileWriter("err.txt", true);
            BufferedWriter err = new BufferedWriter(writer1);

            FileWriter writer2 = new FileWriter("norma.txt", true);
            BufferedWriter norma = new BufferedWriter(writer2);

            FileWriter writer3 = new FileWriter("x.txt", true);
            BufferedWriter xx = new BufferedWriter(writer3);

            for (int i = 0; i <= K; i++) {
                for (int j = 0; j <= N; j++) {
                    if (Math.abs(u[i][j]) > y[i])
                        y[i] = Math.abs(u[i][j]);
                    if (f[i] < Math.abs(u_(dt * i, x0 + dx * j) - u[i][j]))
                        f[i] = Math.abs(u_(dt * i, x0 + dx * j) - u[i][j]);
                }
                file.write("t = " + (dt * i) + "\n");
                file.write("u(t) = " + y[i]+ "\n");
                file.write("e(t) = " + f[i]+ "\n");
                xx.write("t = " + (dt * i) + "\n");
                norma.write("u(t) = " + y[i]+ "\n");
                err.write("e(t) = " + f[i] + "     t" + (dt * i) + "\n");
            }
            file.write("t = " + tn + "\n");

            for (int i = 0; i <= N; i++) {
                file.write("x = " + (x0 + i * dx) + "\n");
                file.write("u(x) = " + u[K][i] + "\n");
                file.write("e(x) = " +( Math.abs(u[K][i] - u_(tn, x0 + i * dx)) )+"\n");

                xx.write("x = " + (x0 + i * dx) + "\n");
                norma.write("u(x) = " + u[K][i] + "\n");
                err.write("e(x) = " + (Math.abs(u[K][i] - u_(tn, x0 + i * dx)))  +"  x   " +(x0 + i * dx)+"\n");
                //err.write((Math.abs(u[K][i] - u_(tn, x0 + i * dx)))+"\n");
                //xx.write((x0 + i * dx) + "\n");
            }
            file.close();
            err.close();
            xx.close();
            norma.close();
        } catch (IOException e) {
            System.out.println(e);
        }
    }

}
public class Main {

    public static void main(String[] args) {
        Hiperbolichna_zadacha hb = new Hiperbolichna_zadacha();
        hb.zadacha();
        }
	// write your code here
    }

