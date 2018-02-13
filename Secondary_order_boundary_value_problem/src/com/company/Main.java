package com.company;

import jdk.nashorn.internal.runtime.OptimisticReturnFilters;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;

class Krayova_zadacha {
    private double y[], x[], vec_b[], vec_c[], vec_a[], function[], error[];
    private double h, f1, f2;
    public int n;
    double a = 2.831;
    double b = -2.557;
    double c = -0.157;
    double d = 1.089;
    double e = 1.821;
    double q = 1;
    double alpha1 = 1;
    double betha2 = 1;
    double betha1 = 2;

    public double y_0(double x) {

        return a * Math.pow(x, 2) + b * x + c + Math.pow((d * x + e), -1);
    }

    public double y_1(double x) {

        return 2 * a * x + b - d * Math.pow((d * x + e), -2);
    }

    public double y_2(double x) {

        return 2 * a + 2 * Math.pow(d, 2) * Math.pow((d * x + e), -3);
    }

    Krayova_zadacha(int an) {

        n = an;
        x = new double[n + 1];
        function = new double[n + 1];
        y = new double[n + 1];
        error = new double[n + 1];
        vec_b = new double[n + 1];
        vec_c = new double[n + 1];
        vec_a = new double[n + 1];
        x[0] = 0.5;
        x[n] = 0.8;
        h = (x[n] - x[0]) / n;
        for (int i = 1; i < n; i++)
            x[i] = x[i - 1] + h;
        f1 = y_0(x[0]) + 2 * y_1(x[0]);
        f2 = y_1(x[n]);
    }

    public void cleanfile (String file){
        //try {
            File res = new File( file + ".txt");
            //BufferedWriter res = new BufferedWriter(writer);
            res.delete();
        //}
//        catch (IOException e) {
//            System.out.println(e);
//        }
    }

    public void TridiagonalMatrixAlgorithmSolving(String Filename )  {
        for (int i = 0; i < n + 1; i++)
            function[i] = y_2(x[i]) + x[i] * y_1(x[i]) + q * y_0(x[i]);

        vec_a[0] = 0;
        vec_b[n] = 0;

        vec_c[0] = 2 * alpha1 / betha1 / h - 2. / h / h - x[0] * alpha1 / betha1 + q;
        vec_b[0] = 2. / h / h;
        function[0] += 2 * f1 / h / betha1 - x[0] * f1 / betha1;

        vec_a[n] = 2. / h / h;
        vec_c[n] = q - 2. / h / h;
        function[n] = function[n] - 2. * f2 / betha2 / h - x[n] * f2 / betha2;


        for (int i = 1; i < n; i++) {
            vec_c[i] = -2 / Math.pow(h, 2) + q;
            vec_b[i] = 1 / Math.pow(h, 2) + x[i] / (2 * h);
            vec_a[i] = 1 / Math.pow(h, 2) - x[i] / (2 * h);
        }

        double m;
        for (int i = 1; i < n + 1; i++) {
            m = vec_a[i] / vec_c[i - 1];
            vec_c[i] -= m * vec_b[i - 1];
            function[i] -= m * function[i - 1];
        }
        y[n] = function[n] / vec_c[n];

        for (int i = n - 1; i >= 0; i--)
            y[i] = (function[i] - vec_b[i] * y[i + 1]) / vec_c[i];

        try {
            FileWriter writer = new FileWriter(Filename + ".txt", true);
            BufferedWriter res = new BufferedWriter(writer);
            res.write("");
            for (int i = 0; i < n + 1; i++)

                error[i] = Math.abs(y[i] - y_0(x[i]));
            res.append("\n" + "x ---------------- = " + "y ----------------  =  " + "yn ---------------  = " + "err--------------  = " + "\n");
            for (int i = 0; i < n + 1; i++) {
                res.append(x[i] + ";");
               // System.out.print(x[i] + ";");
                res.append(y_0(x[i]) + ";");
               // System.out.print(y_0(x[i]) + ";");
                res.append(y[i] + ";");
               // System.out.print(y[i] + ";");
                res.append(error[i] + ";" + "\n");
               // System.out.print(error[i] + ";" + "\n");
            }
            res.close();
        }
        catch (IOException e) {
            System.out.println(e);
        }
    }


    public double Norma() {
        double err = error[0];
        for (int i = 1; i < n + 1; i++)
            if (err < error[i])
                err = error[i];
        return err;
    }

    public double hh(){

        return h;

    }
}



public class Main {

    public static void main(String[] args) {

        Krayova_zadacha obb = new Krayova_zadacha(0);
        String er = "err";
        String res = "res";
        obb.cleanfile(er);
        obb.cleanfile(res);
            for (int i = 2; i < 100; i++) {
                //String res = "res";
                Krayova_zadacha ob = new Krayova_zadacha(i);
                ob.TridiagonalMatrixAlgorithmSolving(res);
                double error1 = ob.Norma();
                double h_1 = ob.hh();
               // System.out.println(" ============================================");
               // System.out.println(h_1 + "   " + error1);
                try {
                    FileWriter writer1 = new FileWriter("err.txt", true);
                    BufferedWriter err = new BufferedWriter(writer1);
                    err.write( h_1 + "               " + error1 + "\n");
                    err.close();
                }catch (IOException e) {
                    System.out.println(e);
                }
            }

       // new Donut();
    }
}
