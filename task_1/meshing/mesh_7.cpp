#include <iostream>
#include <vector>
#include <fstream>
#include <iomanip>
#include <cmath>

using namespace std;

double x1(double ei, double L, double alpha, double beta){

    double A = (beta + 1.0)/(beta - 1.0);
    double exponent = (ei - alpha)/(1.0 - alpha);
    double R = pow(A, exponent);
    double xi=L*((beta + 2.0*alpha)*R - beta + 2.0*alpha)/((1.0 + 2.0*alpha)+ R);
    return xi;
}

double x2(double ei, double L, double beta){

    double A = (beta + 1.0)/(beta - 1.0);
    double exponent = (1-ei);
    double R = pow(A, exponent);
    double xi= L*((beta+1)-(beta-1)*R)/(1+R);
    return xi;
}


int main(){

    double L1, L2;
    double h1, h2;

    double XC, YC;
    double a;

    int Nx1, Nxs, Nx2;
    int Ny1, Nys, Ny2;

    double alpha1, alphas;
    double beta1, betas, beta2;

    string label;

    ifstream input("input.txt");

    while(input >> label)
    {
        if(label == "L1") input >> L1;
        else if(label == "L2") input >> L2;
        else if(label == "h1")input >> h1;
        else if(label == "h2")input >> h2;
        else if(label == "XC")input >> XC;
        else if(label == "YC")input >> YC;
        else if(label == "a")input >> a;
        else if(label == "Nx1")input >> Nx1;
        else if(label == "Nxs")input >> Nxs;
        else if(label == "Nx2")input >> Nx2;
        else if(label == "Ny1")input >> Ny1;
        else if(label == "Nys")input >> Nys;
        else if(label == "Ny2")input >> Ny2;
        else if(label == "alpha1") input >> alpha1;
        else if(label == "alphas") input >> alphas;
        else if(label == "beta1") input >> beta1;
        else if(label == "betas") input >> betas;
        else if(label == "beta2") input >> beta2;

    }

    input.close();

    double XSL = XC - a/2.0;
    double XSR = XC + a/2.0;

    double YSB = YC - a/2.0;
    double YST = YC + a/2.0;

    double XL = XSL - L1;
    double XR = XSR + L2;

    double YB = YSB - h1;
    double YT = YST + h2;

    

    /*
    double dx1 = (XSL-XL)/(Nx1-1);
    double dx2 = (XSR-XSL)/(Nx2-1);
    double dx3 = (XR-XSR)/(Nx3-1);
    */

    double dy1 = (YSB-YB)/(Ny1-1);
    double dys = (YST-YSB)/(Nys-1);
    double dy2 = (YT-YST)/(Ny2-1);
    
    /*
    double r1 = 1.03;
    double r2 = 1.03;

    Nx1 = ceil(1.0+log(1.0+L1*(r1 - 1.0)/dxs)/log(r1));
    Nx2 = ceil(1.0+log(1.0+L1*(r2 - 1.0)/dxs)/log(r2));
    */

    int Nx = Nx1 + Nxs + Nx2 - 2;
    int Ny = Ny1 + Nys + Ny2 - 2;

    vector<double> x(Nx);
    vector<double> y(Ny);

    vector<double> dx(Nx-1);

    for(int j=0; j<Ny1; j++) y[j]=YB + j*dy1;
    for(int j=1; j<Nys; j++) y[Ny1-1+j]=YSB + j*dys;
    for(int j=1; j<Ny2; j++) y[Ny1+Nys-2+j]=YST + j*dy2;
    

    for(int i=0; i<Nx1; i++){
        double ei =double(i)/(Nx1 - 1);
        double xi = x1(ei, L1, alpha1, beta1);
        x[i] = xi;
    }

    for(int i=1; i<Nxs; i++){
        double ei = double (i-1)/(Nxs-1);
        double xi = x1(ei, a, alphas, betas);
        x[Nx1-1+i] = xi;
    }

    for(int i=1; i<Nx2; i++){
        double ei = double (i-1)/(Nx2 - 1);
        double xi = x2(ei, L2, beta2);
        x[Nx1+Nxs-2+i] = xi;
    }

    for(int i=0; i<Nx-1; i++){
        dx[i] = x[i+1] - x[i];
    }
    
    ofstream meshfile("mesh_7.dat");

    meshfile<< fixed << setprecision(8);
    for(int j=0; j<Ny; j++)
    {
        for(int i=0; i<Nx; i++){
            meshfile <<x[i]<<" "<<y[j]<<endl;
        }
    }

    meshfile.close();

    ofstream dxfile("dx_7.dat");

    dxfile<< fixed << setprecision(8);
    for(int i=0; i<Nx; i++){
        dxfile << i <<" "<<dx[i]<<endl;
    }


    return 0;
}