#include <iostream>
#include <vector>
#include <fstream>
#include <iomanip>

using namespace std;

int main(){


    double XL = -10.0;
    double XR = 10.0;

    double YB = -5.0;
    double YT = 5.0;

    double XC = -2.0;
    double YC = 3.0;

    double a = 1.0;


    double XSL = XC - a/2.0;
    double XSR = XC + a/2.0;

    double YSB = YC - a/2.0;
    double YST = YC + a/2.0;

    int Nx1 = 5;
    int Nx2 = 11;
    int Nx3 = 21;

    int Ny1 = 31;
    int Ny2 = 11;
    int Ny3 = 11;

    int Nx = Nx1 + Nx2 + Nx3 - 2;
    int Ny = Ny1 + Ny2 + Ny3 - 2;


    double dx1 = (XSL-XL)/(Nx1-1);
    double dx2 = (XSR-XSL)/(Nx2-1);
    double dx3 = (XR-XSR)/(Nx3-1);

    double dy1 = (YSB-YB)/(Ny1-1);
    double dy2 = (YST-YSB)/(Ny2-1);
    double dy3 = (YT-YST)/(Ny3-1);


    vector<double> x(Nx);
    vector<double> y(Ny);

    for(int i=0; i<Nx1; i++) x[i]=XL + i*dx1;
    for(int i=1; i<Nx2; i++) x[Nx1-1+i]=XSL + i*dx2;
    for(int i=1; i<Nx3; i++) x[Nx1+Nx2-2+i]=XSR + i*dx3;

    for(int j=0; j<Ny1; j++) y[j]=YB + j*dy1;
    for(int j=1; j<Ny2; j++) y[Ny1-1+j]=YSB + j*dy2;
    for(int j=1; j<Ny3; j++) y[Ny1+Ny2-2+j]=YST + j*dy3;

    ofstream meshfile("mesh_4.dat");

    meshfile<< fixed << setprecision(8);
    for(int j=0; j<Ny; j++)
    {
        for(int i=0; i<Nx; i++){
            meshfile <<x[i]<<" "<<y[j]<<endl;
        }
    }

    meshfile.close();


    return 0;
}