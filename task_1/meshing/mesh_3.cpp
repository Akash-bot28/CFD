#include <iostream>
#include <vector>
#include <fstream>
#include <iomanip>

using namespace std;

int main(){

    double XL = -10.0;
    double XM = -5.0;
    double XR = 10.0;

    double YB = -5.0;
    double YM = 3.0;
    double YT = 5.0;

    int Nx1 = 41;
    int Nx2 = 21;
    int Ny1 = 11;
    int Ny2 = 21;

    int Nx = Nx1 + Nx2 - 1;
    int Ny = Ny1 + Ny2 - 1;


    double dx1 = (XM-XL)/(Nx1-1);
    double dx2 = (XR-XM)/(Nx2-1);

    double dy1 = (YM-YB)/(Ny1-1);
    double dy2 = (YT-YM)/(Ny1-1);

    vector<double> x(Nx);
    vector<double> y(Ny);

    for(int i=0; i<Nx1; i++) x[i]=XL + i*dx1 ;
    for(int i=1; i<Nx2; i++) x[Nx1-1+i]=XM + i*dx2 ;
    
    for(int j=0; j<Ny1; j++) y[j]=YB + j*dy1;
    for(int j=0; j<Ny2; j++) y[Ny1-1+j]=YM + j*dy2;

    ofstream meshfile("mesh_3.dat");

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