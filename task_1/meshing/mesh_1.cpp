#include <iostream>
#include <vector>
#include <fstream>
#include <iomanip>

using namespace std;

int main(){

    double XL=-10.0;
    double XR=+10.0;

    double YT=+5.0;
    double YB=-5.0;

    int Nx=61; 
    int Ny=21;




    double dx=(XR-XL)/(Nx-1);
    double dy=(YT-YB)/(Ny-1);

    vector<double> x(Nx);
    vector<double> y(Ny);

    for(int i = 0; i < Nx; i++){
        x[i] = XL + i*dx;
    }
    for(int j = 0; j < Ny; j++){
        y[j] = YB + j*dy;
    }

    ofstream meshfile("mesh_1.dat");

    meshfile << fixed << setprecision(8);

    for(int j = 0; j < Ny; j++)
    {
        for(int i = 0; i < Nx; i++)
        {
            double x = XL + i*dx;
            double y = YB + j*dy;

            meshfile << x << " " << y << endl;
        }
    }

    meshfile.close();



    return 0;
}