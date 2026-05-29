#include <iostream>
#include <vector>
#include <fstream>
#include <iomanip>
#include <string>


using namespace std;

int main(){

    double L1, L2;
    double h1, h2;

    double XC, YC;
    double a;

    int Nx1, Nx2, Nx3;
    int Ny1, Ny2, Ny3;

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
        else if(label == "Nx2")input >> Nx2;
        else if(label == "Nx3")input >> Nx3;
        else if(label == "Ny1")input >> Ny1;
        else if(label == "Ny2")input >> Ny2;
        else if(label == "Ny3")input >> Ny3;
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


    vector<vector<bool>> solid(Nx,vector<bool>(Ny,false));

    for(int j=0; j<Ny; j++){
        for(int i=0; i<Nx; i++){
            if(x[i] > XSL && x[i] < XSR && y[j] > YSB && y[j] < YST){
                solid[i][j] = true;
            }
        }
    }
    





    ofstream meshfile("mesh_5.dat");

    meshfile << Nx << " " << Ny << endl;

    meshfile<< fixed << setprecision(8);
    for(int j=0; j<Ny; j++)
    {
        for(int i=0; i<Nx; i++){
            meshfile <<x[i]<<" "<<y[j]<<" "<<solid[i][j]<<endl;
        }
    }

    meshfile.close();


    return 0;
}