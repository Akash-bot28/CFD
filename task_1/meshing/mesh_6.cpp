#include <iostream>
#include <vector>
#include <fstream>
#include <iomanip>
#include <string>

using namespace std;

struct Geometry
{
    double L1, L2;
    double h1, h2;

    double XC, YC;
    double a;

    double XL, XR;
    double YB, YT;

    double XSL, XSR;
    double YSB, YST;
};

struct MeshParameters
{
    int Nx1, Nx2, Nx3;
    int Ny1, Ny2, Ny3;
};

struct Mesh
{
    int Nx, Ny;

    vector<double> x;
    vector<double> y;

    vector<vector<bool>> solid;
};


void readInput( Geometry& g, MeshParameters& p ){
    
    ifstream input("input.txt");

    string label;

    while(input >> label){

        if(label == "L1") input >> g.L1;
        else if(label == "L2") input >> g.L2;
        else if(label == "h1") input >> g.h1;
        else if(label == "h2") input >> g.h2;
        else if(label == "XC") input >> g.XC;
        else if(label == "YC") input >> g.YC;
        else if(label == "a") input >> g.a;
        else if(label == "Nx1") input >> p.Nx1;
        else if(label == "Nx2") input >> p.Nx2;
        else if(label == "Nx3") input >> p.Nx3;
        else if(label == "Ny1") input >> p.Ny1;
        else if(label == "Ny2") input >> p.Ny2;
        else if(label == "Ny3") input >> p.Ny3;
    }
    input.close();

}

void buildGeometry(Geometry& g){

    g.XSL = g.XC - g.a/2.0;
    g.XSR = g.XC + g.a/2.0;

    g.YSB = g.YC - g.a/2.0;
    g.YST = g.YC + g.a/2.0;

    g.XL = g.XSL - g.L1;
    g.XR = g.XSR + g.L2;

    g.YB = g.YSB - g.h1;
    g.YT = g.YST + g.h2;
}

Mesh generateMesh(const Geometry& g, const MeshParameters& p){

    Mesh mesh;

    mesh.Nx = p.Nx1 + p.Nx2 + p.Nx3 - 2;
    mesh.Ny = p.Ny1 + p.Ny2 + p.Ny3 - 2;

    mesh.x.resize(mesh.Nx);
    mesh.y.resize(mesh.Ny);

    double dx1 = (g.XSL-g.XL)/(p.Nx1-1);
    double dx2 = (g.XSR-g.XSL)/(p.Nx2-1);
    double dx3 = (g.XR-g.XSR)/(p.Nx3-1);

    double dy1 = (g.YSB-g.YB)/(p.Ny1-1);
    double dy2 = (g.YST-g.YSB)/(p.Ny2-1);
    double dy3 = (g.YT-g.YST)/(p.Ny3-1);


    for(int i=0; i<p.Nx1; i++) mesh.x[i]=g.XL + i*dx1;
    for(int i=1; i<p.Nx2; i++) mesh.x[p.Nx1-1+i]=g.XSL + i*dx2;
    for(int i=1; i<p.Nx3; i++) mesh.x[p.Nx1+p.Nx2-2+i]=g.XSR + i*dx3;

    for(int j=0; j<p.Ny1; j++) mesh.y[j]=g.YB + j*dy1;
    for(int j=1; j<p.Ny2; j++) mesh.y[p.Ny1-1+j]=g.YSB + j*dy2;
    for(int j=1; j<p.Ny3; j++) mesh.y[p.Ny1+p.Ny2-2+j]=g.YST + j*dy3;

    return mesh;
}

void solidBody(Mesh& mesh, const Geometry& g){

    mesh.solid.resize(mesh.Nx, vector<bool>(mesh.Ny, false));

    for(int j=0; j<mesh.Ny; j++){
        for(int i=0; i<mesh.Nx; i++){
            if(mesh.x[i] > g.XSL && mesh.x[i] < g.XSR && mesh.y[j] > g.YSB && mesh.y[j] < g.YST){
                mesh.solid[i][j] = true;
            }
        }
    }

}

void exportMesh(const Mesh& mesh, string filename){

    ofstream meshfile(filename);

    meshfile << mesh.Nx << " " << mesh.Ny << endl;
    meshfile<< fixed << setprecision(8);

    for(int j=0; j<mesh.Ny; j++)
    {
        for(int i=0; i<mesh.Nx; i++){
            meshfile <<mesh.x[i]<<" "<<mesh.y[j]<<" "<<mesh.solid[i][j]<<endl;
        }
    }

    meshfile.close();
}

int main()
{
    Geometry g;

    MeshParameters p;

    Mesh mesh;

    readInput(g, p);

    buildGeometry(g);

    mesh = generateMesh(g,p);

    solidBody(mesh, g);

    exportMesh(mesh,"mesh.dat");


    return 0;
}