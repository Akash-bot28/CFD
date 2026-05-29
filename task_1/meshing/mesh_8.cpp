#include <iostream>
#include <vector>
#include <fstream>
#include <iomanip>
#include <string>
#include <cmath>

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
    int Nx1, Nxs, Nx2;
    int Ny1, Nys, Ny2;

    double alpha1, alpha2;
    double betax1, betaxs, betax2;
    double betay1, betays, betay2;
};

struct Mesh
{
    int Nx, Ny;

    vector<double> x;
    vector<double> y;
    vector<double> dx;
    vector<double> dy;

    vector<vector<bool>> solid;
};

void readInput( Geometry& g, MeshParameters& p ){

    string label;

    ifstream input("input.txt");

    while(input >> label)
    {
        if(label == "L1") input >> g.L1;
        else if(label == "L2") input >> g.L2;
        else if(label == "h1") input >> g.h1;
        else if(label == "h2") input >> g.h2;
        else if(label == "XC") input >> g.XC;
        else if(label == "YC") input >> g.YC;
        else if(label == "a") input >> g.a;
        else if(label == "Nx1") input >> p.Nx1;
        else if(label == "Nxs") input >> p.Nxs;
        else if(label == "Nx2") input >> p.Nx2;
        else if(label == "Ny1") input >> p.Ny1;
        else if(label == "Nys") input >> p.Nys;
        else if(label == "Ny2") input >> p.Ny2;
        else if(label == "alpha1") input >> p.alpha1;
        else if(label == "alpha2") input >> p.alpha2;
        else if(label == "betax1") input >> p.betax1;
        else if(label == "betaxs") input >> p.betaxs;
        else if(label == "betax2") input >> p.betax2;
        else if(label == "betay1") input >> p.betay1;
        else if(label == "betays") input >> p.betays;
        else if(label == "betay2") input >> p.betay2;
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

double x1(double zi, double L, double alpha, double beta){

    double A = (beta + 1.0)/(beta - 1.0);
    double exponent = (zi - alpha)/(1.0 - alpha);
    double R = pow(A, exponent);
    double xi=L*((beta + 2.0*alpha)*R - beta + 2.0*alpha)/((1.0 + 2.0*alpha)+ R);
    return xi;
}
double x2(double zi, double L, double beta){

    double A = (beta + 1.0)/(beta - 1.0);
    double exponent = (1-zi);
    double R = pow(A, exponent);
    double xi= L*((beta+1)-(beta-1)*R)/(1+R);
    return xi;
}


Mesh generateMesh(const Geometry& g, const MeshParameters& p){

    Mesh mesh;

    mesh.Nx = p.Nx1 + p.Nxs + p.Nx2 - 2;
    mesh.Ny = p.Ny1 + p.Nys + p.Ny2 - 2;

    mesh.x.resize(mesh.Nx);
    mesh.y.resize(mesh.Ny);
    mesh.dx.resize(mesh.Nx-1);
    mesh.dy.resize(mesh.Ny-1);

    for(int i=0; i<p.Nx1; i++){
        double zi =double(i)/(p.Nx1 - 1);
        double xi = x1(zi, g.L1, p.alpha1, p.betax1);
        mesh.x[i] = g.XL + xi;
    }

    for(int i=1; i<p.Nxs; i++){
        double zi = double (i)/(p.Nxs-1);
        double xi = x1(zi, g.a, p.alpha2, p.betaxs);
        mesh.x[p.Nx1-1+i] = g.XSL + xi;
    }

    for(int i=1; i<p.Nx2; i++){
        double zi = double (i)/(p.Nx2 - 1);
        double xi = x2(zi, g.L2, p.betax2);
        mesh.x[p.Nx1+p.Nxs-2+i] = g.XSR + xi;
    }
    for(int j=0; j<p.Ny1; j++){
        double zj =double(j)/(p.Ny1 - 1);
        double yj = x1(zj, g.h1, p.alpha1, p.betay1);
        mesh.y[j] = g.YB + yj;
    }

    for(int j=1; j<p.Nys; j++){
        double zj = double(j)/(p.Nys-1);
        double yj = x1(zj, g.a, p.alpha2, p.betays);
        mesh.y[p.Ny1-1+j] = g.YSB + yj;
    }

    for(int j=1; j<p.Ny2; j++){
        double zj = double(j)/(p.Ny2 - 1);
        double yj = x2(zj, g.h2, p.betay2);
        mesh.y[p.Ny1+p.Nys-2+j] = g.YST + yj;
    }

    for(int i=0; i<mesh.Nx-1; i++) mesh.dx[i] =mesh.x[i+1]-mesh.x[i];
    for(int j=0; j<mesh.Ny-1; j++) mesh.dy[j] =mesh.y[j+1]-mesh.y[j];

    return mesh;

}


void exportMesh(const Mesh& mesh, string filename){

    ofstream meshfile(filename);

    meshfile<< fixed << setprecision(8);

    for(int j=0; j<mesh.Ny; j++)
    {
        for(int i=0; i<mesh.Nx; i++){
            meshfile <<mesh.x[i]<<" "<<mesh.y[j]<<endl;
        }
    }

    meshfile.close();
}

void exportGridElement(const Mesh& mesh, string xfilename, string yfilename){

    ofstream dxfile(xfilename);
    ofstream dyfile(yfilename);

    dxfile<< fixed<< setprecision(8);
    dyfile<< fixed<< setprecision(8);

    for(int i=0; i<mesh.Nx-1; i++) dxfile<< mesh.x[i]<< " "<< mesh.dx[i]<< endl;
    for(int j=0; j<mesh.Ny-1; j++) dyfile<< mesh.y[j]<< " "<< mesh.dy[j]<< endl;

    dxfile.close();
    dyfile.close();
    
}

int main()
{
    Geometry g;

    MeshParameters p;

    Mesh mesh;

    readInput(g, p);

    buildGeometry(g);

    mesh = generateMesh(g,p);


    exportMesh(mesh,"mesh_8.dat");

    exportGridElement(mesh,"dx_8.dat","dy_8.dat");

    return 0;
}




