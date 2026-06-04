#include <iostream>
#include <vector>
#include <fstream>
#include <iomanip>
#include <string>
#include <cmath>

using namespace std;

struct Geometry{
    double XL, XR;
    double YB, YT;
    double L, H;
};

struct MeshParameters{
    int Nx;
    int Ny;
    double alpha1, alpha2;
    double betax;
    double betay;
};

struct Mesh{
    int Nx, Ny;
    vector<double> x;
    vector<double> y;
    vector<double> dx;
    vector<double> dy;
};

void readInput( Geometry& g, MeshParameters& p ){
    string label;
    ifstream input("input.txt");
    while(input >> label)
    {
        if(label == "XL") input >> g.XL;
        else if(label == "XR") input >> g.XR;
        else if(label == "YB") input >> g.YB;
        else if(label == "YT") input >> g.YT;
        else if(label == "Nx") input >> p.Nx;
        else if(label == "Ny") input >> p.Ny;
        else if(label == "betax") input >> p.betax;
        else if(label == "betay") input >> p.betay;
    }
    input.close();
}

void buildGeometry( Geometry& g){
    g.L=g.XR-g.XL;
    g.H=g.YT-g.YB;
}

double stretchLeft(double zi, double L, double beta){
    double alpha = 0.0;
    double A = (beta + 1.0)/(beta - 1.0);
    double exponent = (zi - alpha)/(1.0 - alpha);
    double R = pow(A, exponent);
    double xi=L*((beta + 2.0*alpha)*R - beta + 2.0*alpha)/((1.0 + 2.0*alpha)*(1+ R));
    return xi;
}

double CompressBoth(double zi, double L, double beta){
    double alpha = 0.5;
    double A = (beta + 1.0)/(beta - 1.0);
    double exponent = (zi - alpha)/(1.0 - alpha);
    double R = pow(A, exponent);
    double xi=L*((beta + 2.0*alpha)*R - beta + 2.0*alpha)/((1.0 + 2.0*alpha)*(1+ R));
    return xi;
}

double stretchRight(double zi, double L, double beta){
    double A = (beta + 1.0)/(beta - 1.0);
    double exponent = (1-zi);
    double R = pow(A, exponent);
    double xi= L*((beta+1)-(beta-1)*R)/(1+R);
    return xi;
}


Mesh generateMesh(const Geometry& g, const MeshParameters& p){

    Mesh mesh;

    mesh.Nx = p.Nx;
    mesh.Ny = p.Ny;

    mesh.x.resize(mesh.Nx);
    mesh.y.resize(mesh.Ny);
    mesh.dx.resize(mesh.Nx-1);
    mesh.dy.resize(mesh.Ny-1);

    for(int i=0; i<p.Nx; i++){
        double zi =double(i)/(p.Nx - 1);
        double xi = CompressBoth(zi, g.L,p.betax);
        mesh.x[i] = g.XL + xi;
    }

    for(int j=0; j<p.Ny; j++){
        double zj =double(j)/(p.Ny - 1);
        double yj = CompressBoth(zj, g.H, p.betay);
        mesh.y[j] = g.YB + yj;
    }

    for(int i=0; i<mesh.Nx-1; i++) mesh.dx[i] =mesh.x[i+1]-mesh.x[i];
    for(int j=0; j<mesh.Ny-1; j++) mesh.dy[j] =mesh.y[j+1]-mesh.y[j];

    return mesh;

}

void exportMesh(const Mesh& mesh, string filename){
    ofstream file(filename);
    file<< fixed << setprecision(8);
    for(int j=0; j<mesh.Ny; j++){
        for(int i=0; i<mesh.Nx; i++){
            file <<mesh.x[i]<<" "<<mesh.y[j]<<endl;
        }
    }
    file.close();
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

    exportMesh(mesh,"mesh.dat");

    exportGridElement(mesh,"dx.dat","dy.dat");

    return 0;
}




