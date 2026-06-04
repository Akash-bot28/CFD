#include <iostream>
#include <vector>
#include <fstream>
#include <iomanip>
#include <string>
#include <cmath>

using namespace std;

struct MeshParameter{
    int Nx, Ny, N;
};

enum class NodeType{
    INTERIOR,   //m.nodeType.assign(p.N, NodeType::INTERIOR);
    OUTER_BOUNDARY, //i==1 i==Nx j==1 j==Ny
    BODY_BOUNDARY,
    BODY_INTERIOR
};

struct Mesh{
    vector<double> x;
    vector<double> y;
    vector<double> dx;
    vector<double> dy;

    vector<NodeType> nodeType;
};

struct Coefficient{
    vector<double> aS;
    vector<double> aW;
    vector<double> aP;
    vector<double> aE;
    vector<double> aN;
    vector<double> bP;
};

struct Field{
    vector<double> phi;
    //vector<double> psi;
};

void readInput(MeshParameter& p, Mesh& m){
    ifstream meshfile("mesh.dat");
    ifstream dxfile("dx.dat");
    ifstream dyfile("dy.dat");

}

void generateCoefficient(const MeshParameter& p,const Mesh& m,Coefficient &c){
    c.aS.assign(p.N,nan(""));
    c.aW.assign(p.N,nan(""));
    c.aP.assign(p.N,nan(""));
    c.aE.assign(p.N,nan(""));
    c.aN.assign(p.N,nan(""));
    c.bP.assign(p.N,nan(""));

    for(j=2; j<=p.Ny-1; j++){
        for(i=2; i<=p.Nx-1; i++){
            l = i + (j)*p.Nx;

            dxW = m.dx[i-1];
            dxE = m.dx[i];
            dyS = m.dy[j-1];
            dyN = m.dy[j];

            // coefficient calculation
            c.aS[l]=dyN/dyS;
            c.aW[l]=(dyN*(dyN+dyS))/(dxW*(dxE+dxW));
            c.aP[l]=-((dyN/dyS)+(dyN*(dyN+dyS))/(dxW*(dxE+dxW))+(dyN*(dyN+dyS))/(dxE*(dxE+dxW))+1);
            c.aE[l]=(dyN*(dyN+dyS))/(dxE*(dxE+dxW));
            c.aN[l]=1;
            c.bp[l]=0;

        }
    }

}

void DirichletBC(const MeshParameter& p,const Mesh& m,Coefficient& c){

    //left boundary bP<---(bP-aW*(phi-not+U*XL) & aW<---0
    //right boundary bP<---(bP-aE*(phi-not+U*XR) & aE<---0
    //bottom boundary bP---<(bP-aS*(phi-not+U*x) {where XL<x<XR) & aS<---0
    //top boundary bP<---(bP-aT*(phi-not+U*x) {where XL<x<XR) & aT--->0

    if(m.nodeType[l] != NodeType::INTERIOR) continue;
}

void initializePhi(onst MeshParameter& p,const Mesh& m,Field& f){
    f.phi.assign(p.N,0.0);
}

//void initializePsi(Field& f){}



int main(){

    MeshParameter p;
    Mesh m;
    Coefficient c;
    Field f;

    readInput(p,m);

    initializePhi(p,m,f);

    generateCoefficient(p,m,c);

    //DirichletBC(p,m,c);

    return 0;
}

