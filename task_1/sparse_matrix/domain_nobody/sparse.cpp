#include <iostream>
#include <vector>
#include <fstream>
#include <iomanip>
#include <string>
#include <cmath>

using namespace std;

enum class NodeType{
    INTERIOR,   //m.nodeType.assign(p.N, NodeType::INTERIOR);
    LEFT_BOUNDARY,       //i==1 i==Nx j==1 j==Ny
    RIGHT_BOUNDARY,
    BOTTOM_BOUNDARY,
    TOP_BOUNDARY
};

struct MeshParameter{
    int Nx, Ny, N;
};

struct Mesh{
    vector<double> x;
    vector<double> y;

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

struct SolverParameter{
    int itermax;
    double tolerance;
    vector<double> error;
};

void readInput(MeshParameter& p, Mesh& m, SolverParameter& s){
    ifstream meshfile("mesh.dat");

    // Read mesh size
    meshfile>>p.Nx>>p.Ny;
    p.N=p.Nx*p.Ny;

    // Allocate memory
    m.x.resize(p.N);
    m.y.resize(p.N);

    m.nodeType.assign(p.N, NodeType::INTERIOR);

    // Read coordinates
    for(int l=0; l<p.N; l++) meshfile >> m.x[l] >> m.y[l];
    
    meshfile.close();

    string label;
    ifstream input("input.txt");
    while(input>> label)
    {
        if(label == "itermax") input >> s.itermax;
        else if(label == "tolerance") input >> s.tolerance;
    }


}

void identifyNodes(const MeshParameter& p, Mesh& m){
    for(int j=0; j<p.Ny; j++)
    {
        for(int i=0; i<p.Nx; i++)
        {
            int l = i + (j)*p.Nx;

            if(i == 0) m.nodeType[l] = NodeType::LEFT_BOUNDARY;
            else if(i == p.Nx-1) m.nodeType[l] = NodeType::RIGHT_BOUNDARY;
            else if(j == 0) m.nodeType[l] = NodeType::BOTTOM_BOUNDARY;
            else if(j == p.Ny-1) m.nodeType[l] = NodeType::TOP_BOUNDARY;
            else m.nodeType[l] = NodeType::INTERIOR;
        }
        
    }    
}



void generateCoefficient(const MeshParameter& p,Mesh& m,Coefficient& c){
    c.aS.assign(p.N,nan(""));
    c.aW.assign(p.N,nan(""));
    c.aP.assign(p.N,nan(""));
    c.aE.assign(p.N,nan(""));
    c.aN.assign(p.N,nan(""));
    c.bP.assign(p.N,nan(""));

    for(int j=1;j<=p.Ny-2;j++)
    {
        for(int i=1;i<=p.Nx-2;i++)
        {
            int l=i+j*p.Nx;
            double dxW = m.x[l]-m.x[l-1];
            double dxE = m.x[l+1]-m.x[l];
            double dyS = m.y[l]-m.y[l-p.Nx];
            double dyN = m.y[l+p.Nx]-m.y[l];

            // coefficient calculation
            c.aS[l]=dyN/dyS;
            c.aW[l]=(dyN*(dyN+dyS))/(dxW*(dxE+dxW));
            c.aP[l]=-((dyN/dyS)+(dyN*(dyN+dyS))/(dxW*(dxE+dxW))+(dyN*(dyN+dyS))/(dxE*(dxE+dxW))+1);
            c.aE[l]=(dyN*(dyN+dyS))/(dxE*(dxE+dxW));
            c.aN[l]=1;
            c.bP[l]=0;
        }
    }

        
}
//void initializePsi(Field& f){}
void initializePhi(const MeshParameter& p,const Mesh& m,Field& f){
    double phi0=0.0;
    double U = 100;

    f.phi.assign(p.N,0.0);
    
    double XL = m.x[0];
    double XR = m.x[p.Nx-1];
    double YB = m.y[0];
    double YT = m.y[(p.Nx-1)*(p.Ny-1)];

    for(int l=0;l<p.N;l+=p.Nx) f.phi[l]=U*XL+phi0;               //inlet
    for(int l=0;l<p.Nx;l++) f.phi[l]=U*m.x[l]+phi0;              //bottom
    for(int l=(p.Ny-1)*p.Nx;l<p.N;l++) f.phi[l]=U*m.x[l]+phi0;   //top
    for(int l=(p.Nx-1);l<p.N;l+=p.Nx) f.phi[l]=U*XR+phi0;        //outlet
}

void BC(const MeshParameter& p,const Mesh& m,Coefficient& c,Field& f){

    for(int l=(1+p.Nx);l<=(1+(p.Ny-2)*p.Nx);l+=p.Nx)               {c.bP[l]=c.bP[l]-c.aW[l]*f.phi[l-1]; c.aW[l]=0;}    //inlet dirichlet

    //for(int l=((p.Nx-2)+p.Nx);l<=((p.Nx-2)+(p.Ny-2)*p.Nx);l+=p.Nx) {c.aP[l]=c.aP[l]+c.aE[l]; c.aE[l]=0;}              //outlet convective bc

    for(int l=((p.Nx-2)+p.Nx);l<=((p.Nx-2)+(p.Ny-2)*p.Nx);l+=p.Nx) {c.bP[l]=c.bP[l]-c.aE[l]*f.phi[l+1]; c.aE[l]=0;} //outlet dirichlet

    for(int l=(1+p.Nx);l<=((p.Nx-2)+p.Nx);l++)                     {c.bP[l]=c.bP[l]-c.aS[l]*f.phi[l-p.Nx]; c.aS[l]=0;} //bottom dirichlet
    for(int l=(1+(p.Ny-2)*p.Nx);l<=((p.Nx-2)+(p.Ny-2)*p.Nx);l++)   {c.bP[l]=c.bP[l]-c.aN[l]*f.phi[l+p.Nx]; c.aN[l]=0;} //top  dirichlett

    //for(int l=(1+p.Nx);l<=((p.Nx-2)+p.Nx);l++)                     {c.aP[l]=c.aP[l]+c.aS[l]; c.aS[l]=0;} //bottom  numen 
    //for(int l=(1+(p.Ny-2)*p.Nx);l<=((p.Nx-2)+(p.Ny-2)*p.Nx);l++)   {c.aP[l]=c.aP[l]+c.aN[l]; c.aN[l]=0;} //top     numen


}

void applyBoundaryValues(const MeshParameter& p,const Mesh& m,Field& f){
    double phi0=0.0;
    double U = 100;

    double XL = m.x[0];
    double XR = m.x[p.Nx-1];
    double YB = m.y[0];
    double YT = m.y[(p.Nx-1)*(p.Ny-1)];

    for(int l=0;l<p.N;l+=p.Nx) f.phi[l]=U*XL+phi0;
    for(int l=0;l<p.Nx;l++) f.phi[l]=U*m.x[l]+phi0;
    for(int l=(p.Ny-1)*p.Nx;l<p.N;l++) f.phi[l]=U*m.x[l]+phi0;
    //for(int l=(p.Nx-1);l<p.N;l+=p.Nx) f.phi[l-1]=f.phi[l];  // numen condition
    for(int l=(p.Nx-1);l<p.N;l+=p.Nx) f.phi[l]=U*XR+phi0;     // dirichlet condition
}




void Gauss_Seidel(SolverParameter& s,const MeshParameter& p,const Mesh& m,const Coefficient& c,Field& f){
    double rms;
    const int Ninterior = (p.Nx - 2) * (p.Ny - 2);
    
    int n;
    for(n =0; n < s.itermax; n++){

        double sumofsquares = 0.0;
        for(int j = 1; j <= p.Ny - 2; j++){
            for(int i = 1; i <= p.Nx - 2; i++){
                int l = i + j * p.Nx;
                // Gauss-Seidel update
                f.phi[l] = (c.bP[l] - c.aS[l]*f.phi[l-p.Nx] - c.aW[l]*f.phi[l-1] - c.aE[l]*f.phi[l+1] - c.aN[l]*f.phi[l+p.Nx]) / c.aP[l];
            }
        }

        applyBoundaryValues(p,m,f);

        for(int j = 1; j <= p.Ny - 2; j++){
            for(int i = 1; i <= p.Nx - 2; i++){
                int l = i + j * p.Nx;
                // Residual
                double r = c.bP[l] - c.aS[l]*f.phi[l-p.Nx] - c.aW[l]*f.phi[l-1] - c.aP[l]*f.phi[l] - c.aE[l]*f.phi[l+1] - c.aN[l]*f.phi[l+p.Nx]; 
                sumofsquares += r*r;
            }
        }
        
        rms = sqrt(sumofsquares / Ninterior);
        s.error.push_back(rms);

        if( rms<= s.tolerance) break;
    }

    cout << "\nSolver iterations = "<< n << endl;
    if(n == s.itermax)
        cout << "Stopped due to itermax." << endl;
    else
        cout << "Converged." << endl;
}

void Jacobi(SolverParameter& s,const MeshParameter& p,const Mesh& m,const Coefficient& c,Field& f){
    double rms;
    const int Ninterior = (p.Nx - 2) * (p.Ny - 2);

    vector<double> phiNew = f.phi;   // dummy storage

    int n;

    for(n =0; n < s.itermax; n++){
        double sumofsquares = 0.0;
        for(int j = 1; j <= p.Ny - 2; j++){
            for(int i = 1; i <= p.Nx - 2; i++){
                int l = i + j * p.Nx;
                // phi update
                phiNew[l] = (c.bP[l] - c.aS[l]*f.phi[l-p.Nx] - c.aW[l]*f.phi[l-1] - c.aE[l]*f.phi[l+1] - c.aN[l]*f.phi[l+p.Nx]) / c.aP[l];
            }
        }

        // Copy new solution
        for(int j=1; j<= p.Ny-2; j++)
        {
            for(int i=1; i<= p.Nx-2; i++)
            {
                int l = i + j*p.Nx;
                f.phi[l] = phiNew[l];
            }
        }

        applyBoundaryValues(p,m,f);

        // Residual calculation
        for(int j = 1; j <= p.Ny - 2; j++){
            for(int i = 1; i <= p.Nx - 2; i++){
                int l = i + j * p.Nx;
                // Residual
                double r = c.bP[l] - c.aS[l]*f.phi[l-p.Nx] - c.aW[l]*f.phi[l-1] - c.aP[l]*f.phi[l] - c.aE[l]*f.phi[l+1] - c.aN[l]*f.phi[l+p.Nx]; 
                sumofsquares += r*r;
            }
        }
        
        rms = sqrt(sumofsquares / Ninterior);
        s.error.push_back(rms);

        if( rms<= s.tolerance) break;

    }

    cout << "\nSolver iterations = "<< n << endl;
    if(n == s.itermax)
        cout << "Stopped due to itermax." << endl;
    else
        cout << "Converged." << endl;

}



void exportPhi(const MeshParameter& p,const Mesh& m,const Field& f,string filename){
    ofstream file(filename);
    file<<fixed<<setprecision(8);

    file<<p.Nx<<"  "<<p.Ny<<endl;
    for(int l=0; l<p.N; l++) file << m.x[l]<<" "<<m.y[l]<<" "<<f.phi[l]<<endl;

    file.close();
}

void exportCoefficient(const MeshParameter& p,Mesh& m,Coefficient& c,string filename){
    ofstream file(filename);
    file<< fixed << setprecision(8);
    for(int l=0; l<p.N; l++)
    {
        file << l<<"  "<< m.x[l]<<" "<<m.y[l]<<" "<< c.aS[l]<<"  "<< c.aW[l]<<"  "<< c.aP[l]<<"  "<< c.aE[l]<<"  "<< c.aN[l]<< endl;
       
    }
}

void exportSolverData(const SolverParameter& s,const MeshParameter& p,const Mesh& m,const Coefficient& c,Field& f,string filename){
    ofstream file(filename);
    file<<fixed<<setprecision(8);
    for(int n=0; n<s.error.size();n++){
        file<<n<<"  "<<s.error[n]<<endl;
    }
}





int main(){

    MeshParameter p;
    Mesh m;
    Coefficient c;
    Field f;
    SolverParameter s;

    readInput(p,m,s);

    identifyNodes(p,m);

    initializePhi(p,m,f);

    exportPhi(p,m,f,"phi field.dat");

    generateCoefficient(p,m,c);

    BC(p,m,c,f);

    exportCoefficient(p,m,c,"banded_sparse_matrix_.dat");

    Gauss_Seidel(s,p,m,c,f);
    //Jacobi(s,p,m,c,f);
    
    exportSolverData(s,p,m,c,f,"error_converge.dat");
 
    exportPhi(p,m,f,"phi_solved.dat");

    return 0;
}

