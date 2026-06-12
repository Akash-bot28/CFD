#include <iostream>
#include <vector>
#include <fstream>
#include <iomanip>
#include <string>
#include <cmath>

using namespace std;

struct Domain{
    double XL,XR,YB,YT;
};
struct Body{
    double XLs,XRs,YBs,YTs;
    int iL,iR,jB,jT;
    vector<vector<double>> mask;
};
struct MeshParameter{
    int Nx, Ny, N;
};
struct Mesh{
    vector<vector<double>> x;
    vector<vector<double>> y;
};
struct Coefficient{
    vector<vector<double>> aS;
    vector<vector<double>> aW;
    vector<vector<double>> aP;
    vector<vector<double>> aE;
    vector<vector<double>> aN;
    vector<vector<double>> bP;
};
struct Field{
    vector<vector<double>> phi;
};
struct SolverParameter{
    int itermax;
    double tolerance;
    double W;
    vector<double> error;
};
void readInput(Domain& d,Body& b,MeshParameter& p, Mesh& m, SolverParameter& s){
    ifstream file("mesh.dat");
    // Read mesh size
    file>>p.Nx>>p.Ny;
    p.N=p.Nx*p.Ny;
    // Allocate memory
    m.x.assign(p.Nx,vector<double>(p.Ny));
    m.y.assign(p.Nx,vector<double>(p.Ny));
    // Read coordinates
    for(int j=0;j<p.Ny;j++){
        for(int i=0;i<p.Nx;i++)
        {
            file >> m.x[i][j]>> m.y[i][j];
        }
    }
    file.close();
    string label;
    ifstream input("input.txt");
    while(input>> label){
        if(label == "XL") input >> d.XL;
        else if(label == "XR") input >> d.XR;
        else if(label == "YB") input >> d.YB;
        else if(label == "YT") input >> d.YT;
        else if(label == "XLs") input >> b.XLs;
        else if(label == "XRs") input >> b.XRs;
        else if(label == "YBs") input >> b.YBs;
        else if(label == "YTs") input >> b.YTs;
        else if(label == "solver_itermax") input >> s.itermax;
        else if(label == "solver_tolerance") input >> s.tolerance;
        else if(label == "W") input >> s.W;
    }
}

void generateCoefficient(const MeshParameter& p,Mesh& m,Coefficient& c){
    c.aS.assign(p.Nx, vector<double>(p.Ny, nan("")));
    c.aW.assign(p.Nx, vector<double>(p.Ny, nan("")));
    c.aP.assign(p.Nx, vector<double>(p.Ny, nan("")));
    c.aE.assign(p.Nx, vector<double>(p.Ny, nan("")));
    c.aN.assign(p.Nx, vector<double>(p.Ny, nan("")));
    c.bP.assign(p.Nx, vector<double>(p.Ny, nan("")));

    for(int j=1;j<=p.Ny-2;j++){
        for(int i=1;i<=p.Nx-2;i++){
            double dxW = m.x[i][j] - m.x[i-1][j];
            double dxE = m.x[i+1][j] - m.x[i][j];
            double dyS = m.y[i][j] - m.y[i][j-1];
            double dyN = m.y[i][j+1] - m.y[i][j];

            // coefficient calculation
            c.aS[i][j]=dyN/dyS;
            c.aW[i][j]=(dyN*(dyN+dyS))/(dxW*(dxE+dxW));
            c.aP[i][j]=-((dyN/dyS)+(dyN*(dyN+dyS))/(dxW*(dxE+dxW))+(dyN*(dyN+dyS))/(dxE*(dxE+dxW))+1);
            c.aE[i][j]=(dyN*(dyN+dyS))/(dxE*(dxE+dxW));
            c.aN[i][j]=1.0;
            c.bP[i][j]=0.0;
        }
    }
}

void initializePhi(const Domain& d, const MeshParameter& p,const Mesh& m,Field& f){
    double phi0 = 0.0;
    double U = 100.0;

    f.phi.assign(p.Nx,vector<double>(p.Ny,0.0));

    //double XL = m.x[0][0];
    //double XR = m.x[p.Nx-1][0];

    // Left boundary (inlet)
    for(int j=0; j<p.Ny; j++) f.phi[0][j] = U*d.XL + phi0;

    // Right boundary (outlet)
    for(int j=0; j<p.Ny; j++) f.phi[p.Nx-1][j] =U*d.XR + phi0;
    
    // Bottom boundary
    for(int i=0; i<p.Nx; i++) f.phi[i][0] = U*m.x[i][0] + phi0;
    
    // Top boundary
    for(int i=0; i<p.Nx; i++) f.phi[i][p.Ny-1] = U*m.x[i][p.Ny-1] + phi0;
}

void BCPhi(const MeshParameter& p,Coefficient& c,Field& f){
    // Adjacent Left boundary (Dirichlet)
    for(int j=1; j<=p.Ny-2; j++) {int i = 1; c.bP[i][j] = c.bP[i][j] - c.aW[i][j]*f.phi[i-1][j]; c.aW[i][j] = 0.0;}
    
    // Adjacent Right boundary (Neumann)(convective boundary condition) DO NOT APPLY 
    //for(int j=1; j<=p.Ny-2; j++) {int i = p.Nx-2; c.aP[i][j] = c.aP[i][j] + c.aE[i][j]; c.aE[i][j] = 0.0;}
    
    // Adjacent Right boundary (Dirichlet)
    for(int j=1; j<=p.Ny-2; j++) {int i = p.Nx-2; c.bP[i][j] = c.bP[i][j] - c.aE[i][j]*f.phi[i+1][j];c.aE[i][j] = 0.0;}
    
    // Adjacent Bottom boundary (Dirichlet)
    for(int i=1; i<=p.Nx-2; i++) {int j = 1; c.bP[i][j] = c.bP[i][j] - c.aS[i][j]*f.phi[i][j-1]; c.aS[i][j] = 0.0;}
   
    // Adjacent Top boundary (Dirichlet)
    for(int i=1; i<=p.Nx-2; i++) {int j = p.Ny-2; c.bP[i][j] = c.bP[i][j] - c.aN[i][j]*f.phi[i][j+1]; c.aN[i][j] = 0.0;}
}

void findBodyIndices(Body& b,const MeshParameter& p,const Mesh& m){
    double tol = 1e-5;
    for(int i=0;i<p.Nx;i++){
        if(abs(m.x[i][0]-b.XLs)<tol) b.iL=i;
        if(abs(m.x[i][0]-b.XRs)<tol) b.iR=i;
    }

    for(int j=0;j<p.Ny;j++){
        if(abs(m.y[0][j]-b.YBs)<tol) b.jB=j;
        if(abs(m.y[0][j]-b.YTs)<tol) b.jT=j;
    }
}

void buildMask(Body& b,const MeshParameter& p){
    b.mask.assign(p.Nx,vector<double>(p.Ny,1.0));

    for(int j=b.jB; j<=b.jT; j++){
        for(int i=b.iL; i<=b.iR; i++){
            b.mask[i][j]=0.0;
        }
    }
}

void BCBodyPhi(const Body& b,const MeshParameter& p,Coefficient& c,Field& f){
    // Left
    for(int j=b.jB ;j<=b.jT ;j++) {int i=(b.iL-1) ; c.aP[i][j] = c.aP[i][j] + c.aE[i][j];;c.aE[i][j] = 0.0;}

    // Right
    for(int j=b.jB ;j<=b.jT ;j++) {int i=(b.iR+1) ; c.aP[i][j] = c.aP[i][j] + c.aW[i][j];; c.aW[i][j] = 0.0;}

    // Bottom
    for(int i=b.iL ;i<=b.iR ;i++) {int j=(b.jB-1); c.aP[i][j] = c.aP[i][j] + c.aN[i][j];; c.aN[i][j] = 0.0;}

    //Top
    for(int i=b.iL ;i<=b.iR ;i++) {int j=(b.jT+1) ; c.aP[i][j] = c.aP[i][j] + c.aS[i][j];; c.aS[i][j] = 0.0;}
}
/*
void deactivateBodyNodes(const Body& b,Coefficient& c){
    for(int j=b.jB; j<=b.jT; j++){
        for(int i=b.iL; i<=b.iR; i++){
            c.aS[i][j]=0.0;
            c.aW[i][j]=0.0;
            c.aE[i][j]=0.0;
            c.aN[i][j]=0.0;
            c.aP[i][j]=1.0;
            c.bP[i][j]=0.0;
        }
    }
}
*/
double rowProduct(const Coefficient& c,const vector<vector<double>>& phi, int i, int j){
    return c.aS[i][j]*phi[i][j-1] + c.aW[i][j]*phi[i-1][j] + c.aP[i][j]*phi[i][j] + c.aE[i][j]*phi[i+1][j] + c.aN[i][j]*phi[i][j+1];
}

void SOR(SolverParameter& s,const Body& b,const MeshParameter& p,const Mesh& m, const Coefficient& c,Field& f){
    
    cout<< "\nsolving Phi..."<<endl;
    double rms;
    const int Ninterior =(p.Nx-2)*(p.Ny-2);

    int n;
    for(n=1; n<s.itermax; n++){
        for(int j=1; j<p.Ny-1; j++){
            for(int i=1; i<p.Nx-1; i++){

                //if(i>=b.iL && i<=b.iR && j>=b.jB && j<=b.jT) continue;

                double Aphi_P = rowProduct(c,f.phi,i,j);
                f.phi[i][j] += b.mask[i][j]*s.W*(c.bP[i][j]- Aphi_P)/c.aP[i][j];
            }
        }

        // Residual calculation
        int  active_nodes=0;
        double sumofsquares = 0.0;
        for(int j=1; j<p.Ny-1; j++){
            for(int i=1; i<p.Nx-1; i++){

                //if(i>=b.iL && i<=b.iR && j>=b.jB && j<=b.jT) continue;

                double Aphi_P = matVecProduct(c,f.phi,i,j);
                double residual = c.bP[i][j]- Aphi_P;
                sumofsquares += b.mask[i][j]*residual*residual; 
                active_nodes += b.mask[i][j];
            }
        }
        rms = sqrt(sumofsquares/active_nodes);
        s.error.push_back(rms);

        if(rms <= s.tolerance) break;
    }
    cout << "\nPHI SOR iterations = "<< n << endl;
    if(n == s.itermax) cout << "Stopped due to itermax." << endl;
    else cout << "Converged."<< endl;
}

void exportPhi(const MeshParameter& p,const Mesh& m,const Field& f,string filename){
    ofstream file(filename);
    file << fixed << setprecision(8);
    file << p.Nx << " " << p.Ny << endl;
    for(int j=0; j<p.Ny; j++){
        for(int i=0; i<p.Nx; i++){
            file<< m.x[i][j]<<" "<< m.y[i][j]<<" "<< f.phi[i][j] << endl;
        }
    }

    file.close();
}

void exportCoefficient(const MeshParameter& p,const Mesh& m,const Coefficient& c,string filename){
    ofstream file(filename);
    file << fixed << setprecision(8);
    file << p.Nx << " " << p.Ny << endl;
    for(int j=0; j<p.Ny; j++){
        for(int i=0; i<p.Nx; i++){
            file
            << i           << " "
            << j           << " "
            << m.x[i][j]   << " "
            << m.y[i][j]   << " "
            << c.aS[i][j]  << " "
            << c.aW[i][j]  << " "
            << c.aP[i][j]  << " "
            << c.aE[i][j]  << " "
            << c.aN[i][j]  << " "
            << c.bP[i][j]
            << endl;
        }
    }

    file.close();
}

void exportSolverData(const SolverParameter& s,const MeshParameter& p,const Mesh& m,const Coefficient& c,Field& f,string filename){
    ofstream file(filename);
    file<<fixed<<setprecision(8);
    for(int n=0; n<s.error.size();n++){
        file<<n<<"  "<<s.error[n]<<endl;
    }
}

int main(){

    Domain d;
    Body b;
    MeshParameter p;
    Mesh m;
    Coefficient c;
    Field f;
    SolverParameter s;

    readInput(d,b,p,m,s);

    initializePhi(d,p,m,f);

    //exportPhi(p,m,f,"phi field.dat");

    generateCoefficient(p,m,c);

    BCPhi(p,c,f);

    findBodyIndices(b,p,m);
    //cout << b.iL <<" "<< b.iR <<" "<< b.jB <<" "<< b.jT << endl;
    BCBodyPhi(b,p,c,f);
    buildMask(b,p);
    //deactivateBodyNodes(b,c);
    //exportCoefficient(p,m,c,"banded_sparse_matrix_.dat");

    //Gauss_Seidel(s,p,m,c,f);
    //Jacobi(s,p,m,c,f);
    SOR(s,b,p,m,c,f);
    
    //exportSolverData(s,p,m,c,f,"phi_error_converge.dat");
 
    exportPhi(p,m,f,"phi_solved.dat");

    return 0;
}