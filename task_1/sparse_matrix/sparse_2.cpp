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
    //vector<vector<double>> psi;
};

struct SolverParameter{
    int itermax;
    double tolerance;
    vector<double> error;
};

void readInput(MeshParameter& p, Mesh& m, SolverParameter& s){
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
    while(input>> label)
    {
        if(label == "itermax") input >> s.itermax;
        else if(label == "tolerance") input >> s.tolerance;
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

void initializePhi(const MeshParameter& p,const Mesh& m,Field& f){
    double phi0 = 0.0;
    double U = 100.0;

    f.phi.assign(p.Nx,vector<double>(p.Ny,0.0));

    double XL = m.x[0][0];
    double XR = m.x[p.Nx-1][0];

    // Left boundary (inlet)
    for(int j=0; j<p.Ny; j++) f.phi[0][j] = U*XL + phi0;

    // Right boundary (outlet)
    for(int j=0; j<p.Ny; j++) f.phi[p.Nx-1][j] =U*XR + phi0;
    
    // Bottom boundary
    for(int i=0; i<p.Nx; i++) f.phi[i][0] = U*m.x[i][0] + phi0;
    
    // Top boundary
    for(int i=0; i<p.Nx; i++) f.phi[i][p.Ny-1] = U*m.x[i][p.Ny-1] + phi0;
}

// void initializePsi(const MeshParameter& p,const Mesh& m,Field& f){
//     double psi0 = 0.0;
//     double U = 100.0;

//     f.psi.assign(p.Nx,vector<double>(p.Ny,0.0));

//     double YB = m.y[0][0];
//     double YT = m.y[0][p.Ny-1];

//     // Left boundary (inlet)
//     for(int j=0; j<p.Ny; j++) f.psi[0][j] = U*m.y[0][j] + psi0;

//     // Right boundary (outlet)
//     for(int j=0; j<p.Ny; j++) f.psi[p.Nx-1][j] = U*m.y[p.Nx-1][j] + psi0;

//     // Bottom boundary
//     for(int i=0; i<p.Nx; i++) f.psi[i][0] = U*YB + psi0;

//     // Top boundary
//     for(int i=0; i<p.Nx; i++) f.psi[i][p.Ny-1] = U*YT + psi0;

// }

void applyBoundaryValuesPhi(const MeshParameter& p,const Mesh& m,Field& f){
    double phi0 = 0.0;
    double U    = 100.0;

    double XL = m.x[0][0];
    double XR = m.x[p.Nx-1][0];

    // Left boundary (inlet)
    for(int j=0; j<p.Ny; j++) f.phi[0][j] = U*XL + phi0;

    // Right boundary (Dirichlet)
    for(int j=0; j<p.Ny; j++) f.phi[p.Nx-1][j] = U*XR + phi0;        //Dirichlet
    //for(int j=0; j<p.Ny; j++) f.phi[p.Nx-1][j] = f.phi[p.Nx-2][j]; //Neumann
    
    // Bottom boundary
    for(int i=0; i<p.Nx; i++) f.phi[i][0] = U*m.x[i][0] + phi0;
    
    // Top boundary
    for(int i=0; i<p.Nx; i++) f.phi[i][p.Ny-1] = U*m.x[i][p.Ny-1] + phi0;
}

// void applyBoundaryValuesPsi(const MeshParameter& p,const Mesh& m,Field& f){
//     double psi0 = 0.0;
//     double U = 100.0;

//     double YB = m.y[0][0];
//     double YT = m.y[0][p.Ny-1];

//     // Left boundary (inlet)
//     for(int j=0; j<p.Ny; j++) f.psi[0][j] = U*m.y[0][j] + psi0;

//     // Right boundary (outlet)
//     for(int j=0; j<p.Ny; j++) f.psi[p.Nx-1][j] = U*m.y[p.Nx-1][j] + psi0;

//     // Bottom boundary
//     for(int i=0; i<p.Nx; i++) f.psi[i][0] = U*YB + psi0;

//     // Top boundary
//     for(int i=0; i<p.Nx; i++) f.psi[i][p.Ny-1] = U*YT + psi0;
// }


void BC(const MeshParameter& p,const Mesh& m,Coefficient& c,Field& f){
    // Left boundary (Dirichlet)
    for(int j=1; j<p.Ny-1; j++) {int i = 1; c.bP[i][j] = c.bP[i][j] - c.aW[i][j]*f.phi[i-1][j]; c.aW[i][j] = 0.0;}
    
    // Right boundary (Neumann)
    for(int j=1; j<p.Ny-1; j++) {int i = p.Nx-2; c.aP[i][j] = c.aP[i][j] + c.aE[i][j]; c.aE[i][j] = 0.0;}
    
    // Right boundary (Dirichlet)
    //for(int j=1; j<p.Ny-1; j++) {int i = p.Nx-2; c.bP[i][j] = c.bP[i][j] - c.aE[i][j]*f.phi[i+1][j];c.aE[i][j] = 0.0;}
    
    // Bottom boundary (Dirichlet)
    for(int i=1; i<p.Nx-1; i++) {int j = 1; c.bP[i][j] = c.bP[i][j] - c.aS[i][j]*f.phi[i][j-1]; c.aS[i][j] = 0.0;}
   
    // Top boundary (Dirichlet)
    for(int i=1; i<p.Nx-1; i++) {int j = p.Ny-2; c.bP[i][j] = c.bP[i][j] - c.aN[i][j]*f.phi[i][j+1]; c.aN[i][j] = 0.0;}
   
    // Bottom boundary (Neumann)
    //for(int i=1; i<p.Nx-1; i++) {int j = 1; c.aP[i][j] = c.aP[i][j] + c.aS[i][j]; c.aS[i][j] = 0.0;}
   
    // Top boundary (Neumann)
    //for(int i=1; i<p.Nx-1; i++) {int j = p.Ny-2; c.aP[i][j] = c.aP[i][j] + c.aN[i][j]; c.aN[i][j] = 0.0;}
   
    
}

// Vector<double> matVecMultiply(const vector<vector<double>>& A,const vector<double>& x){
//     int rows = A.size();
//     int cols = A[0].size();

//     vector<double> y(rows,0.0);
//     for(int i=0; i<rows; i++){
//         for(int j=0; j<cols; j++){
//             y[i]=Yy[i]+A[i][j]*x[j];
//         }
//     }
//     return y;
// }

void Gauss_Seidel(SolverParameter& s,const MeshParameter& p,const Mesh& m,const Coefficient& c,Field& f){
    double rms;
    const int Ninterior =(p.Nx-2)*(p.Ny-2);

    int n;
    for(n=0; n<s.itermax; n++){
        // Gauss-Seidel sweep
        for(int j=1; j<p.Ny-1; j++){
            for(int i=1; i<p.Nx-1; i++){
                f.phi[i][j] = (c.bP[i][j]- c.aS[i][j]*f.phi[i][j-1]- c.aW[i][j]*f.phi[i-1][j]- c.aE[i][j]*f.phi[i+1][j]- c.aN[i][j]*f.phi[i][j+1])/c.aP[i][j];
            }
        }
        
        // Residual calculation
        double sumofsquares = 0.0;
        for(int j=1; j<p.Ny-1; j++){
            for(int i=1; i<p.Nx-1; i++){
                double r = c.bP[i][j] - c.aS[i][j]*f.phi[i][j-1] - c.aW[i][j]*f.phi[i-1][j] - c.aP[i][j]*f.phi[i][j] - c.aE[i][j]*f.phi[i+1][j] - c.aN[i][j]*f.phi[i][j+1];
                sumofsquares += r*r;
            }
        }
        rms =sqrt(sumofsquares/Ninterior);
        s.error.push_back(rms);

        if(rms <= s.tolerance) break;

        applyBoundaryValuesPhi(p,m,f);

    }

    cout << "\nSolver iterations = "<< n << endl;
    if(n == s.itermax) cout << "Stopped due to itermax." << endl;
    else cout << "Converged."<< endl;
}

void Jacobi(SolverParameter& s,const MeshParameter& p,const Mesh& m,const Coefficient& c,Field& f){
    double rms;
    const int Ninterior =(p.Nx-2)*(p.Ny-2);

    vector<vector<double>> phiNew = f.phi;

    int n;
    for(n=0; n<s.itermax; n++){
        // Jacobi update
        for(int j=1; j<p.Ny-1; j++){
            for(int i=1; i<p.Nx-1; i++){
                phiNew[i][j] = (c.bP[i][j]- c.aS[i][j]*f.phi[i][j-1]- c.aW[i][j]*f.phi[i-1][j]- c.aE[i][j]*f.phi[i+1][j]- c.aN[i][j]*f.phi[i][j+1])/c.aP[i][j];
            }
        }

        // Copy new solution

        for(int j=1; j<p.Ny-1; j++){
            for(int i=1; i<p.Nx-1; i++){
                f.phi[i][j] = phiNew[i][j];
            }
        }

       

        // Residual calculation
        double sumofsquares = 0.0;

        for(int j=1; j<p.Ny-1; j++){
            for(int i=1; i<p.Nx-1; i++){
                double residual = c.bP[i][j]- c.aS[i][j]*f.phi[i][j-1]- c.aW[i][j]*f.phi[i-1][j]- c.aP[i][j]*f.phi[i][j]- c.aE[i][j]*f.phi[i+1][j]- c.aN[i][j]*f.phi[i][j+1];

                sumofsquares += residual*residual;
            }
        }

        rms = sqrt(sumofsquares/Ninterior);
        s.error.push_back(rms);

        if(rms <= s.tolerance) break;

        applyBoundaryValuesPhi(p,m,f);

    }

    cout << "\nSolver iterations = "<< n << endl;

    if(n == s.itermax)cout << "Stopped due to itermax."<< endl;
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

void exportPsi(const MeshParameter& p,const Mesh& m,const Field& f,string filename){
    ofstream file(filename);

    file<<fixed<<setprecision(8);
    file<<p.Nx<<" "<<p.Ny<<endl;

    for(int j=0;j<p.Ny;j++){
        for(int i=0;i<p.Nx;i++){
            file<<m.x[i][j]<<" "<<m.y[i][j]<<" "<<f.psi[i][j]<<endl;
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

    MeshParameter p;
    Mesh m;
    Coefficient c;
    Field f;
    SolverParameter s;

    readInput(p,m,s);

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