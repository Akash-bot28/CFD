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

                double Aphi_P = rowProduct(c,f.phi,i,j);
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

// Computes the regular matrix-vector product: A * phi
vector<vector<double>> matVecProduct(const MeshParameter& p,const Body& b,const Coefficient& c,const vector<vector<double>>& x){
    vector<vector<double>> y(p.Nx,vector<double>(p.Ny,0.0));
    for(int j=1; j<p.Ny-1; j++){
        for(int i=1; i<p.Nx-1; i++){
            double row =c.aS[i][j]*x[i][j-1]+ c.aW[i][j]*x[i-1][j]+ c.aP[i][j]*x[i][j]+ c.aE[i][j]*x[i+1][j]+ c.aN[i][j]*x[i][j+1];
            y[i][j] = row*b.mask[i][j];
        }
    }
    return y;
}
// Computes the inner product (dot product) of two 2D grid fields over active nodes
double dotProduct(const Body& b,const MeshParameter& p,const vector<vector<double>>& x,const vector<vector<double>>& y){
    double sum = 0.0;
    for(int j=1; j<p.Ny-1; j++){
        for(int i=1; i<p.Nx-1; i++){
            sum += b.mask[i][j]*x[i][j]*y[i][j];
        }
    }
    return sum;
}

double rmsNorm(const Body& b,const MeshParameter& p,const vector<vector<double>>& r){
    double sum = 0.0;
    double fluid_nodes = 0.0;
    for(int j=1; j<p.Ny-1; j++){
        for(int i=1; i<p.Nx-1; i++){
            sum += b.mask[i][j]*r[i][j]*r[i][j];
            fluid_nodes += b.mask[i][j];
        }
    }
    return sqrt(sum/fluid_nodes);
}
//Computes the transposed matrix-vector product: A^T * phi
vector<vector<double>> matVecTransposeProduct(const Body& b,const MeshParameter& p,const Coefficient& c,const vector<vector<double>>& x){
    vector<vector<double>> ATx(p.Nx,vector<double>(p.Ny,0.0));
    for(int j=1; j<p.Ny-1; j++){
        for(int i=1; i<p.Nx-1; i++){
            //ATx[i][j] =c.aS[i][j+1]*x[i][j+1]+c.aW[i+1][j]*x[i+1][j]+c.aP[i][j]*x[i][j]+c.aE[i-1][j]*x[i-1][j]+c.aN[i][j-1]*x[i][j-1];
            ATx[i][j]   += b.mask[i][j]*c.aP[i][j] * x[i][j];
            ATx[i][j-1] += b.mask[i][j]*c.aS[i][j] * x[i][j];
            ATx[i-1][j] += b.mask[i][j]*c.aW[i][j] * x[i][j];
            ATx[i+1][j] += b.mask[i][j]*c.aE[i][j] * x[i][j];
            ATx[i][j+1] += b.mask[i][j]*c.aN[i][j] * x[i][j];

            // ATx[i][j]   += c.aP[i][j] * x[i][j];
            // ATx[i][j-1] += c.aS[i][j] * x[i][j];
            // ATx[i-1][j] += c.aW[i][j] * x[i][j];
            // ATx[i+1][j] += c.aE[i][j] * x[i][j];
            // ATx[i][j+1] += c.aN[i][j] * x[i][j];
        }
    }
    // Clean up boundary/masked items leakages from distribution
    for (int j = 1; j < p.Ny - 1; j++) {
        for (int i = 1; i < p.Nx - 1; i++) {
            ATx[i][j] *= b.mask[i][j];
        }
    }
    return ATx;
}

void BiCG(SolverParameter& s,const Body& b,const MeshParameter& p,const Mesh& m, const Coefficient& c,Field& f){
    cout<<"\n solving phi"<<endl;
    vector<vector<double>> r(p.Nx, vector<double>(p.Ny,0.0));
    vector<vector<double>> rHat(p.Nx, vector<double>(p.Ny,0.0));
    vector<vector<double>> pVec(p.Nx, vector<double>(p.Ny,0.0));
    vector<vector<double>> pHat(p.Nx, vector<double>(p.Ny,0.0));
    vector<vector<double>> Ax(p.Nx, vector<double>(p.Ny,0.0));
    vector<vector<double>> Ap(p.Nx, vector<double>(p.Ny,0.0));
    vector<vector<double>> ATpHat(p.Nx, vector<double>(p.Ny,0.0));

    // reset direction every 1000 steps
    const int restart_interval = 1000;

    
    Ax = matVecProduct(p,b,c,f.phi);
    for(int j = 1; j < p.Ny - 1; j++){
        for(int i = 1; i < p.Nx - 1; i++){
            if (b.mask[i][j] > 0.0) {
                r[i][j] = c.bP[i][j] - Ax[i][j]; //Compute initial residual r = b - A*x
                rHat[i][j] = r[i][j]; // Choose rHat_0 = r_0
                pVec[i][j] = r[i][j];  // initial searsh direction
                pHat[i][j] = rHat[i][j];  // initial searsh direction
            }
        }
    }

    double rho_old = dotProduct(b,p,rHat,r);
    double rms = 0.0;
    int n;

    for(n=1; n<s.itermax; n++){

        if(n% restart_interval == 0){}
            Ap = matVecProduct(p,b,c,pVec);

            //if (abs(dotProduct(b,p,pHat,Ap)) < 1e-15) {cout << "BiCG Breakdown: Denominator near zero!" << endl;break;}
            double alpha = rho_old/dotProduct(b,p,pHat,Ap); // compute step length

            // Update solution and residual
            for (int j = 1; j < p.Ny - 1; j++) {
                for (int i = 1; i < p.Nx - 1; i++) {
                    if (b.mask[i][j] > 0.0) {
                        f.phi[i][j] += alpha * pVec[i][j];
                        r[i][j]     -= alpha * Ap[i][j];
                    }
                }
            }

            rms = rmsNorm(b, p, r);
            s.error.push_back(rms);
            if (rms <= s.tolerance) break;

            // Update shadow residual  
            ATpHat = matVecTransposeProduct(b,p,c,pHat);
            for (int j = 1; j < p.Ny - 1; j++) {
                for (int i = 1; i < p.Nx - 1; i++) {
                    if (b.mask[i][j] > 0.0) {
                        rHat[i][j] -= alpha * ATpHat[i][j];
                    }
                }
            }

            // Calculate Beta
            double rho_new = dotProduct(b,p,rHat,r);
            //if (abs(rho_old) < 1e-15) { cout << "BiCG Breakdown: rho_old near zero!" << endl; break;}
            double beta = rho_new / rho_old;
            rho_old = rho_new;

            // Update search directions
            for (int j = 1; j < p.Ny - 1; j++) {
                for (int i = 1; i < p.Nx - 1; i++) {
                    if (b.mask[i][j] > 0.0) {
                        pVec[i][j] = r[i][j] + beta * pVec[i][j];
                        pHat[i][j] = rHat[i][j] + beta * pHat[i][j];
                    }
                }
            }

            if(n%100==0) cout<< "iter = " << n<< " rms = " << rms<< endl;
        }
    

    cout << "\nPHI BiCG iterations = " << n << endl;
    if (n == s.itermax) cout << "Stopped due to itermax." << endl;
    else cout << "Converged. Final RMS: " << rms << endl;
    
}

void Bi(SolverParameter& s, const Body& b, const MeshParameter& p, const Mesh& m, const Coefficient& c, Field& f){
    cout << "\nSolving phi using BiCG with Restarts..." << endl;
    
    vector<vector<double>> r(p.Nx, vector<double>(p.Ny, 0.0));
    vector<vector<double>> rHat(p.Nx, vector<double>(p.Ny, 0.0));
    vector<vector<double>> pVec(p.Nx, vector<double>(p.Ny, 0.0));
    vector<vector<double>> pHat(p.Nx, vector<double>(p.Ny, 0.0));
    vector<vector<double>> Ax(p.Nx, vector<double>(p.Ny, 0.0));
    vector<vector<double>> Ap(p.Nx, vector<double>(p.Ny, 0.0));
    vector<vector<double>> ATpHat(p.Nx, vector<double>(p.Ny, 0.0));

    // Define restart frequency (e.g., reset directions every 100 steps)
    const int restart_interval = 500;

    // 1. Initial calculation: r = b - A*x
    Ax = matVecProduct(p, b, c, f.phi);
    for(int j = 1; j < p.Ny - 1; j++){
        for(int i = 1; i < p.Nx - 1; i++){
            if (b.mask[i][j] > 0.0) {
                r[i][j] = c.bP[i][j] - Ax[i][j];
                rHat[i][j] = r[i][j]; 
                pVec[i][j] = r[i][j];
                pHat[i][j] = rHat[i][j];
            }
        }
    }

    double rho_old = dotProduct(b, p, rHat, r);
    double rms = 0.0;
    int n;

    for (n = 1; n < s.itermax; n++) {
        
        // --- TRIGGER RESTART PERIODICALLY ---
        if (n % restart_interval == 0) {
            Ax = matVecProduct(p, b, c, f.phi);
            for(int j = 1; j < p.Ny - 1; j++){
                for(int i = 1; i < p.Nx - 1; i++){
                    // if (b.mask[i][j] > 0.0) {
                    //     r[i][j] = c.bP[i][j] - Ax[i][j];
                    //     rHat[i][j] = r[i][j]; 
                    //     pVec[i][j] = r[i][j];
                    //     pHat[i][j] = rHat[i][j];
                    // }

                    r[i][j] = c.bP[i][j] - Ax[i][j];
                    rHat[i][j] = r[i][j]; 
                    pVec[i][j] = r[i][j];
                    pHat[i][j] = rHat[i][j];
                }
            }
            rho_old = dotProduct(b, p, rHat, r);
        }

        // Ap = A * p
        Ap = matVecProduct(p, b, c, pVec);

        // alpha = (rHat, r) / (pHat, A * p)
        double denom = dotProduct(b, p, pHat, Ap);
        if (abs(denom) < 1e-15) {
            // If breakdown occurs, force an immediate manual restart 
            Ax = matVecProduct(p, b, c, f.phi);
            for(int j = 1; j < p.Ny - 1; j++){
                for(int i = 1; i < p.Nx - 1; i++){
                    // if (b.mask[i][j] > 0.0) {
                    //     r[i][j] = c.bP[i][j] - Ax[i][j];
                    //     rHat[i][j] = r[i][j]; 
                    //     pVec[i][j] = r[i][j];
                    //     pHat[i][j] = rHat[i][j];
                    // }

                    r[i][j] = c.bP[i][j] - Ax[i][j];
                    rHat[i][j] = r[i][j]; 
                    pVec[i][j] = r[i][j];
                    pHat[i][j] = rHat[i][j];
                }
            }
            rho_old = dotProduct(b, p, rHat, r);
            denom = dotProduct(b, p, pHat, Ap);
            if(abs(denom) < 1e-15) {
                cout << "BiCG Hard Breakdown even after forced restart!" << endl;
                break;
            }
        }
        
        double alpha = rho_old / denom;

        // Update solution and residual
        for (int j = 1; j < p.Ny - 1; j++) {
            for (int i = 1; i < p.Nx - 1; i++) {
                    f.phi[i][j] += alpha * pVec[i][j]*b.mask[i][j];
                    r[i][j] -= alpha * Ap[i][j]*b.mask[i][j];
                
            }
        }

        // Convergence Check
        rms = rmsNorm(b, p, r);
        s.error.push_back(rms);
        
        if (n % 1000 == 0) cout << "iter = " << n << " rms = " << rms << endl;
        
        if (rms <= s.tolerance) break;

        // rHat = rHat - alpha * A^T * pHat
        ATpHat = matVecTransposeProduct(b, p, c, pHat);
        for (int j = 1; j < p.Ny - 1; j++) {
            for (int i = 1; i < p.Nx - 1; i++) {
                //if (b.mask[i][j] > 0.0) rHat[i][j] -= alpha * ATpHat[i][j];
                
                rHat[i][j] -= alpha * ATpHat[i][j];
            }
        }

        // Calculate Beta
        double rho_new = dotProduct(b, p, rHat, r);
        double beta = rho_new / rho_old;
        rho_old = rho_new;

        // Update search directions
        for (int j = 1; j < p.Ny - 1; j++) {
            for (int i = 1; i < p.Nx - 1; i++) {
                // if (b.mask[i][j] > 0.0) {
                //     pVec[i][j] = r[i][j] + beta * pVec[i][j];
                //     pHat[i][j] = rHat[i][j] + beta * pHat[i][j];
                // }
                pVec[i][j] = r[i][j] + beta * pVec[i][j];
                pHat[i][j] = rHat[i][j] + beta * pHat[i][j];
            }
        }
    }

    cout << "\nPHI BiCG iterations = " << n << endl;
    if (n == s.itermax) cout << "Stopped due to itermax." << endl;
    else cout << "Converged. Final RMS: " << rms << endl;
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
            << i << " "
            << j << " "
            << m.x[i][j] << " "
            << m.y[i][j] << " "
            << c.aS[i][j] << " "
            << c.aW[i][j] << " "
            << c.aP[i][j] << " "
            << c.aE[i][j] << " "
            << c.aN[i][j] << " "
            << c.bP[i][j]<< endl;
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
    generateCoefficient(p,m,c);
    BCPhi(p,c,f);
    findBodyIndices(b,p,m);
    BCBodyPhi(b,p,c,f);
    buildMask(b,p);
    
    //SOR(s,b,p,m,c,f);
    //BiCG(s,b,p,m,c,f);
    Bi(s,b,p,m,c,f);
    exportSolverData(s,p,m,c,f,"BiCG_phi_error_converge.dat");
 
    exportPhi(p,m,f,"phi_solved.dat");

    return 0;
}