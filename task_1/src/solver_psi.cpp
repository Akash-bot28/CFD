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
    vector<vector<double>> psi;
};
struct SolverParameter{
    int itermax;
    double tolerance;
    double W;
    vector<double> error;
};

void readInput(Domain& d,Body& b,MeshParameter& p, Mesh& m, SolverParameter& s){
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
    input.close();
    
    // Generating mesh from data
    ifstream file("mesh.dat");
    file>>p.Nx>>p.Ny;
    p.N=p.Nx*p.Ny;
    m.x.assign(p.Nx,vector<double>(p.Ny));
    m.y.assign(p.Nx,vector<double>(p.Ny));
    for(int j=0;j<p.Ny;j++){
        for(int i=0;i<p.Nx;i++){
            file >> m.x[i][j]>> m.y[i][j];
        }
    }
    file.close();  
}
void BoundaryValues(const Domain& d, const MeshParameter& p,const Mesh& m,Field& f){
    double psi0 = 0.0;
    double U = 100.0;

    f.psi.assign(p.Nx,vector<double>(p.Ny,0.0));
       
    for(int j=0; j<p.Ny; j++) f.psi[0][j] = U*m.y[0][j] + psi0; // Left boundary (inlet)
    for(int j=0; j<p.Ny; j++) f.psi[p.Nx-1][j] = U*m.y[p.Nx-1][j] + psi0; // Right boundary (outlet)
    for(int i=0; i<p.Nx; i++) f.psi[i][0] = U*d.YB + psi0;  // Bottom boundary
    for(int i=0; i<p.Nx; i++) f.psi[i][p.Ny-1] = U*d.YT + psi0; // Top boundary
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
void CoefficientBC(const MeshParameter& p,Coefficient& c,Field& f){
    // Adjacent Left boundary (Dirichlet)
    for(int j=1; j<=p.Ny-2; j++) {int i = 1; c.bP[i][j] = c.bP[i][j] - c.aW[i][j]*f.psi[i-1][j]; c.aW[i][j] = 0.0;}
    // Adjacent Right boundary (Neumann)(convective)
    //for(int j=1; j<=p.Ny-2; j++) {int i = p.Nx-2; c.aP[i][j] = c.aP[i][j] + c.aE[i][j]; c.aE[i][j] = 0.0;}
    // Adjacent Right boundary (Dirichlet)
    for(int j=1; j<=p.Ny-2; j++) {int i = p.Nx-2; c.bP[i][j] = c.bP[i][j] - c.aE[i][j]*f.psi[i+1][j];c.aE[i][j] = 0.0;}
    // Adjacent Bottom boundary (Dirichlet)
    for(int i=1; i<=p.Nx-2; i++) {int j = 1; c.bP[i][j] = c.bP[i][j] - c.aS[i][j]*f.psi[i][j-1]; c.aS[i][j] = 0.0;}
    // Adjacent Top boundary (Dirichlet)
    for(int i=1; i<=p.Nx-2; i++) {int j = p.Ny-2; c.bP[i][j] = c.bP[i][j] - c.aN[i][j]*f.psi[i][j+1]; c.aN[i][j] = 0.0;}
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
void CoefficientBCBody(const Body& b,const MeshParameter& p,Coefficient& c,Field& f){
    double U=100.0;
    double psi_body=U*(b.YTs+b.YBs)/2.0;
    // Left
    for(int j=b.jB-1 ;j<=b.jT+1;j++) {int i=(b.iL-1) ; c.bP[i][j] = c.bP[i][j] - c.aE[i][j]*psi_body; c.aE[i][j] = 0.0;}
    // Right
    for(int j=b.jB -1;j<=b.jT+1;j++) {int i=(b.iR+1) ; c.bP[i][j] = c.bP[i][j] - c.aW[i][j]*psi_body; c.aW[i][j] = 0.0;}
    // Bottom wall
    for(int i=b.iL-1 ;i<=b.iR+1 ;i++) {int j=(b.jB-1) ; c.bP[i][j] = c.bP[i][j] - c.aN[i][j]*psi_body; c.aN[i][j] = 0.0;}
    //Top wall
    for(int i=b.iL-1 ;i<=b.iR+1 ;i++) {int j=(b.jT+1) ; c.bP[i][j] = c.bP[i][j] - c.aS[i][j]*psi_body; c.aS[i][j] = 0.0;}

    for(int j = b.jB; j <= b.jT; j++) {
        for(int i = b.iL; i <= b.iR; i++) {
            c.aS[i][j] = 0.0;
            c.aW[i][j] = 0.0;
            c.aE[i][j] = 0.0;
            c.aN[i][j] = 0.0;
            c.aP[i][j] = 1.0;   // identity row: phi = 0
            c.bP[i][j] = 0.0;
        }
    }
}

void Jacobi(SolverParameter& s,const Body& b,const MeshParameter& p,const Coefficient& c,Field& f){
    double rms;
    const int Ninterior =(p.Nx-2)*(p.Ny-2);

    vector<vector<double>> psiNew = f.psi;

    int n;
    for(n=1; n<s.itermax; n++){
        for(int j=1; j<p.Ny-1; j++){
            for(int i=1; i<p.Nx-1; i++){
                //psiNew[i][j] = b.mask[i][j]*(c.bP[i][j]- c.aS[i][j]*f.psi[i][j-1]- c.aW[i][j]*f.psi[i-1][j]- c.aE[i][j]*f.psi[i+1][j]- c.aN[i][j]*f.psi[i][j+1])/c.aP[i][j];
                psiNew[i][j] = (c.bP[i][j]- c.aS[i][j]*f.psi[i][j-1]- c.aW[i][j]*f.psi[i-1][j]- c.aE[i][j]*f.psi[i+1][j]- c.aN[i][j]*f.psi[i][j+1])/c.aP[i][j];
            }
        }
        // Copy new solution
        for(int j=1; j<p.Ny-1; j++){
            for(int i=1; i<p.Nx-1; i++){
                f.psi[i][j] = psiNew[i][j];
            }
        }
        // Residual calculation
        int fluid_nodes=0;
        double sumofsquares = 0.0;

        for(int j=1; j<p.Ny-1; j++){
            for(int i=1; i<p.Nx-1; i++){
                double residual = c.bP[i][j]- c.aS[i][j]*f.psi[i][j-1]- c.aW[i][j]*f.psi[i-1][j]- c.aP[i][j]*f.psi[i][j]- c.aE[i][j]*f.psi[i+1][j]- c.aN[i][j]*f.psi[i][j+1];
                sumofsquares += b.mask[i][j]*residual*residual;
                //sumofsquares += residual*residual;
                fluid_nodes += b.mask[i][j];
            }
        }

        rms = sqrt(sumofsquares/fluid_nodes);
        s.error.push_back(rms);

        if(rms <= s.tolerance) break;
        if(n%100==0) cout<<"iter="<<n<<"  rms="<<rms<<endl;
    }

    cout << "\nPSI Jacobi Solver iterations = "<< n << endl;

    if(n == s.itermax)cout << "Stopped due to itermax."<< endl;
    else cout << "Converged."<< endl;
}
double rowProduct(const Coefficient& c,const vector<vector<double>>& x, int i, int j){
    return c.aS[i][j]*x[i][j-1] + c.aW[i][j]*x[i-1][j] + c.aP[i][j]*x[i][j] + c.aE[i][j]*x[i+1][j] + c.aN[i][j]*x[i][j+1];
}
void SOR(SolverParameter& s,const Body& b,const MeshParameter& p,const Coefficient& c,Field& f){
    
    cout<< "\nsolving Psi..."<<endl;
    double rms;
    const int Ninterior =(p.Nx-2)*(p.Ny-2);

    int n;
    for(n=1; n<s.itermax; n++){
        for(int j=1; j<p.Ny-1; j++){
            for(int i=1; i<p.Nx-1; i++){

                //if(i>=b.iL && i<=b.iR && j>=b.jB && j<=b.jT) continue;

                double Apsi_P = rowProduct(c,f.psi,i,j);
                //f.psi[i][j] += b.mask[i][j]*s.W*(c.bP[i][j]- Apsi_P)/c.aP[i][j];
                f.psi[i][j] += s.W*(c.bP[i][j]- Apsi_P)/c.aP[i][j];
            }
        }

        // Residual calculation
        int  active_nodes=0;
        double sumofsquares = 0.0;
        for(int j=1; j<p.Ny-1; j++){
            for(int i=1; i<p.Nx-1; i++){

                //if(i>=b.iL && i<=b.iR && j>=b.jB && j<=b.jT) continue;

                double Apsi_P = rowProduct(c,f.psi,i,j);
                double residual = c.bP[i][j]- Apsi_P;
                //sumofsquares += residual*residual; 
                sumofsquares += b.mask[i][j]*residual*residual;
                active_nodes += b.mask[i][j];
            }
        }
        rms = sqrt(sumofsquares/active_nodes);
        s.error.push_back(rms);

        if(rms <= s.tolerance) break;
        if(n%100==0) cout<<"iter="<<n<<"  rms="<<rms<<endl;
    }
    cout << "\nPSI SOR iterations = "<< n << endl;
    if(n == s.itermax) cout << "Stopped due to itermax." << endl;
    else cout << "Converged."<< endl;
}


vector<vector<double>> matVecProduct(const MeshParameter& p, const Coefficient& c, const vector<vector<double>>& x) {
    vector<vector<double>> Ax(p.Nx, vector<double>(p.Ny, 0.0));
    for (int j = 1; j < p.Ny - 1; j++) {
        for (int i = 1; i < p.Nx - 1; i++) {
            Ax[i][j] = c.aS[i][j] * x[i][j - 1] + 
                      c.aW[i][j] * x[i - 1][j] + 
                      c.aP[i][j] * x[i][j]     + 
                      c.aE[i][j] * x[i + 1][j] + 
                      c.aN[i][j] * x[i][j + 1];
        }
    }
    return Ax;
}
double dotProduct(const MeshParameter& p, const vector<vector<double>>& x, const vector<vector<double>>& y) {
    double sum = 0.0;
    for (int j = 1; j < p.Ny - 1; j++) {
        for (int i = 1; i < p.Nx - 1; i++) {
            sum += x[i][j] * y[i][j];
        }
    }
    return sum;
}
double rmsNorm(const Body& b, const MeshParameter& p, const vector<vector<double>>& r) {
    double sum = 0.0;
    double fluid_nodes = 0.0;
    for (int j = 1; j < p.Ny - 1; j++) {
        for (int i = 1; i < p.Nx - 1; i++) {
            sum += b.mask[i][j] * r[i][j] * r[i][j];
            fluid_nodes += b.mask[i][j];
        }
    }
    if (fluid_nodes == 0.0) return 0.0;
    return sqrt(sum / fluid_nodes);
}

void BiCGSTAB(SolverParameter& s, const Body& b, const MeshParameter& p, const Mesh& m, const Coefficient& c, Field& f) {
    cout << "\nSolving Psi using BiCGSTAB..." << endl;
    s.error.clear();

    // Initialization
    // Compute initial residual: r0 = bP - A*psi
    vector<vector<double>> Ax = matVecProduct(p, c, f.psi);
    vector<vector<double>> r0(p.Nx, vector<double>(p.Ny, 0.0));
    for (int j = 1; j < p.Ny - 1; j++) {
        for (int i = 1; i < p.Nx - 1; i++) {
            r[i][j] = c.bP[i][j] - Ax[i][j];
        }
    }

    // Shadow residual vector r_star (initialized to r0)
    vector<vector<double>> r_star = r;
    vector<vector<double>> p_vec = r;

    double rms = rmsNorm(b, p, r);
    s.error.push_back(rms);

    // Allocate memory for temporary array fields
    vector<vector<double>> v(p.Nx, vector<double>(p.Ny, 0.0));
    vector<vector<double>> s_vec(p.Nx, vector<double>(p.Ny, 0.0));
    vector<vector<double>> t(p.Nx, vector<double>(p.Ny, 0.0));

    double rho_prev = 1.0;
    double alpha = 1.0;
    double omega = 1.0;

    int n;
    for(n=1; n<s.itermax; n++){
        double rho = dotProduct(p, r_star, r);
        //modifie search trajectory p_vec after first step
        if (n > 1) {
            double beta = (rho / rho_prev) * (alpha / omega);
            for (int j = 1; j < p.Ny - 1; j++) {
                for (int i = 1; i < p.Nx - 1; i++) {
                    p_vec[i][j] = r[i][j] + beta * (p_vec[i][j] - omega * v[i][j]);
                }
            }
        }

        // 1. Compute v = A*p
        v = matVecProduct(p, c, p_vec);

        // 2. Compute alpha
        alpha = rho / dotProduct(p, r_star, v);

        // 3. temporary residual s = r - alpha * v
        for (int j = 1; j < p.Ny - 1; j++) {
            for (int i = 1; i < p.Nx - 1; i++) {
                s_vec[i][j] = r[i][j] - alpha * v[i][j];
            }
        }

        //Check Intermediate Residual s
        double s_norm = rmsNorm(b, p, s_vec);
        if (s_norm <= s.tolerance) {
            for (int j = 1; j < p.Ny - 1; j++) {
                for (int i = 1; i < p.Nx - 1; i++) {
                    f.psi[i][j] += alpha * p_vec[i][j];
                }
            }
            rms = s_norm;
            s.error.push_back(rms);
            break;
        }

        // if not converged yet the need to stabalize
        // 4. Compute t = A*s
        t = matVecProduct(p, c, s_vec);

        // 5. Compute omega = (t, s) / (t, t)
        omega = dotProduct(p, t, s_vec) / dotProduct(p, t, t);

        // Update solution: psi = psi + alpha * p + omega * s
        // Update residual: r = s - omega * t
        for (int j = 1; j < p.Ny - 1; j++) {
            for (int i = 1; i < p.Nx - 1; i++) {
                f.psi[i][j] += alpha * p_vec[i][j] + omega * s_vec[i][j];
                r[i][j] = s_vec[i][j] - omega * t[i][j];
            }
        }

        rms = rmsNorm(b, p, r);
        s.error.push_back(rms);

        if (n % 10 == 0) { cout << "iter=" << n << "  rms=" << rms << endl;}

        if (rms <= s.tolerance) break;

        rho_prev = rho;
    }

    cout << "\nPSI BiCGSTAB Solver iterations = " << n << " | Final RMS = " << rms << endl;
    if (n == s.itermax) cout << "Stopped due to itermax." << endl;
    else cout << "Converged." << endl;

}

void exportPsi(const MeshParameter& p,const Mesh& m,const Field& f,string filename){
    ofstream file(filename);
    file << fixed << setprecision(8);
    file << p.Nx << " " << p.Ny << endl;
    for(int j=0; j<p.Ny; j++){
        for(int i=0; i<p.Nx; i++){
            file<< m.x[i][j]<<" "<< m.y[i][j]<<" "<< f.psi[i][j] << endl;
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
    Domain            d;
    Body              b;
    MeshParameter     p;
    Mesh              m;
    Coefficient       c;
    Field             f;
    SolverParameter   s;
    //-------------------------------------------------------------
    readInput(d,b,p,m,s);
    BoundaryValues(d,p,m,f);
    generateCoefficient(p,m,c);
    CoefficientBC(p,c,f);
    findBodyIndices(b,p,m);
    CoefficientBCBody(b,p,c,f);
    //-------------------------------------------------------------
    buildMask(b,p);

    //Jacobi(s,b,p,c,f);
    SOR(s,b,p,c,f);
    //BiCGSTAB(s,b,p,m,c,f);

    exportSolverData(s,p,m,c,f,"psi_SOR_converge.dat");
    exportPsi(p,m,f,"psi_SOR_solved.dat");
    return 0;
}