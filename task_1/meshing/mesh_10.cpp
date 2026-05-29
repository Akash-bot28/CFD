#include <iostream>
#include <vector>
#include <fstream>
#include <iomanip>
#include <string>
#include <cmath>

using namespace std;

struct Axis{
    string name;
    //geometry
    double L1, L2;
    double center;
    double size;
    //boundaries
    double domainMin, domainMax;
    double sMin, sMax;
    //mesh parameters
    int N1, Ns, N2, N;
    //stretching and compressing
    double alpha1, alpha2;
    double beta1, betas, beta2;
    //mesh
    vector<double> coord;
    vector<double> d;
    //grid stretching optimization
    vector<double> error;

    int itermax;
    double tolerance;
    double dbeta;

    vector<double> beta1_hist;
    vector<double> betas_hist;
    vector<double> beta2_hist;

    vector<double> error1_hist;
    vector<double> error2_hist;
};

struct Mesh{
    Axis x;
    Axis y;

};

void readInput(Mesh& m){
    string label;
    ifstream input("input.txt");
    while(input >> label)
    {
        if(label == "L1") input >> m.x.L1;
        else if(label == "L2") input >> m.x.L2;
        else if(label == "XC") input >> m.x.center;
        else if(label == "Nx1") input >> m.x.N1;
        else if(label == "Nxs") input >> m.x.Ns;
        else if(label == "Nx2") input >> m.x.N2;
        else if(label == "betax1") input >> m.x.beta1;
        else if(label == "betaxs") input >> m.x.betas;
        else if(label == "betax2") input >> m.x.beta2;

        else if(label == "h1") input >> m.y.L1;
        else if(label == "h2") input >> m.y.L2;
        else if(label == "YC") input >> m.y.center;
        else if(label == "Ny1") input >> m.y.N1;
        else if(label == "Nys") input >> m.y.Ns;
        else if(label == "Ny2") input >> m.y.N2;
        else if(label == "betay1") input >> m.y.beta1;
        else if(label == "betays") input >> m.y.betas;
        else if(label == "betay2") input >> m.y.beta2;
        

        else if(label == "a"){ input >> m.x.size; m.y.size=m.x.size;}
        else if(label == "alpha1") {input >> m.x.alpha1; m.y.alpha1=m.x.alpha1;}
        else if(label == "alpha2") {input >> m.x.alpha2; m.y.alpha2=m.x.alpha2;}
    
        else if(label == "itermax") {input >> m.x.itermax;m.y.itermax = m.x.itermax;}
        else if(label == "tolerance") {input >> m.x.tolerance;m.y.tolerance = m.x.tolerance;}
        else if(label == "dbeta"){input >> m.x.dbeta;m.y.dbeta = m.x.dbeta;}
    }
    input.close();
}

void buildGeometry(Axis& a){
    a.sMin=a.center-a.size/2.0;
    a.sMax=a.center+a.size/2.0;

    a.domainMin=a.sMin-a.L1;
    a.domainMax=a.sMax+a.L2;
}

double stretchLeft(double zi, double L, double alpha, double beta){
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

void generateAxis(Axis& a){
    a.N=a.N1+a.Ns+a.N2-2;
    a.coord.resize(a.N);
    a.d.resize(a.N-1);

    for(int i=0; i<a.N1; i++){
        double z=double(i)/(a.N1 - 1);
        double s=stretchLeft(z,a.L1,a.alpha1,a.beta1);
        a.coord[i]=a.domainMin + s;
    }

    for(int i=1; i<a.Ns; i++){
        double z=double(i)/(a.Ns - 1);
        double s=stretchLeft(z,a.size,a.alpha2,a.betas);
        a.coord[a.N1-1+i] =a.sMin + s;
    }

    for(int i=1; i<a.N2; i++){
        double z=double(i)/(a.N2 - 1);
        double s=stretchRight(z,a.L2,a.beta2);
        a.coord[a.N1+a.Ns-2+i]=a.sMax + s;
    }

    for(int i=0; i<a.N-1; i++) a.d[i]=a.coord[i+1] - a.coord[i];
}

void betaOptimize(Axis& a){
    a.error.resize(2);
    generateAxis(a);
    a.error[0]= (a.d[a.N1-1]- a.d[a.N1-2])/a.d[a.N1-2];
    a.error[1]= (a.d[a.N1+a.Ns-2]- a.d[a.N1+a.Ns-3])/a.d[a.N1+a.Ns-3];

    int iter = 0;
    while((abs(a.error[0])>a.tolerance || abs(a.error[1])>a.tolerance) && iter < a.itermax){
        a.beta1_hist.push_back(a.beta1);
        a.betas_hist.push_back(a.betas);
        a.beta2_hist.push_back(a.beta2);
        a.error1_hist.push_back(a.error[0]);
        a.error2_hist.push_back(a.error[1]);
        
        if(a.d[a.N1-2]>a.d[a.N1-1]){
            a.beta1=a.beta1-a.dbeta;
            a.betas=a.betas+a.dbeta;
        }
        else{
            a.beta1=a.beta1+a.dbeta;
            a.betas=a.betas-a.dbeta;
            
        }


        if(a.d[a.N1+a.Ns-3]>a.d[a.N1+a.Ns-2]){
            a.betas=a.betas-a.dbeta;
            a.beta2=a.beta2+a.dbeta;
        }
        else{
            a.betas=a.betas+a.dbeta;
            a.beta2=a.beta2-a.dbeta;
            
        }

        generateAxis(a);

        a.error[0]= (a.d[a.N1-1]- a.d[a.N1-2])/a.d[a.N1-2];
        a.error[1]= (a.d[a.N1+a.Ns-2]- a.d[a.N1+a.Ns-3])/a.d[a.N1+a.Ns-3];


        iter++;

    }
    cout << "\nOptimization iterations = "<< iter<< endl;
    cout<<"optimized betas of "<<a.name<<" axis are:"<< endl;
    cout<<"beta1 = "<<a.beta1<<"\n";
    cout<<"betas = "<<a.betas<<"\n";
    cout<<"beta2 = "<<a.beta2<<"\n";
    cout<<"block interface errors are = "<<a.error[0]<<" and "<<a.error[1]<<endl;
    
}

void exportMesh(const Mesh& m, string filename){
    ofstream file(filename);

    file<< fixed << setprecision(8);

    for(int j=0; j<m.y.N; j++)
    {
        for(int i=0; i<m.x.N; i++){
            file <<m.x.coord[i]<<" "<<m.y.coord[j]<<endl;
        }
    }

    file.close();
}

void exportGridSpacing(const Mesh& m,string xfilename, string yfilename){
    ofstream dxfile(xfilename);
    ofstream dyfile(yfilename);
    dxfile<< fixed<< setprecision(8);
    dyfile<< fixed<< setprecision(8);
    for(int i=0; i<m.x.N-1; i++) dxfile<< m.x.coord[i]<< " "<< m.x.d[i]<< endl;
    for(int j=0; j<m.y.N-1; j++) dyfile<< m.y.coord[j]<< " "<< m.y.d[j]<< endl;
    dxfile.close();
    dyfile.close();
    
}

void exportOptimization(const Axis& a){
    ofstream file(a.name + "_optimization.dat");

    for(size_t i=0; i<a.beta1_hist.size(); i++)
    {
        file
        << i << " "
        << a.beta1_hist[i] << " "
        << a.betas_hist[i] << " "
        << a.beta2_hist[i] << " "
        << a.error1_hist[i] << " "
        << a.error2_hist[i]
        << endl;
    }
}


int main(){
    Mesh mesh;
    readInput(mesh);
    mesh.x.name = "X";
    mesh.y.name = "Y";

    buildGeometry(mesh.x);
    buildGeometry(mesh.y);

    betaOptimize(mesh.x);
    betaOptimize(mesh.y);

    generateAxis(mesh.x);
    generateAxis(mesh.y);
    

    exportMesh(mesh,"mesh_10.dat");
    exportGridSpacing(mesh,"dx_10.dat","dy_10.dat");
    exportOptimization(mesh.x);
    exportOptimization(mesh.y);

    return 0;
}









