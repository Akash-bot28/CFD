#include <iostream>
#include <fstream>
#include <vector>
#include <string>

using namespace std;

struct Domain{
    double xmin;
    double xmax;

    double ymin;
    double ymax;
};

struct Body{
    double xc;
    double yc;

    double length;
    double height;

    double xmin;
    double xmax;

    double ymin;
    double ymax;
};

struct Mesh{
    Domain domain;
    vector<Body> bodies;
};

void readInput(Mesh& mesh){
    ifstream input("input.txt");
    //if(!input){cerr << "Cannot open input.txt\n";exit(1);}

    string label;

    while(input >> label){
        if(label[0] == '#'){
            string dummy;
            getline(input,dummy);
            continue;
        }

        if(label == "Xmin") input >> mesh.domain.xmin;
        else if(label == "Xmax") input >> mesh.domain.xmax;      
        else if(label == "Ymin") input >> mesh.domain.ymin;       
        else if(label == "Ymax") input >> mesh.domain.ymax;       
        else if(label == "Bodies"){
            int nBodies;
            input >> nBodies;
            mesh.bodies.resize(nBodies);
            for(int i=0;i<nBodies;i++){
                string dummy;

                input >> dummy;  // xc
                input >> mesh.bodies[i].xc;

                input >> dummy;  // yc
                input >> mesh.bodies[i].yc;

                input >> dummy;  // length
                input >> mesh.bodies[i].length;

                input >> dummy;  // height
                input >> mesh.bodies[i].height;
            }
        }
    }

    input.close();
}

void buildBodyGeometry(Mesh& mesh){
    for(int i=0;i<mesh.bodies.size();i++){
        Body& body = mesh.bodies[i];
        body.xmin =body.xc - body.length/2.0;
        body.xmax =body.xc + body.length/2.0;
        body.ymin =body.yc - body.height/2.0;
        body.ymax =body.yc + body.height/2.0;
    }
}

void exportGeometryPreview(const Mesh& mesh, const string& filename){
    ofstream file(filename);

    //if(!file){cerr << "Cannot create geometry_preview.dat\n";return;}

    file<< mesh.domain.xmin << " "<< mesh.domain.xmax << " "<< mesh.domain.ymin << " "<< mesh.domain.ymax<< endl;
    file << mesh.bodies.size() << endl;
    for(int i=0;i<mesh.bodies.size();i++){
        file<< mesh.bodies[i].xmin << " "<< mesh.bodies[i].xmax << " "<< mesh.bodies[i].ymin << " "<< mesh.bodies[i].ymax<< endl;
    }
    file.close();
}

int main()
{
    Mesh mesh;

    readInput(mesh);

    buildBodyGeometry(mesh);

    exportGeometryPreview(mesh,"preview.dat");

    return 0;
}