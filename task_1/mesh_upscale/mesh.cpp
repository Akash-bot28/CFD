#include <iostream>
#include <vector>
#include <fstream>
#include <iomanip>
#include <string>
#include <cmath>
#include <algorithm>
#include <stdexcept>

using namespace std;

enum class SegmentType{
    Uniform,
    CompressBoth,
    StretchLeft,
    StretchRight
};

struct Segment{
    string name;  // seg0, seg1

    SegmentType type;  // uniform or streched

    double start;
    double end;

    int N;

    int startNode;
    int endNode;

    double beta;
};

struct Interface{
    int leftSegment;
    int rightSegment;
    double error;
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

struct Axis{
    string name;

    double domainMin;
    double domainMax;

    vector<Segment> segments;
    vector<Interface> interfaces;

    vector<double> coord;
    vector<double> d;

    vector<double> errorHistory;

    int itermax;
};

struct Optimization{
    double tolerance;
    double stepSize;
};

struct Mesh{
    Axis x;
    Axis y;

    vector<Body> bodies;

    Optimization opt;

    double betaStretchLeft;
    double betaUniform;
    double betaCompressBoth;
    double betaStretchRight;
};

SegmentType stringToSegmentType(const string& type){
    if(type == "StretchLeft")  return SegmentType::StretchLeft;
    if(type == "Uniform")      return SegmentType::Uniform;
    if(type == "CompressBoth") return SegmentType::CompressBoth;
    if(type == "StretchRight") return SegmentType::StretchRight;

    throw runtime_error("Unknown segment type: " + type);
}

void readInput(Mesh& mesh){
    ifstream input("input.txt");

    if(!input){cerr << "Cannot open input.txt\n";exit(1);}

    string label;

    int nBodies = 0;
    int xSegCount = 0;
    int ySegCount = 0;

    while(input >> label)
    {
        if(!label.empty() && label[0] == '#'){
        string dummy;
        getline(input,dummy);
        continue;
        }
        // Domain
        else if(label == "Xmin") input >> mesh.x.domainMin;
        else if(label == "Xmax") input >> mesh.x.domainMax;
        else if(label == "Ymin") input >> mesh.y.domainMin;
        else if(label == "Ymax") input >> mesh.y.domainMax;

        // Bodies
        else if(label == "Bodies"){
            input >> nBodies;
            mesh.bodies.resize(nBodies);
            for(int i=0; i<nBodies; i++){
                string dummy;

                input >> dummy; // xc
                input >> mesh.bodies[i].xc;

                input >> dummy; // yc
                input >> mesh.bodies[i].yc;

                input >> dummy; // length
                input >> mesh.bodies[i].length;

                input >> dummy; // height
                input >> mesh.bodies[i].height;
            }
        }

        // Segment beta initiaize
        else if(label == "betaStretchLeft") input >> mesh.betaStretchLeft;
        else if(label == "betaUniform") input >> mesh.betaUniform;
        else if(label == "betaCompressBoth") input >> mesh.betaCompressBoth;
        else if(label == "betaStretchRight") input >> mesh.betaStretchRight;

        // X Segments
        else if(label == "XSegments"){
            input >> xSegCount;
            mesh.x.segments.resize(xSegCount);
            for(int i=0; i<xSegCount; i++){
                string segName;
                string segType;
                int N;
                input >> segName >> segType >> N;

                mesh.x.segments[i].name = segName;
                mesh.x.segments[i].type = stringToSegmentType(segType);
                mesh.x.segments[i].N    = N;

                switch(mesh.x.segments[i].type)
                {
                    case SegmentType::StretchLeft:
                        mesh.x.segments[i].beta = mesh.betaStretchLeft;
                        break;

                    case SegmentType::Uniform:
                        mesh.x.segments[i].beta = mesh.betaUniform;
                        break;

                    case SegmentType::CompressBoth:
                        mesh.x.segments[i].beta = mesh.betaCompressBoth;
                        break;

                    case SegmentType::StretchRight:
                        mesh.x.segments[i].beta = mesh.betaStretchRight;
                        break;
                }
            }
        }
        // Y Segments
        else if(label == "YSegments"){
            input >> ySegCount;
            mesh.y.segments.resize(ySegCount);
            for(int i=0; i<ySegCount; i++){
                string segName;
                string segType;
                int N;
                input >> segName >> segType >> N;

                mesh.y.segments[i].name = segName;
                mesh.y.segments[i].type = stringToSegmentType(segType);
                mesh.y.segments[i].N    = N;

                switch(mesh.y.segments[i].type)
                {
                    case SegmentType::StretchLeft:
                        mesh.y.segments[i].beta = mesh.betaStretchLeft;
                        break;

                    case SegmentType::Uniform:
                        mesh.y.segments[i].beta = mesh.betaUniform;
                        break;

                    case SegmentType::CompressBoth:
                        mesh.y.segments[i].beta = mesh.betaCompressBoth;
                        break;

                    case SegmentType::StretchRight:
                        mesh.y.segments[i].beta = mesh.betaStretchRight;
                        break;
                }
            }
        }
        // Optimization
        else if(label == "itermaxX") input >> mesh.x.itermax;
        else if(label == "itermaxY") input >> mesh.y.itermax;
        else if(label == "tolerance") input >> mesh.opt.tolerance;
        else if(label == "step_size") input >> mesh.opt.stepSize;
    }

    input.close();
}

void buildBodyGeometry(Mesh& mesh){
    for(Body& body : mesh.bodies)
    {
        body.xmin = body.xc - body.length/2.0;
        body.xmax = body.xc + body.length/2.0;

        body.ymin = body.yc - body.height/2.0;
        body.ymax = body.yc + body.height/2.0;
    }
}

void buildSegments(Axis& axis, const vector<Body>& bodies){
    vector<double> boundaries;

    boundaries.push_back(axis.domainMin);

    for(const Body& body : bodies)
    {
        if(axis.name == "X")
        {
            boundaries.push_back(body.xmin);
            boundaries.push_back(body.xmax);
        }
        else if(axis.name == "Y")
        {
            boundaries.push_back(body.ymin);
            boundaries.push_back(body.ymax);
        }
        else
        {
            throw runtime_error("Invalid direction in buildSegments()");
        }
    }

    boundaries.push_back(axis.domainMax);

    sort(boundaries.begin(), boundaries.end());
    boundaries.erase(unique(boundaries.begin(), boundaries.end()),boundaries.end());

    int nIntervals = boundaries.size() - 1;

    if(static_cast<int>(axis.segments.size()) != nIntervals)
    {
        throw runtime_error(axis.name +": Number of geometric intervals (" + to_string(nIntervals) +") does not match number of user-defined segments (" + to_string(axis.segments.size()) + ")");
    }

    for(int i = 0; i < nIntervals; i++)
    {
        axis.segments[i].start = boundaries[i];
        axis.segments[i].end = boundaries[i+1];
    }

}

void buildInterfaces(Axis& axis){
    int nSegments = axis.segments.size();
    axis.interfaces.clear();
    for(int i = 0; i < nSegments - 1; i++){
        Interface interface;

        interface.leftSegment  = i;
        interface.rightSegment = i + 1;

        interface.error = 0.0;

        axis.interfaces.push_back(interface);
    }
}

double generateUniform(double zi, double L, double beta){
    double alpha = 0.5;
    double A = (beta + 1.0)/(beta - 1.0);
    double exponent = (zi - alpha)/(1.0 - alpha);
    double R = pow(A, exponent);
    double xi=L*((beta + 2.0*alpha)*R - beta + 2.0*alpha)/((1.0 + 2.0*alpha)*(1+ R));
    return xi;
}
double generateCompressBoth(double zi, double L, double beta){
    double alpha = 0.5;
    double A = (beta + 1.0)/(beta - 1.0);
    double exponent = (zi - alpha)/(1.0 - alpha);
    double R = pow(A, exponent);
    double xi=L*((beta + 2.0*alpha)*R - beta + 2.0*alpha)/((1.0 + 2.0*alpha)*(1+ R));
    return xi;
}
double generateStretchLeft(double zi, double L, double beta){
    double alpha = 0.0;
    double A = (beta + 1.0)/(beta - 1.0);
    double exponent = (zi - alpha)/(1.0 - alpha);
    double R = pow(A, exponent);
    double xi=L*((beta + 2.0*alpha)*R - beta + 2.0*alpha)/((1.0 + 2.0*alpha)*(1+ R));
    return xi;
}
double generateStretchRight(double zi, double L, double beta){
    double A = (beta + 1.0)/(beta - 1.0);
    double exponent = (1-zi);
    double R = pow(A, exponent);
    double xi= L*((beta+1)-(beta-1)*R)/(1+R);
    return xi;
}

vector<double> generateSegment(const Segment& segment){
    vector<double> coord(segment.N);
    double L =segment.end - segment.start;
    for(int i=0; i<segment.N; i++){
        double zi =double(i)/(segment.N - 1);
        double s;
        switch(segment.type){
            case SegmentType::Uniform:
                s = generateUniform(zi,L,segment.beta);
                break;
            case SegmentType::CompressBoth:
                s = generateCompressBoth(zi,L,segment.beta);
                break;
            case SegmentType::StretchLeft:
                s = generateStretchLeft(zi,L,segment.beta);
                break;
            case SegmentType::StretchRight:
                s = generateStretchRight(zi,L,segment.beta);
                break;
        }
        coord[i] =segment.start + s;
    }
    return coord;
}

void generateAxis(Axis& axis){
    axis.coord.clear();
    axis.d.clear();
    for(int i=0; i<axis.segments.size(); i++){
        Segment& seg = axis.segments[i];
        vector<double> segmentCoord = generateSegment(axis.segments[i]);
        if(i == 0){
            seg.startNode = 0;
            seg.endNode = seg.N - 1;
            axis.coord.insert(axis.coord.end(),segmentCoord.begin(),segmentCoord.end());
        }
        else{
            seg.startNode =axis.coord.size() - 1;
            seg.endNode =seg.startNode + seg.N - 1;
            axis.coord.insert(axis.coord.end(),segmentCoord.begin() + 1,segmentCoord.end());
        }
    }
    axis.d.resize(axis.coord.size()-1);
    for(int i=0; i<axis.coord.size()-1; i++) axis.d[i] = axis.coord[i+1] - axis.coord[i];
    
}

void computeInterfaceErrors(Axis& axis){
    for(int i=0; i<axis.interfaces.size(); i++){
        int leftSeg = axis.interfaces[i].leftSegment;
        int interfaceNode = axis.segments[leftSeg].endNode;
        double dLeft = axis.coord[interfaceNode]-axis.coord[interfaceNode-1];
        double dRight = axis.coord[interfaceNode+1]-axis.coord[interfaceNode];
        axis.interfaces[i].error = abs((dRight - dLeft)/dLeft);
    }
}

void updateBetas(Axis& axis, double stepSize){
    for(int i=0; i<axis.interfaces.size(); i++){
        int leftSeg =axis.interfaces[i].leftSegment;
        int rightSeg = axis.interfaces[i].rightSegment;
        int interfaceNode =axis.segments[leftSeg].endNode;
        double dLeft =axis.coord[interfaceNode] - axis.coord[interfaceNode-1];
        double dRight =axis.coord[interfaceNode+1] - axis.coord[interfaceNode];
        double error =axis.interfaces[i].error;
        if(dLeft > dRight){
            axis.segments[leftSeg].beta -= stepSize*error;
            axis.segments[leftSeg].beta = max(axis.segments[leftSeg].beta,1.000001);
            axis.segments[rightSeg].beta += stepSize*error;
            axis.segments[rightSeg].beta = max(axis.segments[rightSeg].beta,1.000001);
        }
        else{
            axis.segments[leftSeg].beta += stepSize*error;
            axis.segments[leftSeg].beta = max(axis.segments[leftSeg].beta,1.000001);
            axis.segments[rightSeg].beta -= stepSize*error;
            axis.segments[rightSeg].beta = max(axis.segments[rightSeg].beta,1.000001);
        }
    }
}

void betaOptimize(Axis& axis, double tolerance, double stepSize){
    axis.errorHistory.clear();
    //generateAxis(axis);
    int iter;
    for(iter=0; iter<axis.itermax; iter++){
        computeInterfaceErrors(axis);
        double maxError = 0.0;
        for(int i=0;i<axis.interfaces.size();i++){
            if(axis.interfaces[i].error > maxError){
                maxError = axis.interfaces[i].error;
            }
        }
        axis.errorHistory.push_back(maxError);
        if(maxError <= tolerance) break;
        updateBetas(axis,stepSize);
        generateAxis(axis);
        
    }
    cout << "\nOptimization iterations = "<< iter << endl;
    if(iter == axis.itermax) cout << "Stopped due to itermax." << endl;
    else cout << "Converged." << endl;
    cout << "\nOptimized betas for "<< axis.name<< " axis:\n";
    for(int i=0; i<axis.segments.size(); i++) cout<< setw(8)<< axis.segments[i].name<< "  "<< setw(15)<< axis.segments[i].beta<< endl;

}

void exportMesh(const Mesh& mesh,const string& filename){
    ofstream file(filename);
    //if(!file) { cerr << "Cannot open "<< filename<< endl; return;}
    file << fixed << setprecision(8);
    for(int j=0;j<mesh.y.coord.size();j++){
        for(int i=0;i<mesh.x.coord.size();i++){
            file<< mesh.x.coord[i]<< " "<< mesh.y.coord[j]<< endl;
        }
    }
    file.close();
}

void exportGridSpacing(const Mesh& mesh,const string& xFilename,const string& yFilename){
    ofstream dxFile(xFilename);
    ofstream dyFile(yFilename);
    //if(!dxFile){cerr << "Cannot open "<< xFilename<< endl;return;}
    //if(!dyFile){cerr << "Cannot open "<< yFilename<< endl;return;}
    dxFile << fixed << setprecision(8);
    dyFile << fixed << setprecision(8);

    for(int i=0;i<mesh.x.d.size();i++) dxFile<< mesh.x.coord[i]<< " "<< mesh.x.d[i]<< endl;
    for(int j=0;j<mesh.y.d.size();j++) dyFile<< mesh.y.coord[j]<< " "<< mesh.y.d[j]<< endl;

    dxFile.close();
    dyFile.close();
}

void exportOptimization(const Axis& axis){
    ofstream file(axis.name +"_optimization.dat");
    //if(!file){cerr << "Cannot open optimization file\n";return;}

    for(int i=0;i<axis.errorHistory.size();i++) file<< i<< " "<< axis.errorHistory[i]<< endl;
    
    file.close();
}

void exportSegmentBoundaries(const Mesh& mesh)
{
    ofstream file("interfaces.dat");

    file << "# X\n";

    for(int i=0; i<mesh.x.segments.size()-1; i++)
    {
        file << mesh.x.segments[i].end << endl;
    }

    file << "# Y\n";

    for(int i=0; i<mesh.y.segments.size()-1; i++)
    {
        file << mesh.y.segments[i].end << endl;
    }

    file.close();
}

int main(){
    Mesh mesh;

    readInput(mesh);

    mesh.x.name = "X";
    mesh.y.name = "Y";

    buildBodyGeometry(mesh);

    buildSegments(mesh.x, mesh.bodies);
    buildSegments(mesh.y, mesh.bodies);

    cout<<"segments along X axis"<<endl;
    for(int i=0;i<mesh.x.segments.size();i++){
        cout<< mesh.x.segments[i].name<< " "<< mesh.x.segments[i].start<< " "<< mesh.x.segments[i].end<< endl;
    }
    cout<<"segments along Y axis"<<endl;
    for(int j=0;j<mesh.y.segments.size();j++){
        cout<< mesh.y.segments[j].name<< " "<< mesh.y.segments[j].start<< " "<< mesh.y.segments[j].end<< endl;
    }

    buildInterfaces(mesh.x);
    buildInterfaces(mesh.y);

    generateAxis(mesh.x);
    generateAxis(mesh.y);

    betaOptimize(mesh.x,mesh.opt.tolerance,mesh.opt.stepSize);
    betaOptimize(mesh.y,mesh.opt.tolerance,mesh.opt.stepSize);

    exportOptimization(mesh.x);
    exportOptimization(mesh.y);

    exportMesh(mesh,"mesh.dat");

    exportGridSpacing(mesh,"dx.dat","dy.dat");

    exportSegmentBoundaries(mesh);

    return 0;

}