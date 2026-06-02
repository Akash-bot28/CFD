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
    int interfaceNode;
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

    //throw runtime_error("Unknown segment type: " + type);
}

void readInput(Mesh& mesh){
    ifstream input("input.txt");

    if(!input)
    {
        cerr << "Cannot open input.txt\n";
        exit(1);
    }

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

generateUniform(double zi, double L, double beta){
    double alpha = 0.5;
    double A = (beta + 1.0)/(beta - 1.0);
    double exponent = (zi - alpha)/(1.0 - alpha);
    double R = pow(A, exponent);
    double xi=L*((beta + 2.0*alpha)*R - beta + 2.0*alpha)/((1.0 + 2.0*alpha)*(1+ R));
    return xi;
}
generateCompressBoth(double zi, double L, double beta){
    double alpha = 0.5;
    double A = (beta + 1.0)/(beta - 1.0);
    double exponent = (zi - alpha)/(1.0 - alpha);
    double R = pow(A, exponent);
    double xi=L*((beta + 2.0*alpha)*R - beta + 2.0*alpha)/((1.0 + 2.0*alpha)*(1+ R));
    return xi;
}
generateStretchLeft(double zi, double L, double beta){
    double alpha = 0.0;
    double A = (beta + 1.0)/(beta - 1.0);
    double exponent = (zi - alpha)/(1.0 - alpha);
    double R = pow(A, exponent);
    double xi=L*((beta + 2.0*alpha)*R - beta + 2.0*alpha)/((1.0 + 2.0*alpha)*(1+ R));
    return xi;
}
generateStretchRight(double zi, double L, double beta){
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
        vector<double> segmentCoord = generateSegment(axis.segments[i]);
        if(i == 0) axis.coord.insert(axis.coord.end(),segmentCoord.begin(),segmentCoord.end());
        else axis.coord.insert(axis.coord.end(),segmentCoord.begin() + 1,segmentCoord.end());
        
    }
    axis.d.resize(axis.coord.size()-1);
    for(int i=0; i<axis.coord.size()-1; i++) axis.d[i] = axis.coord[i+1] - axis.coord[i];
    
}

void computeInterfaceErrors(Axis& axis){
    for(int i=0; i<axis.interfaces.size(); i++){
        int leftSeg =axis.interfaces[i].leftSegment;
        int rightSeg =axis.interfaces[i].rightSegment;
        int interfaceNode =axis.segments[leftSeg].endNode;
        double dLeft =axis.coord[interfaceNode]-axis.coord[interfaceNode-1];
        double dRight =axis.coord[interfaceNode+1]-axis.coord[interfaceNode];
        axis.interfaces[i].error =abs((dRight - dLeft)/dLeft);
    }
}

int main()
{
    Mesh mesh;

    readInput(mesh);

    buildBodyGeometry(mesh);

    buildSegments(mesh.x, mesh.bodies);
    buildSegments(mesh.y, mesh.bodies);

    buildInterfaces(mesh.x);
    buildInterfaces(mesh.y);

    generateAxis(mesh.x);
    generateAxis(mesh.y);

}