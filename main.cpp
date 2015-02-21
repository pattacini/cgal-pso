/*
 * Copyright: (C) 2015 iCub Facility - Istituto Italiano di Tecnologia
 * Authors: Ugo Pattacini
 * CopyPolicy: Released under the terms of the GNU GPL v2.0.
*/

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <deque>
#include <algorithm>
#include <cmath>

#include <yarp/os/all.h>
#include <yarp/sig/all.h>
#include <yarp/math/Math.h>
#include <yarp/math/Rand.h>

#include <iCub/ctrl/math.h>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_face_graph_triangle_primitive.h>

typedef CGAL::Simple_cartesian<double> K;
typedef K::FT FT;
typedef K::Point_3 Point;
typedef K::Segment_3 Segment;
typedef CGAL::Polyhedron_3<K> Polyhedron;
typedef CGAL::AABB_face_graph_triangle_primitive<Polyhedron> Primitive;
typedef CGAL::AABB_traits<K, Primitive> Traits;
typedef CGAL::AABB_tree<Traits> Tree;
typedef Tree::Point_and_primitive_id Point_and_primitive_id;

using namespace std;
using namespace yarp::os;
using namespace yarp::sig;
using namespace yarp::math;
using namespace iCub::ctrl;

/*******************************************************************************/
struct Parameters
{
    int    numParticles;
    int    maxIter;
    double omega;
    double phi_p;
    double phi_g;
    double cost;

    Vector x_lim;
    Vector y_lim;
    Vector z_lim;
    
    Parameters() : x_lim(3,0.0), y_lim(3,0.0), z_lim(3,0.0) { }
};


/*******************************************************************************/
struct Particle
{
    Vector pos;
    Vector vel;
    double cost;
    
    Particle() : pos(6,0.0), vel(6,0.0),
                 cost(numeric_limits<double>::infinity()) { }
};


/*******************************************************************************/
class Swarm
{
    deque<Point>    measurements;
    Polyhedron      model;
    Tree            tree;
    
    Parameters      parameters;
    deque<Particle> x,p;
    Particle        g;
    
    double cost;
    int iter;
    
    /***************************************************************************/
    void randomize()
    {
        for (size_t i=0; i<x.size(); i++)
        {
            Particle &particle=x[i];
            
            particle.pos[0]=Rand::scalar(parameters.x_lim[0],parameters.x_lim[1]);
            particle.pos[1]=Rand::scalar(parameters.y_lim[0],parameters.y_lim[1]);
            particle.pos[2]=Rand::scalar(parameters.z_lim[0],parameters.z_lim[1]);
            particle.pos[3]=Rand::scalar(-M_PI,M_PI);
            particle.pos[4]=Rand::scalar(-M_PI/2.0,M_PI/2.0);
            particle.pos[5]=Rand::scalar(-M_PI,M_PI);
            
            particle.vel[0]=Rand::scalar(-0.01,0.01);
            particle.vel[1]=Rand::scalar(-0.01,0.01);
            particle.vel[2]=Rand::scalar(-0.01,0.01);
            particle.vel[3]=Rand::scalar(-1.0,1.0)*CTRL_DEG2RAD;
            particle.vel[4]=Rand::scalar(-1.0,1.0)*CTRL_DEG2RAD;
            particle.vel[5]=Rand::scalar(-1.0,1.0)*CTRL_DEG2RAD;
        }
    }
    
    /***************************************************************************/
    double evaluate(Particle &particle)
    {
        Matrix H=rpy2dcm(particle.pos.subVector(3,6));
        H(0,3)=particle.pos[0];
        H(1,3)=particle.pos[1];
        H(2,3)=particle.pos[2];
        H=SE3inv(H);
        
        particle.cost=0.0;
        for (size_t i=0; i<measurements.size(); i++)
        {
            Point &m=measurements[i];
            double x=H(0,0)*m.x+H(0,1)*m.y+H(0,2)*m.z+H(0,3);
            double y=H(1,0)*m.x+H(1,1)*m.y+H(1,2)*m.z+H(1,3);
            double z=H(2,0)*m.x+H(2,1)*m.y+H(2,2)*m.z+H(2,3);
            Point p(x,y,z);
            particle.cost+=sqrt(tree.squared_distance(p));
        }
        
        return particle.cost;
    }
    
    /***************************************************************************/
    void print()
    {
        cout<<"iter #"<<iter<<": "
            <<"cost="<<cost<<" ("<<parameters.cost<<")"
            <<endl;
    }
    
public:
    /***************************************************************************/
    deque<Point> &get_measurements() { return measurements; }
    Polyhedron &get_model()          { return model; }
    Parameters &get_parameters()     { return parameters; }
    
    /***************************************************************************/
    void init()
    {
        // constructs AABB tree and computes internal KD-tree 
        // data structure to accelerate distance queries
        tree.insert(faces(model).first,faces(model).second,model);
        tree.accelerate_distance_queries();

        // create particles and init them randomly
        x.assign(parameters.numParticles,Particle());
        randomize();
        p=x;
        
        // evaluate the best particle g before starting
        for (size_t i=0; i<x.size(); i++)
            if (evaluate(p[i])<g.cost)
                g=p[i];
        
        iter=0;
        cost=numeric_limits<double>::infinity();
    }
    
    /***************************************************************************/
    bool step()
    {
        iter++;
        for (size_t i=0; i<x.size(); i++)
        {
            
            x[i].vel=omega*x[i].vel+phi_p*(p[i].pos-x[i].pos)+
                                    phi_g*(g.pos-x[i].pos);
        }
        
        bool term=(iter>=parameters.maxIter) ||
                  (cost<=parameters.cost);
        
        print();
        return term;
    }

    /***************************************************************************/
    void finalize()
    {
    }
};


/*******************************************************************************/
class PSOModule: public RFModule
{
    Swarm swarm;
    
    /***************************************************************************/
    void readMeasurements(ifstream &fin)
    {
        bool readHeader=true;
        char line[255];
        while (!fin.eof())
        {
            fin.getline(line,sizeof(line),'\n');
            Bottle b(line);
            
            if (readHeader)
            {
                if (b.get(0).asString()=="end_header")
                    readHeader=false;
                continue;
            }
                
            if (b.size()>=3)
                swarm.get_measurements().push_back(Point(b.get(0).asDouble(),
                                                   b.get(1).asDouble(),
                                                   b.get(2).asDouble()));
        }
    }
    
    /***************************************************************************/
    bool readLimits(ResourceFinder &rf, const string &tag, Vector &lim)
    {
        if (Bottle *b=rf.find(tag.c_str()).asList())
        {
            if (b->size()>=2)
            {
                lim[0]=b->get(0).asDouble();
                lim[1]=b->get(1).asDouble();
                return true;
            }
        }
        
        return false;
    }
    
public:
    /***************************************************************************/
    bool configure(ResourceFinder &rf)
    {
        if (!rf.check("modelFile"))
        {
            yError()<<"model file not provided!";
            return false;
        }
    
        if (!rf.check("measurementsFile"))
        {
            yError()<<"measurements file not provided!";
            return false;
        }
        
        string modelFileName=rf.find("modelFile").asString().c_str();
        string measurementsFileName=rf.find("measurementsFile").
                                    asString().c_str();

        Parameters &parameters=swarm.get_parameters();
        parameters.numParticles=rf.check("P",Value(20)).asInt();
        parameters.maxIter=rf.check("N",Value(100)).asInt();
        parameters.omega=rf.check("omega",Value(0.8)).asDouble();
        parameters.phi_p=rf.check("phi_p",Value(0.1)).asDouble();
        parameters.phi_g=rf.check("phi_g",Value(0.1)).asDouble();
        parameters.cost=rf.check("cost",Value(0.001)).asDouble();
        readLimits(rf,"x_lim",parameters.x_lim);
        readLimits(rf,"y_lim",parameters.y_lim);
        readLimits(rf,"z_lim",parameters.z_lim);
    
        // read the polyhedron from a .OFF file
        ifstream modelFile(modelFileName.c_str());
        if (!modelFile.is_open())
        {
            yError()<<"problem opening model file!";
            return false;
        }
        
        modelFile>>swarm.get_model();
        
        if (modelFile.bad())
        {
            yError()<<"problem reading model file!";
            modelFile.close();
            return false;
        }
        modelFile.close();
        
        // read the measurements file
        ifstream measurementsFile(measurementsFileName.c_str());
        if (!measurementsFile.is_open())
        {
            yError()<<"problem opening measurements file!";
            modelFile.close();
            return false;
        }
        readMeasurements(measurementsFile);
        measurementsFile.close();
        
        Rand::init();
        swarm.init();
    
        return true;
    }
    
    /***************************************************************************/
    double getPeriod()
    {
        return 0.0;
    }

    /***************************************************************************/
    bool updateModule()
    {
        return !swarm.step();
    }
    
    /***************************************************************************/
    bool close()
    {
        swarm.finalize();
        return true;
    }
};


/*******************************************************************************/
int main(int argc, char *argv[])
{
    ResourceFinder rf;
    rf.configure(argc,argv);
    
    PSOModule module;
    return module.runModule(rf);
}
