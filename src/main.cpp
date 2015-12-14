/*
 * Copyright: (C) 2015 iCub Facility - Istituto Italiano di Tecnologia
 * Authors: Ugo Pattacini
 * CopyPolicy: Released under the terms of the GNU GPL v2.0.
*/

#include <iostream>
#include <fstream>
#include <string>

#include <yarp/os/all.h>
#include <yarp/sig/all.h>
#include <yarp/math/Rand.h>

#include <CGAL/IO/Polyhedron_iostream.h>

#include "pso.h"

using namespace std;
using namespace yarp::os;
using namespace yarp::sig;
using namespace yarp::math;


/*******************************************************************************/
class OptModule: public RFModule
{
    Swarm swarm;
    string outputFileName;
    double maxT,t0;

    /***************************************************************************/
    bool readMeasurements(ifstream &fin)
    {
        int state=0;
        int nPoints;
        char line[255];
        
        while (!fin.eof())
        {
            fin.getline(line,sizeof(line),'\n');
            Bottle b(line);
            Value firstItem=b.get(0);
            bool isNumber=firstItem.isInt() || firstItem.isDouble();
            
            if (state==0)
            {
                string tmp=firstItem.asString().c_str();
                std::transform(tmp.begin(),tmp.end(),tmp.begin(),::toupper);
                if (tmp=="OFF")
                    state++;
            }
            else if (state==1)
            {
                if (isNumber)
                {
                    nPoints=firstItem.asInt();
                    state++;
                }
            }
            else if (state==2)
            {
                if (isNumber && (b.size()>=3))
                {
                    swarm.get_measurements().push_back(Point(b.get(0).asDouble(),
                                                       b.get(1).asDouble(),
                                                       b.get(2).asDouble()));
                    
                    if (--nPoints<=0)
                        return true;
                }
            }
        }
        
        return false;
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
        outputFileName=rf.check("outputFile",Value("output.off")).
                       asString().c_str();

        ParametersPSO &parameters=swarm.get_parameters();
        parameters.numParticles=rf.check("P",Value(20)).asInt();
        parameters.maxIter=rf.check("N",Value(100)).asInt();
        parameters.omega=rf.check("omega",Value(0.8)).asDouble();
        parameters.phi_p=rf.check("phi_p",Value(0.1)).asDouble();
        parameters.phi_g=rf.check("phi_g",Value(0.1)).asDouble();
        parameters.cost=rf.check("cost",Value(0.001)).asDouble();        
        readLimits(rf,"x_lim",parameters.x_lim);
        readLimits(rf,"y_lim",parameters.y_lim);
        readLimits(rf,"z_lim",parameters.z_lim);

        maxT=rf.check("T",Value(10.0)).asDouble();
    
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
        if (!readMeasurements(measurementsFile))
        {
            yError()<<"problem reading measurements file!";
            modelFile.close();
            measurementsFile.close();
            return false;
        }
        measurementsFile.close();
        
        Rand::init();
        swarm.init();
    
        t0=Time::now();
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
        return (!swarm.step() && (Time::now()-t0<maxT));
    }
    
    /***************************************************************************/
    bool close()
    {
        double dt=Time::now()-t0;
        Vector g=swarm.finalize();
        
        cout<<"solution: "<<g.toString(3,3).c_str()<<endl;
        cout<<"found in "<<dt<<" [s]"<<endl;
        
        ofstream fout(outputFileName.c_str());
        fout<<swarm.get_model();
        fout.close();
        
        return true;
    }
};


/*******************************************************************************/
int main(int argc, char *argv[])
{
    ResourceFinder rf;
    rf.configure(argc,argv);
    
    OptModule module;
    return module.runModule(rf);
}
