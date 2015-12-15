/*
 * Copyright: (C) 2015 iCub Facility - Istituto Italiano di Tecnologia
 * Authors: Ugo Pattacini
 * CopyPolicy: Released under the terms of the GNU GPL v2.0.
*/

#include <limits>
#include <deque>
#include <yarp/sig/all.h>

#include "optimizer.h"


/*******************************************************************************/
struct ParametersPSO : public Parameters
{
    int    numParticles;
    int    maxIter;
    double maxT;
    double omega;
    double phi_p;
    double phi_g;
    double cost;

    ParametersPSO() : Parameters(),
                      numParticles(20),
                      maxIter(std::numeric_limits<int>::max()),
                      maxT(std::numeric_limits<double>::infinity()),
                      omega(0.8),
                      phi_p(0.1),
                      phi_g(0.1),
                      cost(0.0) { }
};


/*******************************************************************************/
struct Particle
{
    yarp::sig::Vector pos;
    yarp::sig::Vector vel;
    double cost;    
    
    Particle() : pos(6,0.0), vel(6,0.0),
                 cost(std::numeric_limits<double>::infinity()) { }
};


/*******************************************************************************/
class Swarm : public Optimizer
{
    std::deque<Particle> x,p;
    Particle             g;
    
    yarp::sig::Vector rand_min,rand_max;
    int iter;
    double t,t0;
    
    void randomize();
    double evaluate(Particle &particle);
    void print(const bool randomize_print=false);
    
public:
    Swarm();
    ParametersPSO &get_parameters();
    void init();
    bool step();
    yarp::sig::Vector finalize();
};
