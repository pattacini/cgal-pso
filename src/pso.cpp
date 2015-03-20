/*
 * Copyright: (C) 2015 iCub Facility - Istituto Italiano di Tecnologia
 * Authors: Ugo Pattacini
 * CopyPolicy: Released under the terms of the GNU GPL v2.0.
*/

#include <iostream>
#include <string>
#include <algorithm>
#include <cmath>

#include <yarp/os/all.h>
#include <yarp/math/Math.h>
#include <yarp/math/Rand.h>

#include "pso.h"

#define DEG2RAD     (M_PI/180.0)

typedef CGAL::Aff_transformation_3<K> Affine;

using namespace std;
using namespace yarp::os;
using namespace yarp::sig;
using namespace yarp::math;


/*******************************************************************************/
Swarm::Swarm() : Optimizer()
{
    rand_min.resize(6,0.0);
    rand_max.resize(6,1.0);
    parameters=new ParametersPSO;   // owned by Optimizer
}


/*******************************************************************************/
ParametersPSO &Swarm::get_parameters()
{
    return *static_cast<ParametersPSO*>(parameters);
}


/*******************************************************************************/
void Swarm::randomize()
{
    for (size_t i=0; i<x.size(); i++)
    {
        Particle &particle=x[i];
        ParametersPSO &params=get_parameters();
        
        particle.pos[0]=Rand::scalar(params.x_lim[0],params.x_lim[1]);
        particle.pos[1]=Rand::scalar(params.y_lim[0],params.y_lim[1]);
        particle.pos[2]=Rand::scalar(params.z_lim[0],params.z_lim[1]);
        particle.pos[3]=Rand::scalar(-M_PI,M_PI);
        particle.pos[4]=Rand::scalar(-M_PI/2.0,M_PI/2.0);
        particle.pos[5]=Rand::scalar(-M_PI,M_PI);
        
        particle.vel[0]=Rand::scalar(-0.001,0.001);
        particle.vel[1]=Rand::scalar(-0.001,0.001);
        particle.vel[2]=Rand::scalar(-0.001,0.001);
        particle.vel[3]=Rand::scalar(-1.0,1.0)*DEG2RAD;
        particle.vel[4]=Rand::scalar(-1.0,1.0)*DEG2RAD;
        particle.vel[5]=Rand::scalar(-1.0,1.0)*DEG2RAD;
    }
}


/*******************************************************************************/
double Swarm::evaluate(Particle &particle)
{
    Matrix H=rpy2dcm(particle.pos.subVector(3,5));
    H(0,3)=particle.pos[0];
    H(1,3)=particle.pos[1];
    H(2,3)=particle.pos[2];
    H=SE3inv(H);
    
    particle.cost=0.0;
    for (size_t i=0; i<measurements.size(); i++)
    {
        Point &m=measurements[i];
        double x=H(0,0)*m[0]+H(0,1)*m[1]+H(0,2)*m[2]+H(0,3);
        double y=H(1,0)*m[0]+H(1,1)*m[1]+H(1,2)*m[2]+H(1,3);
        double z=H(2,0)*m[0]+H(2,1)*m[1]+H(2,2)*m[2]+H(2,3);
        particle.cost+=sqrt(tree.squared_distance(Point(x,y,z)));
    }
    
    if (measurements.size()>0)
        particle.cost/=measurements.size();
    
    return particle.cost;
}


/*******************************************************************************/
void Swarm::print(const bool randomize_print)
{
    ParametersPSO &params=get_parameters();
    cout<<"iter #"<<iter<<": "
        <<"cost="<<g.cost<<" ("<<params.cost<<"); ";
    if (randomize_print)
        cout<<"particles scattered away";
    cout<<endl;
}


/*******************************************************************************/
void Swarm::init()
{
    Optimizer::init();
    ParametersPSO &params=get_parameters();
    
    // create particles and init them randomly
    x.assign(params.numParticles,Particle());
    randomize();
    p=x;
    
    // evaluate the best particle g before starting
    for (size_t i=0; i<x.size(); i++)
        if (evaluate(p[i])<g.cost)
            g=p[i];
    
    iter=0;
}


/*******************************************************************************/
bool Swarm::step()
{
    iter++;
    ParametersPSO &params=get_parameters();
    
    for (size_t i=0; i<x.size(); i++)
    {
        Vector r1=Rand::vector(rand_min,rand_max);
        Vector r2=Rand::vector(rand_min,rand_max);
        
        x[i].vel=params.omega*x[i].vel+
                 params.phi_p*r1*(p[i].pos-x[i].pos)+
                 params.phi_g*r2*(g.pos-x[i].pos);
        
        x[i].pos+=x[i].vel;
        x[i].pos[0]=std::min(std::max(x[i].pos[0],params.x_lim[0]),
                             params.x_lim[1]);
        x[i].pos[1]=std::min(std::max(x[i].pos[1],params.y_lim[0]),
                             params.y_lim[1]);
        x[i].pos[2]=std::min(std::max(x[i].pos[2],params.z_lim[0]),
                             params.z_lim[1]);
        x[i].pos[3]=std::min(std::max(x[i].pos[3],-M_PI),M_PI);
        x[i].pos[4]=std::min(std::max(x[i].pos[4],-M_PI/2.0),M_PI/2.0);
        x[i].pos[5]=std::min(std::max(x[i].pos[5],-M_PI),M_PI);
        
        double f=evaluate(x[i]);
        if (f<p[i].cost)
        {
            p[i]=x[i];
            p[i].cost=f;
            if (f<g.cost)
                g=p[i];
        }
    }
    
    bool randomize_print=false;
    if ((iter%100)==0)
    {
        double mean=0.0;
        for (size_t i=0; i<x.size(); i++)
            mean+=norm(g.pos-x[i].pos);
        
        if (x.size()>0)
            mean/=x.size();
        
        if (mean<0.005)
        {
            randomize();
            randomize_print=true;
        }
    }
    
    bool term=(iter>=params.maxIter) || (g.cost<=params.cost);
    
    if ((iter%10)==0)
        print(randomize_print);
    
    return term;
}


/*******************************************************************************/
Vector Swarm::finalize()
{
    print();
    
    Matrix H=rpy2dcm(g.pos.subVector(3,5));
    Affine affine(H(0,0),H(0,1),H(0,2),g.pos[0],
                  H(1,0),H(1,1),H(1,2),g.pos[1],
                  H(2,0),H(2,1),H(2,2),g.pos[2]);
                           
    std::transform(model.points_begin(),model.points_end(),
                   model.points_begin(),affine);

    return g.pos;
}
    