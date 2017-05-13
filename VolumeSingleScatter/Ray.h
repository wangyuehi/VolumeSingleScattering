//
//  Ray.h
//  RayTracer
//
//  Created by WangYue on 16/9/20.
//  Copyright © 2016年 WangYue. All rights reserved.
//

#ifndef Ray_h
#define Ray_h


#include "Vector3.h"

class Ray
{
    Vector3 origin;
    Vector3 direction;
    double time;

public:
     Ray(Vector3 &o,  Vector3 &dir, double time);
    
    Ray(Vector3 &o,  Vector3 &dir);
    Ray();
    
    Vector3 getOrigin(){return origin;};
    Vector3 getDirection(){return direction.norm();};
    void setOrigin(const Vector3 &o){origin = o;};
    void setDirection(const Vector3 &dir){direction = dir;};
    double getTime(){return time;};
    void setTime(double t){time = t;}
    
};

Ray::Ray(){
    origin = Vector3(0,0,0);
    direction = Vector3(1,0,0);
    time = -1;
    
}

Ray::Ray( Vector3 &o,  Vector3 &d){
    origin = o;
    direction = d.norm();
    time = -1;
}

Ray::Ray( Vector3 &o,  Vector3 &d, double t){
    origin = o;
    direction = d.norm();
    time = t;
}


#endif /* Ray_h */


