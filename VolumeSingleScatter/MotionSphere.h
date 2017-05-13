//
//  MotionSphere.h
//  RayTracer
//
//  Created by WangYue on 2016-10-13.
//  Copyright © 2016 WangYue. All rights reserved.
//

#ifndef MotionSphere_h
#define MotionSphere_h

#include <stdio.h>
#include "Vector3.h"
#include "Sphere.h"


class MotionSphere : public Sphere{
 
    Vector3 direction;
    double time;
    Vector3 OrigCenter;
public:
    
    MotionSphere(double r,
                 Vector3 &center,
                 Color &color,
                 Vector3 dir,
                 double time)
    :Sphere(r, center, color),time(time),direction(dir), OrigCenter(center){
    };
    
    void updateSphere(double t){ if(t>0)center = OrigCenter+ direction * t;}
    
    MotionSphere(): time(0.0){};

    virtual Color getObjectColor(){ return Sphere::getObjectColor();}
    
    virtual Vector3 getSphereCenter(){return center;}
    
    virtual double getSphereRadius(){return r;}
    
    virtual int getExponent(){return exponent;}
    
    virtual Vector3 getNormalAt(Vector3 &point){ return Sphere::getNormalAt(point);}
    
    virtual double findIntersection(Ray ray){
       // std::cout<<"previous r "<<r<<std::endl;
        updateSphere(ray.getTime());
        return Sphere::findIntersection(ray);}
    
    virtual bool insideObject(Vector3 lightPos){ return Sphere::insideObject(lightPos);  }


};

#endif /* MotionSphere_h */
