//
//  Plane.h
//  RayTracer
//
//  Created by WangYue on 16/9/21.
//  Copyright © 2016年 WangYue. All rights reserved.
//

#ifndef Plane_h
#define Plane_h

#include <stdio.h>
#include "Vector3.h"
#include "Object.h"

class Plane:  Object{
    
    Vector3 normal;
    double distance;// distance from the origin in 3D space, center of plane to the center of the scene
    Color color;
    
public:
    
    Plane(Vector3 &n, double d, Color &c){
        this->distance = d;
        this->normal = n.norm();
        this->color = color;
    };
    
    Plane(): distance(0.0), normal(Vector3(1.0,0.0,0.0)), color ( Color(0.5,0.5,0.5)){};
    Color getPlaneColor(){return color;}
    Vector3 getPlaneNormal(){return normal;}
    double getPlaneDistance(){return distance;}

    Vector3 getNormalAt(Vector3 &v){
        return normal;
    }
    
    virtual double findIntersection(Ray ray){
        //std::cout<<"L";
        // Return the intersection point as a form of Ray_origin + t * Ray_direction
        
        Vector3 ray_direction = ray.getDirection();
        double a = ray_direction.dot(normal);
        
        if(a == 0){
            // There is no intersection, ray is parellel to the plane.
            return -1;
        }
        else{
            // Ray : P = P0 + t * V;
            // Plane: P * N + d = 0;
            // t = - (d + P * N) / (V * N);
            double b = normal.dot(ray.getOrigin()) + distance;
            return -1 * (b / a);
        }
    }
    
    
};


#endif /* Plane_h */
