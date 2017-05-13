//
//  Triangle.h
//  RayTracer
//
//  Created by WangYue on 16/9/21.
//  Copyright © 2016年 WangYue. All rights reserved.
//

#ifndef Triangle_h
#define Triangle_h

#include <stdio.h>
#include "Vector3.h"
#include "Object.h"

class Triangle:  Object{
    Vector3 PointA, PointB, PointC;
    Vector3 normal;
    double distance;// distance from the origin in 3D space, center of Triangle to the center of the scene
    Color color;
    
    public:
    
    Triangle(Vector3 &a, Vector3 &b,Vector3 &c,Color &color){
        PointA = a;
        PointB = b;
        PointC = c;

        this->color = color;
    };
    
    Triangle(): PointA(Vector3(1,0,0)),PointB(Vector3(0,1,0)),PointC(Vector3(0,0,1)), color ( Color(0.5,0.5,0.5)){};
    
    Color getTriangleColor(){return color;}
   
    Vector3 getTriangleNormal(){
        Vector3 CA(PointA - PointC);
        Vector3 BA(PointA - PointB);
        normal = CA.cross(BA).norm();
        
        return normal;}
    double getTriangleDistance(){
        normal =getTriangleNormal();
        distance = normal.dot(PointA);
        return distance;}
    
    virtual Vector3 getNormalAt(Vector3 &v){
         normal =getTriangleNormal();
        return normal;
    }
    virtual Color getObjectColor(){
        return color;
    };
    virtual double findIntersection(Ray ray){
        //std::cout<<"L";
        // returns the distance between the triangle lies in
        // need to determine whether the intersection is inside the plane or outside the plane.
        
        Vector3 rayOrigin = ray.getOrigin();
        Vector3 rayDirection = ray.getDirection();
       
        double a = rayDirection.dot(normal);
        double distance = getTriangleDistance();
        
        if(a == 0){
            // There is no intersection, ray is parellel to the Triangle.
            return -1;
        }
        else{
            // Ray : P = P0 + t * V;
            // Triangle: P * N + d = 0;
            // t = - (d + P * N) / (V * N);
            double b = normal.dot(ray.getOrigin()) + distance;
            double distance2plane = - 1 * (b / a);
            // CAxQA * n >= 0
            // BAxQA * n >= 0
            // ABxQA * n >= 0
            //  Q is intersection
            Vector3 PointQ = rayDirection*distance2plane+rayOrigin;
            Vector3 CA = PointC - PointA;
            Vector3 QA = PointQ - PointA;
            double test1 = (CA.cross(QA)).dot(normal);
            
            Vector3 BC = PointB - PointC;
            Vector3 QC = PointQ - PointC;
            double test2 = (BC.cross(QC)).dot(normal);
            
            Vector3 AB = PointA - PointB;
            Vector3 QB = PointQ - PointB;
            double test3 = (AB.cross(QB)).dot(normal);
            
            if(test1>=0 && test2 >= 0 && test3>=0){
                return distance2plane;
            }
            else
                return -1;
            
          
        }
    }
    
    
};



#endif /* Triangle_h */
